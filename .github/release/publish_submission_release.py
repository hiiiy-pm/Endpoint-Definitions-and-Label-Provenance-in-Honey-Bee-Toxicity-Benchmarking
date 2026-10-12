"""Publish only hash-verified submission assets using the GitHub Releases API."""
from pathlib import Path
import hashlib
import json
import os
import re
import shutil
import urllib.error
import urllib.parse
import urllib.request

REPOSITORY = "hiiiy-pm/Endpoint-Definitions-and-Label-Provenance-in-Honey-Bee-Toxicity-Benchmarking"
TAG = "ieee-access-submission-final-2026-10-12"
ROOT = Path(__file__).resolve().parents[2]
ASSETS = ROOT / "release_assets"
API = "https://api.github.com/repos/" + REPOSITORY


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def api(method, url, payload=None, content_type="application/json"):
    parsed = urllib.parse.urlparse(url)
    require(parsed.scheme == "https" and parsed.hostname in {"api.github.com", "uploads.github.com"}, "Unexpected API destination")
    data = payload if isinstance(payload, bytes) else (None if payload is None else json.dumps(payload).encode())
    request = urllib.request.Request(url, data=data, method=method, headers={
        "Authorization": "Bearer " + os.environ["GITHUB_TOKEN"],
        "Accept": "application/vnd.github+json",
        "Content-Type": content_type,
        "X-GitHub-Api-Version": "2026-03-10",
        "User-Agent": "IEEE-Access-submission-release",
    })
    with urllib.request.urlopen(request, timeout=180) as response:
        return json.load(response)


def optional_get(url):
    try:
        return api("GET", url)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return None
        raise


def main():
    require(os.environ.get("GITHUB_REPOSITORY") == REPOSITORY, "Unexpected repository")
    require(os.environ.get("GITHUB_REF") == "refs/heads/main", "Only main may publish this snapshot")
    commit = os.environ["GITHUB_SHA"]
    require(re.fullmatch(r"[0-9a-f]{40}", commit) is not None, "Invalid target commit")
    metadata_path = Path(__file__).with_name("submission_zip_metadata.json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    require(metadata["schema_version"] == 2 and len(metadata["archives"]) == 4, "Invalid archive metadata")
    require(isinstance(metadata["source_member_count"], int) and metadata["source_member_count"] > 2, "Invalid source member count")
    report = json.loads((ASSETS / "archive_rebuild_verification.json").read_text(encoding="utf-8"))
    require(report["status"] == "PASS" and report["source_members_verified"] == metadata["source_member_count"], "Archive verification has not passed")
    require(report["research_manifest_entries_verified"] == metadata["research_manifest_entry_count"], "Research manifest verification has not passed")
    require(report["metadata_sha256"] == hashlib.sha256(metadata_path.read_bytes()).hexdigest(), "Verification used different metadata")
    require(len(report["archives"]) == 4, "Incomplete archive set")
    archive_specs = {a["name"]: a for a in metadata["archives"]}
    require(len(archive_specs) == 4 and {a["file"] for a in report["archives"]} == set(archive_specs), "Archive names differ from metadata")
    for item in report["archives"]:
        spec = archive_specs[item["file"]]
        require(item["bytes"] == spec["bytes"] and item["sha256"] == spec["sha256"] and item["byte_identical_to_submission"] is True, "Archive report differs from submission metadata")
        path = ASSETS / item["file"]
        require(path.stat().st_size == item["bytes"] and hashlib.sha256(path.read_bytes()).hexdigest() == item["sha256"], "Archive changed after verification")
    checksums = {}
    for line in (ROOT / "RESEARCH_ARCHIVE_SHA256SUMS.txt").read_text(encoding="utf-8-sig").splitlines():
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        require(match is not None, "Malformed Release checksum line")
        digest, name = match.groups()
        require(name not in checksums, "Duplicate Release checksum entry")
        checksums[name] = digest
    require(checksums == {name: spec["sha256"] for name, spec in archive_specs.items()}, "Release checksum text differs from submission metadata")
    for name in ("RELEASE_README.txt", "RESEARCH_ARCHIVE_SHA256SUMS.txt"):
        shutil.copyfile(ROOT / name, ASSETS / name)
    names = [item["file"] for item in report["archives"]] + ["RELEASE_README.txt", "RESEARCH_ARCHIVE_SHA256SUMS.txt"]
    expected = {name: {"size": (ASSETS / name).stat().st_size, "digest": "sha256:" + hashlib.sha256((ASSETS / name).read_bytes()).hexdigest()} for name in names}
    # Never move an existing tag or replace an asset with different content.
    tag = optional_get(API + "/git/ref/tags/" + TAG)
    if tag is not None:
        require(tag["object"]["type"] == "commit" and tag["object"]["sha"] == commit, "Existing tag points to a different snapshot")
    releases = api("GET", API + "/releases?per_page=100")
    found = [r for r in releases if r["tag_name"] == TAG]
    require(len(found) <= 1, "Ambiguous existing release")
    body = (Path(__file__).with_name("RELEASE_NOTES.txt")).read_text(encoding="utf-8")
    body += "\n\nRepository snapshot: https://github.com/" + REPOSITORY + "/commit/" + commit + "\n"
    if found:
        release = found[0]
        require(release["target_commitish"] == commit, "Existing release target differs")
    else:
        release = api("POST", API + "/releases", {
            "tag_name": TAG, "target_commitish": commit,
            "name": "IEEE Access final submission materials (12 October 2026)",
            "body": body, "draft": True, "prerelease": False,
        })
    release_id = release["id"]
    assets_url = API + f"/releases/{release_id}/assets?per_page=100"
    assets = {a["name"]: a for a in api("GET", assets_url)}
    require(set(assets).issubset(expected), "Unexpected assets in existing release")
    for name in names:
        if name in assets:
            actual = assets[name]
        else:
            require(release["draft"], "Refusing to change an already published release")
            upload = release["upload_url"].split("{", 1)[0] + "?" + urllib.parse.urlencode({"name": name})
            actual = api("POST", upload, (ASSETS / name).read_bytes(), "application/zip" if name.endswith(".zip") else "text/plain")
        require(actual["state"] == "uploaded" and actual["size"] == expected[name]["size"] and actual.get("digest") == expected[name]["digest"], "Remote asset hash/size mismatch: " + name)
        print("Verified release asset:", name, actual["digest"], flush=True)
    final_assets = api("GET", assets_url)
    require({a["name"] for a in final_assets} == set(expected), "Incomplete remote release")
    for asset in final_assets:
        require(asset["state"] == "uploaded" and asset.get("digest") == expected[asset["name"]]["digest"] and asset["size"] == expected[asset["name"]]["size"], "Remote verification failed")
    if release["draft"]:
        release = api("PATCH", API + f"/releases/{release_id}", {"draft": False, "prerelease": False, "make_latest": "true", "body": body})
    require(not release["draft"], "Release remains a draft")
    print("Published:", release["html_url"], flush=True)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a", encoding="utf-8") as handle:
            handle.write("Published verified submission materials: " + release["html_url"] + "\n\nAll four ZIP SHA-256 hashes match the final submission attachments. No experiments were run.\n")


if __name__ == "__main__":
    main()
