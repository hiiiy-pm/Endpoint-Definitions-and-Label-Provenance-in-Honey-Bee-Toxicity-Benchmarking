#!/usr/bin/env python3
"""Publish a fixed, verified research reproducibility release on GitHub."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

from build_assets import ASSET_NAMES, TAG, read_manifest, verify_assets, verify_sources

REPOSITORY = "hiiiy-pm/Endpoint-Definitions-and-Label-Provenance-in-Honey-Bee-Toxicity-Benchmarking"
TITLE = "Research reproducibility materials"
BODY = """Auditing label provenance and endpoint definitions in a honey-bee toxicity benchmark

Authors: Yongye Huang, Runqing Liao, Xiaopeng Yang, and Zhining Wang.

Research.zip contains the analysis and figure-generation code, input data, official splits, source-record evidence, saved models and predictions, statistical results, eleven final figure PDFs, environment specifications, and reproduction instructions. All files are supplied under Research/.
"""


class ApiError(RuntimeError):
    def __init__(self, status: int, method: str, path: str):
        self.status = status
        super().__init__(f"GitHub API {method} {path} returned HTTP {status}")


class SafeRedirectHandler(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        redirected = super().redirect_request(req, fp, code, msg, headers, newurl)
        if redirected is not None:
            if urllib.parse.urlsplit(newurl).scheme != "https":
                raise RuntimeError("Refusing an insecure asset download redirect")
            if urllib.parse.urlsplit(req.full_url).netloc != urllib.parse.urlsplit(newurl).netloc:
                redirected.remove_header("Authorization")
        return redirected


class GitHub:
    def __init__(self, token: str):
        self.token = token
        self.opener = urllib.request.build_opener(SafeRedirectHandler())

    def request(self, method: str, path: str, *, payload=None, raw: bytes | None = None,
                content_type="application/json", accept="application/vnd.github+json"):
        url = path if path.startswith("https://") else "https://api.github.com" + path
        parsed = urllib.parse.urlsplit(url)
        if parsed.scheme != "https" or parsed.netloc not in {"api.github.com", "uploads.github.com"}:
            raise RuntimeError("Refusing an unexpected GitHub API host")
        if payload is not None:
            raw = json.dumps(payload).encode("utf-8")
        headers = {
            "Authorization": f"Bearer {self.token}",
            "Accept": accept,
            "X-GitHub-Api-Version": "2022-11-28",
            "User-Agent": "research-release-publisher",
        }
        if raw is not None:
            headers["Content-Type"] = content_type
        request = urllib.request.Request(url, data=raw, method=method, headers=headers)
        attempts = 3 if method == "GET" else 1
        for attempt in range(attempts):
            try:
                return self.opener.open(request, timeout=180)
            except urllib.error.HTTPError as exc:
                if method == "GET" and exc.code in {429, 500, 502, 503, 504} and attempt + 1 < attempts:
                    exc.close()
                    time.sleep(2 ** (attempt + 1))
                    continue
                status = exc.code
                exc.close()
                raise ApiError(status, method, parsed.path) from None
            except urllib.error.URLError:
                if method == "GET" and attempt + 1 < attempts:
                    time.sleep(2 ** (attempt + 1))
                    continue
                raise RuntimeError(f"GitHub API {method} {parsed.path} could not be reached") from None
        raise RuntimeError("GitHub API retry limit reached")

    def json(self, method: str, path: str, **kwargs):
        with self.request(method, path, **kwargs) as response:
            data = response.read()
        return json.loads(data) if data else None


def release_path() -> str:
    return f"/repos/{REPOSITORY}/releases"


def find_release(api: GitHub):
    matches = []
    # The releases list includes drafts for authenticated callers; lookup by tag
    # is not sufficient for recovering a draft whose tag has not been created.
    for page in range(1, 101):
        batch = api.json("GET", f"{release_path()}?per_page=100&page={page}")
        matches.extend(release for release in batch if release.get("tag_name") == TAG)
        if len(matches) > 1:
            raise RuntimeError("Multiple releases use the intended tag; refusing an ambiguous publication")
        if len(batch) < 100:
            return matches[0] if matches else None
    raise RuntimeError("Unexpectedly large release collection")


def tag_commit(api: GitHub):
    try:
        reference = api.json("GET", f"/repos/{REPOSITORY}/git/ref/tags/{urllib.parse.quote(TAG, safe='')}")
    except ApiError as exc:
        if exc.status == 404:
            return None
        raise
    obj = reference["object"]
    for _ in range(10):
        if obj["type"] == "commit":
            return obj["sha"]
        if obj["type"] != "tag":
            break
        obj = api.json("GET", f"/repos/{REPOSITORY}/git/tags/{obj['sha']}")["object"]
    raise RuntimeError("Release tag does not resolve to a commit")


def verify_target(api: GitHub, release: dict | None, commit: str) -> None:
    tagged = tag_commit(api)
    if tagged is not None and tagged != commit:
        raise RuntimeError("The release tag already points to a different commit; it will not be moved")
    if release is not None:
        if release.get("tag_name") != TAG:
            raise RuntimeError("Unexpected release tag")
        if not release.get("draft") and tagged is None:
            raise RuntimeError("Published release is missing its expected tag")
        if tagged is None:
            target = str(release.get("target_commitish", ""))
            if target != commit:
                resolved = api.json("GET", f"/repos/{REPOSITORY}/commits/{urllib.parse.quote(target, safe='')}")
                if resolved.get("sha") != commit:
                    raise RuntimeError("Draft release targets a different commit")


def list_assets(api: GitHub, release_id: int) -> list[dict]:
    assets = []
    for page in range(1, 101):
        batch = api.json("GET", f"{release_path()}/{release_id}/assets?per_page=100&page={page}")
        assets.extend(batch)
        if len(batch) < 100:
            return assets
    raise RuntimeError("Unexpectedly large release asset collection")


def verify_remote_asset(api: GitHub, asset: dict, expected: dict) -> None:
    name = expected["name"]
    if asset.get("name") != name or asset.get("size") != expected["bytes"] or asset.get("state") != "uploaded":
        raise RuntimeError(f"Remote asset name, size, or upload state differs: {name}")
    digest = asset.get("digest")
    if digest:
        if digest != "sha256:" + expected["sha256"]:
            raise RuntimeError(f"Remote asset SHA-256 differs: {name}")
        return
    hasher = hashlib.sha256()
    length = 0
    with api.request("GET", f"{release_path()}/assets/{asset['id']}", accept="application/octet-stream") as response:
        for chunk in iter(lambda: response.read(1024 * 1024), b""):
            length += len(chunk)
            if length > expected["bytes"]:
                raise RuntimeError(f"Remote asset download exceeds expected byte count: {name}")
            hasher.update(chunk)
    if length != expected["bytes"] or hasher.hexdigest() != expected["sha256"]:
        raise RuntimeError(f"Downloaded remote asset differs from expected SHA-256: {name}")


def inspect_assets(api: GitHub, release: dict, expected: dict[str, dict], complete: bool) -> dict[str, dict]:
    assets = list_assets(api, release["id"])
    names = [asset["name"] for asset in assets]
    if len(names) != len(set(names)) or not set(names).issubset(expected):
        raise RuntimeError("Release contains duplicate or unexpected assets; no existing asset will be deleted")
    if complete and set(names) != set(expected):
        raise RuntimeError("Release does not contain the expected research archive")
    by_name = {asset["name"]: asset for asset in assets}
    for name, asset in by_name.items():
        verify_remote_asset(api, asset, expected[name])
    return by_name


def upload_asset(api: GitHub, release: dict, path: Path, expected: dict) -> None:
    upload_base = release["upload_url"].split("{", 1)[0]
    upload_url = upload_base + "?" + urllib.parse.urlencode({"name": path.name})
    kind = "application/pdf" if path.suffix == ".pdf" else "application/zip"
    for attempt in range(3):
        try:
            asset = api.json("POST", upload_url, raw=path.read_bytes(), content_type=kind)
            verify_remote_asset(api, asset, expected)
            return
        except (ApiError, RuntimeError):
            # An upload can succeed even when the response is lost. Inspect before retrying.
            present = [asset for asset in list_assets(api, release["id"]) if asset["name"] == path.name]
            if len(present) == 1:
                verify_remote_asset(api, present[0], expected)
                return
            if present:
                raise RuntimeError(f"Duplicate asset detected after upload: {path.name}")
            if attempt == 2:
                raise
            time.sleep(2 ** (attempt + 1))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--commit", default=os.environ.get("GITHUB_SHA", ""))
    args = parser.parse_args()
    repository = args.repository_root.resolve()
    manifest = read_manifest(args.manifest or repository / ".github/release/research_metadata.json")
    verify_sources(repository, manifest)
    verify_assets(args.output_dir, manifest)
    if not re.fullmatch(r"[0-9a-f]{40}", args.commit):
        raise ValueError("A full immutable source commit SHA is required")
    configured_repository = os.environ.get("GITHUB_REPOSITORY", REPOSITORY)
    if configured_repository != REPOSITORY:
        raise ValueError("Workflow repository does not match the intended publication repository")
    if os.environ.get("GITHUB_REF") != "refs/heads/main":
        raise ValueError("Publication is permitted only from the main branch")
    token = os.environ.get("GITHUB_TOKEN", "").strip()
    if not token:
        raise ValueError("GITHUB_TOKEN is required")
    api = GitHub(token)
    # Require the intended source commit to exist before creating any release.
    if api.json("GET", f"/repos/{REPOSITORY}/commits/{args.commit}").get("sha") != args.commit:
        raise RuntimeError("Source commit could not be verified")
    release = find_release(api)
    verify_target(api, release, args.commit)
    expected = {record["name"]: record for record in manifest["assets"]}
    if release is not None and not release.get("draft"):
        inspect_assets(api, release, expected, complete=True)
        if release.get("name") != TITLE or release.get("body", "").strip() != BODY.strip():
            raise RuntimeError("Existing published release metadata differs; it will not be overwritten")
        if release.get("prerelease"):
            raise RuntimeError("Existing release is unexpectedly marked as a prerelease")
        if api.json("GET", f"{release_path()}/latest").get("id") != release["id"]:
            raise RuntimeError("The verified release is not GitHub's latest release")
        print(f"Existing published release verified: {release['html_url']}")
        return 0
    if release is None:
        try:
            release = api.json("POST", release_path(), payload={
                "tag_name": TAG, "target_commitish": args.commit,
                "name": TITLE, "body": BODY, "draft": True, "prerelease": False,
            })
        except (ApiError, RuntimeError):
            # Recover the same draft if the create response was lost.
            release = find_release(api)
            if release is None:
                raise
            verify_target(api, release, args.commit)
            if not release.get("draft"):
                raise RuntimeError("Release was published concurrently; rerun verification")
    present = inspect_assets(api, release, expected, complete=False)
    for name in ASSET_NAMES:
        if name not in present:
            upload_asset(api, release, args.output_dir / name, expected[name])
        print(f"Verified release asset: {name}")
    inspect_assets(api, release, expected, complete=True)
    verify_target(api, release, args.commit)
    published = api.json("PATCH", f"{release_path()}/{release['id']}", payload={
        "name": TITLE, "body": BODY, "draft": False, "prerelease": False, "make_latest": "true",
    })
    if published.get("draft") or published.get("prerelease"):
        raise RuntimeError("GitHub did not confirm publication")
    verify_target(api, published, args.commit)
    inspect_assets(api, published, expected, complete=True)
    if api.json("GET", f"{release_path()}/latest").get("id") != published["id"]:
        raise RuntimeError("GitHub did not confirm the new release as latest")
    print(f"Published and verified: {published['html_url']}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        print(f"Release publication stopped: {exc}", file=sys.stderr)
        raise SystemExit(1)
