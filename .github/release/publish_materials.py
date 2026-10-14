"""Publish only the eight verified study assets, using GitHub Actions authentication."""
from pathlib import Path
import json
import os
import re
import urllib.error
import urllib.parse
import urllib.request
from materials_common import ASSET_NAMES, require, verify_blob, load_metadata

REPOSITORY = 'hiiiy-pm/Endpoint-Definitions-and-Label-Provenance-in-Honey-Bee-Toxicity-Benchmarking'
TAG = 'study-materials'
TITLE = 'Study materials'
BODY = ('Materials for "Endpoint Definitions and Label Provenance in Honey-Bee Toxicity Benchmarking".\n\n'
        'The manuscript, Supplementary Information, cover letter and LaTeX sources are accompanied by '
        'the research code, data, model outputs, source records, figures and tables. '
        'The four Research archives extract into one Research directory. '
        'Its README describes the analysis and reproduction commands.')
API = 'https://api.github.com/repos/' + REPOSITORY

def api(method, url, payload=None, content_type='application/json'):
    parsed = urllib.parse.urlparse(url)
    require(parsed.scheme == 'https' and parsed.hostname in {'api.github.com', 'uploads.github.com'},
            'Unexpected GitHub API destination')
    data = payload if isinstance(payload, bytes) else None if payload is None else json.dumps(payload).encode('utf-8')
    request = urllib.request.Request(url, data=data, method=method, headers={
        'Authorization': 'Bearer ' + os.environ['GITHUB_TOKEN'],
        'Accept': 'application/vnd.github+json', 'Content-Type': content_type,
        'X-GitHub-Api-Version': '2026-03-10', 'User-Agent': 'study-materials-publication'})
    with urllib.request.urlopen(request, timeout=240) as response:
        return json.load(response)

def optional_get(url):
    try:
        return api('GET', url)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return None
        raise

def find_release():
    # The tag endpoint documents published releases only. Listing also exposes
    # drafts to contents-write authentication, so interrupted uploads can resume.
    found = []
    page = 1
    while True:
        rows = api('GET', API + '/releases?per_page=100&page=' + str(page))
        require(isinstance(rows, list), 'Unexpected release-list response')
        found.extend(row for row in rows if row['tag_name'] == TAG)
        if len(rows) < 100:
            break
        page += 1
    require(len(found) <= 1, 'Multiple releases use the study-materials tag')
    return found[0] if found else None

def verified_assets(metadata_path, assets_dir):
    metadata, digest = load_metadata(metadata_path)
    report = json.loads((assets_dir / 'materials_verification.json').read_text(encoding='utf-8'))
    require(report['status'] == 'PASS' and report['metadata_sha256'] == digest,
            'The asset verification report is missing or belongs to other metadata')
    require(report['source_members_verified'] == metadata['source_member_count'], 'Source count differs')
    rows = report['assets']
    require(len(rows) == 8 and {r['name'] for r in rows} == set(ASSET_NAMES), 'Incomplete verified asset set')
    report_by_name = {r['name']: r for r in rows}
    for asset in metadata['assets']:
        row = report_by_name[asset['name']]
        require(row['bytes'] == asset['bytes'] and row['sha256'] == asset['sha256'] and
                row['byte_identical_to_submission'] is True, 'Report differs from asset metadata')
        verify_blob(assets_dir / asset['name'], asset['bytes'], asset['sha256'], asset['name'])
    return {a['name']: a for a in metadata['assets']}

def check_remote_asset(actual, expected):
    require(actual['state'] == 'uploaded' and actual['size'] == expected['bytes'] and
            actual.get('digest') == 'sha256:' + expected['sha256'],
            'Remote asset identity differs: ' + expected['name'])

def main():
    require(os.environ.get('GITHUB_REPOSITORY') == REPOSITORY, 'Unexpected repository')
    require(os.environ.get('GITHUB_REF') == 'refs/heads/main', 'Only main may publish study materials')
    commit = os.environ['GITHUB_SHA']
    require(re.fullmatch('[0-9a-f]{40}', commit) is not None, 'Invalid source commit')
    root = Path(__file__).resolve().parents[2]
    assets_dir = root / 'release_assets'
    expected = verified_assets(Path(__file__).with_name('materials_metadata.json'), assets_dir)
    tag = optional_get(API + '/git/ref/tags/' + TAG)
    if tag is not None:
        require(tag['object']['type'] == 'commit' and tag['object']['sha'] == commit,
                'Existing tag identifies a different commit; it will not be moved')
    release = find_release()
    if release is None:
        release = api('POST', API + '/releases', {'tag_name': TAG, 'target_commitish': commit,
                      'name': TITLE, 'body': BODY, 'draft': True, 'prerelease': False})
    else:
        require(release['target_commitish'] == commit, 'Existing release identifies a different commit')
        require(release['name'] == TITLE and release['body'] == BODY, 'Existing release description differs')
    assets_url = API + '/releases/' + str(release['id']) + '/assets?per_page=100'
    assets = {a['name']: a for a in api('GET', assets_url)}
    require(set(assets).issubset(expected), 'Existing release contains unexpected assets')
    for name in ASSET_NAMES:
        if name in assets:
            actual = assets[name]
        else:
            require(release['draft'], 'Published releases cannot be supplemented or replaced')
            upload_url = release['upload_url'].split('{', 1)[0] + '?' + urllib.parse.urlencode({'name': name})
            content_type = 'application/pdf' if name.endswith('.pdf') else 'application/zip'
            data = verify_blob(assets_dir / name, expected[name]['bytes'], expected[name]['sha256'], name)
            actual = api('POST', upload_url, data, content_type)
        check_remote_asset(actual, expected[name])
        print('Verified asset:', name, expected[name]['sha256'], flush=True)
    final_assets = api('GET', assets_url)
    require(len(final_assets) == 8 and {a['name'] for a in final_assets} == set(expected),
            'Remote release must contain exactly the eight deliverables')
    for asset in final_assets:
        check_remote_asset(asset, expected[asset['name']])
    if release['draft']:
        release = api('PATCH', API + '/releases/' + str(release['id']),
                      {'draft': False, 'prerelease': False, 'make_latest': 'true'})
    require(not release['draft'] and not release['prerelease'], 'Release is not public')
    latest = api('GET', API + '/releases/latest')
    require(latest['id'] == release['id'] and latest['tag_name'] == TAG, 'Study materials are not the latest release')
    final_tag = api('GET', API + '/git/ref/tags/' + TAG)
    require(final_tag['object']['type'] == 'commit' and final_tag['object']['sha'] == commit,
            'Published tag does not identify the verified repository commit')
    print('Published verified study materials:', release['html_url'], flush=True)
    if os.environ.get('GITHUB_STEP_SUMMARY'):
        with open(os.environ['GITHUB_STEP_SUMMARY'], 'a', encoding='utf-8') as summary:
            summary.write('Published [Study materials](' + release['html_url'] + ').\n\n'
                          'All eight asset sizes and SHA-256 digests match the submission deliverables.\n')

if __name__ == '__main__':
    main()
