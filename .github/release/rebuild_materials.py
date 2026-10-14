"""Rebuild and verify all eight deliverables. No downloads or scientific experiments."""
from pathlib import Path
import argparse
import json
import os
import zipfile
from materials_common import ZIP_FIELDS, require, sha256, runtime, source_path, verify_blob, load_metadata

def rebuild(repository_root, metadata_path, output_dir):
    root = repository_root.resolve(strict=True)
    output_dir = output_dir.resolve()
    require(output_dir != root and not root.is_relative_to(output_dir), 'Output directory contains source repository')
    metadata, metadata_digest = load_metadata(metadata_path)
    for asset in metadata['assets']:
        if asset['kind'] == 'file':
            verify_blob(source_path(root, asset['source']), asset['bytes'], asset['sha256'], asset['source'])
        else:
            for entry in asset['entries']:
                path = source_path(root, entry['source'])
                require(not path.is_relative_to(output_dir), 'Output overlaps a source file')
                verify_blob(path, entry['file_size'], entry['sha256'], entry['source'])
    output_dir.mkdir(parents=True, exist_ok=True)
    require(not any((output_dir / a['name']).exists() for a in metadata['assets']),
            'Output contains deliverables; choose an empty output directory')
    pending = []
    verified = []
    try:
        for asset in metadata['assets']:
            temp = output_dir / (asset['name'] + '.pending')
            require(not temp.exists(), 'Pending output already exists: ' + temp.name)
            pending.append(temp)
            if asset['kind'] == 'file':
                temp.write_bytes(source_path(root, asset['source']).read_bytes())
            else:
                with zipfile.ZipFile(temp, 'w', zipfile.ZIP_DEFLATED,
                                     compresslevel=metadata['compression_level'], allowZip64=True) as archive:
                    archive.comment = bytes.fromhex(asset['archive_comment_hex'])
                    for entry in asset['entries']:
                        info = zipfile.ZipInfo(entry['name'], tuple(entry['date_time']))
                        for field in ZIP_FIELDS:
                            setattr(info, field, entry[field])
                        info.comment = bytes.fromhex(entry['comment_hex'])
                        info.extra = bytes.fromhex(entry['extra_hex'])
                        info.file_size = entry['file_size']
                        data = verify_blob(source_path(root, entry['source']), entry['file_size'],
                                           entry['sha256'], entry['source'])
                        archive.writestr(info, data, compress_type=entry['compress_type'],
                                         compresslevel=metadata['compression_level'])
                with zipfile.ZipFile(temp) as check:
                    require(check.testzip() is None, 'Rebuilt ZIP CRC error')
                    require(check.namelist() == [e['name'] for e in asset['entries']], 'ZIP member order changed')
            verify_blob(temp, asset['bytes'], asset['sha256'], asset['name'])
            verified.append({'name': asset['name'], 'bytes': asset['bytes'], 'sha256': asset['sha256'],
                             'byte_identical_to_submission': True})
        for temp, asset in zip(pending, metadata['assets'], strict=True):
            os.replace(temp, output_dir / asset['name'])
        report = {'status': 'PASS', 'metadata_sha256': metadata_digest,
                  'source_members_verified': metadata['source_member_count'],
                  'captured_runtime': metadata['captured_runtime'], 'runtime': runtime(),
                  'assets': verified, 'experiments_run': False, 'github_operations': False}
        (output_dir / 'materials_verification.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
        print('PASS: all eight deliverables are byte-identical to the submission files.')
        return report
    finally:
        for temp in pending:
            if temp.exists():
                temp.unlink()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository-root', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, default=Path(__file__).with_name('materials_metadata.json'))
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    rebuild(args.repository_root, args.metadata, args.output_dir)

if __name__ == '__main__':
    main()
