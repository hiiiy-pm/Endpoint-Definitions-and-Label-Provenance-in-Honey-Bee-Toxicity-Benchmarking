"""Capture exact current deliverables after checking every byte against repository sources."""
from pathlib import Path
import argparse
import json
import os
import zipfile
from materials_common import ASSET_NAMES, ZIP_FIELDS, require, sha256, runtime, source_path, member_source, load_metadata

def capture(submission_dir, repository_root, output, compression_level=9):
    submission_dir = submission_dir.resolve(strict=True)
    repository_root = repository_root.resolve(strict=True)
    present = {p.name for p in submission_dir.iterdir()}
    require(present == set(ASSET_NAMES), 'Submission directory must contain exactly the eight deliverables')
    assets = []
    sources = []
    for name in ASSET_NAMES:
        path = submission_dir / name
        blob = path.read_bytes()
        asset = {'name': name, 'bytes': len(blob), 'sha256': sha256(blob)}
        if name.endswith('.pdf'):
            source = 'manuscript/' + name
            require(blob.startswith(b'%PDF-'), 'Invalid PDF header: ' + name)
            require(source_path(repository_root, source).read_bytes() == blob, 'PDF differs from repository: ' + name)
            asset.update(kind='file', source=source)
            sources.append(source)
        else:
            entries = []
            with zipfile.ZipFile(path) as archive:
                require(archive.testzip() is None, 'ZIP CRC error: ' + name)
                for info in archive.infolist():
                    require(not info.is_dir() and not (info.flag_bits & 1), 'Directory/encrypted member is unsupported')
                    source = member_source(name, info.filename)
                    data = archive.read(info)
                    require(source_path(repository_root, source).read_bytes() == data,
                            'ZIP member differs from repository: ' + info.filename)
                    entry = {field: getattr(info, field) for field in ZIP_FIELDS}
                    entry.update(name=info.filename, source=source, date_time=list(info.date_time),
                                 comment_hex=info.comment.hex(), extra_hex=info.extra.hex(),
                                 file_size=info.file_size, sha256=sha256(data))
                    entries.append(entry)
                    sources.append(source)
                asset.update(kind='zip', archive_comment_hex=archive.comment.hex(), entries=entries)
        assets.append(asset)
    metadata = {'schema_version': 3, 'compression_level': compression_level, 'captured_runtime': runtime(),
                'source_member_count': len(sources), 'assets': assets}
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    pending = output.with_name(output.name + '.pending')
    require(not pending.exists(), 'Pending metadata already exists')
    try:
        pending.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + '\n', encoding='utf-8', newline='\n')
        load_metadata(pending)
        os.replace(pending, output)
    finally:
        if pending.exists():
            pending.unlink()
    print('Captured 8 assets and', len(sources), 'verified source members:', output)
    return metadata

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--submission-dir', type=Path, required=True)
    parser.add_argument('--repository-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('materials_metadata.json'))
    parser.add_argument('--compression-level', type=int, default=9, choices=range(10))
    args = parser.parse_args()
    capture(args.submission_dir, args.repository_root, args.output, args.compression_level)

if __name__ == '__main__':
    main()
