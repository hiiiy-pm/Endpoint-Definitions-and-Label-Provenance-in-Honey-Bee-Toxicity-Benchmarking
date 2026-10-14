"""Shared, standard-library validation for eight publication deliverables."""
from pathlib import Path, PurePosixPath
import hashlib
import json
import platform
import re
import sys
import zlib

ASSET_NAMES = (
    'Cover_Letter.pdf', 'Manuscript.pdf', 'Supplementary_Information.pdf',
    'LaTeX_Source.zip', 'Research_Code_and_Data.zip',
    'Research_Figures_and_Tables.zip', 'Research_Model_Outputs.zip',
    'Research_Source_Audits.zip',
)
ZIP_FIELDS = ('create_system', 'create_version', 'extract_version', 'reserved',
              'flag_bits', 'volume', 'internal_attr', 'external_attr', 'compress_type')

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def sha256(data):
    return hashlib.sha256(data).hexdigest()

def runtime():
    return {'python': platform.python_version(), 'zlib_compiled': zlib.ZLIB_VERSION,
            'zlib_runtime': zlib.ZLIB_RUNTIME_VERSION, 'platform': sys.platform}

def safe_relative(name):
    require(isinstance(name, str) and bool(name), 'Empty source/member path')
    path = PurePosixPath(name)
    require(not path.is_absolute() and '\\' not in name and ':' not in name,
            'Unsafe source/member path: ' + name)
    require(all(p not in {'', '.', '..'} for p in name.split('/')),
            'Non-canonical source/member path: ' + name)
    require(path.parts[0] not in {'.git', '.github'}, 'Internal source path: ' + name)
    return path

def source_path(root, name):
    relative = safe_relative(name)
    path = root.joinpath(*relative.parts).resolve(strict=True)
    require(path.is_relative_to(root) and path.is_file(), 'Source outside repository: ' + name)
    return path

def member_source(asset_name, member_name):
    member = safe_relative(member_name)
    prefix = 'LaTeX_Source' if asset_name == 'LaTeX_Source.zip' else 'Research'
    require(len(member.parts) > 1 and member.parts[0] == prefix,
            'Unexpected ZIP prefix: ' + member_name)
    suffix = PurePosixPath(*member.parts[1:])
    return (PurePosixPath('manuscript', 'LaTeX_Source') / suffix).as_posix() if prefix == 'LaTeX_Source' else suffix.as_posix()

def verify_blob(path, size, digest, description):
    require(isinstance(size, int) and size >= 0 and re.fullmatch('[0-9a-f]{64}', digest) is not None,
            'Invalid content identity: ' + description)
    data = path.read_bytes()
    require(len(data) == size and sha256(data) == digest, 'Content mismatch: ' + description)
    return data

def load_metadata(path):
    raw = path.read_bytes()
    metadata = json.loads(raw.decode('utf-8-sig'))
    require(metadata.get('schema_version') == 3, 'Unsupported metadata schema')
    require(metadata.get('compression_level') in range(10), 'Invalid compression level')
    assets = metadata.get('assets', [])
    require(len(assets) == 8 and {a['name'] for a in assets} == set(ASSET_NAMES),
            'Exactly the eight publication assets are required')
    sources = []
    research_sources = []
    for asset in assets:
        name = asset['name']
        require(isinstance(asset['bytes'], int) and asset['bytes'] > 0 and
                re.fullmatch('[0-9a-f]{64}', asset['sha256']) is not None, 'Invalid asset identity')
        if name.endswith('.pdf'):
            require(asset['kind'] == 'file' and asset['source'] == 'manuscript/' + name,
                    'Unexpected PDF source mapping: ' + name)
            sources.append(asset['source'])
        else:
            require(asset['kind'] == 'zip' and len(asset['entries']) > 0, 'Empty/invalid ZIP asset')
            member_names = [e['name'] for e in asset['entries']]
            require(len(member_names) == len(set(member_names)), 'Duplicate ZIP members: ' + name)
            for entry in asset['entries']:
                expected = member_source(name, entry['name'])
                require(entry['source'] == expected, 'ZIP member/source mapping differs')
                require(entry['compress_type'] in {0, 8}, 'Unsupported ZIP compression')
                sources.append(expected)
                if name.startswith('Research_'):
                    research_sources.append(expected)
    require(len(research_sources) == len(set(research_sources)), 'Research ZIPs are not disjoint')
    require(len(sources) == len(set(sources)) == metadata['source_member_count'], 'Source membership mismatch')
    return metadata, sha256(raw)
