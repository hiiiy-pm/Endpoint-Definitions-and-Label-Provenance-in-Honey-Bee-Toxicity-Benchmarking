"""Assemble all supplementary tables in numerical order (first-citation order).

The table generators write their tables into four files. This script splits the
files into single-table blocks, checks that S1-S20 are each present exactly
once, and writes them in numerical order to supplementary_tables_ordered.tex,
which the Supplementary Information source inputs.
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / 'output' / 'latex_tables'
SOURCES = ['supp_tables.tex', 'supp_tables_external.tex', 'supp_tables_source_holdout.tex',
           'supp_tables_tier_boundary.tex']
HEADER = re.compile(r'(?:\\Needspace\{[^}]*\})?(?:\\clearpage)?\\subsection\*\{Supplementary Table S(\d+)([ab]?)\.')

blocks = {}
for name in SOURCES:
    text = (TABLES / name).read_text(encoding='utf-8')
    starts = [m.start() for m in HEADER.finditer(text)]
    if not starts or text[:starts[0]].strip():
        raise SystemExit(f'{name}: content found before the first table header')
    for i, start in enumerate(starts):
        m = HEADER.match(text, start)
        key = (int(m.group(1)), m.group(2))
        if key in blocks:
            raise SystemExit(f'Duplicate supplementary table S{key[0]}{key[1]}')
        end = starts[i + 1] if i + 1 < len(starts) else len(text)
        block = text[start:end].strip()
        # Every table starts on a fresh block with room for its header and first rows.
        block = re.sub(r'^\\clearpage', '', block)
        if not block.startswith(r'\Needspace'):
            block = r'\Needspace{12\baselineskip}' + block
        blocks[key] = block

numbers = sorted({k[0] for k in blocks})
if numbers != list(range(1, 21)):
    raise SystemExit(f'Supplementary tables must be S1-S20 without gaps; found {numbers}')
ordered = [blocks[k] for k in sorted(blocks)]
(TABLES / 'supplementary_tables_ordered.tex').write_text('\n\n'.join(ordered) + '\n', encoding='utf-8')
print('Assembled', len(ordered), 'supplementary table blocks:',
      ', '.join(f'S{n}{s}' for n, s in sorted(blocks)))
