"""Check dependencies and build the manuscript, SI and cover-letter PDFs.

Figures and supplementary tables are first copied from the analysis repository
(default ../Research, or the path given as the first argument) so that the PDFs use
exactly the files generated there.
"""
from pathlib import Path
import argparse
import os
import re
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('repo', nargs='?', type=Path, default=HERE.parent/'Research')
parser.add_argument('--source-only', action='store_true', help='Compile the supplied sources without copying repository assets.')
parser.add_argument('--regenerate', action='store_true', help='Regenerate tables and figures from the saved analysis outputs before compiling.')
parser.add_argument('--check-only', action='store_true', help='Check TeX tools, packages and supplied assets without changing any files.')
args = parser.parse_args()
REPO = args.repo.resolve()
if args.regenerate and (args.source_only or args.check_only):
    parser.error('--regenerate cannot be combined with --source-only or --check-only')

pdflatex = os.environ.get('PDFLATEX') or shutil.which('pdflatex')
bibtex = os.environ.get('BIBTEX') or shutil.which('bibtex')
if not pdflatex or not bibtex:
    raise SystemExit('pdflatex and bibtex are required (PATH or PDFLATEX/BIBTEX).')
tex_bin = Path(pdflatex).resolve().parent
kpse_candidate = tex_bin / ('kpsewhich.exe' if os.name == 'nt' else 'kpsewhich')
kpsewhich = str(kpse_candidate) if kpse_candidate.is_file() else shutil.which('kpsewhich')
if not kpsewhich:
    raise SystemExit('kpsewhich is required for the dependency check; use the same TeX installation as pdflatex.')

def missing_tex(names):
    names = [name for name in names if not (HERE/name).is_file()]
    if not names:
        return []
    try:
        result = subprocess.run([kpsewhich, *names], cwd=HERE, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, encoding='utf-8',
                                errors='replace', timeout=30)
    except subprocess.TimeoutExpired:
        raise SystemExit('TeX package lookup timed out; check the local TeX installation and rerun --check-only.')
    found = {Path(line.strip()).name for line in result.stdout.splitlines() if line.strip()}
    return [name for name in names if name not in found]

# Direct dependencies of all three documents, including the journal class.
packages = set()
for name in ['huang_access.tex', 'supplementary.tex', 'cover_letter.tex',
             'supplementary_style.sty', 'ieeeaccess.cls', 'spotcolor.sty']:
    source = (HERE/name).read_text(encoding='utf-8')
    source = re.sub(r'(?m)(?<!\\)%.*$', '', source)
    if name == 'ieeeaccess.cls':
        # These packages belong to the class's editorial-only branch; this
        # submission uses the ordinary Research Article class configuration.
        source = source.replace(r'\RequirePackage{ltxutil,ltxgrid}', '')
    for group in re.findall(r'\\(?:usepackage|RequirePackage)(?:\[[^\]]*\])?\{([^}]+)\}', source):
        packages.update(x.strip() for x in group.split(','))
# The SI style selects newtx when present and otherwise uses Latin Modern.
packages.difference_update({'newtxtext', 'newtxmath', 'lmodern'})
packages.update({'newtxtext', 'newtxmath'} if not missing_tex(['newtxtext.sty']) else {'lmodern'})
missing = missing_tex([p+'.sty' for p in sorted(packages)]+['IEEEtran.cls', 'article.cls'])
if missing:
    raise SystemExit('Missing TeX packages: '+', '.join(missing)+
                     '. Install these with your TeX distribution package manager, then rerun --check-only.')
print(f'TeX dependency check passed ({len(packages)} direct packages).', flush=True)

def check_supplied_assets():
    assets = ['ieeeaccess.cls', 'IEEEtran.bst', 'references.bib',
              'methods_formalization.tex', 'supplementary_tables_ordered.tex',
              'supplementary_audit_tables.tex', 'logo.png', 'notaglinelogo.png', 'bullet.png']
    assets += [stem+'.pdf' for stem in MAIN_FIGS]
    assets += ['figures/'+stem+'.pdf' for stem in SUPP_FIGS]
    assets += ['t1-formata.map', 't1-times.map', 't1-giovannistd.map',
               't1-formata-regular.tfm', 't1-formata-regular.pfb',
               't1-times.tfm', 't1-times.pfb']
    for source in ['huang_access.tex', 'supplementary.tex']:
        assets += re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{(author_photos/[^}]+)\}',
                             (HERE/source).read_text(encoding='utf-8'))
    absent = [name for name in assets if not (HERE/name).is_file()]
    if absent:
        raise SystemExit('Missing supplied assets: '+', '.join(absent))
    print('Source and figure asset check passed.', flush=True)

if args.regenerate:
    scripts = ['code/09_export_submission_tables.py', 'scripts/make_supp_tables.py',
               'scripts/make_external_tables.py', 'scripts/make_source_holdout_table.py',
               'scripts/make_tier_boundary_tables.py', 'scripts/make_audit_tables.py',
               'scripts/assemble_supplementary_tables.py', 'scripts/make_main_tables.py',
               'scripts/generate_manuscript_figures.py', 'scripts/make_label_source_figure.py',
               'scripts/make_crossover_figure.py', 'scripts/assemble_audit_supplement.py']
    for script in scripts:
        command = [sys.executable, str(REPO/script)]
        if script == 'scripts/generate_manuscript_figures.py':
            command.append('--publication-layout')
        subprocess.run(command, cwd=REPO, check=True)
    subprocess.run([sys.executable, str(REPO/'scripts/generate_manuscript_figures.py'),
                    '--only', 'FigS1_full_performance', 'FigS2_model_decomposition',
                    'FigS3_structural_diagnostics'], cwd=REPO, check=True)
MAIN_FIGS = ['Fig1_study_design', 'Fig2_endpoint_discriminability', 'Fig3_composition_and_models',
             'Fig4_class_count_and_tier_boundary', 'Fig5_label_provenance', 'Fig6_chemical_structures',
             'Fig7_label_crossover']
SUPP_FIGS = ['FigS1_full_performance', 'FigS2_model_decomposition', 'FigS3_structural_diagnostics',
             'FigS4_label_source_contrasts']
if args.check_only:
    check_supplied_assets()
    raise SystemExit(0)
if not args.source_only:
    for stem in MAIN_FIGS:
        shutil.copy2(REPO / 'figures' / f'{stem}.pdf', HERE / f'{stem}.pdf')
    (HERE / 'figures').mkdir(exist_ok=True)
    for stem in SUPP_FIGS:
        shutil.copy2(REPO / 'figures' / f'{stem}.pdf', HERE / 'figures' / f'{stem}.pdf')
    shutil.copy2(REPO / 'output' / 'latex_tables' / 'supplementary_tables_ordered.tex',
                 HERE / 'supplementary_tables_ordered.tex')
    shutil.copy2(REPO / 'output' / 'latex_tables' / 'supplementary_audit_tables.tex',
                 HERE / 'supplementary_audit_tables.tex')

check_supplied_assets()

LOGS = HERE / 'build_logs'
LOGS.mkdir(exist_ok=True)
for stem, use_bib in [('huang_access', True), ('supplementary', False), ('cover_letter', False)]:
    tex = [pdflatex, '-interaction=nonstopmode', '-halt-on-error', stem + '.tex']
    cmds = [tex] + ([[bibtex, stem]] if use_bib else []) + [tex, tex]
    for i, cmd in enumerate(cmds, 1):
        p = subprocess.run(cmd, cwd=HERE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log = p.stdout.decode('utf-8', errors='replace')
        (LOGS / f'{stem}_{i}.txt').write_text(log, encoding='utf-8')
        if p.returncode:
            print(log[-1500:])
            raise SystemExit(f'{stem}: step {i} failed')
    print(stem, 'built')
