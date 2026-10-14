"""Generate supplementary sensitivity tables from scientific outputs without fitting models."""
from pathlib import Path
import csv
import json

ROOT = Path(__file__).resolve().parents[1]
REV = ROOT / 'supplementary_audits'
OUT = ROOT / 'output/latex_tables/supplementary_audit_tables.tex'


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def table(headers, rows, widths=None):
    spec = widths or ('l' + 'r' * (len(headers) - 1))
    lines = [r'\begin{center}', r'\small', r'\begin{tabular}{' + spec + '}', r'\toprule',
             ' & '.join(headers) + r' \\', r'\midrule']
    lines.extend(' & '.join(str(value) for value in row) + r' \\' for row in rows)
    return '\n'.join(lines + [r'\bottomrule', r'\end{tabular}', r'\end{center}'])


def num(value):
    return f'{float(value):.3f}'


def main():
    stat = REV / 'statistics'
    ident = REV / 'identity_sources'
    pieces = [r'\section*{Supplementary sensitivity analyses}',
        r'''These analyses assess label-rule sensitivity, source provenance and identifier consistency. They retain the original
SMILES, labels and official split files. Current identifier responses are a dated
consistency screen; the conditional qualifier replay and additional model
sensitivities are reported separately from the historical reconstruction.''']

    pieces += [r'\subsection*{Supplementary Table S25a. Exclusion-log scope and actual filtering}',
        r'''The upstream log selects all rows for a flagged CAS identity, whereas
the filter removes individual CAS-by-route groups. The additional log entries
do not imply that the executed filter removed those observations.''']
    ex = json.loads((ident / 'exclusion_log_summary.json').read_text())
    pieces.append(table(['Quantity', 'Count'], [
        ['Log rows', ex['log_rows']], ['CAS identities in log', ex['log_CAS']],
        ['Removed CAS-by-route groups', ex['actual_rejected_groups']],
        ['Rows in removed groups', ex['actual_rejected_records']],
        ['Eligible groups also included in log', ex['retained_groups_also_logged']],
        ['Eligible-group rows also included in log', ex['retained_group_rows_also_logged']],
        ['Logged rows in final ECOTOX determining groups', 119]]))
    assert ex['log_rows'] == ex['actual_rejected_records'] + ex['retained_group_rows_also_logged']

    pieces += [r'\subsection*{Supplementary Table S25b. Uniform rule D and unresolved-label handling}',
        r'''All-D and the primary D/C/D rules were applied to all accepted records
and the same upstream-eligible routes at 100/11/1~\ugbee. At 11~\ugbee, CAS
115086-54-9 and 736994-63-1 were unresolved under C and positive under D. Both
had positive distributed labels. With the keep policy, the three 1,035-entry
final label vectors were identical, containing 789/296/170 positives. This is an
input-equivalence check under the same fitting protocol, not a second independent
performance experiment. Under exclude, all-D retained 915 identities instead of
913; the two additions were positive at every cut-off.''']
    members = read_csv(stat / 'all_d_exclusion_membership.csv')
    member_rows = []
    for split in ['Random', 'MaxMin', 'Time']:
        tr = next(row for row in members if row['Split'] == split.lower() and row['Partition'] == 'train')
        te = next(row for row in members if row['Split'] == split.lower() and row['Partition'] == 'test')
        member_rows.append([split, f"{tr['MixedExcludeN']}/{te['MixedExcludeN']}",
                            f"{tr['UniformExcludeN']}/{te['UniformExcludeN']}"])
    vectors = read_csv(stat / 'all_d_vector_comparison.csv')
    assert all(int(row['Length']) == 1035 and int(row['KeepDifferences']) == 0
               and row['MixedKeepUint8SHA256'] == row['UniformKeepUint8SHA256'] for row in vectors)
    pieces.append(table(['Split', 'D/C/D train/test', 'All-D train/test'], member_rows, 'lcc'))

    pieces += [r'\subsection*{Supplementary Table S25c. ECFP crossover after all-D exclusion}',
        r'''Original and provisional labels use identical retained members within
each all-D crossover. Both training-label versions were fitted for Random and
MaxMin; the Time training inputs were unchanged and their fitted scores were
reused. The original mixed-rule exclusion fits were independently reproduced.
Intervals use 2,000 paired scaffold-cluster draws, condition on the fitted models
and retained cohorts, and are not adjusted for the added scenario. This selected
cohort differs from the primary keep cohort. O/O denotes original training and
evaluation labels; R/R denotes provisionally revised labels for both. O/R label-source
notation is distinct from aggregation rules A--D and outcome codes P/N/U.''']
    gaps = read_csv(stat / 'uniform_d_exclude_gap.csv')
    rows = []
    for split in ['Random', 'MaxMin', 'Time']:
        for severe in [11, 1]:
            group = [r for r in gaps if r['Scenario'] == 'uniform_DDD'
                     and r['Split'] == split and int(r['Severe']) == severe]
            def get(train):
                return next(r for r in group if r['TrainLabels'] == train and r['EvalLabels'] == train)
            old, new, change = get('original'), get('provisional'), get('change:provisional-vs-original')
            assert abs(float(new['Gap']) - float(old['Gap']) - float(change['Gap'])) < 1e-12
            rows.append([split, f'{severe}--100', num(old['Gap']), num(new['Gap']),
                         num(change['Gap']), f"[{num(change['CI_low'])}, {num(change['CI_high'])}]"])
    pieces.append(table(['Split', 'Contrast', 'O/O gap', 'R/R gap', 'Change', r'95\% CI'], rows, 'llrrrr'))

    pieces += [r'\Needspace{18\baselineskip}', r'\subsection*{Supplementary Table S25d. Benomyl qualifier-direction replay}',
        r'''The archived oral record for CAS 17804-35-2 reports $<10$~\ugbee;
the cited Rothamsted report's Table 2 (printed p.~56) reports $>10$~\ugbee.
The table below changes only that operator while retaining the numerical bound,
historical conditions and route eligibility. It does not adjudicate the full
experimental identity. The primary keep inputs are unchanged, and the identity
is excluded before and after the replay under the exclude policy. The compound
is Tier2, so the 212-compound Tier1 audit counts are unaffected. Full scenario
input comparisons are in \path{supplementary_audits/statistics/benomyl_scenario_input_impact.csv}.
The original report and database also differ in material and stage metadata;
these fields are retained separately in the source ledger.''']
    pieces.append(table(['Cut-off', 'Reconstructed, $<10$', 'Reconstructed, $>10$', 'Keep, both'], [
        ['100', 'Positive', 'Unresolved', 'Positive'], ['11', 'Positive', 'Unresolved', 'Positive'],
        ['1', 'Unresolved', 'Negative', 'Negative']], 'lccc'))
    pieces.append(r'Original source: \url{https://www.era.rothamsted.ac.uk/eradoc/book/16}. The cited PDF SHA256 and record-to-table locations are supplied with the source audit.')

    pieces += [r'\Needspace{30\baselineskip}', r'\subsection*{Supplementary Table S26a. Current CID--structure consistency screen}',
        r'''PubChem property records were retrieved for all 1,030 distinct stored
CIDs on 3 October 2026. RDKit 2025.09.4 converted the frozen SMILES to InChIKeys.
Complete-key and connectivity-block matches were assessed before comparison of
standardized fragment parents. Differences may reflect salts, stereochemistry,
mixtures or representation changes, as well as incorrect identifiers. They do
not by themselves establish toxicity-label errors. Neither the frozen structures
nor their labels were replaced.''']
    summary = json.loads((ident / 'identity_summary.json').read_text())
    counts = summary['status_counts']
    pieces.append(table(['Comparison to current stored-CID structure', 'Rows'], [
        ['Complete InChIKey match', counts['full_InChIKey_match']],
        ['Connectivity match; complete key differs', counts['connectivity_match_full_key_difference']],
        ['Fragment-parent full-key match only', counts['fragment_parent_match_full_structure_difference']],
        ['Fragment-parent connectivity match only', counts['fragment_parent_connectivity_match']],
        ['No match under these checks; review required', counts['structure_mismatch_requires_adjudication']],
        ['Total', summary['benchmark_rows']]]))
    assert sum(counts.values()) == 1035
    pieces.append(r'''The fixed screening set contains the 33 identities without a
direct connectivity-block match (six fragment-parent cases and 27 remaining
differences). Its membership file and SHA256 were recorded before evaluating
the exclusion sensitivity. The 127 direct connectivity matches were retained;
the ECFP models do not encode chirality. Complete responses and the full 1,035-row
ledger accompany the archive.''')

    pieces += [r'\subsection*{Supplementary Table S26b. Review of five non-unique CID groups}',
        r'''Five stored CIDs are shared by ten benchmark rows. The table separates
current identifier conflicts from alternative entity representations. CAS lookup
and structure agreement corroborate the listed candidate CIDs; they are not
substituted into the frozen benchmark.''']
    pieces.append(table(['Stored CID', 'Case', 'Source-comparison finding'], [
        ['13836', 'Chlormequat / chloride', 'Cation matches 13837; salt matches 13836'],
        ['5284484', 'Metiram / zineb', 'Distinct materials; zineb matches 5284484'],
        ['16116', 'Silthiofam / substituted urea', 'Silthiofam matches 9881821'],
        ['93341', 'Diclomezine / diflumetorim', 'Diflumetorim matches 15390131'],
        ['5483872', 'Metominostrobin / phthalide', 'Phthalide matches 33791']],
        r'lp{0.29\linewidth}p{0.43\linewidth}'))
    pieces.append(r'''All sources obtain their year proxy through CID-based PubChem
queries. ECOTOX use flags also originate from PubChem; PPDB/BPDB use flags were
parsed from provider-page text. The historical per-row literature-query responses
were not retained, and the archived literature endpoint returned HTTP 404 in the
current retrieval. A current creation date cannot replace an unknown historical
literature date. This analysis therefore establishes dependency paths and current
identifier conflicts, but does not adjudicate replacement historical years or
changed Time membership. In the ten-row source comparison, the ECOTOX identifiers match
their frozen structures; the confirmed distinct-structure CID conflicts are in
PPDB rows whose use flags came from provider pages. The ledger records these
source-specific distinctions instead of inferring use errors from every CID flag.''')

    screen = stat / 'si_identifier_screen_insert.tex'
    if screen.exists():
        pieces += [r'\subsection*{Supplementary Table S26c. ECFP sensitivity to the identifier screening set}',
                   screen.read_text(encoding='utf-8')]
    else:
        raise FileNotFoundError('Identifier sensitivity insert is required before final supplement generation')
    # Keep the inferential boundary with its existing output pointer, avoiding
    # a redundant closing paragraph that can strand a few lines on a new page.
    pieces[-1] = pieces[-1].rstrip() + r''' These retrospective sensitivity intervals
exclude uncertainty from historical identifier assignment or study-level adjudication.'''
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text('\n\n'.join(pieces) + '\n', encoding='utf-8')
    print(OUT)


if __name__ == '__main__':
    main()
