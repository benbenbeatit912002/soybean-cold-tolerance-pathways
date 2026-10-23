"""Independently verify fresh scores, full simulations, and manuscript comparisons.

Archived scores and manuscript tables are read ONLY by this verifier. They are
never inputs to recalculation or simulation. Run after both full analyses and
both historical diagnostic runs finish, and after extracting Table S1/S2.
"""
from collections import defaultdict, deque
import csv
from decimal import Decimal
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.stats import false_discovery_control, hypergeom
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"python"))
from ctpathway import read_tsv, write_tsv, canonical_id, read_gmt

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/"results/recalculated_2026-10-21"
SOURCE_ROOT = ROOT.parent


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def diagnostic_limit(a, b):
    """Six-SE screen for two independent 10,000-draw estimates, not equivalence."""
    pooled = (np.asarray(a)+np.asarray(b))/2
    return 6*np.sqrt(pooled*(1-pooled)*2/10000) + 2/10001


def check_scores():
    source = SOURCE_ROOT/"研究法/大豆基因/數據/外部基因校正"
    raw = read_tsv(source/"gene_p_len.txt")
    saved = read_tsv(source/"Yamasaki_Go校正後1213_G_log10.txt")
    groups = defaultdict(deque)
    for row in saved:
        groups[(canonical_id(row["gene_2.0"]), int(row["len"]))].append(row)
    expected = [groups[(canonical_id(r["gene_2.0"]), int(r["len"]))].popleft() for r in raw]
    require(not any(groups.values()), "Unmatched saved records")
    report = {}
    for language in ("python", "r"):
        actual = read_tsv(BASE/f"{language}_scores/recalculated_scores.tsv")
        require(len(actual) == len(raw), "Score record count differs")
        for a, old in zip(actual, raw):
            require(a["gene_id"] == canonical_id(old["gene_2.0"]), "Input order differs")
            for field, original in (("raw_short","p_1"),("raw_middle","p_24"),("gene_length","len"),("historical_snp","snp")):
                require(float(a[field]) == float(old[original]), f"Rebuilt {field} differs")
        report[language] = {}
        for phase, column in (("short","p_1"),("middle","p_24")):
            matches = sum(Decimal(a[f"score_{phase}"]) == Decimal(old[column]) for a,old in zip(actual,expected))
            require(matches == len(actual), f"{language}/{phase}: saved score mismatch")
            require(all(np.isfinite(float(a[f"score_{phase}"])) for a in actual), "Nonfinite score")
            report[language][phase] = dict(records=len(actual), exact_decimal_matches=matches)
            strict = read_tsv(BASE/f"{language}_strict/recalculated_scores.tsv")
            report[language][phase]["strict_exact_matches"] = sum(Decimal(a[f"score_{phase}"]) == Decimal(old[column]) for a,old in zip(strict,expected))
            report[language][phase]["strict_nonfinite"] = sum(not np.isfinite(float(a[f"score_{phase}"])) for a in strict)
        report[language]["exceptions"] = read_tsv(BASE/f"{language}_scores/exception_audit.tsv")
    report["historical_repair_provenance_confirmed"] = False
    report["reference_used_only_in_verifier"] = True
    return report


def verify_runs():
    tables, report = {}, {"runs":[], "language_comparisons":[]}
    pathways=read_gmt(ROOT/'data/input/pathways.gmt')
    # Independent input joins and duplicate means, not only output-to-output checks.
    for language in ('python','r'):
        raw=read_tsv(BASE/f'{language}_scores/recalculated_scores.tsv')
        for phase in ('short','middle'):
            groups=defaultdict(list)
            for row in raw: groups[row['gene_id']].append(float(row[f'score_{phase}']))
            actual=read_tsv(BASE/f'{language}_full/{phase}_analysis_scores.tsv')
            require(set(groups)=={r['gene_id'] for r in actual},'Analysis gene IDs differ from recalculated inputs')
            np.testing.assert_allclose([float(r['analysis_score']) for r in actual],
                [np.mean(groups[r['gene_id']]) for r in actual],rtol=1e-12,atol=1e-14)
    for profile in ("full", "historical"):
        for language in ("python", "r"):
            folder = BASE/f"{language}_{profile}"
            summary = read_tsv(folder/"input_summary.tsv")  # Completion sentinel.
            report["runs"].append(dict(language=language, profile=profile, input_summary=summary, tests=[]))
            for qc in summary:
                phase = qc["phase"]
                raw=read_tsv(BASE/f'{language}_scores/recalculated_scores.tsv')
                universe={r['gene_id'] for r in raw}
                ct={r['gene_id'] for r in read_tsv(ROOT/f'data/input/ct_{phase}.tsv')}
                expected_membership={}
                for pathway in pathways:
                    members=set(pathway['genes']) & universe
                    selected=(6<=len(members)<=1500) if profile=='full' else (5<=len(pathway['genes'])<1500 and bool(set(pathway['genes'])&ct))
                    if selected:
                        expected_membership[pathway['go_id']]=(len(pathway['genes']),len(members),len(members&ct))
                for method in (("hypergeometric","sumstat","gsea","maxmean") if profile=="full" else ("sumstat",)):
                    rows = read_tsv(folder/f"{phase}_{method}.tsv")
                    tables[profile,language,phase,method] = rows
                    expected = int(qc["tested_pathways"] if profile=="full" else qc["candidates"])
                    require(len(rows)==expected and len({r["go_id"] for r in rows})==expected, "Unexpected GO family")
                    require({r['go_id'] for r in rows}==set(expected_membership),'GO family differs from independently joined inputs')
                    if profile=='full':
                        for row in rows:
                            require(tuple(int(row[k]) for k in ('annotation_size','tested_size','ct_count'))==expected_membership[row['go_id']], 'Incorrect pathway membership counts')
                    pcol,qcol = ("hypergeometric_p","hypergeometric_bh") if method=="hypergeometric" else ("p_value","p_bh")
                    p=np.asarray([float(r[pcol]) for r in rows]); q=np.asarray([float(r[qcol]) for r in rows])
                    require(np.isfinite(p).all() and np.isfinite(q).all(), "Nonfinite probabilities")
                    require(np.all((0<=p)&(p<=q+1e-14)&(q<=1)), "Invalid probabilities")
                    np.testing.assert_allclose(q,false_discovery_control(p),rtol=1e-12,atol=1e-14)
                    if method != "hypergeometric":
                        require(all(int(r["permutations"])==10000 and int(r["seed"])==20261021+(phase=="middle") for r in rows), "Incomplete or wrong run")
                        counts=np.asarray([int(r["null_exceedances"]) for r in rows])
                        require(np.all((0<=counts)&(counts<=10000)), "Invalid exceedances")
                        np.testing.assert_allclose(p, (counts+1)/10001 if profile=="full" else counts/10000,rtol=1e-12,atol=1e-15)
                        require(np.isfinite([float(r["statistic"]) for r in rows]).all(), "Invalid statistics")
                    else:
                        exact = hypergeom.sf(np.asarray([int(r['ct_count']) for r in rows])-1,
                            int(qc['universe_size']),np.asarray([int(r['tested_size']) for r in rows]),int(qc['ct_in_universe']))
                        np.testing.assert_allclose(p,exact,rtol=1e-10,atol=0)
                    report["runs"][-1]["tests"].append(dict(phase=phase,method=method,pathways=len(rows),
                        distinct_sizes=len({r["tested_size"] for r in rows}),min_p=float(p.min()),min_bh=float(q.min()),
                        bh_below_0_05=int(np.sum(q<.05)),bh_below_0_005=int(np.sum(q<.005))))
        require(read_tsv(BASE/f"python_{profile}/input_summary.tsv")==read_tsv(BASE/f"r_{profile}/input_summary.tsv"), "Language input counts differ")
        for phase in ("short","middle"):
            for method in (("hypergeometric","sumstat","gsea","maxmean") if profile=="full" else ("sumstat",)):
                left,right=(tables[profile,language,phase,method] for language in ("python","r"))
                for a,b in zip(left,right):
                    require(all(a[k]==b[k] for k in ("go_id","annotation_size","tested_size","ct_count")), "Language memberships differ")
                result=dict(profile=profile,phase=phase,method=method)
                if method=="hypergeometric":
                    for k in ("hypergeometric_p","hypergeometric_bh"):
                        np.testing.assert_allclose([float(r[k]) for r in left],[float(r[k]) for r in right],rtol=1e-10,atol=0)
                    result["exact_probability_agreement"]=True
                else:
                    a=np.asarray([float(r['statistic']) for r in left]); b=np.asarray([float(r['statistic']) for r in right])
                    np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-12)
                    result["max_statistic_difference"]=float(np.max(abs(a-b)))
                    a=np.asarray([float(r['p_value']) for r in left]); b=np.asarray([float(r['p_value']) for r in right])
                    result["max_p_difference"]=float(np.max(abs(a-b)))
                    result["outside_six_SE_screen"]=int(np.sum(abs(a-b)>diagnostic_limit(a,b)))
                report["language_comparisons"].append(result)
    # Verify the intentionally labeled legacy hypergeometric profiles too.
    for phase in ("short","middle"):
        left,right=(read_tsv(BASE/f"{lang}_historical/{phase}_hypergeometric_variants.tsv") for lang in ("python","r"))
        require([(r['profile'],r['go_id']) for r in left]==[(r['profile'],r['go_id']) for r in right],"Diagnostic identities differ")
        for k in ("p_value","p_bh"):
            np.testing.assert_allclose([float(r[k]) for r in left],[float(r[k]) for r in right],rtol=1e-10,atol=0)
    report['simulation_tables']=sum(t['method']!='hypergeometric' for run in report['runs'] for t in run['tests'])
    report['simulation_result_rows']=sum(t['pathways'] for run in report['runs'] for t in run['tests'] if t['method']!='hypergeometric')
    report['null_set_evaluations']=sum(t['distinct_sizes']*10000 for run in report['runs'] for t in run['tests'] if t['method']!='hypergeometric')
    return tables,report


def compare_manuscript(tables):
    manuscript=json.loads((BASE/'manuscript_reference.json').read_text(encoding='utf-8'))
    allrows=[]; summaries=[]
    for reference_name, refs in (("S2",manuscript['records']),("S1",manuscript['selected']['records'])):
        for profile in ('full','historical'):
            for language in ('python','r'):
                for phase in ('short','middle'):
                    for method in (('hypergeometric','sumstat','gsea','maxmean') if profile=='full' else ('sumstat',)):
                        reference=[r for r in refs if r['phase']==phase and r['method']==method]
                        if not reference: continue
                        actual={r['go_id']:r for r in tables[profile,language,phase,method]}
                        pcol,qcol=('hypergeometric_p','hypergeometric_bh') if method=='hypergeometric' else ('p_value','p_bh')
                        differences=[]; exclusions=0; outside=0; sig_agree=0
                        for ref in reference:
                            found=actual.get(ref['go_id'])
                            row=dict(reference=reference_name,profile=profile,language=language,phase=phase,method=method,
                                go_id=ref['go_id'],description=ref['description'],reference_p=ref['p_value'],reference_bh=ref['p_bh'],
                                reference_annotation_size=ref['annotation_size'],reference_ct_count=ref['ct_count'],
                                source_locator=ref.get('source_sheet',str(ref.get('source_table')))+':'+ref.get('source_range',str(ref.get('source_row'))))
                            if found:
                                p,q=float(found[pcol]),float(found[qcol])
                                difference=abs(p-ref['p_value']); differences.append(difference)
                                screen=bool(difference>diagnostic_limit(p,ref['p_value'])) if method!='hypergeometric' else ''
                                outside+=int(screen is True); sig_agree+=int((q<.05)==(ref['p_bh']<.05))
                                row.update(status='compared',tested_size=found['tested_size'],ct_count=found['ct_count'],p_value=p,p_bh=q,
                                    absolute_p_difference=difference,outside_six_SE_screen=screen)
                            else:
                                exclusions+=1
                                row.update(status='excluded_by_measured_size_bounds',tested_size='',ct_count='',p_value='',p_bh='',
                                    absolute_p_difference='',outside_six_SE_screen='')
                            allrows.append(row)
                        summaries.append(dict(reference=reference_name,profile=profile,language=language,phase=phase,method=method,
                            reference_rows=len(reference),compared=len(differences),excluded=exclusions,max_p_difference=max(differences,default=None),
                            outside_six_SE_screen=outside if method!='hypergeometric' else None,bh_0_05_classification_agreement=sig_agree))
    write_tsv(allrows,ROOT/'reports/recalculated_pathway_comparison_2026-10-21.tsv',list(allrows[0]))
    # Test the numeric fingerprint of the matrix/size-column counting mistakes.
    fingerprints=[]
    for language in ('python','r'):
        for phase in ('short','middle'):
            actual={r['go_id']:r for r in read_tsv(BASE/f'{language}_historical/{phase}_hypergeometric_variants.tsv') if r['profile']=='legacy_container_bug'}
            refs=[r for r in manuscript['records'] if r['phase']==phase and r['method']=='hypergeometric']
            errors=[abs(float(actual[r['go_id']]['p_value'])-r['p_value']) for r in refs]
            fingerprints.append(dict(language=language,phase=phase,N=int(next(iter(actual.values()))['universe_size']),
                M=int(next(iter(actual.values()))['ct_population']),records=len(refs),max_absolute_error=max(errors),
                within_5e_8=sum(e<=5e-8 for e in errors),interpretation='Numeric fingerprint; historical execution state not directly recorded'))
    # S1 and S2 are separate source tables, with demonstrably different short HGT values.
    lookup={(r['phase'],r['method'],r['go_id']):r for r in manuscript['records']}
    discrepancies=[]
    for r in manuscript['selected']['records']:
        other=lookup[(r['phase'],r['method'],r['go_id'])]
        if abs(r['p_value']-other['p_value'])>max(5e-3*max(r['p_value'],other['p_value']),5e-16):
            discrepancies.append(dict(phase=r['phase'],method=r['method'],go_id=r['go_id'],S1_p=r['p_value'],S2_p=other['p_value']))
    selected_hypergeometric=[]
    for phase in ('short','middle'):
        qc = next(r for r in read_tsv(BASE/'python_full/input_summary.tsv') if r['phase']==phase)
        diagnostic=read_tsv(BASE/f'python_historical/{phase}_hypergeometric_variants.tsv')
        legacy=next(r for r in diagnostic if r['profile']=='legacy_container_bug')
        parameters={'genome_measured_CT_total':(int(qc['universe_size']),int(qc['ct_in_universe'])),
                    'legacy_container_bug':(int(legacy['universe_size']),int(legacy['ct_population']))}
        for ref in manuscript['selected']['records']:
            if ref['phase']!=phase or ref['method']!='hypergeometric': continue
            predictions={name:float(hypergeom.sf(ref['ct_count']-1,N,ref['annotation_size'],M)) for name,(N,M) in parameters.items()}
            matches=[name for name,p in predictions.items() if abs(p-ref['p_value'])<=max(.005*max(p,ref['p_value']),5e-16)]
            selected_hypergeometric.append(dict(phase=phase,go_id=ref['go_id'],reported_p=ref['p_value'],
                calculated_from_reported_n_and_x=predictions,matching_profiles_at_print_precision=matches))
    archived_sumstat=[]
    archive_root=SOURCE_ROOT/'研究法/大豆基因/篩選後/SUMSTAT'
    for phase in ('short','middle'):
        source=archive_root/f'SUMSTAT_ct_{phase}_pavlue.txt'
        with source.open(encoding='utf-8-sig',newline='') as stream:
            saved={r[0]:float(r[4]) for r in csv.reader(stream,delimiter='\t') if len(r)>=6 and r[0].startswith('GO:')}
        refs=[r for r in manuscript['records'] if r['phase']==phase and r['method']=='sumstat']
        archived_sumstat.append(dict(phase=phase,source=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            reference_rows=len(refs),same_p_within_1e_12=sum(r['go_id'] in saved and abs(saved[r['go_id']]-r['p_value'])<1e-12 for r in refs)))
    return dict(source_workbook=manuscript['source'],source_sha256=manuscript['source_sha256'],
        source_docx=manuscript['selected']['source'],source_docx_sha256=manuscript['selected']['source_sha256'],
        publisher_identity_verified=False, comparisons=summaries,hypergeometric_fingerprints=fingerprints,
        selected_table_discrepancies=discrepancies,selected_hypergeometric_fingerprints=selected_hypergeometric,
        archived_sumstat_reference_matches=archived_sumstat,
        paper_reproduction_verified=False)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,default=SOURCE_ROOT,
                        help='Original project directory containing archived reference files')
    SOURCE_ROOT=parser.parse_args().source_root
    scores=check_scores()
    tables,report=verify_runs()
    report.update(score_verification=scores,manuscript_comparison=compare_manuscript(tables),
        execution_validation='passed',paper_reproduction_verified=False,
        diagnostic_note='Six-standard-error screens are broad diagnostics, not equivalence tests. Monte Carlo draws use native, different RNGs.')
    execution_path=BASE/'execution_record.json'
    if execution_path.exists():
        execution=json.loads(execution_path.read_text(encoding='utf-8'))
        require({r['run'] for r in execution}=={'python_full','r_full','python_historical','r_historical'},'Incomplete execution record')
        require(all(r['observed_exit_code']==0 for r in execution),'A run did not exit successfully')
        report['execution_records']=execution
    else:
        # CLI users can verify outputs without a separately recorded process log.
        report['execution_records']=[]
        report['execution_exit_status_note']='No process-exit record supplied; verification covers saved outputs.'
    files=list(BASE.rglob('*'))
    for pattern in ('python/*.py','R/*.R','scripts/*recalculat*','scripts/*historical_diagnostics*','config.recalculated*'):
        files.extend(ROOT.glob(pattern))
    report['sha256']={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(files)) if p.is_file()}
    manifest=json.loads((ROOT/'reports/source_manifest.json').read_text(encoding='utf-8'))
    for item in manifest:
        require(hashlib.sha256((ROOT/item['staged']).read_bytes()).hexdigest()==item['staged_sha256'],f"Changed staged input: {item['staged']}")
    report['verified_staged_hashes']=len(manifest)
    destination=ROOT/'reports/recalculated_full_verification_2026-10-21.json'
    destination.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(f"PASS: {report['simulation_tables']} simulation tables; {report['null_set_evaluations']:,} null evaluations")
    print('Both phases match all saved score records in both languages. Paper reproduction is NOT verified.')
    print(destination)
