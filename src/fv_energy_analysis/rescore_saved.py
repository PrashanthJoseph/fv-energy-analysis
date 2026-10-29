"""Coordinate-based scoring-only comparison of an already prepared batch.

No preparation, minimization, numbering or geometric detector is invoked.
Historical ligand.wts is used with pre-Talaris defaults in fresh processes.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
from datetime import datetime, timezone
from hashlib import sha256
import json
from importlib import metadata
import math
import multiprocessing
from pathlib import Path
import traceback

from .paths import native_path
from .report import write_csv

_ENGINE = _SCOREFXN = None


def _json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def _save(path, data):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def collect_inputs(batch: Path) -> dict:
    batch = Path(batch)
    with (batch / "batch_summary.csv").open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    candidates, parent_digest, params_digest, parent, params, ligand_chain = [], None, None, None, None, None
    digests = {}
    def digest(path):
        path = Path(path).resolve()
        if path not in digests:
            digests[path] = sha256(path.read_bytes()).hexdigest()
        return digests[path]
    for row in rows:
        if row["status"] != "complete":
            continue
        folder = native_path(row["result_folder"])
        summary = _json(folder / "final/summary.json")
        manifest = _json(folder / "workflow_manifest.json")
        p, m = native_path(summary["parent"]), native_path(summary["mutant"])
        parameter = native_path(manifest["params"])
        for path in (p, m, parameter):
            if not path.is_file():
                raise FileNotFoundError(f"Saved rescore input is missing: {path}")
        pd, qd = digest(p), digest(parameter)
        recorded = manifest.get("scoring_protocol", {}).get("params_sha256")
        if not recorded or recorded != qd:
            raise ValueError(f"Parameter hash mismatch: {folder}")
        chain = summary["ligand_chains"]
        if len(chain) != 1:
            raise ValueError("Scoring-only rescore requires exactly one ligand chain")
        if parent is None:
            parent, params, parent_digest, params_digest, ligand_chain = p, parameter, pd, qd, chain[0]
        elif (pd, qd, chain[0]) != (parent_digest, params_digest, ligand_chain):
            raise ValueError(f"Candidates do not share identical parent, params and ligand chain: {folder}")
        name = Path(row["candidate"].replace("\\", "/")).name
        candidates.append({"candidate": name, "result_folder": folder, "pdb": m, "source_row": row})
    if not candidates:
        raise ValueError("No completed prepared candidates were found")
    if len({c['candidate'] for c in candidates}) != len(candidates):
        raise ValueError("Duplicate candidate names in source batch")
    return {"parent": parent, "params": params, "parent_sha256": parent_digest,
            "params_sha256": params_digest, "ligand_chain": ligand_chain,
            "candidates": candidates, "source_count": len(rows), "excluded_count": len(rows)-len(candidates)}


def _init_worker(params, weights):
    global _ENGINE, _SCOREFXN
    from .pyrosetta_score import _initialize, _legacy_weights, _scorefunction
    _ENGINE = _initialize([params], legacy=_legacy_weights(weights))
    _SCOREFXN = _scorefunction(_ENGINE, weights)


def _read_score(path):
    lines = Path(path).read_text().splitlines()
    return {k: float(v) for k, v in zip(lines[0].split()[1:-1], lines[1].split()[1:-1])}


def _score_one(pdb, output, ligand_chain, protocol):
    from .pyrosetta_score import write_pair_energies, write_interface_metrics
    from .hbonds import write_pose_hbonds
    output, pdb = Path(output), Path(pdb)
    output.mkdir(parents=True, exist_ok=True)
    identity = {"protocol": protocol, "pdb_sha256": sha256(pdb.read_bytes()).hexdigest()}
    cached = output / "result.json"
    if cached.is_file():
        saved = _json(cached)
        if saved.get("identity") != identity:
            raise ValueError(f"Refusing reuse under different inputs/protocol: {output}")
        if all((output / f).is_file() for f in saved.get('artifact_sha256', {})) and saved.get('artifact_sha256'):
            if any(sha256((output/f).read_bytes()).hexdigest()!=digest for f,digest in saved['artifact_sha256'].items()):
                raise ValueError(f"Cached scoring artifact was modified: {output}")
            return saved
    pose = _ENGINE.pose_from_file(str(pdb))
    description = pdb.stem
    write_pair_energies(pose, _SCOREFXN, output / "pair_energies.sc", description)
    write_interface_metrics(pose, _SCOREFXN, ligand_chain, output / "interface.sc", description)
    write_pose_hbonds(pose, output / "hbonds.csv", scorefxn=_SCOREFXN)
    metrics = _read_score(output / "interface.sc")
    # Weighted pair/one-body records must reconstruct the full bound score.
    pair_sum = sum(float(line.split()[-2]) for line in (output / "pair_energies.sc").read_text().splitlines()[1:])
    if not math.isclose(pair_sum, metrics['bound_total'], abs_tol=.01):
        raise ValueError(f"Pair energy reconstruction failed: {pair_sum} vs {metrics['bound_total']}")
    terms = {k:v for k,v in metrics.items() if k.startswith(f'if_{ligand_chain}_')}
    if not math.isclose(sum(terms.values()), metrics[f'interface_delta_{ligand_chain}'], abs_tol=.00005):
        raise ValueError("Interface term decomposition does not reconstruct interface energy")
    with (output / "hbonds.csv").open(newline="") as handle:
        hbonds = list(csv.DictReader(handle))
    weighted_hb = sum(float(h['hbond_weighted_energy_REU']) for h in hbonds)
    from .pyrosetta_score import _active_terms, _weighted_totals
    bound = _weighted_totals(pose, _SCOREFXN, _active_terms(_SCOREFXN))
    expected_hb = sum(v for k,v in bound.items() if k.startswith('hbond_'))
    if not math.isclose(weighted_hb, expected_hb, abs_tol=.0001):
        raise ValueError(f"H-bond energy reconstruction failed: {weighted_hb} vs {expected_hb}")
    ligand_hb = [h for h in hbonds if h['donor_residue'].startswith(ligand_chain+':') != h['acceptor_residue'].startswith(ligand_chain+':')]
    result = {"identity": identity, "pdb": str(pdb), "metrics": metrics,
              "artifact_sha256": {f:sha256((output/f).read_bytes()).hexdigest() for f in ('pair_energies.sc','interface.sc','hbonds.csv')},
              "ligand_hbond_count": len(ligand_hb), "validation": {
                  "pair_reconstruction_error_REU": pair_sum-metrics['bound_total'],
                  "hbond_reconstruction_error_REU": weighted_hb-expected_hb,
                  "pose_residues": pose.total_residue()}, "status": "complete"}
    _save(cached, result)
    return result


def run(source_batch, output, jobs=3, weights="ligand", limit=None):
    from .pyrosetta_score import _legacy_weights
    inputs = collect_inputs(Path(source_batch))
    output = Path(output).resolve()
    if output == Path(source_batch).resolve():
        raise ValueError("Rescoring must use a separate output folder")
    output.mkdir(parents=True, exist_ok=True)
    # Read the exact installed weights; changing the file invalidates reuse.
    import importlib.util
    database = Path(importlib.util.find_spec('pyrosetta').origin).parent / 'database'
    weightfile = Path(weights) if Path(weights).is_file() else database / 'scoring/weights' / (weights.removesuffix('.wts')+'.wts')
    weighttext = weightfile.read_text()
    protocol = {"version": 2, "pyrosetta_version": metadata.version('pyrosetta'),
                "implementation_sha256": {name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                           for name in ('rescore_saved.py','pyrosetta_score.py','hbonds.py')},
                "weights": weights, "weights_sha256": sha256(weightfile.read_bytes()).hexdigest(),
                "legacy_pre_talaris": _legacy_weights(weights), "params_sha256": inputs['params_sha256'],
                "parent_sha256": inputs['parent_sha256'], "minimization": False,
                "charge_model": "unchanged supplied ligand params", "ligand_chain": inputs['ligand_chain']}
    selected = inputs['candidates'][:limit] if limit else inputs['candidates']
    manifest_path = output / 'rescore_manifest.json'
    manifest = {"status":"running", "source_batch":str(Path(source_batch).resolve()), "protocol":protocol,
                "source_count":inputs['source_count'], "selected_count":len(selected),
                "excluded_source_failures":inputs['excluded_count'], "jobs":jobs,
                "started_utc":datetime.now(timezone.utc).isoformat(), "weights_contents":weighttext}
    if manifest_path.is_file() and _json(manifest_path)['protocol'] != protocol:
        raise ValueError("Existing output belongs to a different scoring protocol")
    _save(manifest_path,manifest)
    # A one-process executor scores the parent once, then frees its memory.
    context = multiprocessing.get_context('spawn')
    with concurrent.futures.ProcessPoolExecutor(max_workers=1, mp_context=context,
            initializer=_init_worker, initargs=(str(inputs['params']),weights)) as pool:
        parent = pool.submit(_score_one,str(inputs['parent']),str(output/'parent'),inputs['ligand_chain'],protocol).result()
    (output/'rescore.log').touch(exist_ok=True)
    results, errors = {}, []
    with concurrent.futures.ProcessPoolExecutor(max_workers=jobs,mp_context=context,
            initializer=_init_worker, initargs=(str(inputs['params']),weights)) as pool:
        futures={pool.submit(_score_one,str(item['pdb']),str(output/item['candidate'].removesuffix('.pdb')),
                            inputs['ligand_chain'],protocol):item for item in selected}
        for future in concurrent.futures.as_completed(futures):
            item=futures[future]
            try:
                results[item['candidate']]=future.result()
            except Exception as exc:
                errors.append({'candidate':item['candidate'],'error':str(exc),'traceback':traceback.format_exc()})
            with (output/'rescore.log').open('a') as handle:
                handle.write(json.dumps({'candidate':item['candidate'],'status':'complete' if item['candidate'] in results else 'failed'})+'\n')
    rows=[]
    prefix=f"if_{inputs['ligand_chain']}_"
    for item in selected:
        name=item['candidate']
        if name not in results: continue
        m,p=results[name],parent
        row={'candidate':name,'substitutions':item['source_row'].get('substitutions',''),
             'source_result_folder':str(item['result_folder']),'weights':weights,
             'parent_interface_REU':p['metrics'][f"interface_delta_{inputs['ligand_chain']}"],
             'mutant_interface_REU':m['metrics'][f"interface_delta_{inputs['ligand_chain']}"],
             'interface_delta_REU':m['metrics'][f"interface_delta_{inputs['ligand_chain']}"]-p['metrics'][f"interface_delta_{inputs['ligand_chain']}"],
             'parent_ligand_hbonds':p['ligand_hbond_count'],'mutant_ligand_hbonds':m['ligand_hbond_count'],
             'parent_bound_total_REU':p['metrics']['bound_total'],'mutant_bound_total_REU':m['metrics']['bound_total'],
             'ref2015_parent_interface_REU':item['source_row'].get('interface_parent_REU',''),
             'ref2015_mutant_interface_REU':item['source_row'].get('interface_mutant_REU',''),
             'ref2015_interface_delta_REU':item['source_row'].get('interface_delta_REU','')}
        for key in p['metrics']:
            if key.startswith(prefix):
                term=key.removeprefix(prefix)
                row.update({f'parent_{term}_REU':p['metrics'][key],f'mutant_{term}_REU':m['metrics'][key],
                            f'delta_{term}_REU':m['metrics'][key]-p['metrics'][key]})
        rows.append(row)
    rows.sort(key=lambda r:r['interface_delta_REU'])
    for rank,row in enumerate(rows,1): row['rank_by_interface_delta']=rank
    write_csv(output/'rosettaligand_rescore.csv',rows)
    _save(output/'failures.json',errors)
    manifest.update(status='complete_with_failures' if errors else 'complete',completed_count=len(rows),failure_count=len(errors),
                    finished_utc=datetime.now(timezone.utc).isoformat(),summary_csv=str(output/'rosettaligand_rescore.csv'))
    _save(manifest_path,manifest)
    print(json.dumps({k:manifest[k] for k in ('status','completed_count','failure_count','summary_csv')}),flush=True)
    return 1 if errors else 0


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_batch',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--jobs',type=int,default=3)
    parser.add_argument('--weights',default='ligand')
    parser.add_argument('--limit',type=int,help='Small validation subset; parent is still scored once')
    args=parser.parse_args(argv)
    if args.jobs<1: parser.error('--jobs must be positive')
    try:
        return run(args.source_batch,args.output,args.jobs,args.weights,args.limit)
    except Exception as exc:
        path=args.output/'rescore_manifest.json'
        if path.is_file():
            manifest=_json(path)
            manifest.update(status='failed',error=str(exc),traceback=traceback.format_exc())
            _save(path,manifest)
        raise


if __name__=='__main__':
    raise SystemExit(main())
