"""Run the published protein–ligand Flex ddG XML using installed PyRosetta.

Mutations use input PDB numbering: {"H": "E99F"} or
{"H": ["E99F", "R100D"], "L": "Y91F"}. All entries form ONE mutant.
Each independent replica runs in a fresh process with its own RNG seed.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import sqlite3
import statistics
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

from .validation import validate_ligand_against_params
from .resources import resource_path

SOURCE_XML = resource_path("vendor/flex_ddg/aldeghi_2018/SI_Files/Rosetta/ddG-backrub.xml")
AA3 = dict(zip(
    "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split(),
    "ACDEFGHIKLMNPQRSTVWY",
))
GAM = {
    "fa_sol": (6.940, -6.722), "hbond_sc": (1.902, -1.999),
    "hbond_bb_sc": (0.063, 0.452), "fa_rep": (1.659, -0.836),
    "fa_elec": (0.697, -0.122), "hbond_lr_bb": (2.738, -1.179),
    "fa_atr": (2.313, -1.649),
}


def parse_mutations(mutations: dict) -> list[dict]:
    if not isinstance(mutations, dict) or not mutations:
        raise ValueError("Supply a nonempty mutation dictionary, e.g. {'H': 'E99F'}")
    parsed, seen = [], set()
    for chain, values in mutations.items():
        if not isinstance(chain, str) or not re.fullmatch(r"[A-Za-z0-9]", chain):
            raise ValueError(f"Chain must be one PDB letter/digit: {chain!r}")
        if isinstance(values, str):
            values = [values]
        if not isinstance(values, list) or not values:
            raise ValueError(f"Mutations for {chain} must be a string or nonempty list")
        for value in values:
            match = re.fullmatch(r"([ACDEFGHIKLMNPQRSTVWY])(-?\d+)([A-Za-z]?)([ACDEFGHIKLMNPQRSTVWY])", value) if isinstance(value, str) else None
            if not match:
                raise ValueError(f"Invalid mutation {value!r}; use E99F (or E99AF for insertion code A)")
            old, number, insertion, new = match.groups()
            key = chain, int(number), insertion
            if key in seen:
                raise ValueError(f"Duplicate mutation position: {chain}:{number}{insertion}")
            if old == new:
                raise ValueError(f"Mutation {chain}:{value} does not change the residue")
            seen.add(key)
            parsed.append(dict(chain=chain, old=old, number=int(number), insertion=insertion, new=new))
    return sorted(parsed, key=lambda m: (m['chain'], m['number'], m['insertion']))


def validate_input(pdb: Path, params: Path, ligand_chain: str, mutations: list[dict]) -> dict:
    residues, ligand_residues = {}, set()
    for line in pdb.read_text().splitlines():
        if line.startswith("MODEL"):
            raise ValueError("Supply a single-model PDB without MODEL/ENDMDL records")
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[16:17] not in (" ", ""):
            raise ValueError("Resolve alternate atom locations before Flex ddG")
        key = line[21], int(line[22:26]), line[26].strip()
        if key in residues and residues[key] != line[17:20].strip():
            raise ValueError(f"Conflicting residue names at {key}")
        residues[key] = line[17:20].strip()
        if key[0] == ligand_chain:
            ligand_residues.add(key)
    if len(ligand_residues) != 1:
        raise ValueError("Ligand chain must contain exactly one small-molecule residue")
    for m in mutations:
        key = m['chain'], m['number'], m['insertion']
        found = AA3.get(residues.get(key, ''), residues.get(key, 'missing'))
        if m['chain'] == ligand_chain or found != m['old']:
            raise ValueError(f"{m['chain']}:{m['number']}{m['insertion']} expected {m['old']}, found {found}")
    report = validate_ligand_against_params(pdb, params, ligand_chain)
    if not report['valid']:
        raise ValueError(f"Ligand atom/parameter validation failed: {report}")
    return report


def resfile_text(mutations: list[dict], ligand_chain: str, ligand_number: str) -> str:
    # Keep the ligand out of side-chain packing; its bound pose defines this calculation.
    rows = ["NATAA", "start", f"{ligand_number} {ligand_chain} NATRO"]
    rows += [f"{m['number']}{m['insertion']} {m['chain']} PIKAA {m['new']}" for m in mutations]
    return "\n".join(rows) + "\n"


def make_xml(resfile: Path, weights: str, ligand_chain: str, steps: int,
             max_iter: int, threshold: float, ligand_index: int) -> str:
    tree = ET.fromstring(SOURCE_XML.read_text())
    for score in tree.findall("./SCOREFXNS/ScoreFunction"):
        score.set("weights", weights)
    tree.find("./TASKOPERATIONS/ReadResfile").set("filename", str(resfile))
    ET.SubElement(tree.find('RESIDUE_SELECTORS'), 'Index',
                  {'name':'ligand_residue', 'resnums':str(ligand_index)})
    freeze = ET.SubElement(tree.findall('TASKOPERATIONS')[-1], 'OperateOnResidueSubset',
                          {'name':'freeze_ligand', 'selector':'ligand_residue'})
    ET.SubElement(freeze, 'PreventRepackingRLT')
    movers = tree.find("MOVERS")
    for pack in movers.findall('PackRotamersMover') + movers.findall('BackrubProtocol'):
        pack.set('task_operations', pack.get('task_operations')+',freeze_ligand')
    minimum = movers.find("MinMover")
    minimum.set("max_iter", str(max_iter))
    minimum.set("abs_score_convergence_threshold", str(threshold))
    # Protein backbone/side chains may adjust; freeze ligand torsions and all jumps.
    mm = ET.SubElement(minimum, "MoveMap", {"bb": "1", "chi": "1", "jump": "0"})
    ET.SubElement(mm, "Span", {"begin": str(ligand_index), "end": str(ligand_index), "bb": "0", "chi": "0"})
    movers.find("InterfaceDdGMover").set("chain_name", ligand_chain)
    backrub = movers.find("BackrubProtocol")
    backrub.set("ntrials", str(steps))
    backrub.set("trajectory_stride", str(steps))
    # JD2 structure reporters require a CLI JobDistributor. Save bound PDBs below
    # through PyRosetta instead, leaving the original four-state score reporter.
    movers.remove(movers.find("ReportToDB[@name='structreport']"))
    movers.insert(0, ET.Element('DumpPdb', {'name':'dump_wt', 'fname':'wt_bound.pdb', 'tag_time':'0'}))
    movers.insert(1, ET.Element('DumpPdb', {'name':'dump_mut', 'fname':'mutant_bound.pdb', 'tag_time':'0'}))
    for protocol in movers.findall("ParsedProtocol"):
        for add in list(protocol):
            if add.get("mover_name") == "structreport":
                protocol.remove(add)
        for offset, add in enumerate(list(protocol)):
            if add.get('mover_name') == 'save_wt_bound_pose':
                protocol.insert(offset+1, ET.Element('Add', {'mover_name':'dump_wt'}))
                break
        for offset, add in enumerate(list(protocol)):
            if add.get('mover_name') == 'int_ddG_mover':
                protocol.insert(offset, ET.Element('Add', {'mover_name':'dump_mut'}))
                break
    return ET.tostring(tree, encoding="unicode")


def read_scores(database: Path) -> dict:
    with sqlite3.connect(database) as conn:
        rows = conn.execute('''SELECT batches.name, score_types.score_type_name,
                              structure_scores.score_value
            FROM structure_scores
            JOIN batches ON batches.batch_id=structure_scores.batch_id
            JOIN score_types ON score_types.batch_id=structure_scores.batch_id
              AND score_types.score_type_id=structure_scores.score_type_id''').fetchall()
    states = {}
    for name, term, value in rows:
        name = name.removesuffix("_dbreport")
        if term in states.setdefault(name, {}):
            raise ValueError(f"Multiple checkpoints found for {name}/{term}")
        states[name][term] = value
    required = {'bound_wt', 'unbound_wt', 'bound_mut', 'unbound_mut'}
    if not required.issubset(states):
        raise ValueError(f"Missing four-state scores: {required - states.keys()}")
    terms = states['bound_wt'].keys()
    delta = {term: states['bound_mut'][term] - states['unbound_mut'][term]
             - states['bound_wt'][term] + states['unbound_wt'][term] for term in terms}
    if not all(math.isfinite(v) for state in states.values() for v in state.values()):
        raise ValueError("Non-finite Rosetta scores")
    return dict(states=states, ddg_terms=delta,
                wt_dg=states['bound_wt']['total_score'] - states['unbound_wt']['total_score'],
                mut_dg=states['bound_mut']['total_score'] - states['unbound_mut']['total_score'],
                ddg=delta['total_score'])


def gam_score(terms: dict) -> float:
    return sum(math.exp(a) * math.tanh(terms.get(term, 0.0) * math.exp(b) / 2)
               for term, (a, b) in GAM.items())


def worker(jobfile: Path, replica: int) -> None:
    import pyrosetta
    job = json.loads(jobfile.read_text())
    out = jobfile.parent / f"replica_{replica:03d}"
    seed = job['seed'] + replica
    opts = ['-extra_res_fa', job['params'], '-ex1', '-ex2', '-constant_seed',
            '-jran', str(seed), '-mute', 'all', '-ignore_zero_occupancy', 'false']
    if job['weights'] == 'beta_nov16':
        opts += ['-beta_nov16']
    elif job['weights'] == 'talaris2014':
        opts += ['-restore_talaris_behavior']
    pyrosetta.init(' '.join(shlex.quote(x) for x in opts), silent=True)
    pose = pyrosetta.pose_from_file(job['pdb'])
    for m in job['parsed_mutations']:
        i = pose.pdb_info().pdb2pose(m['chain'], m['number'], m['insertion'] or ' ')
        if not i or pose.residue(i).name1() != m['old']:
            raise ValueError(f"PyRosetta residue mapping failed: {m}")
    ligand = [i for i in range(1, pose.total_residue()+1)
              if pose.pdb_info().chain(i) == job['ligand_chain']]
    if len(ligand) != 1 or not pose.residue(ligand[0]).is_ligand():
        raise ValueError("PyRosetta did not load exactly one ligand")
    if pose.fold_tree().get_jump_that_builds_residue(ligand[0]) <= 0:
        raise ValueError("Place the ligand after the protein chains so it has its own jump")
    ligand_xyz = [pose.residue(ligand[0]).xyz(a) for a in range(1, pose.residue(ligand[0]).natoms()+1)]
    with tempfile.TemporaryDirectory(prefix='fv_flex_ddg_') as work:
        previous = Path.cwd()
        try:
            os.chdir(work)
            xml = make_xml(jobfile.parent/'mutations.resfile', job['weights'], job['ligand_chain'],
                           job['steps'], job['max_iter'], job['threshold'], ligand[0])
            (out/'protocol.xml').write_text(xml)
            objects = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(xml)
            print(f"Replica {replica}: protocol parsed; running {job['steps']} backrub steps", flush=True)
            objects.get_mover('ParsedProtocol').apply(pose)
            result = read_scores(Path('ddG.db3'))
            # Backrub restores its sampled WT pose after scoring the callback.
            for m in job['parsed_mutations']:
                index = pose.pdb_info().pdb2pose(m['chain'], m['number'], m['insertion'] or ' ')
                if pose.residue(index).name1() not in (m['old'], m['new']):
                    raise ValueError(f"Unexpected output residue at {m}")
            pose.dump_pdb(str(out/'backrub_final.pdb'))
            for filename in ('wt_bound.pdb', 'mutant_bound.pdb'):
                shutil.copy2(filename, out/filename)
            mutant_residues = {}
            for line in Path('mutant_bound.pdb').read_text().splitlines():
                if line.startswith('ATOM  '):
                    mutant_residues[line[21], int(line[22:26]), line[26].strip()] = line[17:20]
            for m in job['parsed_mutations']:
                key = m['chain'], m['number'], m['insertion']
                if AA3.get(mutant_residues.get(key,'')) != m['new']:
                    raise ValueError(f'Mutant-bound snapshot did not contain requested mutation: {m}')
            result.update(replica=replica, seed=seed, pyrosetta_version=pyrosetta.version().splitlines()[-1])
            result['ligand_max_displacement_A'] = max(
                (pose.residue(ligand[0]).xyz(a)-xyz).norm()
                for a, xyz in enumerate(ligand_xyz, 1))
            # Frozen jump DOFs can still follow upstream protein motion. Check
            # conformation using internal distances, independent of rigid movement.
            residue = pose.residue(ligand[0])
            result['ligand_max_internal_distance_change_A'] = max(
                abs((residue.xyz(a)-residue.xyz(b)).norm() - (ligand_xyz[a-1]-ligand_xyz[b-1]).norm())
                for a in range(1, residue.natoms()+1) for b in range(1,a))
            if result['ligand_max_internal_distance_change_A'] > 1e-6:
                raise ValueError('Ligand internal conformation changed despite frozen torsions')
            result['gam_ddg'] = gam_score(result['ddg_terms'])
            shutil.copy2('ddG.db3', out/'ddG.db3')
            (out/'result.json').write_text(json.dumps(result, indent=2))
        finally:
            os.chdir(previous)


def run_flex_ddg(pdb, params, mutations, output, *, ligand_chain='X', replicas=35,
                 steps=35000, workers=1, weights='beta_nov16', seed=731,
                 max_iter=5000, threshold=1.0, smoke=False, dry_run=False) -> dict:
    """Run one combined mutant; resume only replicas matching identical inputs/settings."""
    pdb, params, output = (Path(p).expanduser().resolve() for p in (pdb, params, output))
    if weights not in {'beta_nov16', 'ref2015', 'talaris2014'}:
        raise ValueError('Supported weights: beta_nov16, ref2015, talaris2014')
    if not re.fullmatch(r'[A-Za-z0-9]', ligand_chain):
        raise ValueError('Ligand chain must be one PDB letter/digit')
    if smoke:
        replicas, steps, max_iter, threshold = 2, 10, 20, 200.0
    if min(replicas, steps, workers, max_iter) < 1 or threshold <= 0 or not math.isfinite(threshold):
        raise ValueError('Replica, step, worker and minimization settings must be positive')
    if seed < 1 or seed + replicas >= 2147483647:
        raise ValueError('Seed and replica seeds must be positive 32-bit integers')
    parsed = parse_mutations(mutations)
    validation = validate_input(pdb, params, ligand_chain, parsed)
    ligand_row = next(l for l in pdb.read_text().splitlines()
                      if l.startswith('HETATM') and l[21] == ligand_chain)
    job = dict(pdb=str(pdb), params=str(params), mutations=mutations, parsed_mutations=parsed,
               ligand_chain=ligand_chain, replicas=replicas, steps=steps, weights=weights,
               seed=seed, max_iter=max_iter, threshold=threshold, smoke=smoke,
               pdb_sha256=hashlib.sha256(pdb.read_bytes()).hexdigest(),
               params_sha256=hashlib.sha256(params.read_bytes()).hexdigest(),
               protocol_sha256=hashlib.sha256(SOURCE_XML.read_bytes()).hexdigest(),
               implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    try:
        job['pyrosetta_distribution_version'] = importlib.metadata.version('pyrosetta')
    except importlib.metadata.PackageNotFoundError:
        job['pyrosetta_distribution_version'] = 'not installed in this interpreter'
    output.mkdir(parents=True, exist_ok=True)
    jobfile = output/'job.json'
    if jobfile.exists() and json.loads(jobfile.read_text()) != job:
        raise ValueError('Output folder belongs to different inputs/settings; choose a new folder')
    jobfile.write_text(json.dumps(job, indent=2))
    (output/'ligand_validation.json').write_text(json.dumps(validation, indent=2))
    (output/'mutations.resfile').write_text(resfile_text(parsed, ligand_chain, ligand_row[22:27].strip()))
    if dry_run:
        print(f'Validated mutation job: {jobfile}', flush=True)
        return job
    (output/'status.json').write_text(json.dumps({'status':'running'}))
    def launch(replica):
        out = output/f'replica_{replica:03d}'
        out.mkdir(exist_ok=True)
        if (out/'result.json').exists():
            return json.loads((out/'result.json').read_text())
        command = [sys.executable, '-u', '-m', 'fv_energy_analysis.flex_ddg',
                   '--worker-job', str(jobfile), '--replica', str(replica)]
        env = os.environ.copy()
        env['PYTHONPATH'] = str(Path(__file__).resolve().parent.parent) + os.pathsep + env.get('PYTHONPATH','')
        env['OMP_NUM_THREADS'] = '1'
        with (out/'run.log').open('w') as log:
            proc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env)
        if proc.returncode:
            raise RuntimeError(f'Replica {replica} failed; see {out}/run.log')
        return json.loads((out/'result.json').read_text())
    results = []
    try:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(launch, r) for r in range(1, replicas+1)]
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                print(f"{len(results)}/{replicas} replicas complete; replica {result['replica']} ddG = {result['ddg']:.3f} REU", flush=True)
        results.sort(key=lambda r: r['replica'])
        values = [r['ddg'] for r in results]
        average_terms = {t: statistics.mean(r['ddg_terms'][t] for r in results)
                         for t in results[0]['ddg_terms']}
        summary = dict(status='complete', mutations=mutations, weights=weights,
            replicas=replicas, backrub_steps=steps, smoke=smoke, units='Rosetta energy units (REU)',
            wt_dg_mean=statistics.mean(r['wt_dg'] for r in results),
            mut_dg_mean=statistics.mean(r['mut_dg'] for r in results),
            ddg_mean=statistics.mean(values), ddg_median=statistics.median(values),
            ddg_sd=statistics.stdev(values) if replicas > 1 else None,
            ddg_sem=statistics.stdev(values)/math.sqrt(replicas) if replicas > 1 else None,
            gam_ddg_from_mean_terms=gam_score(average_terms),
            ddg_terms_mean=average_terms, pyrosetta_version=results[0]['pyrosetta_version'],
            interpretation='Negative favors binding; positive disfavors binding.',
            limitations=['Smoke settings test execution only, not affinity.' if smoke else 'Sampling spread is not experimental prediction accuracy.',
                        'Fixed internal ligand conformation; no explicit ligand pose or water sampling.',
                        'GAM is the historical nonlinear reweighting, not an estradiol-specific calibration.',
                        'Current PyRosetta build does not reproduce Rosetta 2017.52 exactly.'])
        (output/'summary.json').write_text(json.dumps(summary, indent=2))
        with (output/'replicas.csv').open('w', newline='') as handle:
            keys = ['replica', 'seed', 'wt_dg', 'mut_dg', 'ddg', 'gam_ddg',
                    'ligand_max_displacement_A', 'ligand_max_internal_distance_change_A']
            writer = csv.DictWriter(handle, fieldnames=keys, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(results)
        (output/'status.json').write_text(json.dumps({'status':'complete'}))
        return summary
    except Exception as exc:
        (output/'status.json').write_text(json.dumps({'status':'failed', 'error':str(exc)}))
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pdb', type=Path)
    parser.add_argument('--params', type=Path)
    parser.add_argument('--mutations', help='JSON dictionary, e.g. {"H":"E99F"}')
    parser.add_argument('--mutations-file', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--ligand-chain', default='X')
    parser.add_argument('--replicas', type=int, default=35)
    parser.add_argument('--steps', type=int, default=35000)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--weights', default='beta_nov16')
    parser.add_argument('--seed', type=int, default=731)
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--worker-job', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--replica', type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker_job:
        worker(args.worker_job, args.replica)
        return 0
    if not args.pdb or not args.params or not args.output or bool(args.mutations) == bool(args.mutations_file):
        parser.error('Supply --pdb, --params, --output, and exactly one of --mutations / --mutations-file')
    mutations = json.loads(args.mutations_file.read_text() if args.mutations_file else args.mutations)
    result = run_flex_ddg(args.pdb, args.params, mutations, args.output,
                         ligand_chain=args.ligand_chain, replicas=args.replicas, steps=args.steps,
                         workers=args.workers, weights=args.weights, seed=args.seed,
                         smoke=args.smoke, dry_run=args.dry_run)
    if not args.dry_run:
        print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
