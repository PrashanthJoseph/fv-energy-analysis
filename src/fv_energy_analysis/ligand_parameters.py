"""Prepare Rosetta parameters from an explicit SDF/MOL2 ligand, never a PDB guess."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

from .resources import molfile_to_params_path
from .ligand_chemistry import inspect_params_chemistry


def _mol2_atoms(path: Path):
    """Read atom rows only when the MOL2 declares finite partial charges."""
    lines = path.read_text(encoding='utf-8').splitlines()
    if lines.count('@<TRIPOS>MOLECULE') != 1:
        raise ValueError('MOL2 must contain exactly one ligand')
    header = lines.index('@<TRIPOS>MOLECULE')
    charge_model = lines[header + 4].strip() if len(lines) > header + 4 else ''
    if not charge_model or charge_model.upper() in {'NO_CHARGES', 'NONE'} or charge_model.startswith('@'):
        raise ValueError('MOL2 has no declared partial-charge model; provide a charged MOL2 or an SDF for explicit Gasteiger assignment')
    if '@<TRIPOS>ATOM' not in lines:
        raise ValueError('MOL2 atom table is missing')
    atoms = []
    for line in lines[lines.index('@<TRIPOS>ATOM') + 1:]:
        if line.startswith('@<TRIPOS>'):
            break
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        fields = line.split()
        if len(fields) < 9:
            raise ValueError('Every MOL2 atom must include a partial-charge column')
        try:
            charge = float(fields[8])
        except ValueError as error:
            raise ValueError(f'Invalid MOL2 partial charge: {fields[8]!r}') from error
        if not np.isfinite(charge):
            raise ValueError(f'Non-finite MOL2 partial charge: {fields[8]!r}')
        atoms.append(fields)
    if not atoms:
        raise ValueError('MOL2 atom table is empty')
    return atoms, charge_model


def prepare_parameters(source: Path, output: Path, name='E2L', chain='X', generator=None):
    """Open Babel conversion -> Rosetta generator -> explicit atom/charge trace.

    The parameter generator rounds MOL2 charges to two decimal places. Restore
    the source MOL2 charges by the verified output-PDB atom-order mapping;
    do not invent a redistributed charge model or patch aromatic atom types.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    source, output = source.resolve(), output.resolve()
    if source.suffix.lower() not in ('.sdf','.mol2'):
        raise ValueError('Use an explicit SDF or MOL2 chemical source')
    if output.exists() and any(output.iterdir()):
        raise ValueError('Parameter output must be empty or new; no shared parameters are overwritten')
    output.mkdir(parents=True,exist_ok=True)
    commands=[]
    def run(command):
        commands.append([str(x) for x in command])
        with (output/'parameterization.log').open('a') as log:
            subprocess.run(command,cwd=output,check=True,stdout=log,stderr=subprocess.STDOUT)
    mol2=output/'ligand.mol2'
    if source.suffix.lower()=='.sdf':
        molecules=list(Chem.SDMolSupplier(str(source),removeHs=False))
        if len(molecules)!=1 or molecules[0] is None:
            raise ValueError('SDF must contain exactly one readable ligand')
        original=molecules[0]
        obabel=shutil.which('obabel')
        if not obabel:
            raise RuntimeError('Open Babel obabel is required to convert SDF; existing charged MOL2 input does not need it')
        run([obabel,'-isdf',str(source),'-omol2','-O',str(mol2),'--partialcharge','gasteiger'])
    else:
        shutil.copy2(source,mol2)
    atom_lines, charge_model = _mol2_atoms(mol2)
    if source.suffix.lower()=='.mol2':
        original=Chem.MolFromMol2File(str(source),removeHs=False)
    converted=Chem.MolFromMol2File(str(mol2),removeHs=False)
    if original is None or converted is None:
        raise ValueError('Chemical input or Open Babel output is unreadable')
    if original.GetNumAtoms()!=converted.GetNumAtoms():
        raise ValueError('Conversion changed atom count')
    if [a.GetAtomicNum() for a in original.GetAtoms()] != [a.GetAtomicNum() for a in converted.GetAtoms()]:
        raise ValueError('Conversion reordered atom elements; index-based mapping is unsafe')
    displacement=float(np.max(np.linalg.norm(original.GetConformer().GetPositions()-converted.GetConformer().GetPositions(),axis=1)))
    if displacement>2e-4:
        raise ValueError(f'Conversion changed coordinates by {displacement} A')
    def bonds(molecule):
        return {(min(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),max(b.GetBeginAtomIdx(),b.GetEndAtomIdx())):
                ('ar' if b.GetIsAromatic() else b.GetBondTypeAsDouble()) for b in molecule.GetBonds()}
    if bonds(original)!=bonds(converted):
        raise ValueError('Conversion changed connectivity or bond orders')
    Chem.AssignStereochemistry(original,cleanIt=True,force=True)
    Chem.AssignStereochemistry(converted,cleanIt=True,force=True)
    if Chem.MolToInchiKey(original)!=Chem.MolToInchiKey(converted):
        raise ValueError('Conversion changed molecular stereochemistry')
    if len(atom_lines)!=original.GetNumAtoms():
        raise ValueError('MOL2 atom table is inconsistent')
    run([sys.executable,str(molfile_to_params_path(generator)),'-n',name,'-p',name,'--chain',chain,'ligand.mol2'])
    template=output/f'{name}_0001.pdb'
    pdb_lines=[line for line in template.read_text().splitlines() if line.startswith('HETATM')]
    if len(pdb_lines)!=len(atom_lines):raise ValueError('Rosetta generator changed atom count')
    pdb_xyz=np.array([[float(line[30:38]),float(line[38:46]),float(line[46:54])] for line in pdb_lines])
    differences=pdb_xyz-converted.GetConformer().GetPositions()
    if float(np.max(np.abs(differences-differences.mean(axis=0))))>0.0015:
        raise ValueError('Rosetta output atom order or relative coordinates changed')
    source_charges={line[12:16].strip():float(atom[8]) for line,atom in zip(pdb_lines,atom_lines)}
    if len(source_charges)!=len(atom_lines):raise ValueError('Rosetta atom names are not unique')
    params=output/f'{name}.params'
    lines=params.read_text().splitlines()
    for i,line in enumerate(lines):
        if line.startswith('ATOM '):
            fields=line.split(); prefix,_,_=line.rpartition(fields[-1])
            lines[i]=prefix+f'{source_charges[fields[1]]:.8f}'
    params.write_text('\n'.join(lines)+'\n')
    reference_smiles=Chem.MolToSmiles(Chem.RemoveHs(original),isomericSmiles=True)
    chemistry=inspect_params_chemistry(params,reference_smiles)
    if not chemistry['valid']:
        raise ValueError('Generated chemical parameters failed direct inspection: '+str(chemistry['errors']))
    atom_mapping=[{'source_atom_index':i+1,'source_element':original.GetAtomWithIdx(i).GetSymbol(),
                  'mol2_name':atom[1],'mol2_type':atom[5], 'rosetta_name':line[12:16].strip(),
                  'source_partial_charge':float(atom[8])} for i,(line,atom) in enumerate(zip(pdb_lines,atom_lines))]
    from .report import write_csv,write_json
    write_csv(output/'atom_mapping.csv',atom_mapping)
    mapped_names={item['source_atom_index']-1:item['rosetta_name'] for item in atom_mapping}
    write_csv(output/'bond_mapping.csv',[
        {'source_atom1_index':b.GetBeginAtomIdx()+1,'source_atom2_index':b.GetEndAtomIdx()+1,
         'rosetta_atom1':mapped_names[b.GetBeginAtomIdx()],'rosetta_atom2':mapped_names[b.GetEndAtomIdx()],
         'bond_type':'aromatic' if b.GetIsAromatic() else str(b.GetBondType())}
        for b in original.GetBonds()])
    report={'source':str(source),'source_sha256':sha256(source.read_bytes()).hexdigest(),
            'source_formula':rdMolDescriptors.CalcMolFormula(original),'source_inchikey':Chem.MolToInchiKey(original),
            'reference_smiles':reference_smiles,'atom_count':original.GetNumAtoms(),'bond_count':original.GetNumBonds(),
            'source_to_mol2_max_coordinate_change_A':displacement,'partial_charge_model':charge_model,
            'charge_policy':'Explicit Open Babel Gasteiger assignment' if source.suffix.lower()=='.sdf' else 'Preserve supplied MOL2 charges',
            'charge_handling':'Preserved MOL2 atom charges at eight decimals; no redistribution',
            'params':str(params),'template_pdb':str(template),'chemistry':chemistry,'commands':commands}
    write_json(output/'ligand_definition.json',report)
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--name',default='E2L')
    parser.add_argument('--chain',default='X')
    parser.add_argument('--molfile-to-params')
    args=parser.parse_args(argv)
    report=prepare_parameters(args.input,args.output,args.name,args.chain,args.molfile_to_params)
    print(json.dumps({key:report[key] for key in ('source_formula','source_inchikey','atom_count','bond_count','params')}))


if __name__=='__main__':main()
