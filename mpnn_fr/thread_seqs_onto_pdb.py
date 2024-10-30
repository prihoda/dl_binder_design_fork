#!/usr/bin/env python

import argparse
import glob
import os
print('Importing PyRosetta...')
from pyrosetta import init, pose_from_pdb
from pyrosetta.rosetta import core
import util_protein_mpnn as mpnn_util
from Bio import SeqIO
import re


def thread_mpnn_seq(pdb_path, binder_seq):
    '''
    Thread the binder sequence onto the pose being designed
    '''
    pose = pose_from_pdb(pdb_path)
    rsd_set = pose.residue_type_set_for_pose(core.chemical.FULL_ATOM_t)

    for resi, mut_to in enumerate(binder_seq):
        resi += 1  # 1 indexing
        name3 = mpnn_util.aa_1_3[mut_to]
        new_res = core.conformation.ResidueFactory.create_residue(rsd_set.name_map(name3))
        pose.replace_residue(resi, new_res, True)
    return pose


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument( "-pdbdir", type=str, default="", help='The name of an input directory of pdbs to run through the model' )
    parser.add_argument( "-fastadir", type=str, default="", help='The name of an input directory of FASTAs with sequences' )
    parser.add_argument( "-outpdbdir", type=str, default="outputs", help='The directory to which the output PDB files will be written, used if the -pdbdir arg is active' )

    args = parser.parse_args()

    print('Initializing PyRosetta...')
    # TODO is this even needed if we are not relaxing?
    init("-beta_nov16 -in:file:silent_struct_type binary -mute all" +
             " -use_terminal_residues true -mute basic.io.database core.scoring")

    pdb_paths = glob.glob(os.path.join(args.pdbdir, '*.pdb'))
    print(f'Processing {len(pdb_paths):,} PDB paths...')
    # iterate over PDBs
    for i, pdb_path in enumerate(pdb_paths):
        tag = os.path.basename(pdb_path).rsplit('.', 1)[0]
        records = list(SeqIO.parse(os.path.join(args.fastadir, tag+'.fa'), "fasta"))
        records = records[1:] # first record is input sequence
        print(f'Threading {len(records):,} sequences onto PDB {i+1}/{len(pdb_paths)}: {tag}')
        # iterate over sequences in fasta file
        for record in records:
            assert 'sample=' in record.description, f'Expected ProteinMPNN header with sample number, got: {record.description}'
            sample_number = re.search('sample=(\d+)', record.description).group(1)
            seq = str(record.seq)
            # thread sequence onto pose
            pose = thread_mpnn_seq(pdb_path, seq)
            # save pdb
            if not os.path.exists(args.outpdbdir):
                os.makedirs(args.outpdbdir)
            pose.dump_pdb(os.path.join(args.outpdbdir, tag + f'_sample_{sample_number}' + '.pdb'))
