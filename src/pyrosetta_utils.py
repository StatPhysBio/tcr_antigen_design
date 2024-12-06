

import pyrosetta
from pyrosetta.rosetta import core, protocols
import numpy as np


default_flags = '-ignore_unrecognized_res 1 -include_current -ex1 -ex2 -mute all -include_sugars -ignore_zero_occupancy false -obey_ENDMDL 1'
wet_flags = '-ignore_unrecognized_res 1 -include_current -ex1 -ex2 -mute all -include_sugars -ignore_zero_occupancy false -obey_ENDMDL 1 -ignore_waters 0'
mute_flag = ' -mute all'



def repack_residues(
    scorefxn,
    positions, # 1-indexed pose numbering
    pose,
):
    ''' Repack the sidechains at the residues in "positions" 
    '''

    tf = core.pack.task.TaskFactory()
    tf.push_back(core.pack.task.operation.InitializeFromCommandline()) # use -ex1 and -ex2 rotamers if requested

    # dont allow any design
    op = core.pack.task.operation.RestrictToRepacking()
    tf.push_back(op)

    # freeze residues not in the positions list
    op = core.pack.task.operation.PreventRepacking()
    for i in range(1,pose.size()+1):
        if i not in positions:
            op.include_residue(i)
        else:
            print('repacking at residue', i)
    tf.push_back(op)
    packer = protocols.minimization_packing.PackRotamersMover()
    packer.task_factory(tf)
    packer.score_function(scorefxn)

    # show the packer task
    print(tf.create_task_and_apply_taskoperations(pose))
    
    packer.apply(pose)

def fastrelax_full_pose(pose,
                        scorefxn,
                        relax_backbone = False,
                        nrepeats = 1):
    
    resnums = list(range(1, pose.size()+1))

    fastrelax_positions(scorefxn,
                        resnums if relax_backbone else [],
                        resnums,
                        pose,
                        nrepeats = nrepeats)

def fastrelax_positions(
        scorefxn,
        backbone_flexible_positions,
        sidechain_flexible_positions,
        pose,
        nrepeats = 1,
):
    ''' "Relax" iterates between repacking and gradient-based minimization
    here we are doing "cartesian" relax, which allows bond lengths and angles to change slightly
    (the positions of the atoms are the degrees of freedom, rather than the internal coordinates)
    So the scorefxn should have terms to constrain these near ideal values, eg ref2015_cart.wts
    '''
    # movemap:
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    mm.set_jump(False)

    for i in backbone_flexible_positions:
        mm.set_bb(i, True)

    for i in sidechain_flexible_positions:
        mm.set_chi(i, True)

    fr = protocols.relax.FastRelax(scorefxn_in=scorefxn,
                                   standard_repeats=nrepeats)
    fr.cartesian(True)
    fr.set_movemap(mm)
    fr.set_movemap_disables_packing_of_fixed_chi_positions(True)

    # For non-Cartesian scorefunctions, use "dfpmin_armijo_nonmonotone"
    fr.min_type("lbfgs_armijo_nonmonotone")
    fr.apply(pose)
    

def find_calpha_neighbors(
    core_positions,
    distance_threshold,
    pose
):
    ''' This function finds neighbors of the residues in "core_positions" based on Calpha-Calpha distance
    '''
    # include all the 'core' positions as neighbors (of themselves, e.g.)
    nbr_positions = set(core_positions)

    distance_threshold_squared = distance_threshold**2
    
    for i in range(1, pose.size()+1): # stupid Rosetta 1-indexing
        rsd1 = pose.residue(i)
        try:
            rsd1_CA = rsd1.xyz("CA") # access by string is a *little* slow; could use integer indices
            for j in core_positions:
                rsd2 = pose.residue(j)
                if rsd1_CA.distance_squared(rsd2.xyz("CA")) <= distance_threshold_squared:
                    nbr_positions.add(i)
                    break
        except:
            continue
    return nbr_positions

def get_pdb_residue_info(
    pose,
    resnum,
):
    pi = pose.pdb_info()
    return (pi.chain(resnum), pi.number(resnum), pi.icode(resnum))

def get_pose_residue_number(
    pose,
    chain,
    resnum,
    icode=' ',
):
    return pose.pdb_info().pdb2pose(chain, resnum, icode)

def make_mutations(
    mutations,
    pose,
    verbose=False
):
    ''' Make sequence changes and repack the mutated positions
    
    mutations is a dictionary mapping from pose number to new 1-letter aa
    mutations is 1-indexed
    
    Note that we don't specify the score function! I guess the packer here is
    using a default fullatom scorefunction... Huh
    '''
    oldseq = pose.sequence()

    tf = core.pack.task.TaskFactory()
    #tf.push_back(core.pack.task.operation.InitializeFromCommandline()) # potentially include extra rotamers

    # freeze non-mutated
    op = core.pack.task.operation.PreventRepacking()
    for i in range(1,pose.size()+1):
        if i not in mutations:
            op.include_residue(i)
    tf.push_back(op)

    # force desired sequence at mutations positions
    for i, aa in mutations.items():
        op = core.pack.task.operation.RestrictAbsentCanonicalAAS()
        op.include_residue(i)
        op.keep_aas(aa)
        tf.push_back(op)
        if verbose:
            print('make mutation:', i, oldseq[i-1], '-->', aa)

    packer = protocols.minimization_packing.PackRotamersMover()
    packer.task_factory(tf)

    # show the packer task
    if verbose:
        print(tf.create_task_and_apply_taskoperations(pose))

    packer.apply(pose)
