"""Module for extracting protein pockets"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
from pyrosetta.rosetta.core.pose import Pose

from pyrosetta_utils import (
    find_calpha_neighbors, get_pose_residue_number, get_pdb_residue_info)
from mhc import (
    mhc_pockets, pdb_to_mhc)
from pep import pdb_to_pep
from tcr import (
    pdb_to_tcr_chains, pdb_to_tcr_pocket_template, tcr_pockets)

RADIUS = 10.

def get_pocket(
    pose: Pose
) -> Tuple[np.ndarray]:
    """Get peptide pocket from a pose of a TCR-pMHC"""
    peptide_resnums = [get_pose_residue_number(pose,'C',x) for x in range(1,11)]
    complex_resnums = list(find_calpha_neighbors(peptide_resnums,RADIUS,pose))
    pocket_resnums = np.setdiff1d(complex_resnums,peptide_resnums)
    return (
        np.array(peptide_resnums),
        np.array(pocket_resnums),
        np.array(complex_resnums))

def get_tcr_pmhc_peptide_pocket_resnums(
    pose: Pose,
    mutation_inds: List,
    chain: str    
) -> Tuple[List]:
    """Get TCR pMHC pocket resnums"""
    peptide_resnums = [
        get_pose_residue_number(pose, c, x) for x,c
        in zip(mutation_inds,chain)]
    complex_resnums = find_calpha_neighbors(peptide_resnums, RADIUS, pose)
    pocket_resnums = np.setdiff1d(
        list(complex_resnums), peptide_resnums)
    pocket_size = len(pocket_resnums)
    return peptide_resnums, pocket_resnums, complex_resnums

def make_peptide_res_ids(
    pdb: str,
    chain: str='C',
    icodes: Optional[List]=None
) -> Dict:
    """Look up res ids for the peptide of a given pdb"""
    L = len(pdb_to_pep[pdb])
    return {'peptide':[(chain, i, ' ') for i in range(1,L+1)]}


def get_mhc_pocket_res_ids(
    pdb: str,
    chain: str='A',
    icodes: Optional[List]=None
) -> Dict:
    """
    Look up the res id for the mhc part of the pocket for a
    given pdb
    """
    try:
        mhc = pdb_to_mhc[pdb]
        pockets = mhc_pockets[mhc]
        pocket_res_ids = {}
        for pocket,seqnums in pockets.items():
            pocket_res_ids[f'mhc_{pocket}'] = [(chain, seqnum, ' ') for seqnum in seqnums]
        return pocket_res_ids
    except:
        logging.warning(f"MHC pockets not yet established for {pdb}")
        return {}

def get_tcr_pocket_res_ids(
    pdb: str,
    icodes: Optional[List]=None
) -> Dict:
    """
    Look up the res ids for the tcr part of the pocket for a 
    given pdb
    """
    try:
        template_pdb = pdb_to_tcr_pocket_template[pdb]
        pockets = tcr_pockets[template_pdb]
        pocket_res_ids = {}
        for pocket,seqnums in pockets.items():
            if 'alpha' in pocket:
                chain = pdb_to_tcr_chains[pdb]['alpha']
            else:
                chain = pdb_to_tcr_chains[pdb]['beta']
            pocket_res_ids[f'tcr_{pocket}'] = [(chain, seqnum, ' ') for seqnum in seqnums]
        return pocket_res_ids
    except:
        logging.warning(f"TCR pockets not yet established for {pdb}")
        return {}

def get_pockets(pdb: str, pocket_resnums: List, peptide_resnums: List, pose: Pose, mhc_chain: str = 'A'):
    """Get subpockets for tcr pmhc complexes"""
    regions = {}
    # regions.update(make_peptide_res_ids(pdb))
    regions.update({'peptide': [get_pdb_residue_info(pose, i) for i in peptide_resnums]})

    try:
        regions.update(get_mhc_pocket_res_ids(pdb, chain=mhc_chain))
        regions.update(get_tcr_pocket_res_ids(pdb))
    except Exception as e:
        logging.warning(f"Could not get pockets for {pdb}, printing error message below.")
        logging.warning(e)
    
    if pocket_resnums is not None:
        regions.update(
            {'pocket':[get_pdb_residue_info(pose,i) for i in pocket_resnums]})
    logging.debug(f"Regions for tracking:")
    #logging.debug(
    #json.dumps(pdb_to_tcr_chains, sort_keys=True, indent=2))
    return regions
    