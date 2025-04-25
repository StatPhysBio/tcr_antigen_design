"""Lookup tables to peptide information"""

# Tables:
# - pdb_to_pep

pdb_to_pep = {
    '6am5':'SMLGIGIVPV',
    '6amu':'MMWDRGLGMM',
    '3qdg':'ELAGIGILTV',
    '3qdm':'ELAGIGILTV',
    '3hg1':'ELAGIGILTV',
    '2bnr':'SLLMWITQC',
    '2bnq':'SLLMWITQV',
    '5c0b':'RQFGPDFPTI',
    '5c0c':'RQFGPDWIVA',
    '3uts':'ALWGPDPAAA',
    '7n2n':'TRLALIAPK',
    '7n2o':'LRVMMLAPF',
    '7n2p':'GQVMVVAPR',
    '7n2q':'LRVMMLAPF',
    '7n2r':'TRLALIAPK',
    '7n2s':'TRLALIAPK',
    '8cx4':'LRVMMLAPF',
    '1ao7':'LLFGYPVYV',
    '1qrn':'LLFGYAVYV',
    '1qsf':'LLFGYPVAV',
    '1qse':'LLFGYPRYV',

    '5d2n':'NLVPMVATV',
    '5brz.pdb.human.MH1.A-01.A.C.DE': 'EVDPIGHLY',
    '5bs0.pdb.human.MH1.A-01.A.C.DE': 'ESDPIVAQY',
    '3mv7.pdb.human.MH1.B-35.A.C.DE': 'HPVGEADYFEY',
    '4prp.pdb.human.MH1.B-35.A.C.DE': 'HPVGQADYFEY'
}
pep_to_pdb = {v:k for k,v in pdb_to_pep.items()}

pdb_to_qualitative_pep = {
    '6am5':'MART',
    '6amu':'MART',
    '3qdg':'MART',
    '3qdm':'MART',
    '3hg1':'MART',
    '2bnr':'NY-ESO-1',
    '2bnq':'NY-ESO-1',
    '5c0b':'RQFGPDFPTI',
    '5c0c':'RQFGPDWIVA',
    '3uts':'ALWGPDPAAA',
    '7n2n':'TRLALIAPK',
    '7n2o':'LRVMMLAPF',
    '7n2p':'GQVMVVAPR',
    '7n2q':'LRVMMLAPF',
    '7n2r':'TRLALIAPK',
    '7n2s':'TRLALIAPK',
    '8cx4':'LRVMMLAPF',
    '1ao7':'Tax',
    '1qrn':'Tax',
    '1qsf':'Tax',

    '5d2n': 'NLV',
    '5brz.pdb.human.MH1.A-01.A.C.DE': 'MAGEA3',
    '5bs0.pdb.human.MH1.A-01.A.C.DE': 'Titin'
}