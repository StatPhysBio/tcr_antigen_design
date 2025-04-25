"""Lookup tables for tcr information"""

# Tables:
# - pdb_to_tcr
# - tcr_pockets
# - tcr_pocket_colors
# - pdb_to_tcr_pocket_template
# - pdb_to_tcr_chains

pdb_to_tcr = {
    '1ao7':'A6',
    '1qrn':'A6',
    '1qsf':'A6',
    '1qse':'A6',
    '6am5':'DMF5',
    '6amu':'DMF5',
    '3qdg':'DMF5',
    '3qdm':'DMF4',
    '3hg1':'Mel5',
    '2bnr':'1G4',
    '2bnq':'1G4',
    '5c0c':'1E6',
    '5c0b':'1E6',
    '5c07':'1E6',
    '5c08':'1E6',
    '5c09':'1E6',
    '3uts':'1E6',
    '7n2n':'AS4.2',
    '7n2o':'AS4.2',
    '7n2p':'AS4.3',
    '7n2q':'AS4.3',
    '7n2r':'AS4.3',
    '7n2s':'AS3.1',
    '8cx4':'AS8.4',
    '5d2n':'I-dont-know',
}

# ad hoc pockets made by eye based on distance and cdr3 definitions
# TODO: Implement intersection of overall pocket with cdr regions
tcr_pockets = {
    '7n2n': {
        'beta cdr3': [96,97,98,99,100],
        'alpha cdr1': [30,31,32],
        'alpha cdr2': [52,53,54,55],
        'alpha cdr3': [95,96,97,98,99,100,101,102,103]
    },
    '6am5': {
        'beta cdr3': [94,95,96,97,98,99,100],
        'alpha cdr1': [24,25,26,27,28,29,30,31,32],
        'alpha cdr2': [50,51,52,53,54,55],
        'alpha cdr3': [89,90,91,92,93,94,95,96,97],
    },
    '3qdm': {
        'beta cdr3': [92,93,94,95,96,97,98,99,100,101,102,103,104],
        'alpha cdr1': [24,25,26,27,28,29,30],
        'alpha cdr2': [49,50,51,52,53,54,55],
        'alpha cdr3': [89,90,91,92,93,94,95,96,97,98],
    },
    '5d2n': {'TODO'}
}

tcr_pocket_colors = {
    'alpha cdr3': 'firebrick',
    'alpha cdr2': 'gold',
    'alpha cdr1':'forestgreen',
    'beta cdr3': 'steelblue',
}

pdb_to_tcr_pocket_template = {
    '6am5':'3qdm', # to edit
    '6amu':'3qdm', # to edit
    '3qdg':'3qdm', # to edit
    '3qdm':'3qdm', 
    '3hg1':'3qdm', # to edit
    '5c0c':'3qdm', # to edit
    '3uts':'3qdm', # to edit
    '5c0b':'3qdm', # to edit
    '7n2n':'7n2n',
    '7n2o':'7n2n',
    '7n2p':'7n2n',
    '7n2q':'7n2n',
    '7n2r':'7n2n',
    '7n2s':'7n2n',
    '8cx4':'7n2n',

    '5d2n':'5d2n',
}    #TODO: Ultimately there should be no templates

# long data structure. Next variable for jumping ahead is N/A
pdb_to_tcr_chains = {
    "1ao7": {
        "alpha": "D",
        "beta": "E"
    },
    "1bd2": {
        "alpha": "D",
        "beta": "E"
    },
    "1d9k": {
        "alpha": "A",
        "beta": "B"
    },
    "1fo0": {
        "alpha": "A",
        "beta": "B"
    },
    "1fyt": {
        "alpha": "D",
        "beta": "E"
    },
    "1g6r": {
        "alpha": "A",
        "beta": "B"
    },
    "1j8h": {
        "alpha": "D",
        "beta": "E"
    },
    "1kj2": {
        "alpha": "A",
        "beta": "B"
    },
    "1lp9": {
        "alpha": "E",
        "beta": "F"
    },
    "1mi5": {
        "alpha": "D",
        "beta": "E"
    },
    "1mwa": {
        "alpha": "A",
        "beta": "B"
    },
    "1nam": {
        "alpha": "A",
        "beta": "B"
    },
    "1oga": {
        "alpha": "D",
        "beta": "E"
    },
    "1qrn": {
        "alpha": "D",
        "beta": "E"
    },
    "1qse": {
        "alpha": "D",
        "beta": "E"
    },
    "1qsf": {
        "alpha": "D",
        "beta": "E"
    },
    "1u3h": {
        "alpha": "A",
        "beta": "B"
    },
    "1ymm": {
        "alpha": "D",
        "beta": "E"
    },
    "1zgl": {
        "alpha": "M",
        "beta": "P"
    },
    "2ak4": {
        "alpha": "N",
        "beta": "P"
    },
    "2bnq": {
        "alpha": "D",
        "beta": "E"
    },
    "2bnr": {
        "alpha": "D",
        "beta": "E"
    },
    "2ckb": {
        "alpha": "A",
        "beta": "B"
    },
    "2e7l": {
        "alpha": "A",
        "beta": "D"
    },
    "2esv": {
        "alpha": "D",
        "beta": "E"
    },
    "2f53": {
        "alpha": "D",
        "beta": "E"
    },
    "2f54": {
        "alpha": "D",
        "beta": "E"
    },
    "2gj6": {
        "alpha": "D",
        "beta": "E"
    },
    "2iam": {
        "alpha": "C",
        "beta": "D"
    },
    "2ian": {
        "alpha": "D",
        "beta": "E"
    },
    "2j8u": {
        "alpha": "E",
        "beta": "F"
    },
    "2jcc": {
        "alpha": "E",
        "beta": "F"
    },
    "2oi9": {
        "alpha": "B",
        "beta": "C"
    },
    "2ol3": {
        "alpha": "A",
        "beta": "B"
    },
    "2p5e": {
        "alpha": "D",
        "beta": "E"
    },
    "2p5w": {
        "alpha": "D",
        "beta": "E"
    },
    "2pxy": {
        "alpha": "A",
        "beta": "B"
    },
    "2pye": {
        "alpha": "D",
        "beta": "E"
    },
    "2uwe": {
        "alpha": "E",
        "beta": "F"
    },
    "2vlj": {
        "alpha": "D",
        "beta": "E"
    },
    "2vlk": {
        "alpha": "D",
        "beta": "E"
    },
    "2vlr": {
        "alpha": "D",
        "beta": "E"
    },
    "2wbj": {
        "alpha": "C",
        "beta": "D"
    },
    "2ypl": {
        "alpha": "D",
        "beta": "E"
    },
    "2z31": {
        "alpha": "A",
        "beta": "B"
    },
    "3c5z": {
        "alpha": "A",
        "beta": "B"
    },
    "3c60": {
        "alpha": "A",
        "beta": "B"
    },
    "3c6l": {
        "alpha": "A",
        "beta": "B"
    },
    "3d39": {
        "alpha": "D",
        "beta": "E"
    },
    "3d3v": {
        "alpha": "D",
        "beta": "E"
    },
    "3dxa": {
        "alpha": "D",
        "beta": "E"
    },
    "3e2h": {
        "alpha": "B",
        "beta": "C"
    },
    "3e3q": {
        "alpha": "D",
        "beta": "E"
    },
    "3ffc": {
        "alpha": "D",
        "beta": "E"
    },
    "3gsn": {
        "alpha": "A",
        "beta": "B"
    },
    "3h9s": {
        "alpha": "D",
        "beta": "E"
    },
    "3hg1": {
        "alpha": "D",
        "beta": "E"
    },
    "3kpr": {
        "alpha": "D",
        "beta": "E"
    },
    "3kps": {
        "alpha": "D",
        "beta": "E"
    },
    "3kxf": {
        "alpha": "M",
        "beta": "O"
    },
    "3mbe": {
        "alpha": "C",
        "beta": "D"
    },
    "3mv7": {
        "alpha": "D",
        "beta": "E"
    },
    "3mv8": {
        "alpha": "D",
        "beta": "E"
    },
    "3mv9": {
        "alpha": "D",
        "beta": "E"
    },
    "3o4l": {
        "alpha": "D",
        "beta": "E"
    },
    "3o6f": {
        "alpha": "C",
        "beta": "D"
    },
    "3pl6": {
        "alpha": "C",
        "beta": "D"
    },
    "3pqy": {
        "alpha": "D",
        "beta": "E"
    },
    "3pwp": {
        "alpha": "D",
        "beta": "E"
    },
    "3qdg": {
        "alpha": "D",
        "beta": "E"
    },
    "3qdj": {
        "alpha": "D",
        "beta": "E"
    },
    "3qdm": {
        "alpha": "D",
        "beta": "E"
    },
    "3qeq": {
        "alpha": "D",
        "beta": "E"
    },
    "3qfj": {
        "alpha": "D",
        "beta": "E"
    },
    "3qib": {
        "alpha": "C",
        "beta": "D"
    },
    "3qiu": {
        "alpha": "C",
        "beta": "D"
    },
    "3qiw": {
        "alpha": "C",
        "beta": "D"
    },
    "3rdt": {
        "alpha": "A",
        "beta": "B"
    },
    "3rgv": {
        "alpha": "A",
        "beta": "B"
    },
    "3sjv": {
        "alpha": "D",
        "beta": "E"
    },
    "3t0e": {
        "alpha": "C",
        "beta": "D"
    },
    "3tf7": {
        "alpha": "C",
        "beta": "C"
    },
    "3tfk": {
        "alpha": "C",
        "beta": "D"
    },
    "3tjh": {
        "alpha": "C",
        "beta": "D"
    },
    "3tpu": {
        "alpha": "C",
        "beta": "D"
    },
    "3uts": {
        "alpha": "D",
        "beta": "E"
    },
    "3utt": {
        "alpha": "I",
        "beta": "J"
    },
    "3vxm": {
        "alpha": "D",
        "beta": "E"
    },
    "3vxr": {
        "alpha": "D",
        "beta": "E"
    },
    "3vxs": {
        "alpha": "D",
        "beta": "E"
    },
    "4e41": {
        "alpha": "D",
        "beta": "E"
    },
    "4eup": {
        "alpha": "I",
        "beta": "J"
    },
    "4ftv": {
        "alpha": "D",
        "beta": "E"
    },
    "4g8g": {
        "alpha": "D",
        "beta": "E"
    },
    "4g9f": {
        "alpha": "D",
        "beta": "E"
    },
    "4gg6": {
        "alpha": "G",
        "beta": "H"
    },
    "4grl": {
        "alpha": "C",
        "beta": "D"
    },
    "4h1l": {
        "alpha": "I",
        "beta": "J"
    },
    "4jrx": {
        "alpha": "D",
        "beta": "E"
    },
    "4jry": {
        "alpha": "D",
        "beta": "E"
    },
    "4l3e": {
        "alpha": "D",
        "beta": "E"
    },
    "4may": {
        "alpha": "C",
        "beta": "D"
    },
    "4mji": {
        "alpha": "D",
        "beta": "E"
    },
    "4mnq": {
        "alpha": "D",
        "beta": "E"
    },
    "4ms8": {
        "alpha": "C",
        "beta": "D"
    },
    "4mvb": {
        "alpha": "C",
        "beta": "D"
    },
    "4mxq": {
        "alpha": "C",
        "beta": "D"
    },
    "4n0c": {
        "alpha": "C",
        "beta": "D"
    },
    "4n5e": {
        "alpha": "C",
        "beta": "D"
    },
    "4nhu": {
        "alpha": "C",
        "beta": "D"
    },
    "4ozf": {
        "alpha": "G",
        "beta": "H"
    },
    "4ozg": {
        "alpha": "G",
        "beta": "H"
    },
    "4ozh": {
        "alpha": "G",
        "beta": "H"
    },
    "4ozi": {
        "alpha": "G",
        "beta": "H"
    },
    "4p23": {
        "alpha": "A",
        "beta": "B"
    },
    "4p2o": {
        "alpha": "C",
        "beta": "D"
    },
    "4p2q": {
        "alpha": "D",
        "beta": "E"
    },
    "4p2r": {
        "alpha": "D",
        "beta": "E"
    },
    "4p46": {
        "alpha": "A",
        "beta": "B"
    },
    "4p4k": {
        "alpha": "C",
        "beta": "D"
    },
    "4p5t": {
        "alpha": "A",
        "beta": "B"
    },
    "4prh": {
        "alpha": "D",
        "beta": "E"
    },
    "4pri": {
        "alpha": "D",
        "beta": "E"
    },
    "4prp": {
        "alpha": "D",
        "beta": "E"
    },
    "4qok": {
        "alpha": "D",
        "beta": "E"
    },
    "4qrp": {
        "alpha": "D",
        "beta": "E"
    },
    "4y19": {
        "alpha": "D",
        "beta": "E"
    },
    "4y1a": {
        "alpha": "D",
        "beta": "E"
    },
    "4z7u": {
        "alpha": "G",
        "beta": "H"
    },
    "4z7v": {
        "alpha": "G",
        "beta": "H"
    },
    "4z7w": {
        "alpha": "G",
        "beta": "H"
    },
    "5brz": {
        "alpha": "D",
        "beta": "E"
    },
    "5bs0": {
        "alpha": "D",
        "beta": "E"
    },
    "5c07": {
        "alpha": "D",
        "beta": "E"
    },
    "5c08": {
        "alpha": "D",
        "beta": "E"
    },
    "5c09": {
        "alpha": "D",
        "beta": "E"
    },
    "5c0a": {
        "alpha": "D",
        "beta": "E"
    },
    "5c0b": {
        "alpha": "D",
        "beta": "E"
    },
    "5c0c": {
        "alpha": "I",
        "beta": "J"
    },
    "5d2l": {
        "alpha": "I",
        "beta": "J"
    },
    "5d2n": {
        "alpha": "D",
        "beta": "E"
    },
    "5e6i": {
        "alpha": "A",
        "beta": "B"
    },
    "5e9d": {
        "alpha": "D",
        "beta": "E"
    },
    "5eu6": {
        "alpha": "D",
        "beta": "E"
    },
    "5euo": {
        "alpha": "E",
        "beta": "F"
    },
    "5hhm": {
        "alpha": "D",
        "beta": "E"
    },
    "5hho": {
        "alpha": "D",
        "beta": "E"
    },
    "5hyj": {
        "alpha": "D",
        "beta": "E"
    },
    "5isz": {
        "alpha": "D",
        "beta": "E"
    },
    "5ivx": {
        "alpha": "E",
        "beta": "F"
    },
    "5jhd": {
        "alpha": "D",
        "beta": "E"
    },
    "5jzi": {
        "alpha": "D",
        "beta": "E"
    },
    "5ks9": {
        "alpha": "G",
        "beta": "H"
    },
    "5ksa": {
        "alpha": "C",
        "beta": "D"
    },
    "5ksb": {
        "alpha": "G",
        "beta": "H"
    },
    "5m00": {
        "alpha": "G",
        "beta": "H"
    },
    "5m01": {
        "alpha": "G",
        "beta": "H"
    },
    "5m02": {
        "alpha": "G",
        "beta": "H"
    },
    "5men": {
        "alpha": "D",
        "beta": "E"
    },
    "5nht": {
        "alpha": "A",
        "beta": "B"
    },
    "5nme": {
        "alpha": "D",
        "beta": "E"
    },
    "5nmf": {
        "alpha": "D",
        "beta": "E"
    },
    "5nmg": {
        "alpha": "D",
        "beta": "E"
    },
    "5nqk": {
        "alpha": "A",
        "beta": "B"
    },
    "5sws": {
        "alpha": "D",
        "beta": "E"
    },
    "5swz": {
        "alpha": "D",
        "beta": "E"
    },
    "5tez": {
        "alpha": "I",
        "beta": "J"
    },
    "5til": {
        "alpha": "G",
        "beta": "H"
    },
    "5tje": {
        "alpha": "G",
        "beta": "H"
    },
    "5w1v": {
        "alpha": "D",
        "beta": "E"
    },
    "5w1w": {
        "alpha": "D",
        "beta": "E"
    },
    "5wkf": {
        "alpha": "D",
        "beta": "E"
    },
    "5wkh": {
        "alpha": "D",
        "beta": "E"
    },
    "5wlg": {
        "alpha": "D",
        "beta": "E"
    },
    "5yxn": {
        "alpha": "A",
        "beta": "B"
    },
    "5yxu": {
        "alpha": "A",
        "beta": "B"
    },
    "6am5": {
        "alpha": "D",
        "beta": "E"
    },
    "6amu": {
        "alpha": "D",
        "beta": "E"
    },
    "6avf": {
        "alpha": "A",
        "beta": "B"
    },
    "6avg": {
        "alpha": "B",
        "beta": "E"
    },
    "6bga": {
        "alpha": "C",
        "beta": "D"
    },
    "6bj2": {
        "alpha": "D",
        "beta": "E"
    },
    "6cql": {
        "alpha": "D",
        "beta": "E"
    },
    "6cqn": {
        "alpha": "D",
        "beta": "E"
    },
    "6cqq": {
        "alpha": "D",
        "beta": "E"
    },
    "6cqr": {
        "alpha": "D",
        "beta": "E"
    },
    "6d78": {
        "alpha": "D",
        "beta": "E"
    },
    "6dfs": {
        "alpha": "A",
        "beta": "B"
    },
    "6dfw": {
        "alpha": "G",
        "beta": "H"
    },
    "6dfx": {
        "alpha": "G",
        "beta": "H"
    },
    "6dkp": {
        "alpha": "D",
        "beta": "E"
    },
    "6eqa": {
        "alpha": "D",
        "beta": "E"
    },
    "6g9q": {
        "alpha": "G",
        "beta": "H"
    },
    "6mkd": {
        "alpha": "A",
        "beta": "B"
    },
    "6mkr": {
        "alpha": "A",
        "beta": "B"
    },
    "6mng": {
        "alpha": "A",
        "beta": "B"
    },
    "6mnm": {
        "alpha": "A",
        "beta": "B"
    },
    "6mnn": {
        "alpha": "A",
        "beta": "B"
    },
    "6mno": {
        "alpha": "A",
        "beta": "B"
    },
    "6mtm": {
        "alpha": "D",
        "beta": "E"
    },
    "6q3s": {
        "alpha": "D",
        "beta": "E"
    },
    "7n2n": {
        "alpha": "D",
        "beta": "F"
    },
    "7n2o": {
        "alpha": "D",
        "beta": "F"
    },
    "7n2p": {
        "alpha": "D",
        "beta": "F"
    },
    "7n2q": {
        "alpha": "D",
        "beta": "F"
    },
    "7n2r": {
        "alpha": "D",
        "beta": "F"
    },
    "7n2s": {
        "alpha": "D",
        "beta": "F"
    },
    "8cx4": {
        "alpha": "D",
        "beta": "F"
    },
    
}