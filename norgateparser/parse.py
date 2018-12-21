import functools
import json
import os
import re

import jdna
import pandas as pd
from colorhash import ColorHash

here = os.path.abspath(os.path.dirname(__file__))
data_dir = os.path.join(here, 'data')


def parts_df():
    PARTS_FILEPATH = os.path.join(data_dir, 'NORGateParts.csv')

    parts_df = pd.read_csv(PARTS_FILEPATH)
    parts_df.fillna(value='None', inplace=True)
    parts_df.tail()
    return parts_df


def parts_by_name_and_cat():
    parts_dict = {}
    by_category = {}
    df = parts_df()
    for _, row in df.iterrows():
        seq = row['Sequence']
        if seq != 'None':
            dna = jdna.Sequence(row['Sequence'])
            dna.annotate(None, None, row['Part'], color=ColorHash(row['Part'] + row['Annotation']).hex)
            parts_dict[row['Part']] = dna
            by_category.setdefault(row['Annotation'], {})
            by_category[row['Annotation']][row['Part']] = dna
    return parts_dict, by_category


def wsetdict():
    path = os.path.join(data_dir, 'wset_conversion.json')
    with open(path, 'r') as f:
        wset = json.load(f)
    return wset


def concat_parts(parts):
    """Concatenates a list of Coral.DNA objects or looks for parts by name in the parts dictionary"""
    parts_dict, _ = parts_by_name_and_cat()
    s = jdna.Sequence()
    for p in parts:
        if p is None:
            continue
        if isinstance(p, str):
            p = p.strip()
            if p == '':
                continue
            s += parts_dict[p].copy()
        else:
            s += p
    return s


def generate_gRNA_cassette(prom=None,
                           promlink=None,
                           target=None,
                           handle=None,
                           ribo5=None,
                           ribo3=None,
                           ribo5INS=(None, None),
                           termlink=None,
                           term=None):
    """Generates gRNA cassette DNA sequence using Coral"""
    parts = [prom, promlink, ribo5INS[0], ribo5, ribo5INS[1], target, handle, ribo3, termlink, term]
    seq = concat_parts(parts)
    return seq


# helper function for uninsulated RGR cassette
helper_rgr = functools.partial(
    generate_gRNA_cassette,
    handle='spCas9 gRNA Handle',
    ribo3='HDV Ribozyme',
    ribo5='HH Ribozyme',
    termlink='TP',
    term='tCYC1'
)


# uninsulated RGR cassette
def rgr(target):
    ins = target.copy()[:6].reverse_complement()
    return helper_rgr(target=target.copy(), ribo5INS=(ins, None))


# insulated pADH1 gRNA cassette
padh1_irgr = functools.partial(
    generate_gRNA_cassette,
    handle='spCas9 gRNA Handle',
    ribo3='HDV Ribozyme',
    ribo5='ASBV1 Ribozyme',
    ribo5INS=('ASBV1 pADH1_INS1', 'ASBV1 pADH1_INS2'),
    termlink='TP',
    term='tCYC1'
)

# insulated pGAZL4 gRNA cassette
pgalz4_irgr = functools.partial(
    generate_gRNA_cassette,
    handle='spCas9 gRNA Handle',
    ribo3='HDV Ribozyme',
    ribo5='ASBV1 Ribozyme',
    ribo5INS=('ASBV1 pGAL1_INS1', 'ASBV1 pGAL1_INS2'),
    termlink='TP',
    term='tCYC1'
)

# insulated pGRR gRNA cassette
pGRR_irgr = functools.partial(
    generate_gRNA_cassette,
    handle='spCas9 gRNA Handle',
    ribo3='HDV Ribozyme',
    ribo5='ASBV1 Ribozyme',
    ribo5INS=('ASBV1 pGRR_INS1', 'ASBV1 pGRR_INS2'),
    termlink='TP',
    term='tCYC1'
)


def pGRR(i, j):
    """Generates a new pGRR promoter from inputs 'i' and 'j'"""

    parts_dict, _ = parts_by_name_and_cat()
    i = parts_dict.get(i, None)
    j = parts_dict.get(j, None)
    if i is not None:
        i = i.copy()
    if j is not None:
        j = j.copy()
        j = j.reverse_complement()

    parts = ['pGRR', i, 'pGRR TATA', j, 'pGRR RBS']
    return concat_parts(parts)


def pMOD(homology, cassette, marker=None):
    """Generates new pMOD vector with a new cassette. Cassette should include
    the promoter, gene, and terminator (no assembly linkers)."""

    homology_dict = {
        "URA": [["URA3 Promoter", "URA3", "tADH1"], ["URA3 3'UTR"]],
        "HIS": [["HIS3 Promoter", "HIS3", "tADH1"], ["HIS3 3'UTR"]],
        "TRP": [["TRP1 Promoter", "TRP1", "tADH1"], ["TRP1 3'UTR"]],
        "LTR1": [["LTR1 homology 1", "NATMX"], ["LTR1 homology 2"]],
        "LTR2": [["LTR2 homology 1", "BLEOMX"], ["LTR2 homology 2"]],
        "LTR3": [["LTR3 homology 1", "HYBMX"], ["LTR3 homology 2"]],
        "HO": [["HO homology 1", "KANMX"], ["HO homology 2"]]
    }

    hom1, hom2 = homology_dict[homology]
    if marker is not None:
        hom1[1] = marker

    _parts = ["PmeI", hom1, 'PP2', cassette, "TS", hom2, "PmeI", "AMPR and ORI"]
    _parts = _parts[1:-2]
    parts = []
    for p in _parts:
        if isinstance(p, list):
            parts += p
        else:
            parts.append(p)
    return concat_parts(parts)


def parse_NOR_gate(name):
    # {'marker': '8', 'prom': 'A', 'i': None, 'j': None, 'grna': 'URGR', 'k': 'W10'}
    imatch = "(?P<i>[WF]\d+)"
    jmatch = "(?P<j>[WF]\d+)"

    pgrr_match = '(pGRR(-{i})?(-?{j})?)'.format(i=imatch, j=jmatch)
    promoter_match = 'p\w+|{pgrr_match}|[GA]'.format(pgrr_match=pgrr_match)

    marker_match = '(?P<marker_int>\d)|(?P<homology>ho|HO|\w+)-?(?P<marker>kan|KAN|Kan|\w+)'

    cassette_match = "(?P<grna>[iU]?RGR)-(?P<k>[WF]\d+)|[y]?eGFP"
    pattern = "pMOD-?({marker})-?(?P<prom>{promoter})-(?P<cassette>{cassette})".format(cassette=cassette_match,
                                                                                       marker=marker_match,
                                                                                       promoter=promoter_match,)
    m = re.match(pattern, name, re.IGNORECASE)
    return m


def parse_nor_gate_name_to_sequence(name):
    parts_dict, _ = parts_by_name_and_cat()
    match = parse_NOR_gate(name)
    if match is None:
        print("No match with '{}'".format(name))
        return None

    parsed = match.groupdict()

    homology_dict = {6: "URA", 8: "HIS", 4: "TRP"}
    promoter_dict = {'A': 'pADH1', 'G': 'pGPD'}
    irgr_dict = {
        'pADH1': padh1_irgr,
        'pGALZ4': pgalz4_irgr,
        'pGRR': pGRR_irgr,
    }

    get_target = lambda x: wsetdict().get(x, x)

    i = get_target(parsed['i'])
    j = get_target(parsed['j'])
    k = get_target(parsed['k'])

    if i is None:
        i = 'c3'
    if j is None:
        j = 'c6'

    marker_int = parsed['marker_int']
    if marker_int is not None:
        marker_int = int(marker_int)

    homology = homology_dict.get(marker_int, parsed['homology'])
    try:
        homology = homology_dict[int(homology)]
    except ValueError:
        pass
    marker = parsed['marker']
    prom = promoter_dict.get(parsed['prom'], parsed['prom'])

    grna = None
    if parsed['grna'] is not None:
        if parsed['grna'] == 'RGR':
            grna = rgr(parts_dict[k])
        elif parsed['grna'].lower() in ['irgr', 'urgr']:
            if 'pGRR'.upper() in prom.upper():
                grna = irgr_dict['pGRR'](k)
            else:
                grna = irgr_dict[prom](k)

    if 'pGRR' in prom:
        prom = pGRR(i, j)

    cassette = parsed['cassette']
    if grna is not None:
        cassette = concat_parts([
            prom,
            'PS',
            grna
        ])
    else:
        if 'egfp' in cassette.lower():
            cassette = 'yeGFP'
        cassette = [
            prom, 'PS', cassette, 'TP', 'tCYC1'
        ]

    if marker is not None:
        marker = marker.upper()
        if 'hyg' in marker.lower():
            marker = 'HYGMX'
        if 'nat' in marker.lower():
            marker = 'NATMX'
        if 'kan' in marker.lower():
            marker = 'KANMX'
        if 'bleo' in marker.lower():
            marker = 'BLEOMX'
    return pMOD(homology.upper(), cassette, marker=marker)

# alias for main function
parse_name = parse_nor_gate_name_to_sequence
