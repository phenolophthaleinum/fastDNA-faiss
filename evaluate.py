import argparse
import json
from pathlib import Path

fh = open('X:/edwards2016/runs/run-dim_30-len_125-n_200-epoch_2-genus/rank/rank-k_60.json')
db = json.load(fh)
fh.close()

fh = open('virus.json')
vjson = json.load(fh)
fh.close()

fh = open('host.json')
hjson = json.load(fh)
fh.close()


TAX_RANKS = (
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "superkingdom",   
)


def compare_taxranks(vd, hd, rank):
    vi = vd['host']['lineage_ranks'].index(rank) if rank in vd['host']['lineage_ranks'] else None
    hi = hd['lineage_ranks'].index(rank) if rank in hd['lineage_ranks'] else None
    if vi != None and hi != None:
        if vd['host']['lineage_names'][vi] == hd['lineage_names'][hi]:
            return True
    return False




d = {}
for vid, hids in db.items():
    d[vid] = []
    for taxlevel in TAX_RANKS:
        match = False
        for tup in hids[:1]:
            hid = tup[0]
            match = compare_taxranks(vjson[vid], hjson[hid], taxlevel)
            if match:
                break
        d[vid].append(match)

l = [0 for _i in TAX_RANKS]
for vid in d:
    for i, taxlevel in enumerate(TAX_RANKS):
        if d[vid][i]:
            l[i] += 1


d = {}
for i, tax_rank in enumerate(TAX_RANKS):
    d[tax_rank] = l[i] / len(vjson) * 100

print(d)