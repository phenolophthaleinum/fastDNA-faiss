import glob
import os
from collections import defaultdict
from typing import Dict, Tuple, List, Iterable, TypedDict
import math
from pathlib import Path
import pandas as pd
import utils
import numpy as np
import faiss
import json
import argparse
from joblib import Parallel, delayed
from joblib.externals.loky import set_loky_pickler
from joblib import parallel_backend
from joblib import Parallel, delayed
from joblib import wrap_non_picklable_objects
from timeit import default_timer as timer
from colorama import Fore, init
from multiprocessing import cpu_count
import scoring_func
import re

scoring = tuple([f for f in scoring_func.__dict__.values() if callable(f)])


def make_result(file: str, func: callable) -> dict:
    with open(file) as jf:
        d = json.load(jf)

    ## important; with json fixing
    try:
        with open(file) as jf:
            d = json.load(jf)
    except json.JSONDecodeError:
        print(f"Fixing bad json: {file}")
        # try:
        with open(file) as jf:
            dfix = jf.read()
            dfix = dfix.replace('\n', '')
            dfix = re.sub(r',]', ']', dfix)
            d = json.loads(dfix)
            print(f"Succesfully fixed {file}")
    #     except:
    #         print("Couldn't fix json. Skipping...")
    #         return {}

    df = pd.DataFrame.from_dict(d).fillna(1e-6)  # https://doi.org/10.1371/journal.pbio.3000106 "we estimate that there exist globally between 0.8 and 1.6 million prokaryotic OTUs"

    ######
    # with scoring function selection
    df_func = func(df)
    ######

    try:
        df_func.sort_values(ascending=False, inplace=True)
    except:
        raise ValueError(df_func)
        # data_json = df_means.to_json()
    data_json = df_func.to_json()
    data_parsed = json.loads(data_json)
    ranks = data_parsed.items()
    p = Path(file)
    part = {"_".join(str(p.stem).split("_")[:2]): list(ranks)}
    return part


def batch_exec(files: Tuple[str], scoring_function: str) -> dict:
    local_rank = {}
    for file in files:
        local_rank.update(make_result(file, scoring_function))
    return local_rank


# always uses max available threads - this should be configurable in the future
def partition_queries(files: Iterable[str], partitions: int = cpu_count() - 1) -> np.ndarray:
    partitions = np.array_split(np.array(files), partitions)
    # debug for partitions
    # print(partitions)
    # print([ar.shape for ar in partitions])
    return partitions


def dump_result(output: Path, result_rank: dict):
    # rank_file = Path(rank_str)
    if output.exists():
        os.system(f"rm {str(output)}")
    with output.open("w") as fd:
        json.dump(result_rank, fd, indent=4)


def score_and_rank(input_dir: str, scoring_function: str):
    # colorama
    init()

    # timer
    start = timer()
    ranks = Parallel(verbose=True, n_jobs=-1)(delayed(batch_exec)(batch, scoring_function) for batch in
                                              partition_queries(glob.glob(f"{input_dir}*.json")))
    global_rank = {}
    for result in ranks:
        global_rank.update(result)

    end = timer()
    runtime = end - start
    print(f"{Fore.GREEN} Done in {runtime:.6f} seconds")
    return global_rank


def compare_taxranks(vd, hd, rank):
    vi = vd['host']['lineage_ranks'].index(rank) if rank in vd['host']['lineage_ranks'] else None
    hi = hd['lineage_ranks'].index(rank) if rank in hd['lineage_ranks'] else None
    if vi != None and hi != None:
        if vd['host']['lineage_names'][vi] == hd['lineage_names'][hi]:
            return True
    return False


def evaluate_ranks(result_rank: dict,
                   viruses: dict,
                   hosts: dict,
                   hostname: dict,
                   rank: str = None):
    TAX_RANKS = ("species",
                 "genus",
                 "family",
                 "order",
                 "class",
                 "phylum",
                 "superkingdom")

    d = {}
    d2 = {}
    all_vid = len(list(result_rank.keys()))

    for vid, hids in result_rank.items():
        d[vid] = []
        d2[vid] = []
        for taxlevel in TAX_RANKS:
            match = False
            for tup in hids[:1]:
                hid = tup[0]
                hid = hostname[hid]['ncbi_id']
                match = compare_taxranks(viruses[vid], hosts[hid], taxlevel)
                if match:
                    break
            d[vid].append(match)
            tax_index = hosts[hid]['lineage_ranks'].index(taxlevel)
            d2[vid].append((hosts[hid]['lineage_names'][tax_index], viruses[vid]["host"]["lineage_names"][tax_index], match))

    l = [0 for _i in TAX_RANKS]

    for vid in d:
        for i, taxlevel in enumerate(TAX_RANKS):
            if d[vid][i]:
                l[i] += 1

    d = {}
    for i, tax_rank in enumerate(TAX_RANKS):
        d[tax_rank] = l[i] / all_vid * 100

    if rank:
        return d[rank]

    else:
        return d