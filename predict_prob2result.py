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

scoring_functions = {
    "max": scoring_func.choose_max,
    "poly": scoring_func.polynomial,
    "avg": scoring_func.avg,
    "harmonic": scoring_func.harmonic,
    "rank_row": scoring_func.rank_adjusted_row,
    "rank_column": scoring_func.rank_adjusted_col
}

hostname = utils.get_hostname_data()


# global_rank_tax = {}

def reverse_non_unique_mapping(d):
    dinv = {}
    for k, v in d.items():
        if v in dinv:
            dinv[v].append(k)
        else:
            dinv[v] = [k]
    return dinv


# @delayed
# @wrap_non_picklable_objects
def make_result_divison(file: str, func: str) -> Tuple[dict, dict]:
    with open(file) as jf:
        d = json.load(jf)
    data = []
    for elem in d:
        data.append(dict(sorted(elem.items(), key=lambda x: x[1], reverse=False)))
    df = pd.DataFrame.from_records(data).fillna(0)
    # df_means = df.apply(lambda x: x.mean())

    # # df ranking: max value from column takes 1 and so on
    # df_rank = df.rank(method='max', ascending=False)
    #
    # # df division: original/rank
    # df_div = df.div(df_rank)
    # #print(df_div)

    ######
    # with scoring function selection
    df_func = scoring_functions[func](df)
    ######

    # df_means.sort_values(ascending=False, inplace=True)
    df_func.sort_values(ascending=False, inplace=True)
    # data_json = df_means.to_json()

    # tax (order) extraction
    df_tax = df_func.rename(lambda x: hostname[x]['lineage_names'][3])
    # n best tax (order)
    df_tax = df_tax.nlargest(10)

    data_json = df_func.to_json()
    data_tax_json = df_tax.to_json()
    data_parsed = json.loads(data_json)
    data_tax_parsed = json.loads(data_tax_json)
    ranks = data_parsed.items()
    ranks_tax = data_tax_parsed.items()
    p = Path(file)
    part = {"_".join(str(p.stem).split("_")[:2]): list(ranks)}
    part_tax = {"_".join(str(p.stem).split("_")[:2]): [value[0] for value in ranks_tax]}
    # inversed_part_tax = {value[0]: key for key, list_value in list(part_tax.items()) for value in list_value}
    # print(list(part_tax.items()))
    # print(part_tax)
    # print(inversed_part_tax)
    # global_rank_tax.update(inversed_part_tax)
    return part, part_tax


def make_result(file: str, func: str) -> dict:
    with open(file) as jf:
        d = json.load(jf)
    df = pd.DataFrame.from_dict(d).fillna(0)
    # old way - unnecessary for loop and sorting
    # data = []
    # for elem in d:
    #     data.append(dict(sorted(elem.items(), key=lambda x: x[1], reverse=False)))
    # df = pd.DataFrame.from_records(data).fillna(0)
    ##############


    # df_means = df.apply(lambda x: x.mean())

    # # df ranking: max value from column takes 1 and so on
    # df_rank = df.rank(method='max', ascending=False)
    #
    # # df division: original/rank
    # df_div = df.div(df_rank)
    # #print(df_div)

    ######
    # with scoring function selection
    df_func = scoring_functions[func](df)
    ######

    # df_means.sort_values(ascending=False, inplace=True)
    df_func.sort_values(ascending=False, inplace=True)
    # data_json = df_means.to_json()

    data_json = df_func.to_json()
    data_parsed = json.loads(data_json)
    ranks = data_parsed.items()
    p = Path(file)
    part = {"_".join(str(p.stem).split("_")[:2]): list(ranks)}
    return part


def batch_exec_division(files: Tuple[str], scoring_function: str) -> Tuple[dict, dict]:
    local_rank = {}
    local_rank_tax = {}
    for file in files:
        # local_rank.update(make_result(file, scoring_function))
        res, res_tax = make_result_divison(file, scoring_function)
        local_rank.update(res)
        local_rank_tax.update(res_tax)
    return local_rank, local_rank_tax


def batch_exec(files: Tuple[str], scoring_function: str) -> dict:
    local_rank = {}
    for file in files:
        local_rank.update(make_result(file, scoring_function))
    return local_rank


def partition_queries(files: Iterable[str], partitions: int = cpu_count() - 1) -> np.ndarray:
    partitions = np.array_split(np.array(files), partitions)
    # debug for partitions
    # print(partitions)
    # print([ar.shape for ar in partitions])
    return partitions


def dump_result(output: str, result_rank: dict, scoring_function: str):
    prename = Path(output)
    rank_file = Path(f"{str(prename.with_suffix(''))}_{scoring_function}.json")
    # rank_file = Path(rank_str)
    if rank_file.exists():
        os.system(f"rm {str(rank_file)}")
    with rank_file.open("w") as fd:
        json.dump(result_rank, fd, indent=4)


def run_procedure_division(input_dir: str, scoring_function: str):
    # colorama
    init()

    # timer
    start = timer()
    # set_loky_pickler("pickle")
    # temporarily changed to loky; add later: prefer="threads"
    # ranks = Parallel(verbose=True, n_jobs=-1)(
    #     delayed(do_search)(file, dim, n_samples, k_nearest, index, map_data) for file in glob.glob(f"{input_dir}*.vec"))
    out = Parallel(verbose=True, n_jobs=-1)(delayed(batch_exec_division)(batch, scoring_function) for batch in
                                            partition_queries(glob.glob(f"{input_dir}*.json")))
    ranks, ranks_tax = zip(*out)
    # print(ranks)
    global_rank = {}
    global_rank_tax = {}
    # global global_rank_tax
    for result in ranks:
        global_rank.update(result)
    for result in ranks_tax:
        global_rank_tax.update(result)
        # global_rank_tax.update(result_tax)
    # print(global_rank_tax)
    # prename = Path(output)
    # rank_file = Path(f"{str(prename.with_suffix(''))}_{scoring_function}.json")
    # #rank_file = Path(rank_str)
    # if rank_file.exists():
    #     os.system(f"rm {str(rank_file)}")
    # with rank_file.open("w") as fd:
    #     json.dump(global_rank, fd, indent=4)

    end = timer()
    runtime = end - start
    print(f"{Fore.GREEN} Done in {runtime:.6f} seconds")
    return global_rank, global_rank_tax


def run_procedure(input_dir: str, scoring_function: str):
    # colorama
    init()

    # timer
    start = timer()
    # set_loky_pickler("pickle")
    # temporarily changed to loky; add later: prefer="threads"
    # ranks = Parallel(verbose=True, n_jobs=-1)(
    #     delayed(do_search)(file, dim, n_samples, k_nearest, index, map_data) for file in glob.glob(f"{input_dir}*.vec"))
    ranks = Parallel(verbose=True, n_jobs=-1)(delayed(batch_exec)(batch, scoring_function) for batch in
                                              partition_queries(glob.glob(f"{input_dir}*.json")))
    # print(ranks)
    global_rank = {}
    # global global_rank_tax
    for result in ranks:
        global_rank.update(result)
        # global_rank_tax.update(result_tax)
    # print(global_rank_tax)
    # prename = Path(output)
    # rank_file = Path(f"{str(prename.with_suffix(''))}_{scoring_function}.json")
    # #rank_file = Path(rank_str)
    # if rank_file.exists():
    #     os.system(f"rm {str(rank_file)}")
    # with rank_file.open("w") as fd:
    #     json.dump(global_rank, fd, indent=4)

    end = timer()
    runtime = end - start
    print(f"{Fore.GREEN} Done in {runtime:.6f} seconds")
    return global_rank

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="predict-prob results 2 one big result")
#     parser.add_argument("-i", "--input_dir", required=True,
#                         help="Directory with virus samples of vectors.")
#     parser.add_argument("-o", "--output", required=True,
#                         help="Path to result file with rankings.")
#     parser.add_argument("-s", '--scoring', required=True,
#                         choices=["max", "powmax", "avg", "harmonic"],
#                         help="Scoring function to be applied to the dataset.")
#
#     args = parser.parse_args()
#
#     global_rank = {}
#     run_procedure(args.input_dir, args.scoring)
