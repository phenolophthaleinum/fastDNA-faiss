import glob
import os
from collections import defaultdict
from typing import Dict, Tuple, List, Iterable
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
    "powmax": None,
    "avg": scoring_func.avg,
    "harmonic": scoring_func.harmonic
}


# @delayed
# @wrap_non_picklable_objects
def make_result(file: str, func: str) -> Dict[str, List[Tuple[str, float]]]:

    with open(file) as jf:
        d = json.load(jf)
    data = []
    for elem in d:
        data.append(dict(sorted(elem.items(), key=lambda x: x[1], reverse=False)))
    df = pd.DataFrame.from_records(data).fillna(0)
    # df_means = df.apply(lambda x: x.mean())

    # df ranking: max value from column takes 1 and so on
    df_rank = df.rank(method='max', ascending=False)

    # df division: original/rank
    df_div = df.div(df_rank)
    #print(df_div)

    ######
    # with scoring function selection
    df_func = scoring_functions[func](df_div)
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


def batch_exec(files: Tuple[str], scoring_function: str) -> Dict[str, List[Tuple[str, float]]]:
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


def run_procedure(input_dir: str, output: str, scoring_function: str):
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

    for result in ranks:
        global_rank.update(result)

    prename = Path(output)
    rank_file = Path(f"{str(prename.with_suffix(''))}_{scoring_function}.json")
    #rank_file = Path(rank_str)
    if rank_file.exists():
        os.system(f"rm {str(rank_file)}")
    with rank_file.open("w") as fd:
        json.dump(global_rank, fd, indent=4)

    end = timer()
    runtime = end - start
    print(f"{Fore.GREEN} Done in {runtime:.6f} seconds")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="predict-prob results 2 one big result")
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Directory with virus samples of vectors.")
    parser.add_argument("-o", "--output", required=True,
                        help="Path to result file with rankings.")
    parser.add_argument("-s", '--scoring', required=True,
                        choices=["max", "powmax", "avg", "harmonic"],
                        help="Scoring function to be applied to the dataset.")

    args = parser.parse_args()

    global_rank = {}
    run_procedure(args.input_dir, args.output, args.scoring)
