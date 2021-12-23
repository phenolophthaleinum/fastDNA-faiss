import glob
import os
from collections import defaultdict
from typing import Dict, Tuple, List, Iterable
import math
from pathlib import Path

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


# @delayed
# @wrap_non_picklable_objects
def do_search(file: str, dim: int, n_samples: int, k_nearest: int, index: faiss.IndexFlatL2,
              map_data: Dict[int, str]) -> Dict[str, List[Tuple[str, float]]]:
    """
    Args:
        file (str): path to a binary file with query vectors (virus)
        dim (int): dimensionality of a vector
        n_samples (int): number of samples that query file contains (number of vectors)
        k_nearest (int): k-nearest vectors to be found
        index (IndexFlatL2): faiss index object to be searched
        map_data (dict[int, str]): dictionary describing index - genome_name relations

    Returns: dict: dictionary of viruses with assigned list of the nearest hosts ranked by summed up scores (sum of
    scores 0-1 for each matched query sample) in a descending order
    """
    query = np.fromfile(file, dtype="float32")
    query = query.reshape(n_samples, dim)
    distances, indices = index.search(query, k_nearest)
    classification = np.vectorize(map_data.get)(indices)
    coef = 1 / (2 * math.sqrt(dim))
    pre_rank = defaultdict(list)
    enum = 1
    dict_sample = defaultdict(list)
    # print(f"Distances row number (number of samples): {np.shape(distances)[0]}")
    # print(f"Distances col number (number of k-nearest):{np.shape(distances)[1]}")
    # print(f"Indices row number (number of samples){np.shape(indices)[0]}")
    # print(f"Indices col number (number of k-nearest){np.shape(indices)[1]}")
    # print(f"Max dist: {distances.max()}")
    # print(f"Min dist: {distances.min()}")
    # print(len(np.stack((classification.flatten(), distances.flatten()), axis=1)))

    for idx, distance, in np.stack((classification.flatten(), distances.flatten()), axis=1):
        score = float(distance) * coef
        if score < 0:
            print(
                f"Negative score: {score}; virus, host: {('_'.join(file.split('/')[-1].split('.')[0].split('_')[:2]), idx)}; distance: {distance}\n")
        dict_sample[idx].append(score)
        if enum % dim == 0:
            for host in dict_sample:
                pre_rank[idx].append(max(dict_sample[host]))
            # print(f"Dict_sample size: {len(dict_sample)}")
            dict_sample.clear()
        enum += 1

    sum_rank = {key: sum(values) for key, values in pre_rank.items()}
    ranks = sum_rank.items()
    p = Path(file)
    part = {"_".join(str(p.stem).split("_")[:2]): sorted(ranks, key=lambda elem: elem[1], reverse=True)}
    return part


def batch_search(files: Tuple[str], dim: int, n_samples: int, k_nearest: int, index: faiss.IndexFlatL2,
                 map_data: Dict[int, str]) -> Dict[str, List[Tuple[str, float]]]:
    local_rank = {}
    for file in files:
        local_rank.update(do_search(file, dim, n_samples, k_nearest, index, map_data))
    return local_rank


def partition_queries(files: Iterable[str], partitions: int = cpu_count() - 1) -> np.ndarray:
    partitions = np.array_split(np.array(files), partitions)
    # debug for partitions
    # print(partitions)
    # print([ar.shape for ar in partitions])
    return partitions


def run_procedure(input_dir: str, output: str, k_nearest: int, dim: int, n_samples: int, faiss_index: str, map: str):
    # colorama
    init()

    # timer
    start = timer()
    with open(map, "r") as fh:
        map_data = {int(key): value for key, value in json.load(fh).items()}
    index = faiss.read_index(faiss_index)
    print(type(index))
    # set_loky_pickler("pickle")
    # temporarily changed to loky; add later: prefer="threads"
    # ranks = Parallel(verbose=True, n_jobs=-1)(
    #     delayed(do_search)(file, dim, n_samples, k_nearest, index, map_data) for file in glob.glob(f"{input_dir}*.vec"))
    ranks = Parallel(verbose=True, n_jobs=-1)(
        delayed(batch_search)(batch, dim, n_samples, k_nearest, index, map_data) for batch in
        partition_queries(glob.glob(f"{input_dir}*.vec")))

    for result in ranks:
        global_rank.update(result)

    rank_file = Path(output)
    if rank_file.exists():
        os.system(f"rm {output}")
    with open(output, "w") as fd:
        json.dump(global_rank, fd, indent=4)

    end = timer()
    runtime = end - start
    print(f"{Fore.GREEN} Done in {runtime:.6f} seconds")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="faiss search deployer")
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Directory with virus samples of vectors.")
    parser.add_argument("-o", "--output", required=True,
                        help="Path to result file with rankings.")
    parser.add_argument("-d", "--dim", required=True,
                        help="Dimensionality of vectors")
    parser.add_argument("--n_samples", required=True,
                        help="Number of samples from each host genome")
    parser.add_argument("-k", "--k_nearest", required=True,
                        help="How many (k) nearest samples should be found")
    parser.add_argument("-f", "--faiss_index", required=True,
                        help="Path to faiss index")
    parser.add_argument("-m", "--map", required=True,
                        help="Path to sample map")

    args = parser.parse_args()

    global_rank = {}
    run_procedure(args.input_dir, args.output, int(args.k_nearest), int(args.dim), int(args.n_samples),
                  args.faiss_index, args.map)
