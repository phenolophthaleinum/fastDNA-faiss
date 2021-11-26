import glob
from collections import defaultdict
from typing import Dict, Tuple, List
import math

import numpy as np
import faiss
import json
import argparse
from joblib import Parallel, delayed
from timeit import default_timer as timer
from colorama import Fore, init

# dim = 10
# k_nearest = 10
# query = np.loadtxt("/home/hyperscroll/fastDNA/virus_vectors.txt", dtype="float32")
# index = faiss.read_index("/home/hyperscroll/fastDNA/host_index.index")
# with open("/home/hyperscroll/edwards2016/host/sample_map.json", "r") as fh:
#     map_data = {int(key): value for key, value in json.load(fh).items()}
# distances, indices = index.search(query, k_nearest)
# print(f"Nearest sequences: {indices}")
# classification = np.vectorize(map_data.get)(indices)
# print(f"Classification: {classification}")
#
# pre_rank = defaultdict(list)
# for index, distance in np.stack((classification.flatten(), distances.flatten()), axis=1):
#     pre_rank[index].append(2 - float(distance))
# sum_rank = {key: sum(values) for key, values in pre_rank.items()}
# final_rank = dict(sorted(sum_rank.items(), key=lambda elem: elem[1], reverse=True))
# print(json.dumps(final_rank, indent=4))


def do_search(file: str, dim: int, n_samples: int, k_nearest: str, index: object, map_data: Dict[int, str]) -> Dict[str, List[Tuple[str, float]]]:
    # with txt file;
    # query = np.loadtxt(file, dtype="float32")
    # dim = np.shape(query)[1]
    # with binary:
    query = np.fromfile(file, dtype="float32")
    #n_host = query.size / (n_samples * dim)
    query = query.reshape(n_samples, dim)
    hyper = 2 * math.sqrt(dim)
    #index = faiss.read_index(faiss_index)
    distances, indices = index.search(query, int(k_nearest))
    #distances, indices = faiss_index.search(query, int(k_nearest))
    classification = np.vectorize(map_data.get)(indices)
    #ranks = []
    pre_rank = defaultdict(list)
    # DONE new scoring method (diagonal of a hypercube)
    # TODO filter scores for each sample to only contain scores from unique host
    # DONE raise error to see where negative scores appear - probably no errors (i will be sure if i test on better machine)
    # print(f"Max dist: {distances.max()}")
    # print(f"Min dist: {distances.min()}")

    for index, distance in np.stack((classification.flatten(), distances.flatten()), axis=1):
        #ranks.append((index, 2 - float(distance))) # this is wrong
        # 2 - dist
        #pre_rank[index].append(2 - float(distance))
        # hyper - dist

        score = hyper - float(distance)
        pre_rank[index].append(score)
        # if score >= 0:
        #     pre_rank[index].append(score)
        # else:
        #     raise ValueError(f"Negative score: {score}; index: {index}; distance: {distance}")

    sum_rank = {key: sum(values) for key, values in pre_rank.items()}
    ranks = sum_rank.items()
    part = {"_".join(file.split("/")[-1].split(".")[0].split("_")[:2]): sorted(ranks, key=lambda elem: elem[1], reverse=True)}
    return part
    #global_rank.update(part)


def run_procedure(input_dir, output, k_nearest, dim, n_samples, faiss_index, map):
    # colorama
    init()

    # timer
    start = timer()
    with open(map, "r") as fh:
        map_data = {int(key): value for key, value in json.load(fh).items()}
    index = faiss.read_index(faiss_index)
    ranks = Parallel(verbose=True, n_jobs=-1)(
            delayed(do_search)(file, dim, n_samples, k_nearest, index, map_data) for file in glob.glob(f"{input_dir}*.vec"))
    # txt version
    # ranks = Parallel(verbose=True, n_jobs=-1)(
    #     delayed(do_search)(file, k_nearest, index, map_data) for file in glob.glob(f"{input_dir}*.txt"))
    for result in ranks:
        global_rank.update(result)
    with open(output, "w") as fd:
        json.dump(global_rank, fd, indent=4)

    end = timer()
    runtime = end - start
    print(f"{Fore.GREEN} Done in {runtime:.6f} seconds")

#input: /home/hyperscroll/edwards2016/virus/vectors/
#output: virus_rank.json
#index: /home/hyperscroll/fastDNA/host_index_100.index
#map: /home/hyperscroll/edwards2016/host/sample_map_100.json
# python faiss_search.py -i /home/hyperscroll/edwards2016/virus/vectors/ -o virus_rank.json -k 10 -f /home/hyperscroll/fastDNA/host_index_100.index -m /home/hyperscroll/edwards2016/host/sample_map_100.json


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
    run_procedure(args.input_dir, args.output, args.k_nearest, int(args.dim), int(args.n_samples), args.faiss_index, args.map)
