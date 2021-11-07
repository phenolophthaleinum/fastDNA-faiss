import glob
from collections import defaultdict

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


def do_search(file, k_nearest, faiss_index, map_data):
    query = np.loadtxt(file, dtype="float32")
    index = faiss.read_index(faiss_index)
    distances, indices = index.search(query, int(k_nearest))
    classification = np.vectorize(map_data.get)(indices)
    #ranks = []
    pre_rank = defaultdict(list)
    for index, distance in np.stack((classification.flatten(), distances.flatten()), axis=1):
        #ranks.append((index, 2 - float(distance)))
        pre_rank[index].append(2 - float(distance))
    sum_rank = {key: sum(values) for key, values in pre_rank.items()}
    ranks = sum_rank.items()
    part = {file.split("/")[-1].split(".")[0]: sorted(ranks, key=lambda elem: elem[1], reverse=True)}
    return part
    #global_rank.update(part)


def run_procedure(input_dir, output, k_nearest, faiss_index, map):
    # colorama
    init()

    # timer
    start = timer()
    with open(map, "r") as fh:
        map_data = {int(key): value for key, value in json.load(fh).items()}
    ranks = Parallel(n_jobs=-1)(
        delayed(do_search)(file, k_nearest, faiss_index, map_data) for file in glob.glob(f"{input_dir}*.txt"))
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
    # parser.add_argument("-d", "-dim", required=True,
    #                     help="Dimensionality of vectors")
    parser.add_argument("-k", "--k_nearest", required=True,
                        help="How many (k) nearest samples should be found")
    parser.add_argument("-f", "--faiss_index", required=True,
                        help="Path to faiss index")
    parser.add_argument("-m", "--map", required=True,
                        help="Path to sample map")

    args = parser.parse_args()

    global_rank = {}
    run_procedure(args.input_dir, args.output, args.k_nearest, args.faiss_index, args.map)
