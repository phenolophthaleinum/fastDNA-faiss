import glob
import json
import random
import secrets
from collections import defaultdict
from timeit import default_timer as timer
from typing import List, Union, Tuple

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed
import argparse
import utils
import os
from colorama import Fore, init

all_tax_levels = {
    "superkingdom": 0,
    "phylum": 1,
    "class": 2,
    "order": 3,
    "family": 4,
    "genus": 5,
    "species": 6
}

# importing json data about hosts
host_data = utils.get_host_data()

# read config
config = utils.get_config()


def copy_parallel(src: str, dst: str):
    print(f'copy "{src}" "{dst}"')
    os.system(f'copy "{src}" "{dst}"')
    # with open("X:/edwards2016/host/random_phylum-training_fastDNA.fasta", "a") as w_fh:
    #     SeqIO.write(record, w_fh, "fasta")
    # #records.append(record)
    # with open("X:/edwards2016/host/random_phylum-training_labels.txt", "a") as fh:
    #     fh.write(label + "\n")


def tax_filtering(filter: str, reps: int) -> List[str]:
    # filtering all hosts by a chosen level
    filter_host = defaultdict(list)
    for host in host_data:
        level = host_data[host]["lineage_names"][all_tax_levels[filter]]  # tax level code
        filter_host[level].append(host)

    # optionally show original dataset composition
    # counts = []
    # for key in filter_host:
    #     counts.append((key, len(filter_host[key])))
    # counts.sort(key=lambda x: x[1], reverse=True)
    # print(counts)

    # print(f"{key}: {len(filter_host[key])}")
    # print(filter_host.items())
    # random sampling of a single host from a level
    random_level_host = {
        level: random.sample(filter_host[level], len(filter_host[level]) if len(filter_host[level]) < reps else reps)
        for
        level in filter_host}
    # print(random_level_host)

    # optionally show new dataset composition
    # counts = []
    # for key in random_level_host:
    #     counts.append((key, len(random_level_host[key])))
    # counts.sort(key=lambda x: x[1], reverse=True)
    # print(counts)

    # list of filenames to be used
    filenames = list(random_level_host.values())
    # print(filenames)

    # flatten if needed
    if isinstance(filenames[0], list):
        temp_files = [item for sublist in filenames for item in sublist]
        filenames = temp_files
    print(filenames)

    return filenames


if __name__ == "__main__":
    # colorama
    init()

    # timeit
    start = timer()

    filenames = tax_filtering("species", 1)

    # par = Parallel(n_jobs=-1, verbose=11, pre_dispatch='all', batch_size="auto", backend="loky")(
    #     delayed(copy_parallel)(f"D:\\edwards2016\\host\\fasta\\{file}.fna",
    #                            f"D:\\edwards2016\\host\\fasta-species-rep_1\\{file}.fna") for file in filenames)
    for file in filenames:
        copy_parallel(f"D:\\edwards2016\\host\\fasta\\{file}.fna", f"D:\\edwards2016\\host\\fasta-species-rep_1\\{file}.fna")

    end = timer()
    runtime = end - start
    print(f"{Fore.GREEN}[subsetting dataset] Done in {runtime:.6f}")
