from joblib import Parallel, delayed
import os
import glob
from colorama import Fore, init
from timeit import default_timer as timer
import argparse
import utils
import random_sampling
import index_building
import pathlib
import numpy as np
from typing import Dict, Tuple, List, Iterable
from multiprocessing import cpu_count


def run_predict_prob(file: str, k_best: int, wd: str):
    name = file.split("/")[-1].split(".")[0]
    os.system(
        f"{config['GENERAL']['fastdna_dir']}fastdna predict-prob {config['GENERAL']['active_model']} {file} {k_best} > {wd}virus/preds/{name}_pred.json"
    )


def batch_exec(files: Tuple[str], k_best: int, wd: str):
    for file in files:
        run_predict_prob(file, k_best, wd)


def partition_queries(files: Iterable[str], partitions: int = cpu_count() - 1) -> np.ndarray:
    partitions = np.array_split(np.array(files), partitions)
    # debug for partitions
    # print(partitions)
    # print([ar.shape for ar in partitions])
    return partitions


def main_procedure(wd: str, length: int, n_vir_samples: int, thread: int, n_nuc_threshold: float, k_best: int):
    # colorama
    init()

    # global timer
    total_start = timer()

    # create necessary directories in the main working directory
    dirs = ["samples", "preds"]
    for dir in dirs:
        pathlib.Path(f"{wd}virus/{dir}").mkdir(parents=True, exist_ok=True)
    pathlib.Path(f"{wd}rank").mkdir(parents=True, exist_ok=True)

    # SAMPLING
    random_sampling.main_procedure(wd, False, True, False, config["HOST"]["host_genomes"],
                                   config["VIRUS"]["virus_genomes"], length, n_vir_samples, 0,
                                   n_nuc_threshold)

    # pred = Parallel(n_jobs=thread, verbose=True)(delayed(run_predict_prob)(file, k_best, wd)
    #                                for file in glob.glob(f"{wd}virus/samples/*.fasta"))
    pred = Parallel(n_jobs=thread, verbose=True)(delayed(batch_exec)(batch, k_best, wd)
                                                 for batch in partition_queries(glob.glob(f"{wd}virus/samples/*.fasta")))

    total_end = timer()
    total_runtime = total_end - total_start
    print(f"{Fore.GREEN} Total elapsed time: {total_runtime:.6f} seconds")


if __name__ == "__main__":
    # main_procedure()
    # do not delete, this is actually a main code
    parser = argparse.ArgumentParser(description="fastDNA+faiss virus-host interaction analysis")
    parser.add_argument("-w", "--wd", required=True,
                        help="Working directory, where all files will be deployed")
    # parser.add_argument("-m", "--model", required=True,
    #                     help="fastDNA trained model full path")
    # parser.add_argument("-o", "--output", required=True,
    #                     help="Path to result FASTA file, labels file and model file.")
    parser.add_argument("--length", required=True,
                        help="Length of the samples")
    parser.add_argument("--n_vir", required=True,
                        help="Number of samples to take from a virus genome")
    parser.add_argument("--n_nucleotide_threshold", required=True,
                        help="Ambiguous nucleotide content threshold in sampled sequences (in percents)")
    parser.add_argument("-t", "--thread", required=True,
                        help="Number of threads to use")
    parser.add_argument("-k", "--k_best", required=True,
                        help="Number of best matches")

    config = utils.get_config()
    args = parser.parse_args()

    main_procedure(args.wd, int(args.length), int(args.n_vir), int(args.thread), float(args.n_nucleotide_threshold),
                   int(args.k_best))
