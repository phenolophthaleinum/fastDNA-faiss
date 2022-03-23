from joblib import Parallel, delayed
import os
import glob
from colorama import Fore, init
from timeit import default_timer as timer
import argparse
import utils
import pathlib
import numpy as np
from typing import Dict, Tuple, List, Iterable
from multiprocessing import cpu_count
import predict_prob2result as pred
import taxonomic_discordance2 as td
import json


def run_predict_prob(file: str, k_best: int, wd: str, div: dict):
    """Runs fastdna predict-prob module. This function is run from [`batch_exec` function](#workflow2.batch_exec)

    Args:
        file: Path to FASTA file to be used in prediction
        k_best: Number of best results from prediction
        wd: Path to working directory, which completes path to output file
    """
    name = file.split("/")[-1].split(".")[0]
    pathlib.Path(f"{wd}virus/preds/{name}").mkdir(parents=True, exist_ok=True)
    results_table = {}
    order_result_dict = {}
    for order in div[name]:
        os.system(
            f"{config['GENERAL']['fastdna_dir']}fastdna predict-prob {config['MODELS'][order.lower()]} {file} {k_best} > {wd}virus/preds/{name}/{name}_pred.json"
        )

        for name in pred.scoring_functions.keys():
            preds = pred.run_procedure(input_dir=f"{wd}virus/preds/{name}",
                                       scoring_function=name)
            accordance = td.taxonomic_accordance_sp(dists, preds, virus_dict)
            results_table[name] = {'result_dict': preds,
                                   'accordance': accordance}

        best_score = max(d['accordance'] for d in results_table.values())
        best_func = [name for name, data in results_table.items() if data['accordance'] == best_score][0]
        order_result_dict[order] = {
            "rank": preds,
            "best_score": best_score,
            "best_func": best_func
        }

    best_score_order = max(d['best_score'] for d in order_result_dict.values())
    best_order = [order for order, data in order_result_dict.items() if data['best_score'] == best_score_order][0]

    ## save here
    pathlib.Path(f"{wd}rank/{best_order}").mkdir(parents=True, exist_ok=True)
    with open(f"{wd}rank/{best_order}/{best_order}_{name}_rank_{order_result_dict[best_order]['best_func']}.json") as fh:
        json.dump(order_result_dict[best_order]['rank'], fh, indent=4)

        # pred.dump_result(output=f'{workflow_wd}rank/{search_final_rank}',
        #                  result_rank=results_table[best_func]['result_dict'],
        #                  scoring_function=best_func)


def batch_exec(files: Tuple[str], k_best: int, wd: str, div: dict):
    """Batch fastdna predict-prob module execution

    Args:
        files: Tuple of paths to FASTA files to be used in prediction
        k_best: Number of best results from prediction
        wd: Path to working directory, which completes path to output file
    """
    for file in files:
        run_predict_prob(file, k_best, wd, div)


def partition_queries(files: Iterable[str], partitions: int = cpu_count() - 1) -> np.ndarray:
    """Partitions set of paths to FASTA files with virus genomes to be batched.

    Args:
        files: Iterable of paths FASTA files with virus genomes
        partitions: Number of paths to be collected into a batch.
            !!! info
            Currently, it is a predefined parameter which equals a number of CPU threads.

    Returns:
        ndarray: n-dimensional Numpy array with batches of file paths
    """
    partitions = np.array_split(np.array(files), partitions)
    # debug for partitions
    # print(partitions)
    # print([ar.shape for ar in partitions])
    return partitions


def main_procedure(wd: str, length: int, n_vir_samples: int, thread: int, n_nuc_threshold: float, k_best: int,
                   div: str):
    """Main function of the module that runs:
        - filesystem creation for working directory
        - random sampling of virus genomes
        - fastna predict-prob

    Args:
        div:
        wd: Path to working directory
        length: Length of fragments of the genomes
        n_vir_samples: Number of samples to take from a virus genome
        thread: Number of CPU threads to be used (more threads - faster execution, but higher CPU usage)
        n_nuc_threshold: Maximum allowed ambiguous nucleotide content threshold in sampled sequences (in percents)
        k_best: Number of best results from prediction
    """
    # colorama
    init()

    # global timer
    total_start = timer()

    # create necessary directories in the main working directory
    # dirs = ["samples", "preds"]
    # for dir in dirs:
    #     pathlib.Path(f"{wd}virus/{dir}").mkdir(parents=True, exist_ok=True)
    # pathlib.Path(f"{wd}rank").mkdir(parents=True, exist_ok=True)
    with open(div, 'r') as fh:
        div_dict = json.load(fh)

    # pred = Parallel(n_jobs=thread, verbose=True)(delayed(run_predict_prob)(file, k_best, wd)
    #                                for file in glob.glob(f"{wd}virus/samples/*.fasta"))
    pred = Parallel(n_jobs=thread, verbose=True)(delayed(batch_exec)(batch, k_best, wd, div_dict)
                                                 for batch in
                                                 partition_queries(glob.glob(f"{wd}virus/samples/*.fasta")))

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
    parser.add_argument("--div", required=True,
                        help="Division table")

    config = utils.get_config()
    args = parser.parse_args()
    host_json = pathlib.Path('host.json')
    with host_json.open() as hj:
        host_dict = json.load(hj)

    virus_json = pathlib.Path('virus.json')
    with virus_json.open() as hj:
        virus_dict = json.load(hj)

    dists = td.DistanceMatrix(host_dict)

    main_procedure(args.wd, int(args.length), int(args.n_vir), int(args.thread), float(args.n_nucleotide_threshold),
                   int(args.k_best), str(args.div))
