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


def virus2vector(file: str, wd: str):
    # print(file)
    # fastdna_dir = "/home/hyperscroll/fastDNA/"
    name = file.split("/")[-1].split(".")[0]
    # txt version
    # os.system(f"{config['GENERAL']['fastdna_dir']}fastdna print-word-vectors {config['GENERAL']['active_model']} < {file} > {wd}virus/vectors/{name}_vector.txt")
    os.system(
        f"{config['GENERAL']['fastdna_dir']}fastdna print-word-vectors {config['GENERAL']['active_model']} {wd}virus/vectors/{name}_vector < {file}")


def host2vector(file: str, wd: str):
    # txt version
    # os.system(f"{config['GENERAL']['fastdna_dir']}fastdna print-word-vectors {config['GENERAL']['active_model']} < {file} > {wd}host/vectors/host_vectors.txt")
    os.system(
        f"{config['GENERAL']['fastdna_dir']}fastdna print-word-vectors {config['GENERAL']['active_model']} {wd}host/vectors/host_vectors < {file}")


def main_procedure(wd: str, host: bool, virus: bool, full: bool, length: int, n_vir_samples: int, n_host_samples: int,
                   dim: int, thread: int):
    # colorama
    init()

    # global timer
    total_start = timer()

    dirs = ["samples", "vectors", "maps", "index"]
    for dir in dirs:
        pathlib.Path(f"{wd}host/{dir}").mkdir(parents=True, exist_ok=True)
        if dirs.index(dir) < 2:
            pathlib.Path(f"{wd}virus/{dir}").mkdir(parents=True, exist_ok=True)

    # SAMPLING
    random_sampling.main_procedure(wd, host, virus, full, config["HOST"]["host_genomes"],
                                   config["VIRUS"]["virus_genomes"], length, n_vir_samples, n_host_samples)
    # fastdna_dir = "/home/hyperscroll/fastDNA/"
    # name = "/home/hyperscroll/edwards2016/virus/samples/NC_000866.4_samples.fasta".split("/")[-1].split(".")[0]
    # print(name)
    # os.system(
    #    f"{fastdna_dir}fastdna print-word-vectors {fastdna_dir}edwards_random_model.bin < /home/hyperscroll/edwards2016/virus/samples/NC_000866.4_samples.fasta > /home/hyperscroll/edwards2016/virus/vectors/{name}_vector.txt")

    # VECTORISING
    start = timer()
    if virus:
        vectors = Parallel(verbose=True, n_jobs=-1)(delayed(virus2vector)(file, wd)
                                                    for file in glob.glob(f"{wd}virus/samples/*.fasta"))
    if host:
        host2vector(f"{wd}host/samples/host_samples.fasta", wd)
    if full:
        vectors = Parallel(n_jobs=-1)(delayed(virus2vector)(file, wd)
                                      for file in glob.glob(f"{wd}virus/samples/*.fasta"))
        host2vector(f"{wd}host/samples/host_samples.fasta", wd)
    end = timer()
    runtime = end - start
    print(f"{Fore.GREEN} [fasta2vector] time: {runtime:.6f} seconds")

    # INDEXING
    if host or full:
        start = timer()
        index_building.build_index(f"{wd}host/vectors/host_vectors.vec", dim, n_host_samples,
                                   f"{wd}host/index/host_index.index")
        end = timer()
        runtime = end - start
        print(f"{Fore.GREEN} [index_building] time: {runtime:.6f} seconds")

    total_end = timer()
    total_runtime = total_end - total_start
    print(f"{Fore.GREEN} Total elapsed time: {total_runtime:.6f} seconds")


if __name__ == "__main__":
    # main_procedure()
    # do not delete, this actually a main code
    parser = argparse.ArgumentParser(description="fastDNA+faiss virus-host interaction analysis")
    parser.add_argument("-w", "--wd", required=True,
                        help="Working directory, where all files will be deployed")
    parser.add_argument("--host", required=False, action="store_true",
                        help="Host mode: every available host genome is randomly sampled according to a given criteria,"
                             " then a cloud of host vectors is generated and compiled into a faiss index")
    parser.add_argument("--virus", required=False, action="store_true",
                        help="Virus mode: every available virus genome is randomly sampled according to a given "
                             "criteria, then a cloud of virus vectors is generated which is compared with host cloud "
                             "and results are generated in a form of a rank of virus-host pairs.")
    parser.add_argument("--full", required=False, action="store_true",
                        help="Full mode: Combines host and virus mode in one go.")
    # parser.add_argument("-m", "--model", required=True,
    #                     help="fastDNA trained model full path")
    # parser.add_argument("-o", "--output", required=True,
    #                     help="Path to result FASTA file, labels file and model file.")
    parser.add_argument("--length", required=True,
                        help="Length of the samples")
    parser.add_argument("-d", "--dim", required=True,
                        help="Dimensionality of vectors")
    parser.add_argument("--n_vir", required=True,
                        help="Number of samples to take from a virus genome")
    parser.add_argument("--n_host", required=True,
                        help="Number of samples to take from a host genome")
    parser.add_argument("-t", "--thread", required=True,
                        help="Number of threads to use")

    config = utils.get_config()
    args = parser.parse_args()

    main_procedure(args.wd, args.host, args.virus, args.full, int(args.length), int(args.n_vir), int(args.n_host),
                   int(args.dim), int(args.thread))
