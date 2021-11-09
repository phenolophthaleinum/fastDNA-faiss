from joblib import Parallel, delayed
import os
import glob
from colorama import Fore, init
from timeit import default_timer as timer
import argparse


def fasta2vector(file):
    #print(file)
    fastdna_dir = "/home/hyperscroll/fastDNA/"
    name = file.split("/")[-1].split(".")[0]
    os.system(f"{fastdna_dir}fastdna print-word-vectors {fastdna_dir}edwards_random_family_model.bin < {file} > /home/hyperscroll/edwards2016/virus/sample_100-len_250/vectors/{name}_vector.txt")


def main_procedure():
    # colorama
    init()

    # timer
    start = timer()
    #fastdna_dir = "/home/hyperscroll/fastDNA/"
    #name = "/home/hyperscroll/edwards2016/virus/samples/NC_000866.4_samples.fasta".split("/")[-1].split(".")[0]
    #print(name)
    #os.system(
    #    f"{fastdna_dir}fastdna print-word-vectors {fastdna_dir}edwards_random_model.bin < /home/hyperscroll/edwards2016/virus/samples/NC_000866.4_samples.fasta > /home/hyperscroll/edwards2016/virus/vectors/{name}_vector.txt")
    vectors = Parallel(n_jobs=-1)(delayed(fasta2vector)(file) for file in glob.glob("/home/hyperscroll/edwards2016/virus/sample_100-len_250/samples/*.fasta"))
    final_records = []

    # for host
    # for sublist in new_records:
    #     final_records.extend(sublist)
    # with open("X:/edwards2016/host/random_100_samples-training_fastDNA.fasta", "a") as w_fh:
    #     SeqIO.write(final_records, w_fh, "fasta")

    # for virus, but second slower
    # for sublist in new_records:
    #     with open(f"X:/edwards2016/virus/samples/{sublist[0].id}_samples.fasta", "a") as w_fh:
    #         SeqIO.write(sublist, w_fh, "fasta")

    # mapping samples to nbci ids and dumping them into a file
    # p_records = list(SeqIO.parse("X:/edwards2016/host/random_100_samples-training_fastDNA.fasta", "fasta"))
    # ids_records = [record.id.split(".")[0] for record in p_records]
    # d = {}
    # for id in set(ids_records):
    #     keys = [index for index, value in enumerate(ids_records) if value == id]
    #     for key in keys:
    #         d[key] = id
    # with open("X:/edwards2016/host/sample_map_100.json", "w", encoding='utf-8') as fh:
    #     json.dump(d, fh, indent=4)

    end = timer()
    runtime = end - start
    print(f"{Fore.GREEN} Done in {runtime:.6f} seconds")


if __name__ == "__main__":
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
    # parser.add_argument("-o", "--output", required=True,
    #                     help="Path to result FASTA file, labels file and model file.")
    parser.add_argument("--length", required=True,
                        help="Length of the samples")
    parser.add_argument("-n", "--n_samples", required=True,
                        help="Number of samples to take from a genome")
    parser.add_argument("-t", "--thread", required=True,
                        help="Number of threads to use")

    args = parser.parse_args()

    main_procedure(args.input_dir, args.output, args.filter, args.dim, args.length, args.minn, args.maxn, args.epoch,
                   args.thread)