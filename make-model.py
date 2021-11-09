import json
import secrets
from collections import defaultdict
from timeit import default_timer as timer
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed
import argparse
import utils
import os


all_tax_levels = {
    "superkingdom": 0,
     "phylum": 1,
     "class": 2,
     "order": 3,
     "family": 4,
     "genus": 5,
     "species": 6
}


def fasta_parallel(file):
    record = SeqIO.read(file, "fasta")
    return record
    # with open("X:/edwards2016/host/random_phylum-training_fastDNA.fasta", "a") as w_fh:
    #     SeqIO.write(record, w_fh, "fasta")
    # #records.append(record)
    # with open("X:/edwards2016/host/random_phylum-training_labels.txt", "a") as fh:
    #     fh.write(label + "\n")


def main_procedure(input_dir, out_dir, filter, dim, length, minn, maxn, epoch, thread):
    # timeit
    start = timer()

    # importing json data about hosts
    host_data = utils.get_host_data()

    # reading config
    config = utils.get_config()

    # filtering all hosts by a chosen level
    filter_host = defaultdict(list)
    for host in host_data:
        level = host_data[host]["lineage_names"][all_tax_levels[filter]]  # tax level code
        filter_host[level].append(host)

    # random sampling of a single host from a phylum
    random_level_host = {level: secrets.choice(filter_host[level]) for level in filter_host}
    print(random_level_host)

    # list of filenames to be used
    filenames = list(random_level_host.values())
    print(filenames)

    # get list of taxid (labels in fastDNA)
    labels = [host_data[host]["taxid"] for host in filenames]
    print(len(labels))

    # parse selected files and merge them into single fasta file; create labels file
    #records = [list(SeqIO.parse(f"D:/praktyki2020/edwards2016/host/fasta/{file}.fna", "fasta"))[0] for file in filenames]
    # for file, label in zip(filenames, labels):
    #     print(f"{file} - {label}")

    par = Parallel(n_jobs=-1, verbose=11, pre_dispatch='all', batch_size="auto", backend="loky")(delayed(fasta_parallel)(f"{input_dir}{file}.fna") for file in filenames)
    #print(par)
    # print(len(records))
    # print(records)

    filtered_fasta_file = f"random_train-{filter}-dim_{dim}-len_{length}.fasta"
    with open(f"{out_dir}{filtered_fasta_file}", "w") as w_fh:
        SeqIO.write(par, w_fh, "fasta")

    labels_file = f"random_labels-{filter}-dim_{dim}-len_{length}.txt"
    with open(f"{out_dir}{labels_file}", "w") as fh:
        for label in labels:
            fh.write(label + "\n")

    # run fastDNA training
    model_file = f"random_model-{filter}-dim_{dim}-len_{length}"
    os.system(f"{config['GENERAL']['fastdna_dir']}fastdna supervised -input {out_dir}{filtered_fasta_file} -labels {out_dir}{labels_file} -output {out_dir}{model_file} -dim {dim} -length {length} -minn {minn} -maxn {maxn} -epoch {epoch} -thread {thread}")
    #par = Parallel(n_jobs=-1)(delayed(SeqIO.write)(records, "D:/praktyki2020/edwards2016/host/random_phylum-training_fastDNA.fasta", "fasta") for record in records)
    #SeqIO.write(records, "X:/edwards2016/host/random_phylum-training_fastDNA.fasta", "fasta")
    #print(f"Written {len(records)} records.")

    end = timer()
    runtime = end - start
    print(f"Done in {runtime:.6f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="fastDNA model creation")
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Directory with host genomes.")
    parser.add_argument("-o", "--output", required=True,
                        help="Path to result FASTA file, labels file and model file.")
    parser.add_argument("--filter", required=True,
                        help="Taxonomy level to which genomes should be filtered.")
    parser.add_argument("-d", "-dim", required=True,
                        help="Dimensionality of vectors")
    parser.add_argument("--length", required=True,
                        help="Length of sequences")
    parser.add_argument("--minn", required=True,
                        help="Minn")
    parser.add_argument("--maxn", required=True,
                        help="Maxn")
    parser.add_argument("--epoch", required=True,
                        help="Epoch")
    parser.add_argument("-t", "--thread", required=True,
                        help="Number of threads to use")

    args = parser.parse_args()

    main_procedure(args.input_dir, args.output, args.filter, args.dim, args.length, args.minn, args.maxn, args.epoch,
                   args.thread)
