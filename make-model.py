import glob
import json
import random
import secrets
from collections import defaultdict
from timeit import default_timer as timer
from typing import List, Union, Tuple
from pathlib import Path
import re

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
    "species": 6,
    "hybrid": {
        "superkingdom": 0,
        "phylum": 1,
        "class": 2,
        "order": 3,
        "family": 4,
        "genus": 5,
        "species": 6,
    },
    "debug": {
        "superkingdom": 0,
        "phylum": 1,
        "class": 2,
        "order": 3,
        "family": 4,
        "genus": 5,
        "species": 6,
    }
}

# importing json data about hosts
host_data = utils.get_host_data()

# importing json data about hostvir
hostvir_data = utils.get_hostvir_data()

#importing json data about viruses
virus_data = utils.get_virus_data()

# read config
config = utils.get_config()


def fasta_parallel(file: str) -> SeqRecord:
    """Reads FASTA files and returns parsed data. This function is run in parallel.

    Args:
        file: Path to FASTA file.

    Returns:
        SeqRecord: Parsed data from FASTA file as a SeqRecord object.
    """
    record = SeqIO.read(file, "fasta")
    return record
    # with open("X:/edwards2016/host/random_phylum-training_fastDNA.fasta", "a") as w_fh:
    #     SeqIO.write(record, w_fh, "fasta")
    # #records.append(record)
    # with open("X:/edwards2016/host/random_phylum-training_labels.txt", "a") as fh:
    #     fh.write(label + "\n")


def tax_filtering(filter: str, reps: int) -> Tuple[List[str], List[str]]:
    """Edwards' dataset of host is filtered to certain taxonomy level and specific host genomes that represent
    taxonomy level are randomly chosen.

    Args:
        filter: A taxonomy level by which filtering is conducted
        reps: A number of representatives chosen randomly at selected level (defined by `filter`)

    Returns:
        Tuple[List[str], List[str]]: A list of ncbi_id of selected hosts genomes and a list of taxid of selected hosts.
    """
    # filtering all hosts by a chosen level
    filter_host = defaultdict(list)
    for host in host_data:
        level = host_data[host]["lineage_names"][all_tax_levels[filter]]  # tax level code
        filter_host[level].append(host)

    # random sampling of a single host from a level
    random_level_host = {
        level: random.sample(filter_host[level], len(filter_host[level]) if len(filter_host[level]) < reps else reps)
        for
        level in filter_host}
    print(random_level_host)

    # list of filenames to be used
    filenames = list(random_level_host.values())
    # print(filenames)

    # flatten if needed
    if isinstance(filenames[0], list):
        temp_files = [item for sublist in filenames for item in sublist]
        filenames = temp_files
    print(filenames)

    # get list of taxid (labels in fastDNA)
    labels = ["_".join(re.split(' |\; |\. |\, ', host_data[host]["lineage_names"][-1])) for host in filenames]
    # labels = [host_data[host]["taxid"] for host in filenames]
    print(len(labels))

    return filenames, labels


def no_filtering(input_dir: str) -> Tuple[List[str], List[str]]:
    filenames = [path.split("/")[-1].split(".")[0] for path in glob.glob(f"{input_dir}*.fna")]
    print(filenames)
    # labels = [host_data[host]["taxid"] for host in filenames]
    labels = ["_".join(re.split(' |\; |\. |\, ', host_data[host]["lineage_names"][-1])) for host in filenames]

    return filenames, labels


def hybrid_filtering(filter: str, reps: int) -> Tuple[List[str], List[str]]:
    # filtering all hosts by a chosen level
    filter_host = defaultdict(list)
    for host in host_data:
        level = host_data[host]["lineage_names"][all_tax_levels[filter]['family']]  # tax level code
        filter_host[level].append(host)

    # random sampling of a single host from a level
    random_level_host = {
        level: random.sample(filter_host[level], len(filter_host[level]) if len(filter_host[level]) < reps else reps)
        for
        level in filter_host}
    # print(random_level_host)

    # list with full paths
    random_level_host_paths = {}
    for level in random_level_host:
        random_level_host_paths.update({level: [f"{config['HOST']['host_genomes']}{host}.fna" for host in random_level_host[level]]})

    # add known viruses
    for level in random_level_host:
        to_add = []
        for host in random_level_host[level]:
            try:
                # virus_sample = random.sample(hostvir_data[host]['virus_id'], len(random_level_host[level]) if len(random_level_host[level]) < reps else reps)
                #virus_sample = random.sample(hostvir_data[host]['virus_id'], 1)
                virus_sample = random.sample(hostvir_data[host]['virus_id'], len(hostvir_data[host]['virus_id']))
                to_add.extend(virus_sample)
            except:
                continue
        random_level_host_paths[level].extend([f"{config['VIRUS']['virus_genomes']}{virus}.fna" for virus in to_add])
        random_level_host[level].extend(to_add)
    print(random_level_host)

    # list of filenames to be used
    filenames = list(random_level_host.values())
    # print(filenames)

    # list of paths to be returned
    filepaths = list(random_level_host_paths.values())

    # flatten if needed
    if isinstance(filenames[0], list):
        temp_files = [item for sublist in filenames for item in sublist]
        temp_files_path = [item for sublist in filepaths for item in sublist]
        filenames = temp_files
        filepaths = temp_files_path
    print(filenames)

    # get list of taxid (labels in fastDNA)
    labels = []
    for org in filenames:
        try:
            labels.append("_".join(re.split(' |\; |\. |\, ', host_data[org]["lineage_names"][-1])))
        except KeyError:
            labels.append("_".join(re.split(' |\; |\. |\, ', virus_data[org]["host"]["lineage_names"][-1])))
    #labels = ["_".join(re.split(' |\; |\. |\, ', host_data[host]["lineage_names"][-1])) for host in filenames]
    # labels = [host_data[host]["taxid"] for host in filenames]
    print(len(labels))

    return filepaths, labels


def debug_filtering(filter: str) -> Tuple[List[str], List[str]]:
    # filtering all hosts by a chosen level
    filter_host = defaultdict(list)
    target = ["Clostridiales"]
    filenames = []
    for name in target:
        for host in host_data:
            if host_data[host]["lineage_names"][all_tax_levels[filter]['order']] == name:
                filenames.append(host)

    # list of paths to be returned
    filepaths = [f"{config['HOST']['host_genomes']}{host}.fna" for host in filenames]

    # flatten if needed
    if isinstance(filenames[0], list):
        temp_files = [item for sublist in filenames for item in sublist]
        temp_files_path = [item for sublist in filepaths for item in sublist]
        filenames = temp_files
        filepaths = temp_files_path
    print(filenames)

    # get list of taxid (labels in fastDNA)
    # labels = []
    # for org in filenames:
    #     try:
    #         labels.append("_".join(re.split(' |\; |\. |\, ', host_data[org]["lineage_names"][-1])))
    #     except KeyError:
    #         labels.append("_".join(re.split(' |\; |\. |\, ', virus_data[org]["host"]["lineage_names"][-1])))
    labels = ["_".join(re.split(' |\; |\. |\, ', host_data[host]["lineage_names"][-1])) for host in filenames]
    # labels = [host_data[host]["taxid"] for host in filenames]
    print(len(labels))

    return filepaths, labels


def main_procedure(input_dir: str, out_dir: str, filter: str, dim: int, length: int, minn: int, maxn: int, epoch: int,
                   thread: int, reps: int, rm: bool, save_vec: bool):
    # colorama
    init()

    # timeit
    start = timer()

    out_path = Path(out_dir)
    if not out_path.exists():
        os.system(f"mkdir {str(out_path)}")

    filenames = []
    labels = []

    if filter not in ["none", "hybrid", "debug"]:
        filenames, labels = tax_filtering(filter, reps)
    if filter == "none":
        filenames, labels = no_filtering(input_dir)
    if filter == "hybrid":
        filenames, labels = hybrid_filtering(filter, reps)
    if filter == "debug":
        filenames, labels = debug_filtering(filter)

    # parse selected files and merge them into single fasta file; create labels file
    # records = [list(SeqIO.parse(f"D:/praktyki2020/edwards2016/host/fasta/{file}.fna", "fasta"))[0] for file in filenames]
    # for file, label in zip(filenames, labels):
    #     print(f"{file} - {label}")
    if filter not in ["hybrid", "debug"]:
        par = Parallel(n_jobs=-1, verbose=11, pre_dispatch='all', batch_size="auto", backend="loky")(
            delayed(fasta_parallel)(f"{input_dir}{file}.fna") for file in filenames)
    if filter in ['hybrid', "debug"]:
        par = Parallel(n_jobs=-1, verbose=11, pre_dispatch='all', batch_size="auto", backend="loky")(
            delayed(fasta_parallel)(file) for file in filenames)
    # print(par)
    # print(len(records))
    # print(records)

    filtered_fasta_file = f"random_train-{filter}-dim_{dim}-len_{length}.fasta"
    with open(f"{out_dir}{filtered_fasta_file}", "w+") as w_fh:
        SeqIO.write(par, w_fh, "fasta")

    labels_file = f"random_labels-{filter}-dim_{dim}-len_{length}.txt"
    with open(f"{out_dir}{labels_file}", "w") as fh:
        for label in labels:
            fh.write(label + "\n")

    # run fastDNA training
    model_file = f"random_model-{filter}-dim_{dim}-len_{length}-epoch{epoch}"
    if save_vec:
        os.system(
            f"{config['GENERAL']['fastdna_dir']}fastdna supervised -input {out_dir}{filtered_fasta_file} -labels {out_dir}{labels_file} -output {out_dir}{model_file} -dim {dim} -length {length} -minn {minn} -maxn {maxn} -epoch {epoch} -thread {thread} -saveVec")
    else:
        os.system(
            f"{config['GENERAL']['fastdna_dir']}fastdna supervised -input {out_dir}{filtered_fasta_file} -labels {out_dir}{labels_file} -output {out_dir}{model_file} -dim {dim} -length {length} -minn {minn} -maxn {maxn} -epoch {epoch} -thread {thread}")
    # par = Parallel(n_jobs=-1)(delayed(SeqIO.write)(records, "D:/praktyki2020/edwards2016/host/random_phylum-training_fastDNA.fasta", "fasta") for record in records)
    # SeqIO.write(records, "X:/edwards2016/host/random_phylum-training_fastDNA.fasta", "fasta")
    # print(f"Written {len(records)} records.")

    if rm:
        print("removing")
        os.system(f"rm {out_dir}{filtered_fasta_file}")
        os.system(f"rm {out_dir}{labels_file}")

    end = timer()
    runtime = end - start
    print(f"{Fore.GREEN}[make-model] Done in {runtime:.6f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="fastDNA model creation")
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Directory with host genomes.")
    parser.add_argument("-o", "--output", required=True,
                        help="Path to result FASTA file, labels file and model file.")
    parser.add_argument("-f", "--filter", required=True,
                        choices=["phylum", "class", "order", "family", "genus", "species", "none", 'hybrid', 'debug'],
                        help="Taxonomy level to which genomes should be filtered. Choosing 'none' implies no taxonomy filtering.")
    parser.add_argument("-r", "--reps", required=False, default=1,
                        help="Maximum number of representatives from the filtered group. Default value is 1.")
    parser.add_argument("-d", "--dim", required=True,
                        help="Dimensionality of vectors")
    parser.add_argument("-l", "--length", required=True,
                        help="Length of sequences")
    parser.add_argument("--minn", required=True,
                        help="Minimum k-mer size")
    parser.add_argument("--maxn", required=True,
                        help="Maximum k-mer size (max k=15, otherwise fastDNA fails)")
    parser.add_argument("-e", "--epoch", required=True,
                        help="Number of epochs (each added epoch increases runtime significantly)")
    parser.add_argument("-t", "--thread", required=True,
                        help="Number of threads to use")
    parser.add_argument("--rm", required=False, action="store_true", default=False,
                        help="Remove potentially redundant files after model creation. Default is 'false'.")
    parser.add_argument("--saveVec", required=False, action="store_true", default=False,
                        help="Enables saving of a readable model file (.vec). Enabling this may significantly increase execution time. Default is 'false'.")

    args = parser.parse_args()

    main_procedure(args.input_dir, args.output, args.filter, int(args.dim), int(args.length), int(args.minn),
                   int(args.maxn), int(args.epoch),
                   int(args.thread), int(args.reps), args.rm, args.saveVec)
