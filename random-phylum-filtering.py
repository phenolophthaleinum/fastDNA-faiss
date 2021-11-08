import json
import secrets
from collections import defaultdict
from timeit import default_timer as timer
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed

#records = []


def fasta_parallel(file):
    record = SeqIO.read(file, "fasta")
    return record
    # with open("X:/edwards2016/host/random_phylum-training_fastDNA.fasta", "a") as w_fh:
    #     SeqIO.write(record, w_fh, "fasta")
    # #records.append(record)
    # with open("X:/edwards2016/host/random_phylum-training_labels.txt", "a") as fh:
    #     fh.write(label + "\n")


def main():
    # timeit
    start = timer()

    # importing json data about hosts
    with open("host.json", "r") as fh:
        host_data = json.load(fh)

    # filtering all hosts by phylum
    phylum_host = defaultdict(list)
    for host in host_data:
        phylum = host_data[host]["lineage_names"][4]  # 1 - phylum, 4 - family
        phylum_host[phylum].append(host)

    # random sampling of a single host from a phylum
    random_phylum_host = {phylum: secrets.choice(phylum_host[phylum]) for phylum in phylum_host}
    print(random_phylum_host)

    # list of filenames to be used
    filenames = list(random_phylum_host.values())
    print(filenames)

    # get list of taxid (labels in fastDNA)
    labels = [host_data[host]["taxid"] for host in filenames]
    print(len(labels))

    # parse selected files and merge them into single fasta file; create labels file
    #records = [list(SeqIO.parse(f"D:/praktyki2020/edwards2016/host/fasta/{file}.fna", "fasta"))[0] for file in filenames]
    # for file, label in zip(filenames, labels):
    #     print(f"{file} - {label}")

    #par = Parallel(n_jobs=-1, verbose=11, pre_dispatch='all', batch_size="auto", backend="loky")(delayed(fasta_parallel)(f"D:/edwards2016/host/fasta/{file}.fna") for file in filenames)
    #print(par)
    # print(len(records))
    # print(records)

    # with open("D:/edwards2016/host/random_family-training_fastDNA.fasta", "a") as w_fh:
    #     SeqIO.write(par, w_fh, "fasta")
    #
    with open("D:/edwards2016/host/random_family-training_labels_test.txt", "w") as fh:
        for label in labels:
            fh.write(label + "\n")

    #par = Parallel(n_jobs=-1)(delayed(SeqIO.write)(records, "D:/praktyki2020/edwards2016/host/random_phylum-training_fastDNA.fasta", "fasta") for record in records)
    #SeqIO.write(records, "X:/edwards2016/host/random_phylum-training_fastDNA.fasta", "fasta")
    #print(f"Written {len(records)} records.")

    end = timer()
    runtime = end - start
    print(f"Done in {runtime:.6f}")


if __name__ == "__main__":
    main()
