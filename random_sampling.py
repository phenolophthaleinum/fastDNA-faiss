import json
import secrets
import glob
import utils
from timeit import default_timer as timer
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed
from typing import List
from colorama import Fore, init


def sample_sequences(file: str, length: int, n: int, wd: str, virus: bool) -> List[SeqRecord]:
    """Function reads each given fasta file, randomly samples n subsequences of a defined length, creates new records with those subsequences and puts these new records into a list (which is finally returned).

    Args:
        file (str): Full path to a fasta file
        wd (str): Working directory path
        length (int): Length of a sampled subsequence
        n (int): Number of subsequences to be sampled
        virus (bool): Flag which enables or disables virus genomes sampling

    Returns:
        list[SeqRecord]: Returns list of newly created biopython `SeqRecord`, each contains sampled subsequences.
    """
    new_records = []
    record = SeqIO.read(file, "fasta")
    for i in range(n):
        pick = secrets.choice(range(0, len(record.seq) - length))
        desc = f'{record.description} sample_{i} gen_pos:({pick}:{pick + length})'
        new_record = SeqRecord(
            id=record.id,
            seq=Seq(str(record.seq[pick:pick + length])),
            description=desc
        )
        new_records.append(new_record)

    if virus:
        # if virus, this is possible
        with open(f"{wd}virus/samples/{new_records[0].id}_samples.fasta", "a") as w_fh:
            SeqIO.write(new_records, w_fh, "fasta")

    return new_records


def main_procedure(wd, host, virus, full, host_dir, virus_dir, length, n_samples):
    # colorama
    init()

    # timer
    start = timer()

    # importing json data about hosts
    utils.get_host_data()

    final_records = []
    # parallel sampling records of all files and dumping them into one file
    if host:
        new_records = Parallel(n_jobs=-1, verbose=True)(delayed(sample_sequences)(file, length, n_samples, wd, virus) for file in glob.glob(f"{host_dir}*.fna"))

        # for host
        for sublist in new_records:
            final_records.extend(sublist)
        with open(f"{wd}host/samples/host_samples.fasta", "a") as w_fh:
            SeqIO.write(final_records, w_fh, "fasta")

        # mapping samples to nbci ids and dumping them into a file; edit 10.11.21 - much faster
        p_records = list(
            SeqIO.parse(f"{wd}host/samples/host_samples.fasta", "fasta"))
        ids_records = [record.id.split(".")[0] for record in p_records]
        d = {}
        for index, value in enumerate(ids_records):
            d[index] = value
        # for id in set(ids_records):
        #     keys = [index for index, value in enumerate(ids_records) if value == id]
        #     for key in keys:
        #         d[key] = id
        with open(f"{wd}host/maps/sample_map.json", "w", encoding='utf-8') as fh:
            json.dump(d, fh, indent=4)

    if virus:
        new_records = Parallel(n_jobs=-1, verbose=True)(
            delayed(sample_sequences)(file, length, n_samples, wd, virus) for file in glob.glob(f"{virus_dir}*.fna"))

    if full:
        new_records = Parallel(n_jobs=-1, verbose=True)(
            delayed(sample_sequences)(file, length, n_samples, wd, virus=False) for file in glob.glob(f"{host_dir}*.fna"))

        # for host
        for sublist in new_records:
            final_records.extend(sublist)
        with open(f"{wd}host/samples/host_samples.fasta", "a") as w_fh:
            SeqIO.write(final_records, w_fh, "fasta")

        # mapping samples to nbci ids and dumping them into a file; edit 10.11.21 - much faster
        p_records = list(SeqIO.parse(f"{wd}host/samples/host_samples.fasta", "fasta"))
        ids_records = [record.id.split(".")[0] for record in p_records]
        d = {}
        for index, value in enumerate(ids_records):
            d[index] = value
        # for id in set(ids_records):
        #     keys = [index for index, value in enumerate(ids_records) if value == id]
        #     for key in keys:
        #         d[key] = id
        with open(f"{wd}host/maps/sample_map.json", "w", encoding='utf-8') as fh:
            json.dump(d, fh, indent=4)

        new_records = Parallel(n_jobs=-1, verbose=True)(
            delayed(sample_sequences)(file, length, n_samples, wd, virus=True) for file in glob.glob(f"{virus_dir}*.fna"))

    # for virus, but second slower
    # for sublist in new_records:
    #     with open(f"X:/edwards2016/virus/samples/{sublist[0].id}_samples.fasta", "a") as w_fh:
    #         SeqIO.write(sublist, w_fh, "fasta")

    end = timer()
    runtime = end - start
    print(f"{Fore.GREEN}[random-sampling] Done in {runtime:.6f} seconds")


# if __name__ == "__main__":
#     main()
