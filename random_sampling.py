import json
import os
import secrets
import glob
from pathlib import Path
import utils
from timeit import default_timer as timer
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed
from typing import List
from colorama import Fore, init


def sample_sequences(file: str, length: int, n: int, wd: str, n_nuc_threshold: float, virus: bool, ) -> List[SeqRecord]:
    """Function reads each given fasta file, randomly samples n subsequences of a defined length, creates new records with those subsequences and puts these new records into a list (which is finally returned).

    Args:
        file (str): Full path to a fasta file
        wd (str): Working directory path
        length (int): Length of a sampled subsequence
        n (int): Number of subsequences to be sampled
        n_nuc_threshold (float): Ambiguous nucleotide content threshold in sampled sequences (in percents)
        virus (bool): Flag which enables or disables virus genomes sampling

    Returns:
        list[SeqRecord]: Returns list of newly created biopython `SeqRecord`, each contains sampled subsequences.
    """
    new_records = []
    record = SeqIO.read(file, "fasta")
    for i in range(n):
        new_seq = None
        valid = False
        pick = None
        while not valid:
            pick = secrets.choice(range(0, len(record.seq) - length))
            new_seq = Seq(str(record.seq[pick:pick + length]))
            n_content = new_seq.count("N")
            if n_content / length < n_nuc_threshold:
                valid = True
            else:
                print(f"Too high N content ({n_content / length}%). Sampling again")

        desc = f'{record.description} sample_{i} gen_pos:({pick}:{pick + length})'
        new_record = SeqRecord(
            id=record.id,
            seq=new_seq,
            description=desc
        )
        new_records.append(new_record)

    if virus:
        # if virus, this is possible
        with open(f"{wd}virus/samples/{new_records[0].id}_samples.fasta", "w") as w_fh:
            SeqIO.write(new_records, w_fh, "fasta")

    return new_records


# DONE differentiate number of samples from virus and hosts (more from hosts most likely)
# DONE new sample sources - not only from all edwards host genomes but from each species representative
def main_procedure(wd: str, host: bool, virus: bool, full: bool, host_dir: str, virus_dir: str, length: int,
                   n_vir_samples: int, n_host_samples: int, n_nuc_threshold: float):
    """Main function of the module that samples genomes randomly:

    Args:
        wd: Working directory path
        host: Flag, which allows sampling only host genomes
            !!! p-obsolete "Obsolete"
                In current vision of the whole process, this parameter is obsolete. Using it might be undesirable and cause program to fail. Might be removed in later version.
        virus: Flag, which allows sampling only virus genomes.
        full: Flag, which allows sampling both virus and host genomes.
            !!! p-obsolete "Obsolete"
                In current vision of the whole process, this parameter is obsolete. Using it might be undesirable and cause program to fail. Might be removed in later version.
        host_dir: Path to directory where sampled host fragments will be saved
        virus_dir: Path to directory where sampled virus fragments will be saved
        length: Length of sampled fragments
        n_vir_samples: Number of samples to take from a virus genome
        n_host_samples: Number of samples to take from a host genome
        n_nuc_threshold: Maximum allowed ambiguous nucleotide content threshold in sampled sequences (in percents)
    """
    # colorama
    init()

    # timer
    start = timer()

    # importing json data about hosts
    utils.get_host_data()

    final_records = []
    # parallel sampling records of all files and dumping them into one file
    if host:
        new_records = Parallel(n_jobs=-1, verbose=True)(
            delayed(sample_sequences)(file, length, n_host_samples, wd, virus, n_nuc_threshold) for file in
            glob.glob(f"{host_dir}*.fna"))

        # for host
        for sublist in new_records:
            final_records.extend(sublist)

        # host_file = Path(f"{wd}host/samples/host_samples.fasta")
        # if host_file.exists():
        #     os.system(f"{wd}host/samples/host_samples.fasta")
        with open(f"{wd}host/samples/host_samples.fasta", "w") as w_fh:
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
            delayed(sample_sequences)(file, length, n_vir_samples, wd, virus, n_nuc_threshold) for file in
            glob.glob(f"{virus_dir}*.fna"))

    if full:
        new_records = Parallel(n_jobs=-1, verbose=True)(
            delayed(sample_sequences)(file, length, n_host_samples, wd, n_nuc_threshold, virus=False) for file in
            glob.glob(f"{host_dir}*.fna"))

        # for host
        for sublist in new_records:
            final_records.extend(sublist)
        with open(f"{wd}host/samples/host_samples.fasta", "w") as w_fh:
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
            delayed(sample_sequences)(file, length, n_vir_samples, wd, n_nuc_threshold, virus=True) for file in
            glob.glob(f"{virus_dir}*.fna"))

    # for virus, but second slower
    # for sublist in new_records:
    #     with open(f"X:/edwards2016/virus/samples/{sublist[0].id}_samples.fasta", "a") as w_fh:
    #         SeqIO.write(sublist, w_fh, "fasta")

    end = timer()
    runtime = end - start
    print(f"{Fore.GREEN}[random-sampling] Done in {runtime:.6f} seconds")

# if __name__ == "__main__":
#     main()
