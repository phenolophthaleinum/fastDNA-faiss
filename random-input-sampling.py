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


def sample_sequences(file: str, length: int, n: int) -> List[SeqRecord]:
    """Function reads each given fasta file, randomly samples n subsequences of a defined length, creates new records with those subsequences and puts these new records into a list (which is finally returned).

    Args:
        file (str): Full path to a fasta file
        length (int): Length of a sampled subsequence
        n (int): Number of subsequences to be sampled

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
            seq=Seq(record.seq[pick:pick + length]),
            description=desc
        )
        new_records.append(new_record)
    return new_records


def main():
    # colorama
    init()

    # timer
    start = timer()

    # importing json data about hosts
    utils.get_host_data()

    # parallel sampling records of all files and dumping them into one file
    new_records = Parallel(n_jobs=-1, verbose=11, pre_dispatch='all', batch_size="auto", backend="loky")(
        delayed(sample_sequences)(file, 200, 10) for file in glob.glob("X:/edwards2016/host/fasta/*.fna"))
    final_records = []
    for sublist in new_records:
        final_records.extend(sublist)
    with open("X:/edwards2016/host/random_samples-training_fastDNA.fasta", "a") as w_fh:
        SeqIO.write(final_records, w_fh, "fasta")

    # mapping samples to nbci ids and dumping them into a file
    p_records = list(SeqIO.parse("X:/edwards2016/host/random_samples-training_fastDNA.fasta", "fasta"))
    ids_records = [record.id.split(".")[0] for record in p_records]
    d = {}
    for id in set(ids_records):
        keys = [index for index, value in enumerate(ids_records) if value == id]
        for key in keys:
            d[key] = id
    with open("X:/edwards2016/host/sample_map.json", "w", encoding='utf-8') as fh:
        json.dump(d, fh, indent=4)

    end = timer()
    runtime = end - start
    print(f"{Fore.GREEN} Done in {runtime:.6f} seconds")


if __name__ == "__main__":
    main()
