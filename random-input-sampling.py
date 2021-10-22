import json
import secrets
import glob
from timeit import default_timer as timer
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed
from contextlib import ExitStack
import multiprocessing
import itertools


def sample_sequences(file, length, n):
    #print(file)
    record = SeqIO.read(file, "fasta")
    #print(record)
    # if record.seq is None:
    #     print(f"{record.id} has none seq")
    #print(f"{record.id} : {len(record.seq)}")
    #samples = []
    #with open("X:/edwards2016/host/random_samples-training_fastDNA.fasta", "a") as w_fh:
    with open("X:/edwards2016/virus/random_samples_virus.fasta", "a") as w_fh:
        for i in range(n):
            pick = secrets.choice(range(0, len(record.seq) - length))
            #new_seq = record.seq[pick:pick + length]
        #     #samples.append(str(record.seq[pick:pick + length]))
            desc = f'{record.description} sample_{i} gen_pos:({pick}:{pick + length})'
            new_record = SeqRecord(
                id=record.id,
                seq=Seq(record.seq[pick:pick + length]),
                description=desc
            )
            SeqIO.write(new_record, w_fh, "fasta")
    #print(f"{record.id}: {samples}")


def main():
    start = timer()

    # importing json data about hosts
    with open("host.json", "r") as fh:
        host_data = json.load(fh)

    # debugging --------------------------------------------------------------------------------------------------------
    #print(os.scandir("X:/edwards2016/host/fasta/"))
    #print(len(glob.glob("X:/edwards2016/host/fasta/*.fna")))
    #print(glob.glob("X:/edwards2016/host/fasta/*.fna"))
    # debugging --------------------------------------------------------------------------------------------------------

    # parallel sampling records of all files and dumping them into one file
    ## X:/edwards2016/host/fasta/*.fna
    #par = Parallel(n_jobs=-2, pre_dispatch="all", batch_size="auto", backend="threading")(delayed(sample_sequences)(file, 200, 10) for file in glob.glob("X:/edwards2016/host/fasta/NC_001900.fna"))
    par = Parallel(n_jobs=-2, pre_dispatch="all", batch_size="auto", backend="threading")(
        delayed(sample_sequences)(file, 200, 10) for file in ["X:/edwards2016/virus/fasta/NC_001900.fna"])
    # with multiprocessing.Pool(processes=8) as pool:
    #     pool.starmap(sample_sequences, zip(glob.glob("X:/edwards2016/host/fasta/*.fna"), itertools.repeat(200), itertools.repeat(10)))


    # mapping samples to taxid and dumping them into a file
    p_records = list(SeqIO.parse("X:/edwards2016/host/random_samples-training_fastDNA.fasta", "fasta"))
    ids_records = [record.id.split(".")[0] for record in p_records]
    #print(ids_records)
    d = {}
    # for id in ids_records:
    #     index = ids_records.index(id)
    #     d[index] = id
    for id in set(ids_records):
        keys = [index for index, value in enumerate(ids_records) if value == id]
        for key in keys:
            d[key] = id
    #print(d)
    with open("X:/edwards2016/host/sample_map.json", "w", encoding='utf-8') as fh:
        json.dump(d, fh, indent=4)

    # potential map reading
    # with open("X:/edwards2016/host/sample_map.json", "r") as fh:
    #     map_data = {int(key): value for key, value in json.load(fh).items()}
    # map_data[10000]

    # debugging (can be skipped)---------------------------------------------------------------------------------------------

    #print(json.dumps(d, sort_keys=True))
    # debugging purposes
    from collections import Counter
    print(json.dumps(Counter(ids_records), indent=4))

    # with multiprocessing.Pool(processes=8) as pool:
    #     pool.starmap(sample_sequences, zip(glob.glob("X:/edwards2016/host/fasta/*.fna"), itertools.repeat(200), itertools.repeat(10)))

    # pool = multiprocessing.Pool(8)
    # pool.starmap(sample_sequences, zip(glob.glob("X:/edwards2016/host/fasta/*.fna"), itertools.repeat(200), itertools.repeat(10)))
    # pool.close()
    # pool.join()


    # for file in glob.glob("X:/edwards2016/host/fasta/*.fna"):
    #     record = SeqIO.read(file, "fasta")
    #     print(record.id)

    # debugging --------------------------------------------------------------------------------------------------------

    end = timer()
    runtime = end - start
    print(f"Done in {runtime:.6f}")


if __name__ == "__main__":
    main()
