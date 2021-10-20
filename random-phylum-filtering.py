import json
import secrets
from collections import defaultdict
from timeit import default_timer as timer
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed

# timeit
start = timer()

# importing json data about hosts
with open("host.json", "r") as fh:
    host_data = json.load(fh)

# filtering all hosts by phylum
phylum_host = defaultdict(list)
for host in host_data:
    phylum = host_data[host]["lineage_names"][1]
    phylum_host[phylum].append(host)

# random sampling of a single host from a phylum
random_phylum_host = {phylum: secrets.choice(phylum_host[phylum]) for phylum in phylum_host}
print(random_phylum_host)

# list of filenames to be used
filenames = list(random_phylum_host.values())
print(len(filenames))

# get list of taxid (labels in fastDNA)
labels = [host_data[host]["taxid"] for host in filenames]
print(labels)

# parse selected files and merge them into single fasta file; create labels file
records = [list(SeqIO.parse(f"D:/praktyki2020/edwards2016/host/fasta/{file}.fna", "fasta"))[0] for file in filenames]
print(len(records))
print(records)

with open("D:/praktyki2020/edwards2016/host/random_phylum-training_labels.txt", "w") as fh:
    for label in labels:
        fh.write(label + "\n")

#par = Parallel(n_jobs=-1)(delayed(SeqIO.write)(records, "D:/praktyki2020/edwards2016/host/random_phylum-training_fastDNA.fasta", "fasta") for record in records)
SeqIO.write(records, "D:/praktyki2020/edwards2016/host/random_phylum-training_fastDNA.fasta", "fasta")
print(f"Written {len(records)} records.")

end = timer()
runtime = end - start
print(f"Done in {runtime:.6f}")
