from Bio import SeqIO
from Bio import SeqRecord
import json
import base64
import dill
from timeit import default_timer as timer
import gzip
from PIL import Image
# treepoem._read_file("NC_000866.png")
# with open('NC_000866.png', 'rb') as f:
#     o = f.read()
# # j_obj = json.loads(o.decode('utf-8'))
# j_obj = base64.encodebytes(o).decode('utf-8')
# print(j_obj)

# print(json.load(o))
#

# record, = decode(Image.open('NC_000866.png'))
# print()

# tests on 10 mb NC_017186 - Amycolatopsis mediterranei

# second/first fastest to load, file about 70% smaller in size, quite fast to make (1.6 s)
total_start = timer()
r_img = Image.open("NC_017186.png")
obj = dill.loads(r_img.tobytes())
total_end = timer()
total_runtime = total_end - total_start
print(f"Total elapsed time (SeqRecord read from png file): {total_runtime:.6f} seconds")

# fastest to load, file about 70% smaller in size, slow to make (because of 'ultra' compression - about 10 s)
total_start = timer()
with gzip.open('NC_dilled.dill.gz', 'rb') as f:
    obj2 = dill.load(f)
total_end = timer()
total_runtime = total_end - total_start
print(f"Total elapsed time (SeqRecord read from compressed dill file): {total_runtime:.6f} seconds")

# second slowest to load, no space saving
total_start = timer()
seq_obj, = SeqIO.parse("NC_017186.fna", 'fasta')
total_end = timer()
total_runtime = total_end - total_start
print(f"Total elapsed time (biopython parse from normal fasta): {total_runtime:.6f} seconds")

# slowest to load, file about 70% smaller in size, slow to make (because of 'ultra' compression - about 10 s)
total_start = timer()
with gzip.open('NC_017186.fna.gz', 'rt') as f:
    seq_obj2, = SeqIO.parse(f, "fasta")
total_end = timer()
total_runtime = total_end - total_start
print(f"Total elapsed time (biopython parse from compressed fasta): {total_runtime:.6f} seconds")

print(obj)
print(obj2)
print(seq_obj)
print(seq_obj2)
# print(obj)
# seq_obj, = SeqIO.parse("X:/edwards2016/host/fasta/NC_017186.fna", 'fasta')
# print(seq_obj.seq == obj.seq)
