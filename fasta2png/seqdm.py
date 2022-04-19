#import treepoem
import math

from Bio import SeqIO
import json
from pylibdmtx.pylibdmtx import encode
from PIL import Image
import dill
import glob
import os
from timeit import default_timer as timer

# seq_obj, = SeqIO.parse("X:/edwards2016/virus/fasta/NC_000866.fna", 'fasta')
# record = SeqIO.to_dict(seq_obj)
# #record['seq'] = str(record['seq'])
# b_record = json.dumps(record, default=vars).encode('utf-8')
#p_seq = dill.dumps(seq_obj)
# img = treepoem.generate_barcode(
#     barcode_type='datamatrix',
#     data=b_record
# )
# img.convert("1").save("NC_000866.png")
# dm = encode(b_record)
# img = Image.frombytes('RGB', (dm.width, dm.height), dm.pixels)
# img.save('NC_000866.png')

# NC_017186.fna - bac
# NC_000866.fna - vir

# loop version
# for file in glob.glob("X:/edwards2016/virus/fasta/*.fna"):
#     filename = os.path.basename(file)
#     seq_obj, = SeqIO.parse(file, 'fasta')
# # seq_obj, = SeqIO.parse("X:/edwards2016/caenorhabditis_elegans.fna", 'fasta')
#     p_seq = dill.dumps(seq_obj)
# # with open("NC_dilled.dill", 'wb') as f:
# #     dill.dump(seq_obj, f, protocol=dill.HIGHEST_PROTOCOL)
# # dim = int(math.sqrt(len(p_seq)))
#     img = Image.frombytes('P', (len(p_seq), 1), p_seq)
#     img.save(f"X:/edwards2016/virus/png-fasta/{filename}.png")

# nx1 size
# total_start = timer()
# seq_obj, = SeqIO.parse("X:/edwards2016/host/fasta/NC_017186.fna", 'fasta')
# # seq_obj, = SeqIO.parse("X:/edwards2016/caenorhabditis_elegans.fna", 'fasta')
# p_seq = dill.dumps(seq_obj)
# # with open("NC_dilled.dill", 'wb') as f:
# #     dill.dump(seq_obj, f, protocol=dill.HIGHEST_PROTOCOL)
# # dim = int(math.sqrt(len(p_seq)))
# img = Image.frombytes('P', (len(p_seq), 1), p_seq)
# img.save(f"NC_017186.png")
# total_end = timer()
# total_runtime = total_end - total_start
# print(f"Total elapsed time (SeqRecord save to png file): {total_runtime:.6f} seconds")

# square sized img
total_start = timer()
# seq_obj, = SeqIO.parse("X:/edwards2016/virus/fasta/NC_000866.fna", 'fasta')
seq_obj, = SeqIO.parse("X:/edwards2016/caenorhabditis_elegans.fna", 'fasta')
p_seq = dill.dumps(seq_obj)
# with open("NC_dilled.dill", 'wb') as f:
#     dill.dump(seq_obj, f, protocol=dill.HIGHEST_PROTOCOL)
# dim = int(math.sqrt(len(p_seq)))
comp = math.floor(math.sqrt(len(p_seq))) + 1
target = comp * comp
diff = target - len(p_seq)
n_seq = p_seq + bytes(diff)
dim = int(math.sqrt(len(n_seq)))
img = Image.frombytes('P', (dim, dim), n_seq)
img.save(f"caenorhabditis_elegans_sq.png")
total_end = timer()
total_runtime = total_end - total_start
print(f"Total elapsed time (SeqRecord save to png file): {total_runtime:.6f} seconds")