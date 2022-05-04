import math
import dill
from Bio import SeqIO
from PIL import Image
from timeit import default_timer as timer

total_start = timer()

# assuming single record for now
seq_obj, = SeqIO.parse("NC_017186.fna", 'fasta')
# object dump
p_seq = dill.dumps(seq_obj)
# closest square root
dim = math.floor(math.sqrt(len(p_seq))) + 1
# final size
target = dim * dim
# how much needed to fill up the target
diff = target - len(p_seq)
# final data
n_seq = p_seq + bytes(diff)
# save image from bytes with P image mode - 8-bit pixels, mapped to any other mode using a color palette - most space saved without disrupting the data
img = Image.frombytes('P', (dim, dim), n_seq)
img.save(f"NC_017186_seqrec2png.png")

# timing
total_end = timer()
total_runtime = total_end - total_start
print(f"Total elapsed time (SeqRecord save to png file): {total_runtime:.6f} seconds")