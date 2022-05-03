from Bio import SeqIO
import itertools as it
import numpy as np
import dill
from PIL import Image
import math


perm12 = it.product("ACTG", repeat=12)
perm12_l = list(perm12)
base12 = [''.join(item) for item in perm12_l]
rgb = [(x,y,z) for x in range(256) for y in  range(256) for z in range(256)]
d = {c: kmer for c, kmer in zip(rgb, base12)}
d_rev = {kmer: c for c, kmer in zip(rgb, base12)}
record, = SeqIO.parse("NC_022116.fna", 'fasta')
s = str(record.seq)
n = len(s)
comp = math.floor(math.sqrt(n)) + 1
target = comp * comp
diff = target - n
iternum = n - (n % 12)
final_diff = diff + (n % 12)
s_coded = []
for i in range(0, iternum, 12):
    s_coded.append(d_rev[s[i:i + 12]])
s_coded2 = [item for tup in s_coded for item in tup]
s_coded_b = bytes(s_coded2)
# coded_final = s_coded_b + bytes(final_diff)
# dim = int(math.sqrt(len(coded_final)))
## from divisors its 1131x2265
print(len(s_coded_b))
# s_arr = np.array(s_coded)
# img = Image.frombytes('P', (1131, 2265), s_coded_b)
# img_conv = img.convert('RGB')
# img_arr = np.array(img_conv)
# img2 = Image.fromarray(s_arr)
# with open("NC_022116.sqz", 'wb') as f:
#     dill.dump(s_coded_b, f, protocol=dill.HIGHEST_PROTOCOL)
img3 = Image.new('RGB', (1131, 755))
img3.putdata(s_coded)
img3.save(f"NC_022116_sqz_c.png")
