import numpy as np
import faiss

dim = 10
database = np.loadtxt("/home/hyperscroll/fastDNA/edwards_vectors_100_sample_100_len.txt", dtype="float32")
index = faiss.IndexFlatL2(dim)
index.add(database)
faiss.write_index(index, "/home/hyperscroll/fastDNA/host_index_100_sample_100_len.index")

