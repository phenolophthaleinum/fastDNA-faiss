import numpy as np
import faiss

dim = 10
database = np.loadtxt("/home/hyperscroll/fastDNA/edwards_vectors.txt", dtype="float32")
#n_database = database.shape[0]
# query = np.loadtxt("virus_vectors.txt", dtype="float32")
# n_query = query.shape[0]
index = faiss.IndexFlatL2(dim)
index.add(database)
faiss.write_index(index, "/home/hyperscroll/fastDNA/host_index.index")

