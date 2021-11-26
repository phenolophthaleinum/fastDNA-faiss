import numpy as np
import faiss


def build_index(vectors_path, dim, n_samples, index_path):
    #dim = 10
    #database = np.loadtxt("/home/hyperscroll/fastDNA/edwards_vectors_100_sample_250_len.txt", dtype="float32")
    #database = np.loadtxt(vectors_path, dtype="float32")
    database = np.fromfile(vectors_path, dtype="float32")
    print(database.size)
    print(n_samples)
    print(dim)
    n_host = database.size / (n_samples * dim)
    print(n_host)
    database = database.reshape(n_samples * int(n_host), dim)
    index = faiss.IndexFlatL2(np.shape(database)[1])
    index.add(database)
    #faiss.write_index(index, "/home/hyperscroll/fastDNA/host_index_100_sample_250_len.index")
    faiss.write_index(index, index_path)

