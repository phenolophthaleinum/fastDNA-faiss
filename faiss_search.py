import numpy as np
import faiss
import json

dim = 10
k_nearest = 10
query = np.loadtxt("/home/hyperscroll/fastDNA/virus_vectors.txt", dtype="float32")
index = faiss.read_index("/home/hyperscroll/fastDNA/host_index.index")
with open("/home/hyperscroll/edwards2016/host/sample_map.json", "r") as fh:
    map_data = {int(key): value for key, value in json.load(fh).items()}
distances, indices = index.search(query, k_nearest)
print(f"Nearest sequences: {indices}")
print(f"Classification: {np.vectorize(map_data.get)(indices)}")
