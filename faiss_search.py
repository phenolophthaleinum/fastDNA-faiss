from collections import defaultdict

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
classification = np.vectorize(map_data.get)(indices)
print(f"Classification: {classification}")

pre_rank = defaultdict(list)
for index, distance in np.stack((classification.flatten(), distances.flatten()), axis=1):
    pre_rank[index].append(2 - float(distance))
sum_rank = {key: sum(values) for key, values in pre_rank.items()}
final_rank = dict(sorted(sum_rank.items(), key=lambda elem: elem[1], reverse=True))
print(json.dumps(final_rank, indent=4))
