import secrets
from collections import defaultdict

import numpy as np
import pandas as pd
import json
from sklearn.decomposition import PCA
from plotly import express
from timeit import default_timer as timer
import umap


import utils

start = timer()

database = np.loadtxt("edwards_vectors.txt", dtype="float32")
df = pd.DataFrame(data=database, columns=[f"dim{i}" for i in range(10)])

#pca = PCA(n_components=3)
reducer = umap.UMAP(n_components=3)
principals= reducer.fit_transform(df)
#principals = pca.fit_transform(df)
principals_df = pd.DataFrame(data=principals, columns=["principal_x", "principal_y", "principal_z"])
with open("sample_map.json", "r") as fh:
    map_data = {int(key): value for key, value in json.load(fh).items()}
principals_df["ncbi_id"] = principals_df.index.map(map_data)
principals_df_grouped = principals_df.groupby("ncbi_id")
principals_df["sample_name"] = principals_df["ncbi_id"] + " sample_" + principals_df_grouped.cumcount().map(str)

host_data = utils.get_host_data()
for key in host_data:
    principals_df.loc[principals_df["ncbi_id"] == key, "full_name"] = host_data[key]["organism_name"]

# filtering to reduce amout of records: take random genome from each phylum
phylum_host = defaultdict(list)
for host in host_data:
    phylum = host_data[host]["lineage_names"][1]
    phylum_host[phylum].append(host)
random_phylum_host = [secrets.choice(phylum_host[phylum]) for phylum in phylum_host]
print(random_phylum_host)

principals_df_filtered_phylum = principals_df[principals_df["ncbi_id"].isin(random_phylum_host)]
print(principals_df_filtered_phylum.index)

#print(principals_df)

# principals_df_sampled = principals_df.sample(frac=0.05)
# principals_df_sampled_grouped = principals_df_sampled.groupby(["ncbi_id"]).apply(lambda x: x.sort_values(["ncbi_id"])).reset_index(drop=True)
# print(principals_df_sampled_grouped)

end = timer()
runtime = end - start
print(f"Done in {runtime:.6f} seconds")
#print(principals_df["sample_name"])
#print(principals_df.index.tolist())

# saving to csv
#principals_df.to_csv("edwards_vectors.csv", sep=",", float_format="%.6f")
principals_df_filtered_phylum.to_csv("edwards_vectors_phylum_filtered.csv", sep=",", float_format="%.6f")


# fig = express.scatter_3d(principals_df,
#                          x="principal_x", y="principal_y", z="principal_z",
#                          color="ncbi_id",
#                          size=principals_df.index.tolist(),
#                          size_max=5,
#                          hover_name="sample_name",
#                          hover_data=["principal_x",
#                                      "principal_y",
#                                      "principal_z",
#                                      "full_name",
#                                      "ncbi_id"])
# fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
# fig.write_html("edwards_vectors_sampled.html")
