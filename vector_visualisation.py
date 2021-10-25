import numpy as np
import pandas as pd
import json
from sklearn.decomposition import PCA
from plotly import express

database = np.loadtxt("edwards_vectors.txt", dtype="float32")
df = pd.DataFrame(data=database, columns=[f"dim{i}" for i in range(10)])

pca = PCA(n_components=3)
principals = pca.fit_transform(df)
principals_df = pd.DataFrame(data=principals, columns=["principal_x", "principal_y", "principal_z"])
with open("sample_map.json", "r") as fh:
    map_data = {int(key): value for key, value in json.load(fh).items()}
principals_df["ncbi_id"] = principals_df.index.map(map_data)
principals_df_grouped = principals_df.groupby("ncbi_id")
principals_df["sample_name"] = principals_df["ncbi_id"] + " sample_" + principals_df_grouped.cumcount().map(str)

print(principals_df)
print(principals_df["sample_name"])
#print(principals_df.index.tolist())
#principals_df.to_csv("edwards_vectors.csv", sep=",")

fig = express.scatter_3d(principals_df,
                         x="principal_x", y="principal_y", z="principal_z",
                         color="ncbi_id",
                         size=principals_df.index.tolist(),
                         size_max=100,
                         hover_name="sample_name",
                         hover_data=["principal_x",
                                     "principal_y",
                                     "principal_z",
                                     "ncbi_id"])
fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
fig.write_html("edwards_vectors.html")
