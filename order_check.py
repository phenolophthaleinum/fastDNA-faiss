import json
import utils
import pandas as pd
from IPython.display import display
import plotly.express as px

virus_data = utils.get_virus_data()
with open("testdebugF_general_max.json", 'r') as fh:
    preds = json.load(fh)
res = {}
for virus in preds:
    target = virus_data[virus]["host"]["lineage_names"][3]
    try:
        search = [(*sublist, preds[virus].index(sublist) + 1) for sublist in preds[virus] if target in sublist][0]
    except:
        search = [None, None, None]
    # print(f"{virus}: {search}")
    res[virus] = search

# print(res['NC_024392'][2])
df = pd.DataFrame.from_dict(res)
display(df)
df_topX = df.iloc[[2]]
counts = df_topX.apply(pd.Series.value_counts, axis=1)
counts_percent = counts.apply(lambda x: (int(x) / 820) * 100)
fig = px.bar(counts, labels=dict(x="position", y="counts", color="position"))
fig.show()