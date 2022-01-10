import pandas as pd
import json
import utils

tax_data = utils.get_tax_data()
target_dir = '/home/hyperscroll/fastDNA-faiss/predict-prob_test/'
filename = 'predict_prob_out_none_60dim_epoch2_EColi'
level = ""

all_tax_levels = {
    "superkingdom": 0,
    "phylum": 1,
    "class": 2,
    "order": 3,
    "family": 4,
    "genus": 5,
    "species": 6
}


def tax2name(taxid):
    name = tax_data[str(taxid)]['organism_name']
    global level
    level = "species"
    return name


def tax2genus(taxid):
    name = tax_data[str(taxid)]['lineage_names'][all_tax_levels["genus"]]
    global level
    level = "genus"
    return name


def tax2family(taxid):
    name = tax_data[str(taxid)]['lineage_names'][all_tax_levels["family"]]
    global level
    level = "family"
    return name


def tax2order(taxid):
    name = tax_data[str(taxid)]['lineage_names'][all_tax_levels["order"]]
    global level
    level = "order"
    return name


def tax2class(taxid):
    name = tax_data[str(taxid)]['lineage_names'][all_tax_levels["class"]]
    global level
    level = "class"
    return name


with open("fastDNA/predict_prob_out_none_60dim_epoch2_n5.json") as jf:
    d = json.load(jf)
data = []
for elem in d:
    data.append(dict(sorted(elem.items(), key=lambda x: x[1], reverse=False)))
df = pd.DataFrame.from_records(data).fillna(0)
# df.rename(columns=lambda x: tax2name(x), inplace=True)
df.to_csv("predict_prob_rearranged_names.csv", sep=",")
compressions = ['infer', 'gzip', 'bz2', 'zip', 'xz', None]
for c in compressions:
    if c is None:
        df.to_pickle(f"predict_prob_rearranged.pkl", None)
        continue
    df.to_pickle(f"predict_prob_rearranged.{c}", c)
df_means = df.apply(lambda x: x.mean())
df_means.sort_values(ascending=False, inplace=True)
print(df_means)
ranks = df_means.to_json()
print(type(ranks))
ranks_parsed = json.loads(ranks)
print(type(ranks_parsed))
print(ranks_parsed.items())
json.dumps(ranks_parsed, indent=4)
df_means.to_json("predict_prob_rearranged_means.json", indent=4)


# data = pd.read_csv(f"{target_dir}{filename}.txt", header=None, sep=" ")
# tax_col = [col for col in data.columns if col % 2 == 0]
# for col in tax_col:
#     data[col] = data[col].apply(lambda x: tax2name(x))
# data_reduced = data[tax_col]
# print(data_reduced)
# # data_counts = data_reduced[0].value_counts(normalize=True)
# #data_reduced.reset_index(drop=True)
# data_counts = data_reduced.apply(lambda x: x.value_counts(normalize=True), axis=1).fillna(0)
# data_counts.to_csv(f"{target_dir}{filename}_{level}_counts.csv")
# data_means = data_counts.apply(lambda x: x.mean())
# data_means.to_csv(f"{target_dir}{filename}_{level}_means.csv")
# data.to_csv(f"{target_dir}{filename}_{level}.csv")

#data[tax_col] = data[tax_col].apply(lambda x: tax2name(x))
