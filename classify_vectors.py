import utils
import pandas as pd


hostvir = utils.get_hostvir_data()

viruses = []
for host in hostvir:
    viruses.extend(hostvir[host]["virus_id"])
virus_set = set(viruses)
df = pd.DataFrame(columns=list(virus_set))
for host in hostvir:
    if not hostvir[host]["virus_id"]:
        df.loc[host] = 0
    for id in hostvir[host]["virus_id"]:
        df.loc[host, id] = 1
df = df.fillna(0)
print(df)
df.to_pickle("v-h_classify.pkl")
