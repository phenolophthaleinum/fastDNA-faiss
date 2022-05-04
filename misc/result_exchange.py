import utils
import json
import pandas
import pandas as pd


with open("bayes_test_avg.json", 'r') as fh:
    res_json = json.load(fh)
virus_data = utils.get_virus_data()
hostname_data = utils.get_hostname_data()
new_d = {(virus_data[k]["organism_name"], k): (hostname_data[v[0][0]]["lineage_names"], v[0][1]) for k, v in res_json.items()}
df = pd.DataFrame.from_dict(new_d, orient="index")
df.to_excel("bayes_test_avg_exchanged.xlsx")

