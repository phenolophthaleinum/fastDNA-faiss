import json


json_data = {}
with open("virus_rank-sample_100-len_250.json", "r") as fd:
    json_data = json.load(fd)

for virus in json_data:
    for res_tuple in json_data[virus]:
        if float(res_tuple[1]) > 50 and float(res_tuple[1]) < 100:
            print(f"{virus}: {res_tuple}")

