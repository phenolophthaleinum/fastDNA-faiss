import json


filename = 'preds_tax_test.json'
order = "Lactobacillales"
fh = open(filename)
db = json.load(fh)
fh.close()

score_dict = {}
for vir in db:
    if order in db[vir]:
        score_dict[vir] = True
    else:
        score_dict[vir] = False
score = sum(score_dict.values()) / len(score_dict)
print(score)