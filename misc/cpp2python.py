from subprocess import Popen, PIPE
import json

output = Popen(["/home/hyperscroll/fastDNA-faiss/fastDNA/fastdna", "predict-prob", "/home/hyperscroll/edwards2016/models/debugF_test_bacillales_premodel/random_model-order-dim_10-len_125-epoch1.bin", "NC_001884.1_samples.fasta", "10"], stdout=PIPE)
out_by = output.communicate()[0]
out_js = json.loads(out_by)
print(json.dumps(out_js[0], indent=4))
