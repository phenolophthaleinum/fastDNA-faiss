import utils

host_data = utils.get_host_data()

with open("D:/edwards2016/host/random_family-training_labels_test.txt", "r") as fd:
    labels = fd.read().splitlines()

with open("D:/edwards2016/host/random_family-training_ncbi_test.txt", "w") as f:
    for label in labels:
        print(label)
        for host in host_data:
            if host_data[host]["taxid"] == label:
                print(host_data[host]["version"])
                f.write(host_data[host]["version"].split(".")[0] + " " + host_data[host]["organism_name"] + "\n")
