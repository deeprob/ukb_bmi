filename = "/data5/deepro/sw/annovar/humandb/hg38_dbnsfp42c.txt"

with open(filename, "r") as f:
    data = f.readlines()

# change column names
data[0] = "\t".join([val.replace("-", "_").replace("+", "") for val in data[0].split()])
data[0] += "\n"

with open(filename, "w") as f:
    f.writelines(data)
