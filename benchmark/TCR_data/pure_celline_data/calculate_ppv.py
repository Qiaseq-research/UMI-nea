import pandas as pd
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

df = pd.read_csv(infile,sep="\t")
df = df.groupby(by=["read set","is_truth","input level","cellline"]).sum().reset_index()
df["cellline"] = df["cellline"].astype(int).astype(str)
df1 = pd.pivot_table(df, values="number of UMI", index=["read set","input level","cellline"], columns=["is_truth"]).reset_index()
df1["Positive predictive value"] = df1[True]/(df1[True]+df1[False])
df1 = df1.round(3)
cols = ["read set", "input level", "cellline", "False postive", "True positive", "Positive predictive value"]
df1.columns = cols
df1 = df1[["cellline", "input level", "False postive", "True positive", "Positive predictive value"]]
df1["input level"] = df1["input level"].replace(["high","medium","low"], [1,2,3])
df1 = df1.sort_values(by=['cellline', 'input level'])
df1["input level"] = df1["input level"].replace([1,2,3], ["high","medium","low"])
df1.to_csv(outfile, index=False, header=True, sep="\t")
