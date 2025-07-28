import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

infile = "umi_count_7_clonotypes.csv"
df = pd.read_csv(infile)
df['feature'] = (
    df['platform'] + '_' + 
    df['condition'].astype(str) + '_' + 
    df['replicate'] + '_' +
    df['clonotype'].astype(str)
)
tools_ord = ["Calib","UMIc-seq","UMI-tools","UMI-nea", "MiXCR"]
df['tools'] = pd.Categorical(df['tools'],categories=tools_ord,ordered=True)
df = df.sort_values(by=["tools"])
tool_matrix = df.pivot(index='tools', columns='feature', values='number of UMI')
corr_matrix = tool_matrix.T.corr()
plt.figure(figsize=(10, 8))
g = sns.clustermap(
    corr_matrix,
    cmap='coolwarm',
    annot=True,
    fmt=".3f",
    annot_kws={"size": 15},
    vmin=0.84,
    vmax=1,
    linewidths=1,
    figsize=(10, 9)
)
g.ax_heatmap.set_xlabel("Tools", fontsize=15)
g.ax_heatmap.set_ylabel("Tools", fontsize=15)
g.ax_heatmap.set_xticklabels(
    g.ax_heatmap.get_xticklabels(), 
    fontsize=14)
g.ax_heatmap.set_yticklabels(
    g.ax_heatmap.get_yticklabels(), 
    fontsize=14)
cb = g.cax
cb.tick_params(labelsize=14)
plt.savefig("clonotype_umi_count_linkage.png", bbox_inches="tight", dpi=350)