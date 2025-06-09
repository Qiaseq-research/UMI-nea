import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
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
pca = PCA(n_components=2)
pca_results = pca.fit_transform(tool_matrix)
colors = ["#DDA0DD","#E6C290","#9FD88E","#4775BF","#FF8888"]
plt.figure(figsize=(10, 6))
for i, tool in enumerate(tool_matrix.index):
    plt.scatter(pca_results[i, 0], pca_results[i, 1], s=100, label=tool, c=colors[i])
    
plt.xlabel("PC1 " + str((pca.explained_variance_ratio_[0]*100).round(2))+"%")
plt.ylabel("PC2 " + str((pca.explained_variance_ratio_[1]*100).round(2))+"%")
plt.legend()
plt.savefig("clonotype_umi_count_PCA.png", bbox_inches="tight")