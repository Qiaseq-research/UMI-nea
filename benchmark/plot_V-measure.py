import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.lines as lines
import math
import sys

infile = "performance.csv"
df = pd.read_csv(infile,sep=",")
read_types = {12:"short_reads", 18:"short_reads", 25:"long_reads", 50:"long_reads"}
df["read_type"] = df["umi_len"].apply(lambda x: read_types[x])
nf = {1000:"1,000",10000:"10,000"}
df["num_founder"] = df["num_founder"].apply(lambda x: nf[x])
df_b = df[["num_founder","read_type","tool","V-measure"]]

tools_ord = ["UMI-nea","umi-tools","UMIc-seq","calib"]
df_b_gp = df_b.groupby(["num_founder","read_type","tool"])["V-measure"].mean()
df_b_gp = df_b_gp.reset_index()
df_b_gp['tool'] = pd.Categorical(df_b_gp['tool'],categories=tools_ord,ordered=True)
df_b_gp = df_b_gp.sort_values(by=["num_founder","read_type","tool"])
df_b_gp = df_b_gp.round(4)

bar_xcoord = [[-0.29,-0.08,0.11,0.31],[0.71,0.91,1.11,1.31]]
bar_ycoord = [[y for y in df_b_gp.loc[(df_b_gp["num_founder"]=="1,000") & (df_b_gp["read_type"]=="short_reads"),"V-measure"].tolist()],
             [y for y in df_b_gp.loc[(df_b_gp["num_founder"]=="10,000") & (df_b_gp["read_type"]=="short_reads"),"V-measure"].tolist()],
             [y for y in df_b_gp.loc[(df_b_gp["num_founder"]=="1,000") & (df_b_gp["read_type"]=="long_reads"),"V-measure"].tolist()],
             [y for y in df_b_gp.loc[(df_b_gp["num_founder"]=="10,000") & (df_b_gp["read_type"]=="long_reads"),"V-measure"].tolist()]]

sns.set_style('whitegrid')
colors = sns.color_palette(["#4775BF","#9FD88E","#E6C290","#DDA0DD"])
rtypes = ["short_reads","long_reads"]
fig, axes = plt.subplots(1,2, figsize = (25,10), frameon=True)
for i in range(len(rtypes)):
    rt = rtypes[i]
    df1 = df_b[df_b["read_type"]==rt]
    sns.barplot(ax=axes[i],data=df1, x="num_founder", y="V-measure", 
                hue="tool",hue_order=tools_ord,
                ci="sd", palette=colors)
    for j in range(len(bar_xcoord)):
        for k in range(len(bar_xcoord[j])):
            x_loc = bar_xcoord[j][k]
            y_loc = bar_ycoord[i*2+j][k]
            t = str(bar_ycoord[i*2+j][k])
            axes[i].text(x_loc, y_loc, t, fontsize=18)
    if i < len(rtypes)-1:
        lines, labels = axes[i].get_legend_handles_labels()
    if i > 0:
        axes[i].set_ylabel("")
    else:
        axes[i].set_ylabel("V-measure", fontsize=25)
    axes[i].get_legend().remove()
    axes[i].set(ylim=(0.86,1.005))
    axes[i].set_xlabel(rt,fontsize=25)
    axes[i].tick_params(axis='both', which='major', labelsize=20)
fig.tight_layout()
fig.legend(lines, labels, loc='upper right', bbox_to_anchor=(1.09,1), fontsize=20)
fig.savefig("simulation_Vmeasure_4_tools.png", bbox_inches="tight", dpi=500)