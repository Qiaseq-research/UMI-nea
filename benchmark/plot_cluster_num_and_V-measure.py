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
df_cluster = df[["num_founder","read_type","tool","dedup_umi_cluster"]]
df_est = df.loc[df["tool"]=="UMI-nea", ["num_founder","read_type","tool","estimated_molecule"]]
df_est["tool"] = "UMI-nea_est_mol"
df_est["dedup_umi_cluster"] = df_est["estimated_molecule"]
df_est = df_est[["num_founder","read_type","tool","dedup_umi_cluster"]]
df_a = pd.concat([df_cluster, df_est], axis=0)
df_a_gp = df_a.groupby(["num_founder","read_type","tool"])["dedup_umi_cluster"].mean()
df_a_gp = df_a_gp.reset_index()
df_a_gp["dedup_umi_cluster"] = df_a_gp["dedup_umi_cluster"].astype(int)
df_a_gp = df_a_gp.sort_values(by=["num_founder","read_type","tool"])

bar_xcoord = [[-0.29,-0.08,0.11,0.31],[0.71,0.91,1.11,1.31]]
bar_ycoord = [[y for y in df_a_gp.loc[(df_a_gp["num_founder"]==1000) & (df_a_gp["read_type"]=="short_reads"),"dedup_umi_cluster"].tolist()],
             [y for y in df_a_gp.loc[(df_a_gp["num_founder"]==1000) & (df_a_gp["read_type"]=="long_reads"),"dedup_umi_cluster"].tolist()],
             [y for y in df_a_gp.loc[(df_a_gp["num_founder"]==10000) & (df_a_gp["read_type"]=="short_reads"),"dedup_umi_cluster"].tolist()],
             [y for y in df_a_gp.loc[(df_a_gp["num_founder"]==10000) & (df_a_gp["read_type"]=="long_reads"),"dedup_umi_cluster"].tolist()]]
fs=[1000, 10000]

sns.set_style('whitegrid')
colors = sns.color_palette(["#4775BF", "#719DE5","#9FD88E","#E6C290"])
fig, axes = plt.subplots(1,2, figsize = (25,10), frameon=True)
for i in range(len(fs)):
    num_f = fs[i]
    df1 = df_a[df_a["num_founder"]==num_f]
    sns.barplot(ax=axes[i],data=df1, x="read_type", y="dedup_umi_cluster", 
                hue="tool", hue_order=["UMI-nea", "UMI-nea_est_mol", "UMIC-seq", "umi-tools"],
                ci=85, palette=colors)
    for j in range(len(bar_xcoord)):
        for k in range(len(bar_xcoord[j])):
            x_loc = bar_xcoord[j][k]
            y_loc = min(bar_ycoord[i*2+j][k]+num_f*0.005,num_f+num_f*0.092)
            t = str(bar_ycoord[i*2+j][k])
            axes[i].text(x_loc, y_loc, t, fontsize=18)
    if i < len(fs)-1:
        lines, labels = axes[i].get_legend_handles_labels()
    if i > 0:
        axes[i].set_ylabel("")
    else:
        axes[i].set_ylabel("Number of clusters", fontsize=25)
    axes[i].get_legend().remove()
    axes[i].set(ylim=(num_f-num_f*0.1, num_f+num_f*0.1))
    axes[i].set_xlabel("read types",fontsize=25)
    axes[i].tick_params(axis='both', which='major', labelsize=20)
fig.tight_layout()
fig.legend(lines, labels, loc='upper right', bbox_to_anchor=(1.135,1), fontsize=20)
fig.savefig("simulation_cluster_num.png", bbox_inches="tight")


df_b = df[["num_founder","read_type","tool","V-measure"]]
df_b_gp = df_b.groupby(["num_founder","read_type","tool"])["V-measure"].mean()
df_b_gp = df_b_gp.reset_index()
df_b_gp = df_b_gp.sort_values(by=["num_founder","read_type","tool"])
df_b_gp = df_b_gp.round(4)
bar_xcoord = [[-0.25,0.01,0.28],[0.74,1.01,1.28]]
bar_ycoord = [[y for y in df_b_gp.loc[(df_b_gp["num_founder"]==1000) & (df_b_gp["read_type"]=="short_reads"),"V-measure"].tolist()],
             [y for y in df_b_gp.loc[(df_b_gp["num_founder"]==1000) & (df_b_gp["read_type"]=="long_reads"),"V-measure"].tolist()],
             [y for y in df_b_gp.loc[(df_b_gp["num_founder"]==10000) & (df_b_gp["read_type"]=="short_reads"),"V-measure"].tolist()],
             [y for y in df_b_gp.loc[(df_b_gp["num_founder"]==10000) & (df_b_gp["read_type"]=="long_reads"),"V-measure"].tolist()]]
fs=[1000, 10000]

colors = sns.color_palette(["#4775BF","#9FD88E","#E6C290"])
fig, axes = plt.subplots(1,2, figsize = (25,10), frameon=True)
for i in range(len(fs)):
    num_f = fs[i]
    df1 = df_b[df_b["num_founder"]==num_f]
    sns.barplot(ax=axes[i],data=df1, x="read_type", y="V-measure", 
                hue="tool",hue_order=["UMI-nea", "UMIC-seq", "umi-tools"],
                ci=85, palette=colors)
    for j in range(len(bar_xcoord)):
        for k in range(len(bar_xcoord[j])):
            x_loc = bar_xcoord[j][k]
            y_loc = bar_ycoord[i*2+j][k]
            t = str(bar_ycoord[i*2+j][k])
            axes[i].text(x_loc, y_loc, t, fontsize=18)
    if i < len(fs)-1:
        lines, labels = axes[i].get_legend_handles_labels()
    if i > 0:
        axes[i].set_ylabel("")
    else:
        axes[i].set_ylabel("V-measure", fontsize=25)
    axes[i].get_legend().remove()
    axes[i].set(ylim=(0.93,1.005))
    axes[i].set_xlabel("read types",fontsize=25)
    axes[i].tick_params(axis='both', which='major', labelsize=20)
fig.tight_layout()
fig.legend(lines, labels, loc='upper right', bbox_to_anchor=(1.09,1), fontsize=20)
fig.savefig("simulation_Vmeasure.png", bbox_inches="tight")