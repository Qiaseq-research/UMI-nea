import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.lines as lines
# match V*:CDR3:J
dfc = pd.read_csv("umi_count_7_clonotypes.txt",sep="\t")
dfc = dfc[dfc["number of UMI"] >= 1]
clono0 = np.unique(dfc["clonotype"].tolist())
clono_ordered = [5,6,3,7,4,2,1]
clono = []
for i in clono_ordered:
    clono.append(clono0[i-1])
plfm = ["NextSeq2000", "Miseq"]
dfc['platform'] = dfc['platform'].apply(lambda x: "NextSeq2000" if x == '230117_VH01211_6_AAC7YLMM5' else "Miseq")
cond = {'0.01% RNA':0,'0.1% RNA':1,'1% RNA':2,'10% RNA':3}
dfc["clonotype_full"] = dfc["clonotype"]
dfc["clonotype"] = dfc["clonotype"].replace(clono, [i for i in range(1,len(clono)+1)])
dfc["log10(number of UMI)"] = np.log10(dfc["number of UMI"])
dfc.to_csv("clonotype_umi_count.csv", index=False)
dfc["condition"] = dfc["condition"].apply(lambda x: cond[x])
conditions = list(cond.keys())
conditions.sort()
fig, axes = plt.subplots(2,len(clono), figsize = (5*len(clono),9), frameon=True)
sns.set_style('whitegrid')
for i in range(len(clono)):
    cl = i+1
    for j in range(len(plfm)):
        p = plfm[j]
        dfc1 = dfc[(dfc["clonotype"]==cl) & (dfc["platform"]==p)]
        sns.lineplot(ax=axes[j,i],x=dfc1["condition"],y=dfc1["log10(number of UMI)"],
                     hue=dfc1["replicate"],marker='o', palette=sns.color_palette("mako_r", 2))
        if i == len(clono)-1:
            lines, labels = axes[j,i].get_legend_handles_labels()
        cd1 = dfc1["condition"].tolist()
        lg1 = dfc1["log10(number of UMI)"].tolist()
        numi = dfc1["number of UMI"].tolist()
        r1 = dfc1["replicate"].tolist()
        for k in range(len(cd1)):
            if r1[k] == "rep1":
                v = 'top'
                h = 'left'
            else:
                v = 'bottom'
                h = 'right'
            axes[j,i].text(cd1[k], lg1[k], numi[k], va=v, ha=h, fontsize=16)
        axes[j,i].get_legend().remove()
        axes[j,i].set_ylim(-0.1,4)
        axes[j,i].set_yticks([0,1,2,3,4])
        axes[j,i].set_yticklabels([1,10,100,1000,10000], fontsize=12)
        axes[j,i].set_ylabel("number of UMI", fontsize=14)
        axes[j,i].set_xlim(-0.1,3.1)
        axes[j,i].set_xticks([0,1,2,3])
        axes[j,i].set_xticklabels(conditions, fontsize=14)
        axes[j,i].set_xlabel("condition", fontsize=14)
        axes[j,i].set_title(p+" clonotype "+str(cl), fontsize=16)
        if j==0:
            axes[j,i].set_xlabel("")
        if i>0:
            axes[j,i].set_ylabel("")
fig.tight_layout()
fig.legend(lines, labels, loc='upper right', bbox_to_anchor=(1.13,1), fontsize=14)
fig.suptitle("Reproducibility of UMI count for 7 spiked-in clonotypes",fontsize=25)
fig.subplots_adjust(top=0.85)
fig.savefig("clonotype_umi_count.png", bbox_inches="tight")