import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.lines as lines
import math

infile = "./summary.sort.thread12-48.txt"
df = pd.read_csv(infile,sep="\t",header=None)
df.columns = ["replicates","founder","progeny","UMI length","thread","error rate","runtime"]
outpng = "./runtime.benchmark.png"
df = df.round()
err_list = ["0.5%","5%"]
len_list = [12,50]
m_list = [50,2000]
colors = [['b'],['r']]
se_list = [(200,650),(5000,22000)]
fig, axes = plt.subplots(1,2, figsize = (11,8), frameon=True)
sns.set_style('whitegrid')
for i in range(len(err_list)):
    err = err_list[i]
    l = len_list[i]
    df1 = df[df["UMI length"]==l]
    m = m_list[i]
    #sns.lineplot(ax=axes[i],x=df1["thread"],y=df1["runtime"],err_style="bars")
    #sns.scatterplot(ax=axes[i],x=df1["thread"],y=df1["runtime"],color='k')
    sns.lineplot(ax=axes[i], data=df1, x=df1["thread"],y=df1["runtime"], errorbar=(("se", 3) ), hue="error rate", palette=colors[i], legend=None )
    axes[i].legend(title="UMI length="+str(l) + "bp\n\nError rate:", loc=1, labels=[err_list[i]],  handlelength=5, fontsize='x-large', title_fontsize='x-large' )
    ys = se_list[i][0]
    ye = se_list[i][1]
    yt = list(range(ys,ye,m_list[i] ) )
    print (yt)
    axes[i].set_ylim(ys,ye)
    axes[i].set_yticks(yt)
    axes[i].set_yticklabels(list(yt), fontsize=12)
    axes[i].set_xticks(list(range(12,49,12)))
    axes[i].set_xticklabels(list(range(12,49,12)),fontsize=14)
    #axes[i].set_xlabel("# of Threads"+"\n"+"UMI="+str(l)+"bp, "+"Error rate="+str(err*100)+"%",fontsize=18)
    axes[i].set_xlabel("# of Threads", fontsize=14)
    if i == 1:
        axes[i].yaxis.tick_right()
        axes[i].set_ylabel("")
    else:
        axes[i].set_ylabel("run-time in sec",fontsize=18)
fig.tight_layout(pad=-0.55)
fig.savefig(outpng,bbox_inches="tight")
