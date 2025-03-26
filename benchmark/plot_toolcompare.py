import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.lines as lines
import math
import sys

infile = "./supplemental.table.1.txt"
df_b = pd.read_csv(infile,sep="\t")

#print (df_b[df_b.groupby(['num_founder','num_progeny', 'UMI_length', 'error_rate' ])[['V-measure']].idxmax().to_list()]["tool"] )

maxidx=df_b.groupby(['num_founder','num_progeny', 'UMI_length', 'error_rate' ])['V-measure'].idxmax()

print (df_b.loc[maxidx, ['tool', 'num_founder', 'num_progeny', 'UMI_length', 'error_rate', 'V-measure' ] ] )
#print (df_b.loc[[maxidx]['V-measure']].loc[:, ['tool', 'num_founder','num_progeny', 'UMI_length', 'error_rate', 'V-measure']] )

df_b_UMI_nea=df_b.loc[df_b["tool"]=="UMI-nea"]
print (abs(df_b_UMI_nea['perc_error']).mean())
print (abs(df_b_UMI_nea['perc_error']).quantile(0.25))
print (abs(df_b_UMI_nea['perc_error']).quantile(0.75))

#sys.exit()

print ("below is UMI-nea v-measure")
print (df_b.loc[df_b["tool"]=="UMI-nea"]['V-measure'].mean())
print (df_b.loc[df_b["tool"]=="UMI-nea"]['V-measure'].std())
print ("below is umi-tools v-measure")
print (df_b.loc[df_b["tool"]=="umi-tools"]['V-measure'].mean())
print (df_b.loc[df_b["tool"]=="umi-tools"]['V-measure'].std())
print ("below is UMIC-seq v-measure")
print (df_b.loc[df_b["tool"]=="UMIC-seq"]['V-measure'].mean())
print (df_b.loc[df_b["tool"]=="UMIC-seq"]['V-measure'].std())

def plot_bar(df_b,founder, metric):
	df_f=df_b[df_b['num_founder']==founder]
	mean_by_error=df_f.groupby(['tool','error_rate'])[metric].mean()
	print (mean_by_error)
	print(mean_by_error.loc[('UMI-nea')])
	std_by_error=df_f.groupby(['tool','error_rate'])[metric].std()
	print (std_by_error)
	q1_by_error=df_f.groupby(['tool','error_rate'])[metric].quantile(0.25)
	q3_by_error=df_f.groupby(['tool','error_rate'])[metric].quantile(0.75)
	iqr_by_error=q3_by_error-q1_by_error

	mean_by_umi=df_f.groupby(['tool','UMI_length'])[metric].mean()
	print (mean_by_umi)
	std_by_umi=df_f.groupby(['tool','UMI_length'])[metric].std()
	print (std_by_umi)
	q1_by_umi=df_f.groupby(['tool','UMI_length'])[metric].quantile(0.25)
	q3_by_umi=df_f.groupby(['tool','UMI_length'])[metric].quantile(0.75)
	iqr_by_umi=q3_by_umi-q1_by_umi

	df = pd.DataFrame()

	grp1 = 'UMI-tools'
	grp2 = 'UMIC-seq'
	grp3 = 'UMI-nea'

	df['label'] = ['error=0.1%','error=0.5%','error=3%', 'UMI=12bp','UMI=24bp','UMI=48bp']
	df[grp1] = pd.concat( [  mean_by_error.loc['umi-tools'],   mean_by_umi.loc['umi-tools'] ]  ).values.reshape(-1,).tolist()
	df[grp1+'_SD'] = pd.concat( [  std_by_error.loc['umi-tools'],   std_by_umi.loc['umi-tools'] ]  ).values.reshape(-1,).tolist()
	df[grp2] = pd.concat( [  mean_by_error.loc['UMIC-seq'],   mean_by_umi.loc['UMIC-seq'] ]  ).values.reshape(-1,).tolist()
	df[grp2+'_SD'] = pd.concat( [  std_by_error.loc['UMIC-seq'],   std_by_umi.loc['UMIC-seq'] ]  ).values.reshape(-1,).tolist()
	df[grp3] = pd.concat( [  mean_by_error.loc['UMI-nea'],   mean_by_umi.loc['UMI-nea'] ]  ).values.reshape(-1,).tolist()
	df[grp3+'_SD'] =  pd.concat( [  std_by_error.loc['UMI-nea'],   std_by_umi.loc['UMI-nea'] ]  ).values.reshape(-1,).tolist()

	ax = df.plot.bar(x='label',
                y=[grp1,grp2,grp3],
                yerr=df[[grp1+'_SD',grp2+'_SD',grp3+'_SD']].T.values,
		color={"UMI-tools": "red","UMIC-seq": "green", "UMI-nea": "blue"})
	if metric=='V-measure':
		plt.yticks([0.90, 0.92, 0.94, 0.96,0.98,1.00], fontsize=6)
		plt.ylim(0.9,1.09)
		ax.set_ylabel("V-measurement", fontsize=8)
	else:
		ax.set_ylabel("Time(sec)", fontsize=8)
		#plt.subplots_adjust(upper=0.5)
		plt.legend(fontsize=6,bbox_to_anchor=(1.15, 1), loc='upper center')
	plt.xticks(fontsize=6)
	ax.set_xlabel("error-rate/UMI-length", fontsize=8)
	plt.subplots_adjust(left=0.2, right=0.8, top=0.6, bottom=0.2)
	plt.savefig('supplemental.fig1.'+str(founder)+"."+metric+"."+'png', dpi=900, bbox_inches='tight')
	plt.show()

plot_bar(df_b,10000, 'V-measure')
plot_bar(df_b,100000, 'V-measure')
plot_bar(df_b,200000, 'V-measure')

plot_bar(df_b,10000, 'time')
plot_bar(df_b,100000, 'time')
plot_bar(df_b,200000, 'time')
