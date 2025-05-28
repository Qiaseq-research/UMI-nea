import sys
import re
from random import randint, random, shuffle
import numpy as np
from skbio.alignment import StripedSmithWaterman
from multiprocessing import Pool
from itertools import repeat, chain
import pandas as pd
from scipy.stats import poisson

def sample_umi(l):
    umi = ""
    bases = "GATC"
    for _ in range(l):
        umi += bases[randint(0,3)]
    return umi

def generate_Mut(c, do_indel, ratio):
    # ratio = [(ins),n1(del),n2(sub)]
    if do_indel:
        bases = ["G"+c, "A"+c, "T"+c, "C"+c] * 3 * ratio[0]
        bases += [''] * 12 * ratio[1]
        for b in re.findall('[A-Z]', "GATC".replace(c, '')):
            bases += [b] * 4 * ratio[2]
        #print(len(bases),bases)
    else:
        bases = re.findall('[A-Z]', "GATC".replace(c, ''))
    return bases[randint(0,len(bases)-1)]

def generate_Mut_old(c, do_indel):
    if do_indel:
        bases = "G_A_T_C".replace(c, '').split("_")
        bases += ["G"+c, "A"+c, "T"+c, "C"+c]
    else:
        bases = re.findall('[A-Z]', "GATC".replace(c, ''))
    return bases[randint(0,len(bases)-1)]

def find_nb_para_k(mu, k):
    var = mu + (mu**2)/k
    print(var)
    p = mu/var
    n = (mu*p)/(1-p)
    return n,p

def read_offN(file):
    df = pd.read_csv(file)
    return df["offN"].tolist()

def find_nb_para(mu, var):
    #var = sig**2
    #print(var)
    p = mu/var
    n = (mu*p)/(1-p)
    return n,p

def find_nb_para_p(mu, p):
    n = (mu*p)/(1-p)
    return n,p

def sim_UMI(simf, umi_len, acc, pN, oN, do_indel, mut_ratio, p1):
    umi_lengthses = [int(umi_len)]
    off_accuracieses = [float(acc)]
    error_rates = [1-e for e in off_accuracieses]
    tripses = [1]

    for umi_length in umi_lengthses:
        print("\n\n~~~~~~~~~~~~~~~~~~~~~~")
        print("NOW WORKING ON UMI LENGTH %d" % umi_length)
        for trip in tripses:
            for off_accuracy in off_accuracieses:
                #Set parameters
                par_N = int(pN)   #Number of parent sequences
                off_N = int(oN)      #Average number of offspring sequences (mean of normal distribution)
                off_N_scale = 1   #Standard deviation for normal distrubution for offspring number
                mu,var = off_N, p1*off_N
                p = mu/var
                r = (mu*mu)/(var-mu)
                name = "%s_%s_ul%d_acc%.3f_%d" % (par_N, off_N, umi_length, off_accuracy, trip) #For plotting
                f = open(simf+".out", "w")
                f1 = open(simf+".truth.labels", "w")
                print("\nWorking %s" % name)
                #N_par random UMI DNA sequences of length ul
                parents_1 = set()
                for i in range(int(par_N*2)):
                    parents_1.add(sample_umi(umi_length))
                    if len(parents_1) >= par_N:
                        break
                parents = list(parents_1)
                #parents = []
                #for i in range(par_N):
                #    parents.append(sample_umi(umi_length))
                #Generate offspring/"reads"
                off_randombase = 1 - off_accuracy
                reads = []    #List of sequenced UMIs/offspring
                reads_err = [] #sanity check
                read_parent_id = []  #List of sequence parents (true clusters)
                read_parent_label = []
                num_sig = 0
                for i in range(par_N):
                    rn_off = np.random.negative_binomial(r,p)
                    if rn_off == 0:
                        num_sig +=1
                        reads_err.append("sig")
                    else:
                        reads_err.append("")
                    reads.append(parents[i])
                    read_parent_id.append(parents[i])
                    read_parent_label.append(i)
                    for j in range(rn_off):
                        kid = ""
                        kid_err = ""
                        k_idx = 0
                        for c in parents[i]:
                            if random() < off_randombase:
                                c1 = generate_Mut(c, do_indel, mut_ratio)
                                kid += c1
                                if kid_err == "":
                                    kid_err += (str(k_idx)+":"+c+"/"+c1)
                                else:
                                    kid_err += ("_"+str(k_idx)+":"+c+"/"+c1)
                                #print(kid_err)
                            else:
                                kid += c
                            k_idx += 1
                        reads.append(kid)        #"Sequenced" UMIs: off_N times par_N read with error rate as defined
                        reads_err.append(kid_err)
                        read_parent_id.append(parents[i]) #Parent of each of the reads
                        read_parent_label.append(i)
                print("Number of singleton is", num_sig)
                print("%d total sequence reads generated" % len(reads))
                for idx in range(len(reads)):
                    f.write(reads[idx]+"\t"+str(read_parent_id[idx])+"\t"+reads_err[idx]+"\n")
                    f1.write(reads[idx]+" "+str(read_parent_label[idx])+"\n")
                f.close()
                f1.close()

simf = sys.argv[1]
umi_len = sys.argv[2]
acc = sys.argv[3]
pN = sys.argv[4]
oN = sys.argv[5]
do_indel = bool(int(sys.argv[6]))
mut_ratio = [int(x) for x in sys.argv[7].split("-")]
p1 = float(sys.argv[8])
#generate_Mut("A", do_indel, mut_ratio)
sim_UMI(simf, umi_len, acc, pN, oN, do_indel, mut_ratio, p1)
