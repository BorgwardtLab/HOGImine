import sys
import os
import numpy as np


def countLines(filen):
    cnt = 0
    with open(filen) as file:
        for line in file:
            cnt+=1
    return cnt


def getSiniTime(filen):
    with open(filen) as file:
        lines = list(file)
        line = lines[0]
        tokens = line.rstrip().split()
        return float(tokens[2][0:-6])

def getOurTime(filen):
    with open(filen) as file:
        lines = list(file)
        line = lines[1]
        tokens = line.rstrip().split()
        return float(tokens[1])

n_sampless = [100, 300, 1000, 3000, 10000, 30000, 100000]
n_nodes = 75
n_edges = 100
n_snp_max = 10
seed = 0
ass = 0.7
injectClique = 0

gen_data = "../utils/generate_subgraph_data_012.py"
gen_subgr = "../utils/get_subgraphs_from_network.py"
ours = "../../HOGImine-binary/hogimine_binary"


ours_genes = []
ours_edges = []
ours_sub3 = []
ours_sub4 = []


for n_samples in n_sampless:
    name_base = "sim_n"+str(n_samples)+"_g"+str(n_nodes)+"e"+str(n_edges)+"s"+str(n_snp_max)+"a"+str(ass)+"_seed"+str(seed)
    print("Generating data", n_samples)
    os.system("python3 "+gen_data+" "+str(n_samples)+" "+str(n_nodes)+" "+str(n_edges)+" "+str(n_snp_max)+" "+str(ass)+" "+name_base+" "+str(seed)+" "+str(injectClique))
    os.system("python3 "+gen_subgr+" "+name_base+"_edges.txt "+name_base+"_genes.txt 1 0 > tmp.txt")
    os.system("python3 "+gen_subgr+" "+name_base+"_edges.txt "+name_base+"_sub3.txt 3 0 > tmp.txt")
    os.system("python3 "+gen_subgr+" "+name_base+"_edges.txt "+name_base+"_sub4.txt 4 0 > tmp.txt")

    print("Ours genes")
    os.system(ours+" -i "+name_base+"_X.txt -l "+name_base+"_Y.txt -c "+name_base+"_C.txt -s "+name_base+"_snpID.txt -m "+name_base+"_snp_map.txt -e "+name_base+"_genes.txt -f 0.05 -o out_"+name_base+"_ours_genes > out_"+name_base+"_ours_genes_log.txt")
    ours_genes.append(getOurTime("out_"+name_base+"_ours_genes_log.txt"))
    print("Ours edges")
    os.system(ours+" -i "+name_base+"_X.txt -l "+name_base+"_Y.txt -c "+name_base+"_C.txt -s "+name_base+"_snpID.txt -m "+name_base+"_snp_map.txt -e "+name_base+"_edges.txt -f 0.05 -o out_"+name_base+"_ours_edges > out_"+name_base+"_ours_edges_log.txt")
    ours_edges.append(getOurTime("out_"+name_base+"_ours_edges_log.txt"))
    print(ours_edges[-1])
    print("Ours subgr3")
    os.system(ours+" -i "+name_base+"_X.txt -l "+name_base+"_Y.txt -c "+name_base+"_C.txt -s "+name_base+"_snpID.txt -m "+name_base+"_snp_map.txt -e "+name_base+"_sub3.txt -f 0.05 -o out_"+name_base+"_ours_sub3 > out_"+name_base+"_ours_sub3_log.txt")
    ours_sub3.append(getOurTime("out_"+name_base+"_ours_sub3_log.txt"))
    print(ours_sub3[-1])
    if n_samples > 40000:
        ours_sub4.append(0)
    else:
        print("Ours subgr4")
        os.system(ours+" -i "+name_base+"_X.txt -l "+name_base+"_Y.txt -c "+name_base+"_C.txt -s "+name_base+"_snpID.txt -m "+name_base+"_snp_map.txt -e "+name_base+"_sub4.txt -f 0.05 -o out_"+name_base+"_ours_sub4 > out_"+name_base+"_ours_sub4_log.txt")
        ours_sub4.append(getOurTime("out_"+name_base+"_ours_sub4_log.txt"))
        print(ours_sub4[-1])

    os.system("rm "+name_base+"*")
    os.system("rm out*")
    os.system("rm tmp.txt")


print(ours_genes)
print(ours_edges)
print(ours_sub3)
print(ours_sub4)
