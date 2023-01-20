import sys
import os
import numpy as np


def countLines(filen):
    cnt = 0
    with open(filen) as file:
        for line in file:
            cnt+=1
    return cnt

def getMinPval(filen):
    sol = 1.0
    with open(filen) as file:
        for line in file:
            lne = line.rstrip()
            tokens = line.split(",")
            if(tokens[0] == "p-value"):
                continue
            if float(tokens[0]) < sol:
                sol = float(tokens[0])
    return sol


n_samples = 3000
n_nodes = 75
n_edges = 100
n_snp_max = 10
seed = 0
ass_strengths = [0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
injectClique = 0

gen_data = "../utils/generate_subgraph_data_012.py"
gen_subgr = "../utils/get_subgraphs_from_network.py"
ours = "../../HOGImine-binary/hogimine_binary"


ours_genes = []
ours_edges = []
ours_sub3 = []
ours_sub4 = []

ours_genes_pv = []
ours_edges_pv = []
ours_sub3_pv = []
ours_sub4_pv = []

for ass in ass_strengths:
    name_base = "sim_n"+str(n_samples)+"_g"+str(n_nodes)+"e"+str(n_edges)+"s"+str(n_snp_max)+"a"+str(ass)+"_seed"+str(seed)
    print("Generating data", ass)
    os.system("python3 "+gen_data+" "+str(n_samples)+" "+str(n_nodes)+" "+str(n_edges)+" "+str(n_snp_max)+" "+str(ass)+" "+name_base+" "+str(seed)+" "+str(injectClique))
    os.system("python3 "+gen_subgr+" "+name_base+"_edges.txt "+name_base+"_genes.txt 1 0 > tmp.txt")
    os.system("python3 "+gen_subgr+" "+name_base+"_edges.txt "+name_base+"_sub3.txt 3 0 > tmp.txt")
    os.system("python3 "+gen_subgr+" "+name_base+"_edges.txt "+name_base+"_sub4.txt 4 0 > tmp.txt")


    print("Ours genes")
    os.system(ours+" -i "+name_base+"_X.txt -l "+name_base+"_Y.txt -c "+name_base+"_C.txt -s "+name_base+"_snpID.txt -m "+name_base+"_snp_map.txt -e "+name_base+"_genes.txt -f 0.05 -o out_"+name_base+"_ours_genes > out_"+name_base+"_ours_genes_log.txt")
    ours_genes.append(countLines("out_"+name_base+"_ours_genes_sign.txt"))
    ours_genes_pv.append(getMinPval("out_"+name_base+"_ours_genes_sign.txt"))
    print("Ours edges")
    os.system(ours+" -i "+name_base+"_X.txt -l "+name_base+"_Y.txt -c "+name_base+"_C.txt -s "+name_base+"_snpID.txt -m "+name_base+"_snp_map.txt -e "+name_base+"_edges.txt -f 0.05 -o out_"+name_base+"_ours_edges > out_"+name_base+"_ours_edges_log.txt")
    ours_edges.append(countLines("out_"+name_base+"_ours_edges_sign.txt"))
    ours_edges_pv.append(getMinPval("out_"+name_base+"_ours_edges_sign.txt"))
    print("Ours subgr3")
    os.system(ours+" -i "+name_base+"_X.txt -l "+name_base+"_Y.txt -c "+name_base+"_C.txt -s "+name_base+"_snpID.txt -m "+name_base+"_snp_map.txt -e "+name_base+"_sub3.txt -f 0.05 -o out_"+name_base+"_ours_sub3 > out_"+name_base+"_ours_sub3_log.txt")
    ours_sub3.append(countLines("out_"+name_base+"_ours_sub3_sign.txt"))
    ours_sub3_pv.append(getMinPval("out_"+name_base+"_ours_sub3_sign.txt"))
    print("Ours subgr4")
    os.system(ours+" -i "+name_base+"_X.txt -l "+name_base+"_Y.txt -c "+name_base+"_C.txt -s "+name_base+"_snpID.txt -m "+name_base+"_snp_map.txt -e "+name_base+"_sub4.txt -f 0.05 -o out_"+name_base+"_ours_sub4 > out_"+name_base+"_ours_sub4_log.txt")
    ours_sub4.append(countLines("out_"+name_base+"_ours_sub4_sign.txt"))
    ours_sub4_pv.append(getMinPval("out_"+name_base+"_ours_sub4_sign.txt"))


    os.system("rm "+name_base+"*")
    os.system("rm out_"+name_base+"*")
    os.system("rm tmp.txt")

print(n_samples, n_nodes, n_edges, n_snp_max)
print(ours_genes, ours_genes_pv)
print(ours_edges, ours_edges_pv)
print(ours_sub3, ours_sub3_pv)
print(ours_sub4, ours_sub4_pv)
