from collections import defaultdict
import logging
import random
import os
import sys

import numpy as np
import pandas as pd

#from genetics import network, genes
#from xwelt_utils import simulation, tests, utils


logging.basicConfig(level='INFO', format='.. %(message)s')


def get_rand_subgraph(genes, edges, n_hops):
    subgr = set()
    adj = {}
    for ge in genes:
        adj[ge] = []
    for ed in edges:
        a, b = ed
        adj[a].append(b)
        adj[b].append(a)
    # do a random walk of n_hops hops
    u = random.choice(genes)
    while len(adj[u])==0:
      u = random.choice(genes)
    subgr.add(u)
    for i in range(n_hops):
        u = random.choice(adj[u])
        subgr.add(u)

    return sorted(list(subgr))

def get_intervals_from_subgraph(subgr, min_leng, max_leng, gene2snp_map):
    inter = []
    for gene in subgr:
        leng = np.random.randint(min_leng, max_leng+1)
        st = np.random.randint(0, len(gene2snp_map[gene])-leng+1)
        for i in range(leng):
            inter.append(gene2snp_map[gene][st+i])
    return sorted(inter)

def generate_network(n_edges, n_genes):
  genes = [f'gene{x}' for x in range(n_genes)]
  edges = []

  while len(edges) != n_edges:
    geneA = random.choice(genes)
    geneB = random.choice(genes)
    if geneA == geneB:
      continue
    if geneA < geneB:
      if not (geneA, geneB) in edges:
        edges.append((geneA, geneB))
    else:
      if not (geneB, geneA) in edges:
        edges.append((geneB, geneA))

  return edges, genes

def add_clique(subgr, edges):
    for i in range(len(subgr)):
        for j in range(i+1, len(subgr)):
            geneA = subgr[i]
            geneB = subgr[j]

            if geneA < geneB:
              if not (geneA, geneB) in edges:
                edges.append((geneA, geneB))
            else:
              if not (geneB, geneA) in edges:
                edges.append((geneB, geneA))

def generate_genes_snps(gene_lst, min_n_snps_per_gene, max_n_snps_per_gene, max_overlap):
    snp2gene_map = []
    gene2snp_map = {}
    snp_ids = []
    id = 0
    for gene in gene_lst:
        gene2snp_map[gene] = []
        lenn = random.choice(range(min_n_snps_per_gene, max_n_snps_per_gene))
        for i in range(lenn):
            id2 = id+i
            if len(snp_ids) == id2:
                snp_ids.append("snp"+str(id2))
                snp2gene_map.append([])
            snp2gene_map[id2].append(gene)
            gene2snp_map[gene].append(id2)
        id = id + lenn - random.choice(range(max_overlap+1))

    return snp_ids, snp2gene_map, gene2snp_map









# takes n_samples n_genes n_edges snp_per_gene outfile
def main():
  n_samples = int(sys.argv[1])

  N_GENES = int(sys.argv[2])
  N_EDGES = int(sys.argv[3])


  min_interval_size = 2
  max_interval_size = 3

  min_n_snps_per_gene = 3
  max_n_snps_per_gene = int(sys.argv[4])

  ass_strength = float(sys.argv[5])

  out_pref = sys.argv[6]

  seed = int(sys.argv[7])

  injectClique = int(sys.argv[8])

  random.seed(seed)
  np.random.seed(seed)

  N_SIMULATIONS = 1




  # read and translate the network that should be used.
  edge_lst, gene_lst = generate_network(N_EDGES, N_GENES)
  for _ in range(3): # add some other random cliques
      add_clique(get_rand_subgraph(gene_lst, edge_lst, 3), edge_lst)

  # start the simulation.

  for sim_id in range(0, N_SIMULATIONS):
      # Create the mapping from SNPs to genes with overlaps.
      snp_ids, snp2gene_map, gene2snp_map = generate_genes_snps(gene_lst, min_n_snps_per_gene, max_n_snps_per_gene, 2)
      #print(snp_ids)
      #print(snp2gene_map)
      #print(gene2snp_map)

      n_features = len(snp_ids)
      data = np.random.binomial(1, 0.01, (n_samples, n_features))
      #print(data)

      rand_subgr_true = get_rand_subgraph(gene_lst, edge_lst, 2)
      if injectClique:
          add_clique(rand_subgr_true, edge_lst)

      #print(rand_subgr_true)
      intervs_true = get_intervals_from_subgraph(rand_subgr_true, min_interval_size, max_interval_size, gene2snp_map);
      #print(intervs_true)

      while True:
          rand_subgr_conf = get_rand_subgraph(gene_lst, edge_lst, 2)
          if len([_ for _ in rand_subgr_conf if _ in rand_subgr_true]) == 0: # we want disjoint graphs
              break
      #print(rand_subgr_conf)
      intervs_conf = get_intervals_from_subgraph(rand_subgr_conf, min_interval_size, max_interval_size, gene2snp_map);
      #print(intervs_conf)

      #print(np.sum(np.any(data[:, intervs_true], axis = 1)))


      # assign the encoding and the labels.
      encoding_sig = np.any(data[:, intervs_true], axis = 1)
      encoding_conf = np.any(data[:, intervs_conf], axis = 1)

      labels = np.random.binomial(1, 0.5, (n_samples))
      labels[encoding_sig] = np.random.binomial(1, ass_strength, (np.sum(encoding_sig)))
      labels[encoding_conf] = np.random.binomial(1, ass_strength, (np.sum(encoding_conf)))

      #print(np.sum(labels)/len(labels))
      #print(np.sum(labels[encoding_sig])/np.sum(encoding_sig))
      #print(np.sum(labels[encoding_conf])/np.sum(encoding_conf))

      # flip the covariate with probability of 0.005 
      covariate = []
      for x in encoding_conf:
        p = np.random.uniform(0, 1, 1)
        if p <= 0.005: covariate.append(int(not x))
        else: covariate.append(int(x))
      covariate = np.asarray(covariate)


      # ----------------------------------------------------------------------
      # write description of simulation study.
      with open(f'{out_pref}_parameters.txt', 'w') as fo:
        fo.write(f'# n. samples: {n_samples} \n')
        fo.write(f'# genes: {N_GENES} \n')
        fo.write(f'# edges: {N_EDGES} \n')
        fo.write(f'Min. SNPs per gene: {min_n_snps_per_gene} \n')
        fo.write(f'Max. SNPs per gene: {max_n_snps_per_gene} \n')
        fo.write(f'Min. interval size: {min_interval_size} \n')
        fo.write(f'Max. interval size: {max_interval_size} \n')
        fo.write(f'Assoc. strength: {ass_strength} \n')
        fo.write(f'Seed: {seed} \n')
        fo.write("Significant subgraph: "+str(rand_subgr_true)+" "+str(intervs_true)+"\n")
        fo.write("Confounding subgraph: "+str(rand_subgr_conf)+" "+str(intervs_conf)+"\n")


      np.savetxt(f'{out_pref}_X.txt', np.transpose(data), fmt='%d')
      np.savetxt(f'{out_pref}_Y.txt', labels, fmt='%d')
      np.savetxt(f'{out_pref}_C.txt', covariate, fmt='%d')
      np.savetxt(f'{out_pref}_snpID.txt', snp_ids, fmt='%s')

      # write the map.
      with open(f'{out_pref}_snp_map.txt', 'w') as fout:
          for key, val in gene2snp_map.items():
              fout.write(f'{key} ')
              for snp in val:
                  fout.write(f'{snp_ids[snp]} ')
                  pass
              fout.write('\n')

      # write the edges.
      with open(f'{out_pref}_edges.txt', 'w') as fout:
        for edge in edge_lst:
          fout.write(f'{edge[0]}\t{edge[1]}\n')



  pass


if __name__ == '__main__':
  main()
