# HOGImine
HOGImine (Higher-Order Genetic Interaction miner), is a pattern-mining-based algorithm for finding genetic meta-markers, i.e. combinations of genetic markers, that show a statistical association with a phenotype.
Compared to the state-of-the-art miners, it expands the class of discoverable genetic meta-markers 
by considering higher-order interactions of genes and by allowing multiple encodings for the
genetic variants. Moreover, our method can exploit prior biological knowledge on
gene interactions, such as protein-protein interaction networks, genetic pathways and protein complexes,
to restrict its search space.

HOGImine accepts both binary encodings and additive encodings for the markers.

### Compilation
Compilation uses make. Currently the code is compiled with ```gcc-12```, but other compilers might work as well.
```
cd HOGImine-binary
make
```

### Options
The options for HOGImine are:
- ```-i file```: marker (binary or additive encoding) file
- ```-l file```: labels file
- ```-s file```: snp names file
- ```-c file```: covariate file
- ```-m file```: snp map file
- ```-e file```: edge file
- ```-f level```: target fwer
- ```-o or -O file```: output file name, with -o in compressed format and with -O in verbose format
- ```-p p```: number $p$ of permutations (default $0$), if $p > 0$ it runs a permutation testing procedure
- ```-d d```: maximum interval length
- ```-v```: outputs all the testable patterns



### Usage example
```
./src/hogimine_additive -i data/athaliana/interactome_0kb/avrRpm1/avrRpm1_X.txt -l data/athaliana/interactome_0kb/avrRpm1/avrRpm1_Y.txt -c data/athaliana/covar_snps/avrRpm1/avrRpm1_covar_n2.txt -s data/athaliana/interactome_0kb/avrRpm1/avrRpm1_snpID.txt -m data/athaliana/interactome_0kb/avrRpm1/avrRpm1_snp_map.txt -e data/athaliana/athal_ppi/interactome/AI_interactions_genes.txt -f 0.05 -O out_athaliana
```

