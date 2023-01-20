# HOGImine
HOGImine (Higher-Order Genetic Interaction miner), is a pattern-mining-based algorithm for finding genetic meta-markers, i.e. combinations of genetic markers, that show a statistical association with a phenotype.
Compared to the state-of-the-art miners, it expands the class of discoverable genetic meta-markers 
by considering higher-order interactions of genes and by allowing multiple encodings for the
genetic variants. Moreover, our method can exploit prior biological knowledge on
gene interactions, such as protein-protein interaction networks, genetic pathways and protein complexes,
to restrict its search space.

We provide two versions, HOGImine-binary, that accepts binary encodings for the markers, and HOGImine-additive, that accepts additive encodings.

### Compilation
Compilation uses make. Currently the code is compiled with ```gcc-12```, but other compilers might work as well.
```
cd HOGImine-binary
make
```

### Options
The options for HOGImine-binary are:
- ```-i file```: marker file
- ```-l file```: labels file
- ```-s file```: snp names file
- ```-c file```: covariate file
- ```-m file```: snp map file
- ```-e file```: edge file
- ```-f level```: target fwer
- ```-o file```: output file name
- ```-p p```: number $p$ of permutations (default $0$), if $p > 0$ it runs a permutation testing procedure
- ```-d d```: maximum interval length
- ```-v```: outputs all the testable patterns

The options for HOGImine-additive are the same, with the following changes:
- ```-i file```: marker (with dominant encoding) file
- ```-h file```: marker (with recessive encoding) file


### Usage example
```
./HOGImine-binary/hogimine_binary -i data/athaliana/interactome_0kb/avrRpm1/avrRpm1_X.txt -l data/athaliana/interactome_0kb/avrRpm1/avrRpm1_Y.txt -c data/athaliana/covar_snps/avrRpm1/avrRpm1_covar_n2.txt -s data/athaliana/interactome_0kb/avrRpm1/avrRpm1_snpID.txt -m data/athaliana/interactome_0kb/avrRpm1/avrRpm1_snp_map.txt -e data/athaliana/athal_ppi/interactome/AI_interactions_genes.txt -f 0.05 -o out_athaliana
```
