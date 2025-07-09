#!/bin/bash

#Get required files for PK-Graphs from DSSR
for f in ../bpRNA_refine/PDB/DSSR_results/graphs_canonical/bpRNA*_segment_graphs.txt
do echo -e -n "$f\n"; done | less > DSSR_canonical_graph_paths.txt

for f in ../bpRNA_refine/PDB/DSSR_results/graphs_nc/bpRNA*_segment_graphs.txt
do echo -e -n "$f\n"; done | less > DSSR_noncanonical_graph_paths.txt
