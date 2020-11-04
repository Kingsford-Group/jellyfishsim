# Overview

This repository contains the program to compute the similarity matrix of k-mer counts vectors for a list of datasets, and also contains the programs to perform the seeded-chunking and to compute the Hausdorff distances used in the Hierarchical Representative Set Selection. The Hierarchical Representative Set Selection can be found in this repository: https://github.com/Kingsford-Group/hierrepsetselection
 

# Installation

Download the source code from this repository. And then use the following commands to compile:
```
    autoreconf -fi
    ./configure
    make
```

# Usage

## Compute the similarity matrix of k-mer counts vectors

This program can be used for any applications. The cosine similarity is the similarity measure. Use the following command to get the similarity matrix of k-mer counts vectors for a list of datasets:
```
    sortsim <klen> <list_datasets_file>
```
where `<klen>` is the k-mer size (e.g. 17), and `<list_datasets_file>` is a file containing the names of all the datasets' k-mer counts files (full-path), one filename per line. Note that these k-mer counts files are gzipped.


## Perform the seeded-chunking\n(used in the Hierarchical Representative Set Selection)

The seeded-chunking method produces chunks as separately as possible from the original set of datasets. Use the following command to perform the seeded-chunking:
```
    chunking <klen> <full_set_datasets_file> <rand_set_datasets_file> <num_chunks> <chunk_size>
```

where `<klen>` is the k-mer size (e.g. 17), `<full_set_datasets_file>` is a file containing the names of all the datasets' k-mer counts files (full-path) in the original full set to be chunked (one filename per line), `<rand_set_datasets_file>` is a file containing the names of randomly selected datasets' k-mer counts files (these datasets are selected randomly from the original full set and will be used for selecting seeds), `<num_chunks>` is the number of chunks, and `<chunk_size>` is the size of each chunk. All these these k-mer counts files are gzipped.


## Compute the Hausdorff distances\n(used in the Hierarchical Representative Set Selection)

Both the classical Hausdorff distance and the partial Hausdorff distance are computed.  Use the following command to compute the Hausdorff distances:
```
    hausdorff <klen> <full_set_datasets_file> <rep_set_datasets_file> <q>
```

where `<klen>` is the k-mer size (e.g. 17), `<full_set_datasets_file>` is a file containing the names of all the datasets' k-mer counts files (full-path) in the original full set, `<rep_set_datasets_file>` is a file containing the names of the selected representative datasets' k-mer counts files, and `<q>` is a parameter used in the partial Hausdorff distance: `q = 1 – K / |X|` where `|X|` is the size of the original full set, and `K` is for using the Kth largest value (counting from the minimum) as the partial Hausdorff distance. All these k-mer counts files are gzipped.


