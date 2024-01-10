# KD-Club: An Efficient Exact Algorithm with New Coloring-based Upper Bound for the Maximum k-Defective Clique Problem

This repository implements the KD-Club algorithm proposed in our paper. 

## Compile the code

```sh
$ make clean
$ make
```

## Run the code

```sh
$ ./kDefective -g {path_to_graph} -k {k_value} 
```

## Data format

The format of the first line is "p edge ", followed by the number of vertices and the number of edges. Each subsequent line represents each edge, starting with e, and represented by two endpoints.