# KDBB: An Exact Algorithm with New Upper Bounds for the Maximum k-Defective Clique Problem in Massive Sparse Graphs

Since the KDBB source code is not publicly available, this is the version we implemented in our experiments.

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