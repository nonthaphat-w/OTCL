# OTCL
optimal gene tree completion using RF-distance in linear time

We implement a linear time algorithm of optimal tree completion problem under the Robinson-Foulds distance (RF distance) .
OTCL takes the input as a reference tree T (which is a tree of leaf set S) and a gene tree t (which is a tree of leaf set R where R is a subset of S). It adds a leaf set S\R into t to minimize the RF distance compared to tree T.

Require library from https://github.com/jamolnng/argparse

gcc version 7.2.0

Current version: 1.0

compile:	

	g++ OCTL.cpp -o OTCL

run

	./OCTL -b <reference tree> -i <gene tree> --seed <random seed> -n <size of leaf set S> -m <size of leaf set R> -o <output file>