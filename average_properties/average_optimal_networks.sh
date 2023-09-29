#!/bin/bash
N=64
L=0.1
p=0.1
alpha_MIN=0.
alpha_MAX=1.
Delta_alpha=0.01
M=10
BOOLE=0

./average_optimal_networks $N $L $p $alpha_MIN $alpha_MAX $Delta_alpha $M $BOOLE

