#!/bin/bash
input_file_name=SAMPLE.dat
N=40
lambda_0=0.1
p=0.2
Delta_alpha=0.01

./optimize_QE $input_file_name $N $lambda_0 $p $Delta_alpha
