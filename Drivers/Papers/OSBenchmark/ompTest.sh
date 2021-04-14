#!/bin/bash
set -x
# script to run all simulations for the OS paper
time ./OSSilo -test -omp 8 > ompTest8.out
time ./OSSilo -test -omp 10 > ompTest10.out
time ./OSSilo -test -omp 12 > ompTest12.out
time ./OSSilo -test -omp 14 > ompTest14.out
time ./OSSilo -test -omp 16 > ompTest16.out
time ./OSSilo -test -omp 18 > ompTest18.out
time ./OSSilo -test -omp 20 > ompTest20.out
