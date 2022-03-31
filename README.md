# R-RANSAC with SPRT
A MATLAB implementation of R-RANSAC with the SPRT test (WALDSAC)
based on Matas and Chum's 2008 paper [link](https://cmp.felk.cvut.cz/~matas/papers/chum-waldsac-iccv05.pdf). This repo was motivated by the scarcity of implementations available.

The SPRT test was pioneered by Abraham Wald in his seminal [Sequential Tests of Statistical Hypotheses](https://projecteuclid.org/journals/annals-of-mathematical-statistics/volume-16/issue-2/Sequential-Tests-of-Statistical-Hypotheses/10.1214/aoms/1177731118.full) 1945 paper. The paper is very clear and gives an intuitive presentation of sequential tests, i.e. very read-worthy.

## Overview
Both the T(d,d) test and the SPRT test are implemented, see `rransac_tdd.m` and `rransac_sprt.m` respectively. 

## TODO
  - Fix stability issues for SPRT termination criteria when inlier ratio is 1.
  - Improve stability of numerical solver for `h`.
