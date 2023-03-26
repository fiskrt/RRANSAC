# R-RANSAC with Sequential Probability Ratio Test
A MATLAB implementation of R-RANSAC leveraging the sequential probability ratio test (SPRT) test (WALDSAC)
based on Matas and Chum's 2008 paper [link](https://cmp.felk.cvut.cz/~matas/papers/chum-waldsac-iccv05.pdf). This repo was motivated by the scarcity of implementations available.

The SPRT test was pioneered by Abraham Wald in his seminal [Sequential Tests of Statistical Hypotheses](https://projecteuclid.org/journals/annals-of-mathematical-statistics/volume-16/issue-2/Sequential-Tests-of-Statistical-Hypotheses/10.1214/aoms/1177731118.full) 1945 paper. The paper is very clear and gives an intuitive presentation of sequential tests, i.e. very read-worthy.

## Overview
Both the $T(d,d)$ test and the SPRT test are implemented, see `rransac_tdd.m` and `rransac_sprt.m` respectively. 

## Results
![image](https://user-images.githubusercontent.com/43207511/227797535-86cc21e5-c929-41a3-bfe0-b7e16a6cf1a6.png)
![image](https://user-images.githubusercontent.com/43207511/227797602-464dcf21-ef7a-4401-aa60-f4382486ee2b.png)

## TODO
  - Fix stability issues for SPRT termination criteria when inlier ratio is 1.
  - Improve stability of numerical solver for `h`.
