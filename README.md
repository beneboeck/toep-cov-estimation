# Gohberg-Semencul Covariance Estimation via Autoregressive Parameters

Welcome to the repository for the paper "Gohberg-Semencul Covariance Estimation via Autoregressive Parameters"! This repository contains the code to reproduce the main results of our research work.
## Abstract

In our work, we introduce a class of Toeplitz Covariance Matrix Estimators and their inverses based on the so-called Gohberg-Semencul decomposition, which is closely related to autoregressive parameters.

## Instructions

The script "main_example.m" contains an exemplary matlab file, in which you can comment in your desired (inverse) covariance estimator. Our code comprises several estimators for Toeplitz structured covariances comprising our proposed estimators "PGD" and "PLS" estimator as well as the baselines "EM", "Circ", "Avg", "Band", "Tape", "TSL", "ShU", and "ShB" (see our paper for the acronyms). The script applies any estimator to N P-dimensional samples generated from an AR(3) process with adjustable parameters. The directory `our_estimators` contains our proposed estimators. The directory `cov_generators` contains scripts to generate Toeplitz covariance matrices. The `baselines` directory contains the baselines, and the `utils` directory stores some auxiliary scripts.

## Citation
If you are using this code for your research, please cite

```bibtex
@ARTICLE{boeckToep2025,
  author={Böck, Benedikt and Semmler, Dominik and Fesl, Benedikt and Baur, Michael and Utschick, Wolfgang},
  journal={IEEE Transactions on Signal Processing}, 
  title={Gohberg-Semencul Toeplitz Covariance Estimation via Autoregressive Parameters}, 
  year={2025},
  volume={73},
  number={},
  pages={858-875},
  keywords={Estimation;Covariance matrices;Tuning;Matrix decomposition;Vectors;Array signal processing;Standards;Parallel processing;Optimization;Hands;Covariance estimation;autoregressive processes;Gohberg-Semencul;Toeplitz;likelihood estimation},
  doi={10.1109/TSP.2025.3536101}}

```
## Licence of Contributions
This code is covered by the BSD 3-Clause License:

> BSD 3-Clause License
>
> Copyright (c) 2023 Benedikt Böck.
> All rights reserved.
>
> Redistribution and use in source and binary forms, with or without
>modification, are permitted provided that the following conditions are met:
>
> * Redistributions of source code must retain the above copyright notice, this
>  list of conditions and the following disclaimer.
>
> * Redistributions in binary form must reproduce the above copyright notice,
>  this list of conditions and the following disclaimer in the documentation
>  and/or other materials provided with the distribution.
>
> * Neither the name of the copyright holder nor the names of its
>  contributors may be used to endorse or promote products derived from
>  this software without specific prior written permission.
>
> THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
> AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
> IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
> DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
> FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
> DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
> SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
> CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
> OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
> OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
