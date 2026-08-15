<div align="justify">


# RARE

## Brief Introduction

**RARE** (MVMR incorporating **R**are variants **A**ccounting for multiple **R**isk factors and shared horizontal pl**E**iotropy) is a multivariable Mendelian randomization (MVMR) method designed to identify the causal relationship between the main exposure ($`X_1`$) and disease outcome (_Y_), while further elucidating the modifying effects of rare genetic variants on the causal relationship. **Fig.1** shows the structure of RARE.

<p align="center">
  <img src="https://github.com/Hide-in-lab/RARE/blob/Supplementary/Github_RARE/Figure1.jpg?raw=true" width="40%" />
  <br>
  <b>Fig.1 </b>
</p>

## Motivation

When investigating the causal effect of common trait on disease outcome, the **estimated effect** ($\color{red}\beta_1$) **may be biased by the potential effects of other correlated traits**. For example, when evaluating the causal effect of high-density lipoprotein (HDL) on cardiovascular disease, it is important to account for the potential influences of low-density lipoprotein (LDL) and triglycerides on the effect estimate ($\color{red}\beta_i$). However, such scenarios cannot be adequately addressed by conventional univariable Mendelian randomization (UVMR), motivating the extension of MR to a multivariable framework.

In addition, in the conventional MR analysis framework, IVs are typically selected based on their minor allele frequencies (MAFs), which may lead to the exclusion of rare variants (**MAF < 0.01**, $\color{red}G_r$) and consequently **overlook the potential contributions of rare variants to causal effect estimation**. To further reduce estimation bias, RARE incorporates multiple exposures while simultaneously accounting for rare variants within the model.

Compared with common variants ($\color{red}G_c$), rare variants often exhibit greater trait specificity. However, **the individual effect sizes of rare variants are generally smaller and their estimates are associated with larger variances** in GWAS data. Therefore, RARE aggregates rare variants into a **polygenic risk score** ($\color{red}PRS$) to jointly capture their effects and enhance their contribution to causal effect estimation.

RARE is the only method (2024-07-31) that accounts for the impact of rare variants in causal inference while simultaneously considers UHP and CHP.

# Data Input

Suppose we have _**K**_ exposures of interest, denoted as $X_1, X_2, \cdots, X_K$, and one outcome _**Y**_. We let $X_1$ to be the $\color{red}{\textbf{primary exposure}}$, and $X_i$ (_i_ = 2, ... ,_K_) are considered as the $\color{red}{\textbf{secondary exposures}}$. Our primary objective is to estimate the causal effect between the primary exposure $X_1$ and the outcome _**Y**_. Suppose there are _**p**_ shared SNPs, among which $n_c$ are common variants and $n_r$ are rare variants, that are simultaneously associated with these exposures, with _**P**_ values below a given threshold (e.g., 5*e-8) from GWAS studies. What we need are:

**1)** _**Estimated IV-to-exposure effect size**_ $`\hat{\boldsymbol{\gamma}}_{ij}`$ and its _**corresponding standard error**_ $`\hat{\mathbf{S}}_{\gamma_{ij}}`$, where

i ($i = 1, 2, \cdots, K$) represents the i-th tissue, and j ($p = 1, 2, \cdots, p$) represents the j-th IV.


**2)** _**Estimated IV-to-outcome effect size**_ $`\hat{\boldsymbol{\Gamma}}_{ij}`$ and its _**corresponding standard error**_ $`\hat{\mathbf{S}}_{\Gamma_{ij}}`$




Let  and  be the true marginal associations of the j-th IV with the i-th exposure  and outcome , respectively. The estimated effects and standard errors of the associations between the j-th SNP and the i-th exposure and outcome from GWAS are denoted as  and , respectively.



# Installation
Install this tool by use of the 'devtools' package. Note that RARE partly depends on the C++ languange, thus you should appropriately set 
Rtools and X code for Windows, Mac OS/X, and Linux, respectively.
```
install.packages( 'devtools' )  
library( devtools )  
install_github( 'Hide-in-lab/RARE@main', force  = T )
```


# Reference

$\color{red}{\textbf{Please kindly cite the following paper if you use the}}$  `RARE` $\color{red}{\textbf{package:}}$

Yu Cheng<sup>+</sup>, Xinjia Ruan<sup>+</sup>, Xiaofan Lu, Yuqing Yang, Yuhang Wang, Shangjin Yan, Yuzhe Sun, Fangrong Yan<sup> #</sup>, Liyun Jiang<sup> #</sup>, Tiantian Liu<sup> #</sup>, **Accounting for the impact of rare variants on causal inference with RARE: a novel multivariable Mendelian randomization method**, Briefings in Bioinformatics, Volume 26, Issue 3, May 2025, bbaf214, https://doi.org/10.1093/bib/bbaf214

# Development

This package is developed and maintained by **Yu Cheng**.\
**Contact e-mail**: yucheng.cpu@foxmail.com


</div>
