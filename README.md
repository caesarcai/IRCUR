# IRCUR for Rapid Robust PCA
This is the Matlab code repo for a rapid non-convex Robust Principal Component Analysis (RPCA) algorithm, coined (Iterative Robust CUR) IRCUR [1].

###### To display math symbols properly, one may have to install a MathJax plugin. For example, [MathJax Plugin for Github](https://chrome.google.com/webstore/detail/mathjax-plugin-for-github/ioemnmodlmafdkllaclgeombjnmnbima?hl=en).

## Robust Principal Component Analysis
In this project, we focus on RPCA problem under fully observed setting, that is about separating a low rank matrix $L\in \mathbb{R}^{m\times n}$ and a sparse outlier matrix $S\in \mathbb{R}^{m\times n}$ from their sum $D = L + S$.

## Key Idea for Acceleration




## Syntex
Using all default parameters:
```
[C1, pinv_U1, R1, ircur_timer, ircur_err] = IRCUR( D, r, '');
```

Using custom parameters:
```
para.beta_init = 2*max(max(abs(L_true)));
para.beta      = para.beta_init;
para.tol       = 1e-5;
para.con       = 3;
para.resample  = true;
[C1, pinv_U1, R1, ircur_timer, ircur_err] = IRCUR( D, r, para);
```

## Input Description
1. D : Observed matrix. Sum of underlying low rank matrix and underlying sparse matrix. 
1. r : Target rank of underlying low rank matrix.
1. params : parameters for the algorithm
   * .max_iter : Maximum number of iterations. (default 200)
   * .tol : Desired Frobenius norm error. (default 1e-6)
   * .beta_init : Parameter for thresholding at initialization. (default 4\*beta)
   * .beta : Parameter for thresholding. (default 1/(2*nthroot(m\*n,4)))
   * .gamma : Parameter for desired convergence rate. Value should between 0 and 1. Turn this parameter bigger will slow the convergence speed but tolerate harder problem, such as higher p, r or mu. (default 0.7)   
   * .mu : Incoherence of underlying low rank matrix. Input can be in format of .mu = mu_max, or .mu = [mu_U, mu_V]. (default 5) 
   * .con : constant for row/column samples. con*r*log(n) rows and columns will be sampled. (default 4)
   * .resample : Whether the program resamples the rows and columns every iteration. (default true)

## Output Description
1. $C$， $pinv\_U$， $R$ : CUR decomposition of $D = C\*pinv\_U\*R$, pinv_U is the seudo inverse of U.
1. timer : time consumed in each iteration.
1. err: relative error of each iteration.


## Reference
[1] HanQin Cai, Keaton Hamm, Longxiu Huang, Jiaqi Li, and Tao Wang. Rapid Robust Principal Component Analysis: CUR Accelerated Inexact Low Rank Estimation, IEEE Signal Processing Letters, 28 (2021): 116-120.
