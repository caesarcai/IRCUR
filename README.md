# IRCUR for Rapid Robust PCA
This is Matlab repo for a rapid non-convex Robust Principal Component Analysis (RPCA) algorithm, coined Iterative Robust CUR (IRCUR):

[1] HanQin Cai, Keaton Hamm, Longxiu Huang, Jiaqi Li, and Tao Wang. <a href=https://doi.org/10.1109/LSP.2020.3044130>Rapid Robust Principal Component Analysis: CUR Accelerated Inexact Low Rank Estimation</a>, *IEEE Signal Processing Letters*, 28 (2021): 116-120.

###### To display math symbols properly, one may have to install a MathJax plugin. For example, [MathJax Plugin for Github](https://chrome.google.com/webstore/detail/mathjax-plugin-for-github/ioemnmodlmafdkllaclgeombjnmnbima?hl=en).

## Robust Principal Component Analysis
In this project, we focus on RPCA problem under fully observed setting, that is about separating a low rank matrix $L\in \mathbb{R}^{m\times n}$ and a sparse outlier matrix $S\in \mathbb{R}^{m\times n}$ from their sum $D = L + S$.

## Key Idea for Acceleration
We used fast CUR decomposition in the place of low rank approximation and redesign all the costly steps in the classic alternating projections framework to reduce the computational complexity to $O(\max \lbrace m,n \rbrace r^2 \log(m) \log(n))$ flops. More details can be found in our paper [1].


## Syntex
Using all default parameters:
```
[C, pinv_U, R, ircur_timer, ircur_err] = IRCUR( D, r, '');
```

Using custom parameters:
```
para.beta_init = 2*max(max(abs(L_true)));
para.beta      = para.beta_init;
para.tol       = 1e-5;
para.con       = 3;
para.resample  = true;
[C, pinv_U, R, ircur_timer, ircur_err] = IRCUR( D, r, para);
```

## Input Description
1. D : Observed matrix. Sum of underlying low rank matrix and underlying sparse matrix. 
1. r : Target rank of underlying low rank matrix.
1. params : parameters for the algorithm
   * .max_iter : Maximum number of iterations. (default 200)
   * .tol : Desired Frobenius norm error. (default 1e-6)
   * .beta_init : Parameter for thresholding at initialization. (default 4\*beta)
   * .beta : Parameter for thresholding. (default 1/(2*nthroot(m\*n,4)))
   * .gamma : Parameter for desired convergence rate. Value should between 0 and 1. Turn this parameter bigger will slow the convergence speed but tolerate harder problem, such as higher $\alpha$, $r$ or $\mu$. (default 0.7)   
   * .mu : Incoherence of underlying low rank matrix. Input can be in format of .mu = mu_max, or .mu = [mu_U, mu_V]. (default 5) 
   * .con : constant for row/column samples. $con\*r\*\log(m)$ rows and $con\*r\*\log(n)$ columns will be sampled. (default 4)
   * .resample : Whether the program resamples the rows and columns every iteration. (default true)

## Output Description
1. C，pinv_U，R : CUR decomposition of $D = C U^\dagger R$, where $U^\dagger$ is the pseudo-inverse of $U$.
1. timer : time consumed in each iteration.
1. err: relative error of each iteration.

## Demo
Clone the codes and run the demo file *test_IRCUR.m*. It contains 2 demos, one for IRCUR-R (resample rows/columns every iterations), another for IRCUR-F (no iterative resample). We also plot the relative err vs time per iteration for you, which should show the linear convergence of IRCUR.

## Video Demo
Shoppingmall: https://youtu.be/05CdKIy2KEA

Restaurant: https://youtu.be/rfkRkeOGJHM


