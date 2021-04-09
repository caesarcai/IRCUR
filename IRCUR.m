function [ C,pinv_U,R, timer, err] = IRCUR( D, r, para )
% 
% Inputs:
% D : Observed matrix. Sum of underlying low rank matrix and underlying
%     sparse matrix. 
% r : Target rank of underlying low rank matrix.
% params : parameters for the algorithm
%   .max_iter : Maximum number of iterations. (default 200)
%   .tol : Desired Frobenius norm error. (default 1e-6)
%   .beta_init : Parameter for thresholding at initialization. (default
%                4*beta)
%   .beta : Parameter for thresholding. (default 1/(2*nthroot(m*n,4)))
%   .gamma : Parameter for desired convergence rate. Value should between 0
%            and 1. Turn this parameter bigger will slow the convergence
%            speed but tolerate harder problem, such as higher p, r or mu. 
%            (default 0.7)   
%   .mu : Incoherence of underlying low rank matrix. Input can be in format
%         of .mu = mu_max, or .mu = [mu_U, mu_V]. (default 5) 
%   .con : constant for row/column samples. con*r*log(n) rows and columns
%          will be sampled. (default 4)
%   .resample : Whether the program resamples the rows and columns every
%               iteration. (default true)
%
% Outputs:
% C， pinv_U， R : CUR decomposition of D, pinv_U is the seudo inverse of U.
% timer : time consumed in each iteration.
% err: relative error of each iteration.
%
% Please cite our paper "Rapid Robust Principal Component Analysis: 
% CUR Accelerated Inexact Low Rank Estimation" if you use this code
%
%
%
%

if exist('PROPACK', 'dir')==7
    addpath PROPACK;
    propack_exist = true;
else
    propack_exist = false;
end

[m,n]     = size(D);

%% Default/Inputed parameters
max_iter  = 200;
tol       = 1e-6;
beta      = 1/(2*nthroot(m*n,4));
beta_init = 4*beta;
gamma     = 0.7;    
mu        = 5;     
con       = 4;
resample  = true;

%% parameter setting
if isfield(para,'beta_init') 
    beta_init = para.beta_init; 
    fprintf('beta_init = %f set.\n', beta_init);
else
    fprintf('using default beta_init = %f.\n', beta_init);
end

if isfield(para,'beta') 
    beta = para.beta; 
    fprintf('beta = %f set.\n', beta);
else
    fprintf('using default beta = %f.\n', beta);
end

if isfield(para,'gamma') 
    gamma = para.gamma; 
    fprintf('gamma = %f set.\n', tol);
else
    fprintf('using default gamma = %f.\n', gamma);
end

if isfield(para,'mu') 
    mu = para.mu; 
    fprintf('mu = [%f,%f] set.\n', mu(1), mu(end));
else
    fprintf('using default mu = [%f,%f].\n', mu, mu);
end

if isfield(para,'max_iter')   
    max_iter = para.max_iter; 
    fprintf('max_iter = %d set.\n', max_iter);
else
    fprintf('using default max_iter = %d.\n', max_iter);
end

if isfield(para,'tol')        
    tol= para.tol; 
    fprintf('tol = %e set.\n', tol);
else
    fprintf('using default tol = %e.\n', tol);
end 

if isfield(para,'con')
    con = para.con;
    fprintf('sample const = %d set. \n',con);
else
    fprintf('using default sample const = %d. \n',con);
end

if isfield(para,'resample')
    resample = para.resample;
    fprintf('resample = %d set. \n',resample);
else
    fprintf('using default resample = %d. \n',resample);
end

err    = -1*ones(max_iter,1);
timer  = zeros(max_iter,1);



tic
[m,n] = size(D);

siz_col = ceil(con*r*log(n));
siz_row = ceil(con*r*log(m));

if ~resample
    rows = randi(m,1,siz_row);
    cols = randi(n,1,siz_col);
    rows = unique(rows);
    cols = unique(cols);
    D_cols = D(:,cols);
    D_rows = D(rows,:);
    norm_of_D = (norm(D_rows, 'fro')+ norm(D_cols, 'fro'));
end
init_timer = toc;

%% main algorithm
for t = 1 : max_iter
    tic;
    
    %% Resample
    if resample
        rows = randi(m,1,siz_row);
        cols = randi(n,1,siz_col);
        rows = unique(rows);
        cols = unique(cols);
        D_cols = D(:,cols);
        D_rows = D(rows,:);
        norm_of_D = (norm(D_rows, 'fro')+ norm(D_cols, 'fro'));
    end
    
    %% update S
    
    if t == 1
        zeta = beta_init;
        S_cols = wthresh( D_cols,'h',zeta);
        S_rows = wthresh( D_rows,'h',zeta);
    else
        zeta = gamma * zeta;
        L_cols = C*pinv_U*(R(:,cols));
        L_rows = (C(rows,:))*pinv_U*R;
        S_rows = wthresh( D_rows-L_rows,'h',zeta);
        S_cols = wthresh( D_cols-L_cols,'h',zeta);
    end

    
    %% Update C pinv_U R
    C = D_cols-S_cols;
    R = D_rows-S_rows;
    MU = C(rows,:);
    [Uu,Su,Vu] = svd(MU);
    d = diag(Su);
    Su = diag(1./d(1:r));
    pinv_U = Vu(:,1:r)*Su*(Uu(:,1:r))';

    %% calculate L to get error
    L_temp=C*pinv_U*R;
    L_temp_cols = L_temp(:,cols);
    L_temp_rows = L_temp(rows,:);
    
    %% Stop Condition
    err(t) = (norm(D_rows-L_temp_rows-S_rows, 'fro') + norm(D_cols-L_temp_cols-S_cols, 'fro'))/norm_of_D;
    timer(t) = toc;
    if err(t) < tol  
        fprintf('Total %d iteration, final error: %e, total time: %f  \n', t, err(t), sum(timer(timer>0)));
        err(1) = err(1) + init_timer;
        err = err(1:t);
        timer = timer(1:t);
        return;
    else
        fprintf('Iteration %d: error: %e, timer: %f \n', t, err(t), timer(t));
    end
    
end
fprintf('Maximum iterations reached, final error: %e.\n======================================\n', err(t));
err(1) = err(1) + init_timer;
end

