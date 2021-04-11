clear;
close all;
n = 5000;     % Problem dimension m-by-n
m = n;
r = 5;        % rank
alpha = 0.3;  % Percentage of outliers 
c = 5;        % Parameter controls the size of outliers

%% Generate a RPCA problem
disp('Generating a RPCA problem.')
A_generater = randn(m,r);
B_generater = randn(r,n);
L_true = A_generater * B_generater;
norm_of_L_true = norm(L_true,'fro');

S_supp_idx = randsample(m*n, round(alpha*m*n), false);
S_range = c*mean(mean(abs(L_true)));
S_temp = 2*S_range*rand(m,n)-S_range; 
S_true = zeros(m, n);
S_true(S_supp_idx) = S_temp(S_supp_idx);                       
norm_of_S_true = norm(S_true,'fro');

D = L_true + S_true;


%% IRCUR-R    % Resample row/column version
disp('Running IRCUR-R now.')
para.beta_init = 2*max(max(abs(L_true)));
para.beta      = para.beta_init;
para.tol       = 1e-5;
para.con       = 3;
para.resample  = true;
[C1, pinv_U1, R1, ircur_r_timer, ircur_r_err] = IRCUR( D, r, para);

recover_err_ircur_r = norm(L_true - C1 * pinv_U1 * R1, 'fro') / norm(L_true,'fro')



%% IRCUR-F    % Fix row/column version
disp('Running IRCUR-F now.')
para2.beta_init = 2*max(max(abs(L_true)));
para2.beta      = para.beta_init;
para2.tol       = 1e-5;
para2.con       = 3;
para2.resample  = false;
[C2, pinv_U2, R2, ircur_f_timer, ircur_f_err] = IRCUR( D, r, para2);

recover_err_ircur_f = norm(L_true - C2 * pinv_U2 * R2, 'fro') / norm(L_true,'fro')


%% Plot the converegence
figure;
plot(cumsum(ircur_r_timer),ircur_r_err,'bD-',cumsum(ircur_f_timer),ircur_f_err,'r+-');
legend('ICUR-R','ICUR-F');
title('Relative Error vs Runtime');
ylabel('Relative Error');
xlabel('Time(secs)');
set(gca,'YScale', 'log')
