clear all,
close all,

addpath(genpath('./Support'))

fprintf('Truncated Multiterm LSQR: \n')
n_vec = 20^2;
epsilon = 1e-1;
n = n_vec;

%h=1/(n-1);
h = 2 / (n - 1);

% Construct the data
% [OP{1,1},OP{1,3},OP{1,4},PREC{1,1},OP{2,2},OP{2,4},...
%     OP{2,3},PREC{1,2},C{1},C{2},C{3}]=build1(h,epsilon,n);
% 
% OP{1,2}=speye(n);
% OP{2,1}=speye(n);
% p=length(OP);
OP{1, 1} = randn(n);
OP{1, 3} = randn(n);
OP{1, 4} = randn(n);
OP{2, 2} = randn(n);
OP{2, 4} = randn(n);
OP{2, 3} = randn(n); 
OP{1, 2} = randn(n);
OP{2, 1} = randn(n);
p = length(OP);

its_max = 300;
Params.r = 10000;
Params.tol = 1e-10;
Params.tol_tr = 1e-14;
Params.imax = 100;
Params.imax = 10000;
m = 10;
rng(5)
R1 = randn(n, m); R1 = R1 / norm(R1); 
R2 = randn(n, m); R2 = R2 / norm(R2);

KK = sparse(n*n, m*m);
for k=1:4
  A{k,1}=OP{1,k}*R1; 
  A{k,2}=OP{2,k}*R2;
  AG{1, k} = OP{1, k} * R1; 
  AG{2, k} = OP{2, k} * R2;
  KK = KK + kron(A{k, 2}, A{k, 1});
end
mm = 4;
C{1} = randn(n, mm);
C{2} = randn(n, mm);
rhs = C{1} * C{2}'; rhs = rhs(:);

% Computing exact solution
truesol = KK \ rhs;

fprintf('Running built-in LSQR on the vectorized problem:\n')
tic;
[xx,dd,rr]=lsqr(KK,rhs,1e-10,300);
t_vec = toc;
res_rel_vec = norm(KK * xx - rhs) / norm(rhs);
res_rel_vec_ne = norm(KK' * (KK * xx - rhs)) / norm(KK' * rhs);

fprintf('Time: %.4e, Rel Res: %.4e, Rel Res Normal Eq.: %.4e\n', t_vec, res_rel_vec, res_rel_vec_ne)

fprintf('Running Multi-term truncated matrix LSQR:\n')
tic;
[X_1,X_2, r_res,estim_re,rks,DD] = Truncated_LSQR_multiterm1(AG, C, Params);
t_mat = toc;
X_sol = X_1 * X_2';
x_sol = X_sol(:);
res_rel_mat = norm(KK * x_sol - rhs) / norm(rhs);
res_rel_mat_ne = norm(KK' * (KK * x_sol - rhs)) / norm(KK' * rhs);

fprintf('Time: %.4e, Rel Res: %.4e, Rel Res Normal Eq.: %.4e\n', t_mat, res_rel_mat, res_rel_mat_ne)

errorLSmultivslsqr = norm(reshape(xx, m, m) - X_1 * X_2', 'fro') / norm(xx);
errorLS_lsqr = norm(reshape(xx, m, m) - reshape(truesol, m, m), 'fro') / norm(xx);
errorLS_multi = norm(reshape(truesol, m, m) - X_1 * X_2', 'fro') / norm(xx);

fprintf('Error Multi-term truncated matrix LSQR vs Vec LSQR: %.4e\n', errorLSmultivslsqr)
fprintf('Error Backslash vs Vec LSQR: %.4e\n', errorLS_lsqr)
fprintf('Error Multi-term truncated matrix LSQR vs Backslash: %.4e\n', errorLS_multi)
