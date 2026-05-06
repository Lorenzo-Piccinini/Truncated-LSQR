function[X, r_res] = TRUNC_LSQR_MULTITERM(A, rhs, Params)
% function[X, r_res] = TRUNC_LSQR_MULTITERM(A, rhs, Params)
%
% We want to solve the multiterm generalyzed Sylvester equation
%
%      A{1,1}*X*A{1,2}' + A{2,1}*X*A{2,2}' + ... + A{p,1}*X*A{p,2} = rhs{1}*rhs{2}'
%
% using a truncated LSQR method. The method is based on the same principles as the one described in TRUNC_LSQR.m, but it is adapted to the multiterm case.
%
% INPUT: 
% - A: structure containing the coefficient matrices as above.
% - rhs: structure containing the two factors of the right-hand side as above.
% - Params.tol: tolerance chosen for the stopping criteria.
% - Params.imax: maximum number of iterations allowed.
% - Params.tol_tr: truncation tolerance.
% - Params.r: maximum rank allowed when doing the truncation.
% 
% OUTPUT:
% - X_1,X_2: approximated solution (it will be computed as X=X1*X2').
% - r_res: vector of real relative residuals computed at each iteration.
% - estimated_res: vector of the estimates of the realtive residual at each iteration.
% - rks: vector of ranks of the approximated solution ar each iteration.
% - DD: array containing all the direction matrices D^(i).
%
% If you use this code please cite the following paper:
% L. Piccinini, and V. Simoncini. Truncated LSQR for matrix least squares problems, 
% Computational Optimization and Applications 91 (2), 905-932.
% DOI: https://doi.org/10.1007/s10589-024-00629-w
%
% If you note any bug or you have any suggestion please contact the authors.


% Computing the norm of the right-hand side
[Q1, R1] = qr(rhs{1}, 0);
[Q2, R2] = qr(rhs{2}, 0);
% Decomposing the right-hand side in low rank factors L1*D*L2^T
C.L1 = sparse(Q1); C.D = sparse(R1*R2'); C.L2 = sparse(Q2);
beta = norm(C.D,'fro');
beta = sqrt(trace(rhs{2}'*rhs{2}*rhs{1}'*rhs{1}));
sq_beta = sqrt(beta);

% Storing the right-hand side of the normal equation (needed for the
% stopping criteria)
C_nrm_eq = OpL_T(A, C, Params.r, Params.tol_tr);

res0_2 = norm(C_nrm_eq.D, 'fro');

U.L1 = C.L1;
U.L2 = C.L2;
U.D = C.D/beta;

res = beta;
totres = res;
truenormres = res;

maximum_rank = Params.r;

V = OpL_T(A, U, Params.r, Params.tol_tr);

maxrank = 0;
maxrank = max([maxrank, length(V.D)]);
% V.D = V.D / norm(V.D, 'fro');
alfa = norm(V.D, 'fro');
V.D = V.D/alfa;

phi = beta;
rho = alfa;
i = 0;

X.L1 = eye(size(V.L1));
X.L2 = eye(size(V.L2));
X.D = zeros(size(V.D));

D.L1 = V.L1;
D.L2 = V.L2;
D.D = V.D;

DD = 0;

res0 = beta;
res_old = res0;
r_res = [];
r_res = [r_res; res0 / res0];
estimated_res = [];
estimated_res = [estimated_res; res0];
res_2 = res0_2;

rks = [];

while i < Params.imax

    i = i + 1;
    Utmp = OpL(A, V, Params.r, Params.tol_tr);
    QU1tmp = Utmp.L1; RU1tmp = sqrt(Utmp.D);
    QU2tmp = Utmp.L2; RU2tmp = sqrt(Utmp.D);
    [U.L1, U.L2, U.D] = ThreeTermsUpdate(U.L1, -alfa * diag(U.D), QU1tmp, RU1tmp, U.L2, eye(size(U.L2, 2)), QU2tmp, RU2tmp, Params.r, Params.tol_tr);
    % U.D = diag(U.D);
    

    maxrank = max([maxrank, size(U.L2, 2)]);
    % beta = norm(U.D, 'fro');
    beta = sqrt(trace((U.L2'*U.L2)*diag(U.D)'*(U.L1'*U.L1)*diag(U.D)));
    
    U.D = U.D / beta; 

    Vtmp = OpL_T(A, U, Params.r, Params.tol_tr);
    QV1tmp = Vtmp.L1; RV1tmp = sqrt(Vtmp.D);
    QV2tmp = Vtmp.L2; RV2tmp = sqrt(Vtmp.D);
    % keyboard
    [V.L1, V.L2, V.D] = ThreeTermsUpdate(V.L1, -beta * diag(V.D), QV1tmp, diag(RV1tmp), V.L2, eye(size(V.L2, 2)), QV2tmp, diag(RV2tmp), Params.r, Params.tol_tr);
    % V.D = diag(V.D);
    maxrank = max([maxrank, size(V.L2, 2)]);
    % alfa = norm(V.D, 'fro');
    alfa = sqrt(trace((V.L2'*V.L2)*diag(V.D)'*(V.L1'*V.L1)*diag(V.D)));
    V.D = V.D / alfa;

    rho1 = (rho^2 + beta^2)^0.5;
    c = rho / rho1;
    s = beta / rho1;
    theta = s * alfa;
    rho = -c * alfa;
    phi1 = c * phi;
    phi = s * phi;

    coef = phi1 / rho1;
    
    [X.L1, X.L2, X.D] = ThreeTermsUpdate(X.L1, diag(X.D), D.L1, coef * diag(D.D), X.L2, eye(size(X.L2, 2)), D.L2, eye(size(D.L2, 2)), Params.r, Params.tol_tr);
    maxrank = max([maxrank, size(X.L2, 2)]);
    % X.D = diag(X.D);

    rks = [rks; size(X.D, 2)];
    res_tmp = OpL(A, X, Params.r, Params.tol_tr);
    res_old = truenormres;
    [res_tmp.L1, res_tmp.L2, res_tmp.D] = ThreeTermsUpdate(C.L1, diag(C.D), res_tmp.L1, -diag(res_tmp.D), C.L2, eye(size(C.L2, 2)), res_tmp.L2, eye(size(res_tmp.L2, 2)), Params.r, Params.tol_tr);
    % truenormres = norm(res_temp.D, 'fro');
    
    truenormres = sqrt(trace( (res_tmp.L2'*res_tmp.L2)*diag(res_tmp.D)'*(res_tmp.L1'*res_tmp.L1) *  diag(res_tmp.D)));
    totres = [totres, truenormres];

    coef = theta / rho1;
    [D.L1, D.L2, D.D] = ThreeTermsUpdate(V.L1, diag(V.D), D.L1, -coef * diag(D.D), V.L2, eye(size(V.L2, 2)), D.L2, eye(size(D.L2, 2)), Params.r, Params.tol_tr);
    maxrank = max([maxrank, size(D.L2, 2)]);

    res = phi;
    res_2 = phi * alfa * abs(c);
    r_res = [r_res; truenormres / res0];
    estimated_res = [estimated_res; res / res0];

    if truenormres / res0 <= Params.tol || abs(truenormres - res_old) / truenormres <= (Params.tol / 10) 
        break, end
    

end
disp([i, Params.imax, Params.tol, truenormres/res0, res/res0])
end

function Y = OpL(A, X, r, tol_tr)
% Using the ThreeTermsUpdate.m
% X = X.L1*X.D*X.L2';
QY1 = A{1,1} * X.L1; RY1 = sqrt(X.D);
QY2 = A{1,2} * X.L2; RY2 = sqrt(X.D);
for i = 2:size(A,1)



    Qtmp1 = A{i,1} * X.L1; Rtmp1 = sqrt(X.D);
    Qtmp2 = A{i,2} * X.L2; Rtmp2 = sqrt(X.D);
    [Y.L1, Y.L2, Y.D] = ThreeTermsUpdate(QY1, RY1, Qtmp1, Rtmp1, QY2, RY2, Qtmp2, Rtmp2, r, tol_tr);
    if i < size(A,1)
        QY1 = Y.L1; RY1 = sqrt(Y.D);
        QY2 = Y.L2; RY2 = sqrt(Y.D);

    end
end
end

function Y = OpL_T(A, X, r, tol_tr)
% Using the ThreeTermsUpdate.m
% X = X.L1*X.D*X.L2';
QY1 = A{1,1}' * X.L1; RY1 = sqrt(X.D);
QY2 = A{1,2}' * X.L2; RY2 = sqrt(X.D);
for i = 2:size(A,1)



    Qtmp1 = A{i,1}' * X.L1; Rtmp1 = sqrt(X.D);
    Qtmp2 = A{i,2}' * X.L2; Rtmp2 = sqrt(X.D);
    [Y.L1, Y.L2, Y.D] = ThreeTermsUpdate(QY1, RY1, Qtmp1, Rtmp1, QY2, RY2, Qtmp2, Rtmp2, r, tol_tr);
    if i < size(A,1)
        QY1 = Y.L1; RY1 = sqrt(Y.D);
        QY2 = Y.L2; RY2 = sqrt(Y.D);

    end
end
end