function[X_1,X_2,r_res,estimated_res,rks,DD] = Truncated_LSQR_multiterm(A,C,Params)
% function[X, r_res] = Truncated_LSQR_multiterm(A, rhs, Params)
%
% We want to solve the multiterm generalyzed Sylvester equation
%
%      A{1,1}*X*A{1,2}' + A{2,1}*X*A{2,2}' + ... + A{p,1}*X*A{p,2}' = C{1}*C{2}'
%      A{1,1}*X*A{2,1}' + A{1,2}*X*A{2,2}' + ... + A{1,p}*X*A{2,p}' = C{1}*C{2}'
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
C1 = C{1}; 
C2 = C{2};
beta = sqrt(trace((C1' * C1) * (C2' * C2)));
sq_beta = sqrt(beta);

CC1 = []; CC2 = [];
nterms = size(A, 2);
for k = 1:nterms 
         CC1 = [CC1, A{1,k}' * C1];
         CC2 = [CC2, A{2,k}' * C2];
end

res0_2 = sqrt(trace((CC1' * CC1) * (CC2' * CC2)));

U_1 = C1; % / norm(C1);
U_2 = C2; % / norm(C2);
res = beta;
totres = res;
truenormres = res;

maximum_rank = Params.r;

% Computing the Three-Terms truncation for the first factor
[QU_1, RU_1] = qr(U_1, 0);      
[QU_2, RU_2] = qr(U_2, 0);

[uu, S, vv] = svd(full(RU_1 * RU_2'), 0);
SS = diag(S);
RU = SS./norm(SS) ;
QU_1 = QU_1 * uu;                 
QU_2 = QU_2 * vv;   
%norm(QU_1*diag(RU)*QU_2','fro')
%[beta,norm(SS)],pause

maxrank = 0;

[V_1, V_2] = L_T(A, QU_1, QU_2, (RU));
[QV_1, QV_2, RV] = ThreeTermsTrunc(V_1, V_2, Params.r, Params.tol_tr);

maxrank = max([maxrank, length(RV)]);
totrank(1) = maxrank;

alfa = norm(RV);
RV = RV / alfa; 

% Initializing parameters
phi = beta;
rho = alfa;
i = 0;

rks = [];

% Initializing the zero solution in the three-terms truncation
X_1 = zeros(size(QV_1));       
X_2 = zeros(size(QV_2));

%QX_1 = eye(size(X_1));         
%QX_2 = eye(size(X_2));
%RX = zeros(size(RV,2),1);

rks = [rks;size(X_1,2)];
QD_1 = QV_1;                   
QD_2 = QV_2;
RD = RV;

DD = 0;

res0 = beta;
res_old = res0;
r_res = [];
r_res = [r_res; res0 / res0];
estimated_res = [];
estimated_res = [estimated_res; res0];
res_2 = res0_2;

while ( i < Params.imax )

    i = i+1;

    % Updating the matrix U

    [Ut_1, Ut_2] = L(A, QV_1, QV_2, RV);

    Iu = eye(length(RU));     
    Iu1 = eye(size(Ut_1, 2));
     [QU_1, QU_2, RU] = ThreeTermsUpdate(QU_1, -alfa * diag(RU), Ut_1, Iu1, QU_2, Iu, Ut_2, Iu1, Params.r, Params.tol_tr);


    maxrank = max([maxrank,size(QU_2,2)]);

    beta = norm(RU);
    RU = RU / beta; 

    % Updating the matrix V

    [Vt_1, Vt_2] = L_T(A, QU_1, QU_2, RU);

     [QV_1, QV_2, RV] = ThreeTermsUpdate(QV_1, -beta * diag(RV), Vt_1, eye(size(Vt_1, 2)), QV_2, eye(length(RV)), Vt_2, eye(size(Vt_2, 2)), Params.r, Params.tol_tr);


    maxrank = max([maxrank, length(RV)]);

    alfa = norm(RV);
    RV = RV / alfa;

    rho1 = (rho^2+beta^2)^(0.5);
    c = rho/rho1;
    s = beta/rho1;
    theta = s*alfa;
    rho = -c*alfa; 
    phi1 = c*phi; 
    phi = s*phi; 

    % Updating the matrix X

    coef = phi1 / rho1;
    if i == 1 
    [QX_1, QX_2, RX] = ThreeTermsTrunc(QD_1 * (coef * diag(RD)), QD_2, Params.r, Params.tol_tr);
coefX=1;
    else
    [QX_1, QX_2, RX] = ThreeTermsUpdate(QX_1, diag(RX), QD_1, coef * diag(RD), QX_2, eye(length(RX)), QD_2, eye(length(RD)), Params.r, Params.tol_tr);
 coefX=coef*diag(RD);

    end

    maxrank = max([maxrank, length(RX)]);
    
    X_1 = QX_1 * diag((RX)); X_2 = QX_2;  
    rks = [rks; length(RX)];

    [wrk1, wrk2] = L(A, QX_1, QX_2, RX);
    res_old = truenormres;
    % Computing the true residual
    ResLS1 = [C1, -wrk1]; ResLS2 = [C2, wrk2];
    truenormres = sqrt(trace( (ResLS2' * ResLS2) * (ResLS1' * ResLS1) ));

    totres = [totres, truenormres];

    % Updating the matrix D

    coef = theta / rho1;
     [QD_1, QD_2, RD] = ThreeTermsUpdate(QV_1, diag(RV), QD_1, -coef * diag(RD), QV_2, eye(length(RV)), QD_2, eye(length(RD)), Params.r, Params.tol_tr);

    maxrank = max([maxrank, length(RD)]);
    
    res = phi;
    res_2 = phi * alfa * abs(c);
    r_res = [r_res; truenormres / res0];
    estimated_res = [estimated_res; res / res0]; 

    % Stopping criteria
    disp([truenormres/res0  norm(coefX)/norm(RX)])
    if truenormres/res0 <= Params.tol || norm(coefX) / norm(RX) < (Params.tol / 10)
   %if truenormres/res0 <= Params.tol || abs(truenormres - res_old) / truenormres < (Params.tol / 10)
        break, end

    totrank = [totrank, maxrank];
end
%disp([i, truenormres/res0])
end

function[Y1,Y2] = L(A, QX1, QX2, RX)
% Operator Y1*Y2'=L(X1*X2^T) 
% RX is a vector
RX = sqrt(RX);
RX = diag(RX);
QX1 = QX1 * RX;
QX2 = QX2 * RX;
Y1 = []; Y2 = [];
for k = 1:size(A, 2)
         Y1 = [Y1, A{1,k} * QX1];
         Y2 = [Y2, A{2,k} * QX2];
end
end



function[Y1, Y2] = L_T(A, QZ1, QZ2, RZ)
% Operator Y1*Y2'=L^T(Z1*Z2^T) 
% RZ is a vector
RZ = diag(sqrt(RZ));
QZ1 = QZ1 * RZ;
QZ2 = QZ2 * RZ;
Y1 = []; Y2 = [];
for k = 1:size(A,2)
         Y1 = [Y1, A{1,k}' * QZ1];
         Y2 = [Y2, A{2,k}' * QZ2];
end
end
