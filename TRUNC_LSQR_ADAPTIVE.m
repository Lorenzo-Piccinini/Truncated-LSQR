function[X_1,X_2,r_res,estimated_res,rks,DD] = TRUNC_LSQR_ADAPTIVE(A,B,C1,C2,F,E,Params)
%
% function[X,D,r_res,a_res,rks]=lsqr_gen_trunc_V4(A,B,C1,C2,F,E,tol,imax,tol_tr,r)
%
% Truncated version of the LSQR algorithm (from the paper of Page and Saunders), in the sense that it uses low
% rank approximations instead of full matrices. This implementation solves
% the Generalized Sylvester leasts squares problem 
% 
%                   A*X*F+E*X*B=C1*C2^T,
%
% but it can be
% generalized also for other version of the Sylvester equation.
%
% See the preprint by Valeria Simoncini, Lorenzo Piccinini 
% https://hal.science/hal-04437719
%
% INPUT: 
% - A,B,F,E: coefficient matrices.
% - C1,C2: low rank rhs.
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


% Computing norm of the operator \cc{A} equivalent to L (needed for
% stopping criteria, see Paige and Saunders.
% n_A = norm(A,'fro')*norm(F,'fro')+norm(E,'fro')*norm(B,'fro');

% Uncomment the next line to check the loss of orthogonality
% global sigmadrop

% Computing the norm of the right-hand side
beta = sqrt(trace((C1'*C1)*(C2'*C2)));
sq_beta = sqrt(beta);

% Storing the right-hand side of the normal equation (needed for the
% stopping criteria)
CC1 = [A'*C1, E'*C1];
CC2 = [F*C2, B*C2];
% AtC = CC1*CC2';

res0_2 = sqrt(trace((CC1'*CC1)*(CC2'*CC2)));

% Updatingss
U_1 = C1/sq_beta;             
U_2 = C2/sq_beta;
res = beta;
totres = res;
truenormres = res;

maximum_rank = Params.r

% COmputing the Three-Terms truncation for the first factor
[QU_1,RU_1] = qr(U_1,0);      
[QU_2,RU_2] = qr(U_2,0);

[uu,S,vv] = svd(full(RU_1*RU_2'));
SS = diag(S);
RU = SS;
QU_1 = QU_1*uu;                 
QU_2 = QU_2*vv;   

% check again for p>1
%norm(U_1*U_2' - QU_1*S*QU_2',1),pause

maxrank = 0;

[V_1, V_2] = L_T(A,B,F,E,U_1,U_2,diag(RU));
[QV_1,QV_2,RV] = ThreeTermsTrunc(V_1,V_2,Params.r,Params.tol_tr);

maxrank = max([maxrank,length(RV)]);
totrank(1) = maxrank;

alfa = norm(RV);
RV = RV/alfa;

% Initializing parameters.
phi = beta;
rho = alfa;
i = 0;

rks = [];

% Initializing the zero solution in the three-terms truncation
X_1 = zeros(size(QV_1));       
X_2 = zeros(size(QV_2));

QX_1 = eye(size(X_1));         
QX_2 = eye(size(X_2));
RX = zeros(size(RV,2),1);

rks = [rks;size(X_1,2)];
QD_1 = QV_1;                   
QD_2 = QV_2;
RD = RV;

DD = 0;

res0 = beta;
res_old = res0;
r_res = [];
r_res = [r_res;res0/res0];
estimated_res = [];
estimated_res = [estimated_res;res0];
res_2 = res0_2;

% Uncomment the next 2 lines to chek loss of orthogonality
% sigmatot=[];
% sigmatot1=[];


while ( i < Params.imax )

    i = i+1;
   
    % Uncomment the next line to chek the loss of orthogonality
    % uold=QU_1*diag(RU)*QU_2';uold=uold(:);


    % Updating the matrix U

    [Ut_1, Ut_2] = L(A,B,F,E,QV_1,QV_2,RV);
    Iu = eye(length(RU));     
    Iu1 = eye(size(Ut_1,2));
    [QU_1,QU_2,RU] = ThreeTermsUpdate(QU_1,-alfa*diag(RU),Ut_1,Iu1,QU_2,Iu,Ut_2,Iu1,Params.r,Params.tol_tr);

    % Uncomment the next line to chek loss of orthogonality
    % sigmatot=[sigmatot,sigmadrop];
    
    maxrank = max([maxrank,size(QU_2,2)]);

    beta = norm(RU);
    RU = RU/beta; 

    % Uncomment the next 4 lines to check the loss of orthogonality
    % unew=QU_1*diag(RU)*QU_2';unew=unew(:);
    % uorth(i)=unew'*uold;
    % [i,unew'*uold]
    % vold=QV_1*diag(RV)*QV_2';vold=vold(:);

    % Updating the matrix V

    [Vt_1, Vt_2] = L_T(A,B,F,E,QU_1,QU_2,RU);
    [QV_1,QV_2,RV] = ThreeTermsUpdate(QV_1,-beta*diag(RV),Vt_1,eye(size(Vt_1,2)),QV_2,eye(length(RV)),Vt_2,eye(size(Vt_2,2)),Params.r,Params.tol_tr);
    
    % Uncomment the next line to chek loss of orthogonality
    % sigmatot1=[sigmatot1,sigmadrop];

    maxrank = max([maxrank,length(RV)]);

    alfa = norm(RV);
    RV = RV/alfa;

    % Uncomment the next 3 lines to chek loss of orthogonality
    % vnew=QV_1*diag(RV)*QV_2';vnew=vnew(:);
    % vorth(i)=vnew'*vold;
    % [i,vnew'*vold]

    rho1 = (rho^2+beta^2)^(0.5);

    c = rho/rho1;
    s = beta/rho1;

    theta = s*alfa;
    rho = -c*alfa;
    phi1 = c*phi;
    phi = s*phi;

    % Updating the matrix X

    coef = phi1/rho1;
    [QX_1,QX_2,RX] = ThreeTermsUpdate(QX_1,diag(RX),QD_1,coef*diag(RD),QX_2,eye(length(RX)),QD_2,eye(length(RD)),Params.r,Params.tol_tr);

    maxrank = max([maxrank,length(RX)]);
    % Uncomment the next line to chek loss of orthogonality
    % sigmadrop
  
    rks = [rks;length(RX)];

    [wrk1, wrk2] = L(A,B,F,E,QX_1,QX_2,RX);
    res_old = truenormres;
    % Computing the true residual
    ResLS1 = [C1, -wrk1]; ResLS2 = [C2, wrk2];
    truenormres = sqrt(trace( (ResLS2'*ResLS2)*(ResLS1'*ResLS1) ));

    totres = [totres,truenormres];

    % Updating the matrix D

    coef=theta/rho1;
    [QD_1,QD_2,RD] = ThreeTermsUpdate(QV_1,diag(RV),QD_1,-coef*diag(RD),QV_2,eye(length(RV)),QD_2,eye(length(RD)),Params.r,Params.tol_tr);
    maxrank = max([maxrank,length(RD)]);

    % Uncomment the next line to chek loss of orthogonality
    % sigmadrop
    
    res = phi;
    res_2 = phi*alfa*abs(c);
    r_res = [r_res; truenormres/res0];
    estimated_res = [estimated_res; res/res0]; 

    % disp([i, res/res0, truenormres/res0, abs(truenormres-res_old)/truenormres, maxrank])
    % disp([i, res/res0, res_2/res0_2, abs(res-res_old)/res, maxrank])

    % Stopping criteria
    %if abs(truenormres-res_old)/truenormres < Params.tol || abs(truenormres-res_old) < Params.tol
    if abs(truenormres-res_old) < Params.tol
        Params.r = Params.r+10;
    end

    if truenormres/res0 <= Params.tol || abs(truenormres-res_old)/truenormres < Params.tol/50
        break, end

    totrank = [totrank,maxrank];
end
% Uncomment the next line to chek loss of orthogonality
% sigmatot

X_1 = QX_1*diag(sqrt(RX)); X_2 = QX_2*diag(sqrt(RX));
disp([i, Params.imax, Params.tol, res_2/res0_2, res/res0])

% Uncomment to see the plot about the loss of orthogonality
%{
figure(202)
semilogy(totres/totres(1),'d-','linewidth',4)
%semilogy(r_res,'linewidth',4)
hold on
semilogy(abs(uorth),'o','linewidth',4)
semilogy(abs(vorth),'x','linewidth',4)
%semilogy(totrank/totrank(end),'x','linewidth',4)
rf=find(totrank/totrank(end)==1,1,'first');
semilogy([rf,rf],[1e-15,10],'linewidth',4)
hold off
legend('true res normal eqn','orth U','orth V','max rank')
xlabel('number of iterations')
ylabel('loss of optimality')
 axis([0,150,1e-16,10]);
hold off
figure(203)
semilogy(sigmatot,'o','linewidth',4)
hold on
semilogy(sigmatot1,'x','linewidth',4)
hold off 
 axis([0,150,1e-16,10]);
xlabel('number of iterations')
ylabel('magnitude of dropped sing.values')
%}

end

function[Y1,Y2] = L(A,B,D,E,QX1,QX2,RX)
% Operator L(X1*X2^T) := A*X1*X2^T*D+E*X1*X2^T*B = Y1*Y2^T.

% RX is a vector
RX = sqrt(RX);
RX = diag(RX);
QX1 = QX1*RX;
QX2 = QX2*RX;
Y1 = [A*QX1, E*QX1];
Y2 = [D'*QX2, B'*QX2];
end



function[Y1, Y2] = L_T(A,B,D,E,QZ1,QZ2,RZ)
% Operator L^T(Z1*Z2^T) := A^T*Z1*Z2^T*D^T+E^T*Z1*Z2^T*B^T = Y1*Y2^T.

% RZ is a vector
RZ = diag(sqrt(RZ));
QZ1 = QZ1*RZ;
QZ2 = QZ2*RZ;
Y1 = [A'*QZ1, E'*QZ1];
Y2 = [D*QZ2, B*QZ2];
end