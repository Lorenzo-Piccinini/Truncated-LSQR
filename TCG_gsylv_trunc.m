function[X1,X2,r_res,a_res,rks,Param,Param2,PP]=TCG_gsylv_trunc(A,B,C1,C2,D,E,tol,imax,X1,X2,tol_tr,r)
% function[X,r_res,a_res,rks,Param,Param2]=TCG_gsylv_trunc(A,B,C1,C2,D,E,tol,imax)
%
% Truncate Conjugate Gradient implemented in order to solve:
%
%                 A*X*D+E*X*B=C1*C2^T
% 
% It can be generalized to other linear matrix equation (e.g. T-Sylvester
% quation).
%
% INPUT:
% - A,B,D,E: coefficient matrices.
% - C1,C2: right-hand side C=C1*C2^T.
% - tol: tolerance chosen for stopping criteria.
% - imax: maximum number of iterations allowed.
% 
% OUTPUT:
% - X: approximated solution X=X1*X2^T.
% - r_res: vector of relative residuals.
% - a_res: vector of absolute residuals.
% - rks: vector of ranks of X_k.
% - Param: vector of parameters alpha_k.
% - Param2: vector of parameters beta_k.


% Initializing the approximation.
%X1 = zeros(size(A,2),1);
%X2 = zeros(size(E,2),1);

% Setting of the truncation function.
%tol_tr = 1e-6; 
%r = 200; 
flag = 1;
N=size(X1,1);
s=200;
maxrank=r;
%hsA = setup_sketching_handle(N,s);
%hsB = setup_sketching_handle(N,s);
%tol_tr=tol_tr;

% Computing the residual matrix:
% R = A'*C*D'+E'*C*B'-A'*(A*X*D+E*X*B)*D'-E'*(A*X*D+E*X*B)*B';

R1 = [A'*C1, E'*C1, -A'*(A*X1), -A'*(E*X1), -E'*(A*X1), -E'*(E*X1)];
R2 = [D*C2, B*C2, D*(D'*X2), D*(B'*X2), B*(D'*X2), B*(B'*X2)];
%[R1,R2] = trunc(R1,R2,r,tol_tr,flag);

%P = R;
P1 = R1;
P2 = R2;
[P1,P2] = trunc(P1,P2,r,tol_tr,flag);
%[P1,P2] = sketch_trunc(P1,P2,hs,r,tol_tr);
%PP{1} = P1*P2';
PP{1} = 1;
rank(PP{1});
%Q = A'*(A*P1*P2'*D+E*P1*P2'*B)*D'+E'*(A*P1*P2'*D+E*P1*P2'*B)*B';
Q1 = [A'*(A*P1), A'*(E*P1), E'*(A*P1), E'*(E*P1)];
Q2 = [D*(D'*P2), D*(B'*P2), B*(D'*P2), B*(B'*P2)];
size(Q1)
[Q1,Q2] = trunc(Q1,Q2,r,tol_tr,flag);

%e = trace(P'*Q);
e = trace(((Q2)'*(P2))*((P1)'*(Q1)));
k = 0;
r_res = [];
a_res = [];
rks = [];

%res0 = norm(R,'fro');
res0 = sqrt(trace((R1'*R1)*(R2'*R2)));
trueres0 = sqrt(trace((C1'*C1)*(C2'*C2)));
res = res0;
a_res = [a_res; 1];
r_res = [r_res; res/res0];

% Computing matrices we are going to use at each iteration (avoiding to
% compute them every time).
%ATC = A'*C1; ETC = E'*C1; ATA = A'*A;
%ATE = A'*E; ETA = ATE'; ETE = E'*E;
%DC = D*C2; BC = B*C2; DDT = D*D';
%DBT = D*B'; BDT = DBT'; BBT = B*B';
AtC1=A'*C1; EtC1=E'*C1; DC2=D*C2; BC2=B*C2;
trueres=1;trueres_old=0;

Param = []; Param1 = []; Param2 = []; Param3 = [];

while (res/res0 > tol && k < imax)  % && abs(trueres-trueres_old)/trueres>tol)

    %alfa = trace(R'*P)/e;
    alfa = trace((P2'*R2)*(R1'*P1))/e;
    %alfa = trace((hsB(P2)'*hsB(R2))*(hsA(R1)'*hsA(P1)))/e;
    %alfa1 = a_res(end)^2/e;
    Param = [Param; (alfa)];
    %Param1 = [Param1; norm(alfa1)];
    sq_alfa = sqrt(alfa);

    %X = X+alfa*P;
    X1 = [X1, sq_alfa*P1];
    X2 = [X2, sq_alfa*P2];
    % X1 = [X1, sqrt(alfa)*P1];
    % X2 = [X2, P2];
    [X1,X2] = trunc(X1,X2,r,tol_tr,flag);

    rks = [rks; size(X1,2)];

    %R = A'*C*D'+E'*C*B'-A'*(A*X*D+E*X*B)*D'-E'*(A*X*D+E*X*B)*B';
    %R1 = [A'*C1, E'*C1, -A'*A*X1, -A'*E*X1, -E'*A*X1, -E'*E*X1];
    %R1 = [ATC, ETC, -ATA*X1, -ATE*X1, -ETA*X1,-ETE*X1];
    wrk1=A*X1;wrk2=E*X1;
    R1 = [AtC1, EtC1, -A'*wrk1, -A'*wrk2, -E'*wrk1, -E'*wrk2];
    %R1 = [A'*C1, E'*C1, -A'*(A*X1), -A'*(E*X1), -E'*(A*X1), -E'*(E*X1)];
    wrk11=D'*X2;wrk12=B'*X2;
    R2 = [DC2, BC2, D*wrk11, D*wrk12, B*wrk11, B*wrk12];
    %R2 = [D*C2, B*C2, D*(D'*X2), D*(B'*X2), B*(D'*X2), B*(B'*X2)];
 trueres_old=trueres;
%                 A*X*D+E*X*B=C1*C2^T
  ResLS1=[C1, -wrk1, -wrk2]; ResLS2=[C2, wrk11, wrk12];
  trueres=sqrt(trace( (ResLS2'*ResLS2)*(ResLS1'*ResLS1) ));


    %R2 = [D*C2, B*C2, D*D'*X2, D*B'*X2, B*D'*X2, B*B'*X2];
    %R2 = [DC, BC, DDT*X2, DBT*X2, BDT*X2, BBT*X2];
    %disp('R1')
    %[R1,R2] = trunc(R1,R2,r,tol_tr,flag);
    %size(R1)

    %res = norm(R,'fro');
    res = sqrt((trace((R1'*R1)*(R2'*R2))));
if imag(res)~=0,res,res=real(res);end

    %if (k>20 && res/a_res(end)>0.95), break, end

    r_res = [r_res; real(res)/res0];
    a_res = [a_res; real(trueres)/trueres0];
    %a_res = [a_res; res];
    %beta = -trace((R1*R2')'*(Q1*Q2'))/e
    beta = -trace((Q2'*R2)*(R1'*Q1))/e;
    %beta = -trace((hsB(Q2)'*hsB(R2))*(hsA(R1)'*hsA(Q1)))/e;
    Param2 = [Param2; (beta)];

    if beta<0,fprintf('Beta negative'), break, end
    if r_res(end) < tol, fprintf('small r_res\n'), r_res(end),break, end
    if abs(trueres-trueres_old)/trueres<tol/10,break,end
    %pause
    %beta1 = res^2/(a_res(end-1)^2);
    %Param2 = [Param2; (beta)];
    %Param3 = [Param3; norm(beta1)];
    %if (k>1), res^2/(r_res(end-1)*res0)^2,end
    %pause
    sq_beta = sqrt(beta);
    %P = R+beta*P;
    P1 = [R1, sq_beta*P1];
    P2 = [R2, sq_beta*P2];
    % P1 = [R1, sqrt(beta)*P1];
    % P2 = [R2, P2];
    %P=P1*P2';
    %Q=Q1*Q2';
    %trace(P'*Q)
    %disp('P')
    [P1,P2] = trunc(P1,P2,r,tol_tr,flag);
    %[P1,P2] = sketch_trunc(P1,P2,hs,r,tol_tr);
    %[norm(P1,1),norm(P2,1),rank(P1,1e-12)],rank(P2,1e-12)
    %svd(P1)
    %pause
    %size(P1)
    %Q = A'*(A*P*D+E*P*B)*D'+E'*(A*P*D+E*P*B)*B';
    %Q1 = [ATA*P1, ATE*P1, ETA*P1, ETE*P1];
    %Q2 = [DDT*P2, DBT*P2, BDT*P2, BBT*P2];   
    wrk1=A*P1;wrk2=E*P1;
    Q1 = [A'*wrk1, A'*wrk2, E'*wrk1, E'*wrk2];
    wrk1=D'*P2;wrk2=B'*P2;
    Q2 = [D*wrk1, D*wrk2, B*wrk1, B*wrk2];
%size(Q1)
%   Q1 = [A'*(A*P1), A'*(E*P1), E'*(A*P1), E'*(E*P1)];
%   Q2 = [D*(D'*P2), D*(B'*P2), B*(D'*P2), B*(B'*P2)];
    %P=P1*P2';
    %Q=Q1*Q2';
    %trace(P'*Q)
    %disp('Q1')
    [Q1,Q2] = trunc(Q1,Q2,r,tol_tr,flag);
    %size(Q1)
    %pause
    %e = trace(P'*Q);
    %e = trace((hsB(Q2)'*hsB(P2))*(hsA(P1)'*hsA(Q1))); %p_k^TAp_k
    e = trace((Q2'*P2)*(P1'*Q1));
   %disp(real([k, r_res(end), a_res(end), alfa, beta]))
    k = k+1;
    %PP{k+1} = P1*P2';
    PP{k+1}=1;
    %pause
end
%X = X1*X2';
disp([k,res/res0, tol])
%pause

% figure(45)
% semilogy(Param,'*-')
% hold on
% semilogy(Param1,'o-')
%legend('\alpha','\beta')

%norm(Param-Param1)
%norm(Param2-Param3)
% disp([Param,Param1,Param2,Param3])
% pause
