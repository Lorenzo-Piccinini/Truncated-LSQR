format short e
format compact

%rng(1)
% Setting the rank
p = 1; 

% Uncomment if you want to directly run this code instead of run.m
%
% k = 10;
% Setting 1
% n = (k+1)*p;
% m = k*p;
% Setting 2
% n = 2000
% m = 900;


% tot_timecg=[]; tot_timelsqr4=[];

% Creation of data

switch pb

case 1
nn = 120;
A = gen5pt(nn,nn,1);
b = rand(nn*nn,1);b=b/norm(b);
[V,H] = arnoldi(A,m,b);
fprintf('Conditioning of the problem: %e\n', condest(H'*H))
n = m+1;
I = speye(n,m);
H1 = H; H2 = H; I1 = I; I2 = I;
rhs1 = eye(n,1);
rhs2 = rhs1;


%ex 2  pde
case 2

E = ones(m,1);     I = speye(n,m);
H = spdiags(E,1,n,m)-spdiags(E,-1,n,m);
H = (m+1)/2*H;
H(m,m-1:m) = (m+1)*[ -1 1]; 
H(1,1:2) = (m+1)*[-1 1];
H = [[1,sparse(1,m-1)];H(1:m,1:m)];
H1 = H;
xx = ones(n,1);
xx = linspace(0,1,n)';
H2 = H1; H2(2:n,1) = zeros(n-1,1);
I1 = I;
I2 = I;
rhs1 = exp(xx);rhs2=exp(xx);
rhs = rhs1*rhs2'+(-exp(xx)+xx)*ones(1,n);

%{
[x,y]=meshgrid(linspace(0,1,n));
u=@(x,y)(exp(x+y)-exp(0)+x);
y=linspace(0,1,n);
rhs(n,:)=rhs(n,:)+u(1,y);
x=linspace(0,1,n);
rhs(:,n)=rhs(n,:)+u(x,1);
%}

%ex 3 dictionary learning example,
%y \approx  D x 
%D larga
%D = kron(D2',D1)+kron(D4',D3)    =>    D1 X D2 + D3 X D4 \approx Y

case 3
    
    addpath('./Data')
    load mnist_all

    n_sample = 200;
    n_pix = 2;
    n_pix_rhs = 1;
    dim_patch = 4;
    digit = 2;
    [H1,H2,I1,I2,rhs1,rhs2] = dict_setup(n_sample,dim_patch,digit);
    rhs = rhs1*rhs2';

%ex 4  toeplitz matrices
case 4

    H1=sparse(toeplitz([3,-1,-1/2,zeros(1,n-3)],[3,1,-0,0,zeros(1,m-4)]));
%good  H
    H2=sparse(toeplitz([-1,3,zeros(1,n-2)],[-1,1/2,-1, 0/2 0 zeros(1,m-5)]));
%bad  H
 %   H2=sparse(toeplitz([-1,3,zeros(1,n-2)],[-1,3/2,-1, 0/2 0 zeros(1,m-5)]));
%very bad  H
%    H2=sparse(toeplitz([-1,3,zeros(1,n-2)],[-1,2,-1, 0/2 0 zeros(1,m-5)]));
%   I1=speye(n,m);
    I1=H1;
%   I2=speye(n,m);
    I2=H2;
    rhs1=ones(n,1);
    rhs2=ones(n,1);
    rhs=rhs1*rhs2';

end


% direct soln
%%{
%M=kron(I1,H1)+kron(H2,I2);
%size(H1)
%condest(M'*M),pause
%rhs=rhs1*rhs2';
%sol=zeros(m,m);
%sol=M\rhs(:);
%ssol=reshape(sol,m,m);

%mesh(ssol)
%{
u=@(x,y)(exp(x+y)-exp(x)+x);
[x,y]=meshgrid(linspace(0,1,n));
subplot(1,2,1);mesh(u(x,y))
subplot(1,2,2);mesh(ssol)
%}
%conv tol (normal eqn res)
tol = 1e-6; %diff in res tol taken as  tol/10
%max # iter
imax = 550;
imax = 350;
%imax=500;

%truncation tol
tol_tr = 1e-12;

%max iterate rank
% r = 1000;
r = 100;
% r = 50;


fprintf('Truncated CG: \n')
tic;
X01 = zeros(size(H1,2),1);
X02 = zeros(size(H2,2),1);
[X1,X2,r_res1,a_res1,rks1,p1,p2,PP1] = TCG_gsylv_trunc(H1,H2',rhs1,rhs2,I1',I2,tol,imax,X01,X02,tol_tr,r);
t_tcg = toc;


fprintf('Truncated LSQR: \n')
tic;
its_max = imax;
Params.r = 10;
Params.tol = tol;
Params.tol_tr = tol_tr;
Params.imax = imax;
[ZZs1,ZZs2,r_res4,a_res4,rks4,DD2]=TRUNC_LSQR(H1,H2',rhs1,rhs2,I1',I2,Params);
t_tlsqr = toc;


fprintf('Truncated Adaptive LSQR: \n')
its_max = imax;
Params.r = 10;
Params.tol = tol;
Params.tol_tr = tol_tr;
Params.imax = imax;
tic;
[ZZ5,ZZ6,r_res6,a_res6,rks6,DD6] = TRUNC_LSQR_ADAPTIVE(H1,H2',rhs1,rhs2,I1',I2,Params);
t_lsqrAd = toc;

R1=[rhs1, -H1*X1, -I2*X1]; R2=[rhs2, I1*X2, H2*X2];
trueres_TCG = sqrt(trace( (R2'*R2)*(R1'*R1) ))/(norm(rhs1)*norm(rhs2))

R1=[rhs1, -H1*ZZs1, -I2*ZZs1]; R2=[rhs2, I1*ZZs2, H2*ZZs2];
trueres_lsqr = sqrt(trace( (R2'*R2)*(R1'*R1) ))/(norm(rhs1)*norm(rhs2))

R1=[rhs1, -H1*ZZ5, -I2*ZZ5]; R2=[rhs2, I1*ZZ6, H2*ZZ6];
trueres_lsqrAdpt = sqrt(trace( (R2'*R2)*(R1'*R1) ))/(norm(rhs1)*norm(rhs2))

fprintf('Time TCG: %e, Time TLSQR: %e, Time Adaptive TLSQR: %e\n', t_tcg,  t_tlsqr,  t_lsqrAd)



figure(1)
semilogy(0:length(r_res1)-1,r_res1,'bo-','linewidth',4)
hold on
semilogy(0:length(r_res4)-1,r_res4,'md-','linewidth',4)
semilogy(0:length(r_res6)-1,r_res6,'g+-','linewidth',4)
hold off
title('Relative normal eqn Residuals')
legend({'TCG','TLSQR','Adaptive LSQR'})
xlabel('Iterations')
ylabel('Norm of Relative normal eqn Residual')
hold off

figure(11)
semilogy(0:length(a_res1)-1,a_res1,'b+-','linewidth',4)
hold on
semilogy(0:length(r_res4)-1,r_res4,'md-','linewidth',4)
semilogy(0:length(r_res6)-1,r_res6,'g+-','linewidth',4)
hold off
title('Relative LS Residuals')
legend({'TCG','TLSQR_V6','LSQR_NEW'})
xlabel('Iterations')
ylabel('Norm of Relative LS Residual')
hold off

figure(2)
semilogy(rks1,'+b-')
hold on
semilogy(rks4,'md-')
semilogy(rks6,'g+-','MarkerSize',10)
title('Ranks plot')
legend({'TCG','TLSQR_V6','LSQR_NEW'})
xlabel('Iterations')
ylabel('Rank')
hold off

% figure(3)
% semilogy(0:length(a_res1)-1,a_res1,'b+-','MarkerSize',5)
% hold on
% semilogy(0:length(a_res2)-1,a_res2,'ro-','MarkerSize',5)
% hold on
% semilogy(0:length(a_res3)-1,a_res3,'g*-','MarkerSize',5)
% hold on
% semilogy(0:length(a_res4)-1,a_res4,'kd-','MarkerSize',5)
% hold on
% semilogy(0:length(a_res5)-1,a_res5,'b+-','MarkerSize',10)
% hold off
% title('Absolute Residuals')
% legend({'TCG','Sketch-TCG','TLSQR','Sketch-TLSQR','Sketch-TLSQRV2'})
% %legend({'Sketch-TCG','TLSQR','Sketch-TLSQR','Sketch-TLSQRV2'})
% xlabel('Iterations')
% ylabel('Norm of Absolute Residual')
% hold off


% fprintf('Checking the computed directions: \n')
% p1sum = []; p2sum = [];
% for j=1:min([length(PP1),length(PP2),5])
% 
%     P1 = PP1{j}; P2 = PP2{j};
%     P1 = P1(:); P2 = P2(:);
%     p1sum = [p1sum, P1];
%     p2sum = [p2sum, P2];
%     sp1 = min(svd(p1sum));
%     sp2 = min(svd(p2sum));
%     disp([j,sp1,sp2,(P1'*P2)/(norm(P1)*norm(P2))])
% 
% end