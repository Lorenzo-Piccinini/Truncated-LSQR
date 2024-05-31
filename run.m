clear all
clc
close all

tot_timecg=[];
tot_timelsqr2=[];
tot_timelsqr4=[];
pb=4;
n = 30000
for m=26000:200:26000,
m
  TCG_vs_STCG,
%pause

end

%{
figure(101)

plot( [tot_timecg' tot_timelsqr2' tot_timelsqr4'])
legend('tcg','tlsqr2','tlsqr4')

xx=X1*X2';
zz0=ZZs1*ZZs2';
%zz=ZZ3*ZZ4';
zzp=ZZ5*ZZ6';

% trueerror_tcg=norm(ssol-xx)/norm(ssol)
% trueerror_tlsqr2=norm(ssol-zz0)/norm(ssol)
% %trueerror_tlsqr4=norm(ssol-zz)/norm(ssol)
% trueerror_tlsqr5=norm(ssol-zzp)/norm(ssol)

figure
%plot(ssol(:,1))
hold on
plot(xx(:,1))
plot(zzp(:,1))
hold off
%}

% figure(301)
% spy(abs(zzp)>2e-4)
% 
% S=[]; for k=1:10, S=[S, norm(zzp(200*(k-1)+1:200*k,200*(k-1)+1:200*k),'fro')]; end
% S