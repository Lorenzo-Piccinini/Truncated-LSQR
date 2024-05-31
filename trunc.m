function [X1,X2]=trunc(A,B,r,tol,type_trunc)
% function [X1,X2]=trunc(A,B,r,tol,type_trunc)
% 
% function used to compute low rank factors X1*X2' that approximates A*B'.
% INPUT:
% - A,B: factors we want to approximate.
% - r: maximum rank after truncation.
% - tol: truncation tolerance.
% - type_trunc: how truncation is done.

[Q1,R1]=qr(A,0);
[Q2,R2]=qr(B,0);
[U,S,V]=svd(full(R1*R2'),0);
SS = diag(S);

%svd(A)
%svd(B)
%diag(S)

% figure(15)
% semilogy(diag(S),'*')
% hold on
% semilogy(svd(A),'o')
% semilogy(svd(B),'d')
% hold off

if type_trunc==1
  t=find(cumsum(diag(S))./sum(diag(S))>1-tol,1);
else
  t=find(sqrt(cumsum(diag(S).^2))./norm(diag(S))>1-tol,1);
end
t=max([t,1]);
k = round(min([r t]));
%pause
%X1=Q1*(U(:,1:k)*sqrt(S(1:k,1:k)));
%X2=Q2*(V(:,1:k)*sqrt(S(1:k,1:k)));
X1 = Q1*(U(:,1:k)*diag(sqrt(SS(1:k))));
X2 = Q2*(V(:,1:k)*diag(sqrt(SS(1:k))));

end
