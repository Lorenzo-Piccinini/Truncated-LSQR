function [QX1,QX2,RX] = ThreeTermsTrunc(A,B,r,tol)
% function [QX1,QX2,RX] = ThreeTermsTrunc(A,B,r,tol)
% 
% function used to compute low rank factors X1*X2' that approximates A*B'.
% INPUT:
% - A,B: factors we want to approximate.
% - r: maximum rank after truncation.
% - tol: truncation tolerance.

[Q1,R1] = qr(A,0);
[Q2,R2] = qr(B,0);
[U,S,V] = svd(full(R1*R2'),0);
SS = diag(S);
SS = SS(:);

t = find(cumsum(SS)./sum(SS)>1-tol,1);

t = max([t,1]);
k = round(min([r t]));

QX1 = Q1*U(:,1:k); 
QX2 = Q2*V(:,1:k);
RX = SS(1:k);

% We store RX1 as a diagonal matrix!
end