function [QX1,QX2,RX] = ThreeTermsUpdate(QA1,RA1,QB1,RB1,QA2,RA2,QB2,RB2,r,tol)
% function [QX1,QX2,RX] = ThreeTermsUpdate(QA1,RA1,QB1,RB1,QA2,RA2,QB2,RB2,r,tol)
%  
% Updating the three terms truncation


% Gram-Schmidt-like for X1 

% Uncomment the next line to chek the loss of orthogonality
% global sigmadrop

R12 = QA1'*QB1;
W1 = QB1 - QA1*R12;
t=  QA1'*W1;
R12 = R12+ t;
W1 = W1 - QA1*t;
[Q2,R22] = qr(W1,0);
Q1 = [QA1,Q2];

% Gram-Schmidt-like for X2

P12 = QA2'*QB2;
U1 = QB2-QA2*P12;
t = QA2'*U1;
P12 = P12 + t;
U1 = U1-QA2*t;
[QQ2,P22] = qr(U1,0);
Q2 = [QA2, QQ2]; 

wrk1 = RB1*RB2';
wrk2 = wrk1*P22';
wrk3 = wrk1*P12';

R = [RA1*RA2'+R12*wrk3, R12*wrk2;
     R22*wrk3, R22*wrk2];

[U,S,V] = svd(full(R),0);
SS = diag(S);SS=SS(:);

t = find(cumsum(SS)./sum(SS)>1-tol,1);

t = max([t,1]);
k = round(min([r t]));

QX1 = Q1*U(:,1:k);  
QX2 = Q2*V(:,1:k);
RX = SS(1:k);

% Uncomment the next lines to chekc the loss of orthogonality
% if k+1 < length(SS)
%     sigmadrop = SS(k+1);
% else
%     sigmadrop = 0;
% end
end
