function [V,H]=arnoldi(AM,m,res);
%function [V,H]=arnoldi(AM,m,res);
% creates a basis V for the krylov 
% subspace [ v, AMv, AM^2 v ,..AM^(m-1)v]
% using Gram-Schimidt orth. with re-orth. 
% H is the (m+1)xm matrix of coeff.
%      
tet = 10.; hinorm=0.;
[n,n]=size(AM);
V=zeros(n,m+1); H=zeros(m+1,m);
V(1:n,1) = res(1:n)/norm(res);
for i=1:m;
  i1=i+1; it = 0; t=0.;
  V(:,i1) = AM*V(:,i);
  while (  it < 2 )
  %while ( t*tet <= hinorm & it < 2 )
    hinorm=0.; it = it+1;
    for j=1:i
      t = V(1:n,j)'*V(1:n,i1);
      hinorm = hinorm + abs(t^2); H(j,i)=H(j,i)+t;
      V(1:n,i1)=V(1:n,i1)-t*V(1:n,j);
    end;
    t = norm(V(1:n,i1));
  end;
  H(i1,i)=t;
  if (t ~= 0.)
    t=1.0/t; V(1:n,i1)=V(1:n,i1)*t;
  end;
end;
