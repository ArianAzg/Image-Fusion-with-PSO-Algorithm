function [F]= expEdge(Pan, lamda,eps)
[FX,FY] = gradient(Pan);
u=((FX.*FX)+(FY.*FY)).^(1/2);
%a = max(max(abs(u)))
[m,n]=size(u); 
F= exp( -((lamda*(ones(m,n)))./((abs(u).*abs(u).*abs(u).*abs(u))+ eps)));