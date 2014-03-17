clear all
close all
load('hw4_y1.mat');
N = length(y);
M = round(sqrt(N));
R = zeros(M);
for i = 1:N-M
   yv = y(i:M+i-1);
   R = R + yv*yv';
end
R = R/(N-M+1);
[V,D] = eig(R);
K = 3;
Qn = V(:,1:M-K);
res = 1000;
w = linspace(-pi,pi,res);
Pm = zeros(1,res);
m = (0:(M-1))';
for i = 1:res
    aw = exp(1j*w(i)*m);
    Pm(i) = abs(1/((aw')*Qn*(Qn')*aw));
end