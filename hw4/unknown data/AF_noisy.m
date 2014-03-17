clear all
close all
load('hw4_y2.mat');
K = 3;
L = K;
N = length(y);
Y = toeplitz(y(L+1:N),y(L+1:-1:1));
Yt = Y;
a = zeros(1,size(Yt,1));
b = zeros(1,size(Y,2));
for iter = 1:100
    [U,S,V] = svd(Yt);
    Uh = U;
    Sh = S(:,1:K);
    Vh = V(:,1:K);
    Yh = Uh*Sh*Vh';
    for kc = 1:size(Yh,1)
        a(kc) = mean(diag(Yh,-kc+1));
    end
    for kr = 1:size(Yh,2)
        b(kr) = mean(diag(Yh,kr-1));
    end
    Yt = toeplitz(a,b);
end

[Ut,St,Vt] = svd(Yt);
h = Vt(:,size(Vt,2));
H = h(end:-1:1);
r = roots(H);
freqs = -angle(r);

