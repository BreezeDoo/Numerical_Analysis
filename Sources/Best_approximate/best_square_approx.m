%用Legendre正交多项式求Runge的最佳平方逼近

clear; clc;
p0 = 1;
p1 = [1,0];
N = zeros(1,9);
N(9) = 0.5*Interg(p0);
N(8) = 0;
NN = [zeros(1,8),N(9)];
for k=2:8
    p0 = [[0,0],p0];
    pk = (2*k-1)/k .*conv([1,0],p1)-(k-1)/k .*p0;
    if mod(k,2)==0
        N(9-k) = (2*k+1)/2*Interg(pk);
    end
    NN = NN+N(9-k) .*[zeros(1,8-k),pk];
    p0 = p1;
    p1 = pk;
end
xx = -1:0.1:1;
Nx = 0;
for r=1:9
    Nx = Nx+NN(r) .*xx .^(9-r);
end
err = 1 ./(1+25 .*xx .^2)-Nx;
plot(xx,err)
grid on
title('最佳平方逼近误差函数图像');
E4 = max(abs(err));
        