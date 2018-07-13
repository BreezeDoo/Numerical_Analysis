%对Runge函数的一元多项式插值（Newton法）
clear; clc;
x = -1:0.25:1;
f = 1 ./(1+25 .*x .^2);
h = 0.25;
N = [zeros(1,8),f(1)];
for k=1:8
    chafen = 0;
    for i=0:k
        z = factorial(k)/(factorial(i)*factorial(k-i));
        chafen = chafen+(-1)^i*z*f(1+k-i);
    end
    n = 1;
    for j=1:k
        n = conv(n,[1,-x(j)]);
    end
    n = [zeros(1,8-k),n];
    N = N+(chafen/(h^k*factorial(k))) .*n;
end
disp(N)
xx = -1:0.1:1;
Nx = 0;
for r=1:9
    Nx = Nx+N(r) .*xx .^(9-r);
end
err = 1 ./(1+25 .*xx .^2)-Nx;
plot(xx,err)
grid on
title('一元Newton插值多项式误差图');
abserr = abs(err);
E1 = max(abserr);




