%Newton
clear all;
x=-pi:2*pi/29:pi;
for i=1:30
    f0(i)=sin(2*x(i));
end
h=2*pi/29;
N=[zeros(1,29),f0(1)];
for k=1:29
    chafen=0;
    for i=0:k
        z=factorial(k)/(factorial(i)*factorial(k-i));
        chafen=chafen+(-1)^i*z*f0(1+k-i);
    end
    n=1;
    for j=1:k
        n=conv(n,[1,-x(j)]);
    end
    n=[zeros(1,29-k),n];
    N=N+(chafen/(h^k*factorial(k))) .*n;
end
disp(N)
Nx=0;
for r=1:30
    Nx=Nx+N(r) .*x .^(30-r);
end
plot(x,Nx,'g-.')
grid on
title('f(x)的Newton插值多项式图像')