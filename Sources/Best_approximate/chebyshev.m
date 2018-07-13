%��Chebyshev����ʽ��������һԪ����ʽ��ֵ

clear; clc;
x = zeros(1,9);
for i=1:9
    x(i) = cos((2*i-1)*pi/18);
end
f = 1 ./(1+25 .*x .^2);
N = [zeros(1,8),f(1)];
for k=1:8
    n = 1;
    for i=1:k
        n = conv(n,[1,-x(i)]);
    end
    n = [zeros(1,8-k),n];
    cs = Diff_Quot(x,f,k);
    N = N+cs(1) .*n;
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
title('��9��Chebyshev����ʽ���Ĳ�ֵ���ͼ');
abserr = abs(err);
E2 = max(abserr);