%用追赶法求解三对角线性方程组
function x=crout(a,c,d,b,n)
%a为对角线向量，c为对角线上面带宽为1的向量，d为对角线下面带宽为1的向量
%n为三对角矩阵的维数
p = zeros(n,1);
q = zeros(n,1);
y = zeros(n,1);
x = zeros(n,1);
p(1) = a(1);
for i=1:n-1
    q(i) = c(i)/p(i);
    p(i+1) = a(i+1)-d(i)*q(i);
end
y(1) = b(1)/p(1);
for i=2:n
    y(i) = (b(i)-d(i-1)*y(i-1))/p(i);
end
x(n) = y(n);
for i=n-1:-1:1
    x(i) = y(i)-q(i)*x(i+1);
end