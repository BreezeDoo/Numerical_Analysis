%��׷�Ϸ�������Խ����Է�����
function x=crout(a,c,d,b,n)
%aΪ�Խ���������cΪ�Խ����������Ϊ1��������dΪ�Խ����������Ϊ1������
%nΪ���ԽǾ����ά��
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