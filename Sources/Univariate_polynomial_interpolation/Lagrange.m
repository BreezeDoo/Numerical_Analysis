%Lagrange

clear all;
x=-pi:2*pi/29:pi;
for i=1:30
    f0(i)=sin(2*x(i));
end
P=zeros(1,30);
for i=1:30
    m=1;
    n=1;
    for j=1:30
        if j ~= i
            m=conv(m,[1,-x(j)]);
            n=n*(x(i)-x(j));
        else
            continue;
        end
    end
    P=P+m ./ n .*f0(i);
end
disp(P);
Px=0;
for r=1:30
    Px=Px+P(r) .* x.^(30-r);
end
plot(x,Px,'r:')
grid on
title('f(x)的Lagrange插值多项式图像')
