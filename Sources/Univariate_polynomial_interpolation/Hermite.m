%Hermite插值
clear all;
x=-pi:pi/7:pi;
for i=1:15
    f0(i)=sin(2*x(i));
    f1(i)=2*cos(2*x(i));
end
H=zeros(1,30);
for k=1:15
    L=1;
    LD=0;
    for i=1:15
        if i~=k
            L=conv(L,[1,-x(i)])./(x(k)-x(i));
            LD=LD+1/(x(k)-x(i));
        end
    end
    L2=conv(L,L);
    H=H+conv([-2*LD,2*LD*x(k)+1],L2).*f0(k)+conv([1,-x(k)],L2).*f1(k);
end
disp(H)
Hx=0;
for r=1:30
    Hx=Hx+H(r) .*x .^(30-r);
end
plot(x,Hx)
grid on
title('f(x)的Hermite插值多项式图像')

