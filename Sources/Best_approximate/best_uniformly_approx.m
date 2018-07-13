%用Remez算法实现最佳一致逼近

clear; clc;
x = zeros(10,1);
for i=1:10
    x(i) = cos((i-1)*pi/9);
end
while 1
    f = 1 ./(1+25 .*x .^2);
    A = zeros(10,10);
    for i=1:10
        A(i,:) = [x(i)^8,x(i)^7,x(i)^6,x(i)^5,x(i)^4,x(i)^3,x(i)^2,x(i),1,(-1)^(i-1)];
    end
    G = A\f;
    P = G(1:9,1);
    u = G(10);
    x0 = -1:0.1:1;
    Px = 0;
    for r=1:9
        Px = Px+P(r) .*x0 .^(9-r);
    end
    err = 1 ./(1+25 .*x0 .^2)-Px;
    if max(abs(err))-abs(u)<=1.0e-6
        break;
    end
    xx = Newton(P,x(3));
    if xx<-1 || xx>1
        disp('error')
    elseif xx>=-1 && xx<x(1)
        if (f(1)-pp(P,x(1)))*(runge(xx)-pp(P,xx))>=0
            x(1) = xx;
        else
            for i=2:10
                x(i) = x(i-1);
            end
            x(1) = xx;
        end
    elseif xx>x(10) && xx<=1
        if (f(10)-pp(P,x(10)))*(runge(xx)-pp(P,xx))>=0
            x(10) = xx;
        else
            for i=1:9
                x(i) = x(i+1);
            end
            x(10) = xx;
        end
    else
        for i=1:9
            if xx>x(i) && xx<x(i+1)
                if (f(i)-pp(P,x(i)))*(runge(xx)-pp(P,xx))>=0
                    x(i) = xx;
                else
                    x(i+1) = xx;
                end
            end
        end
    end
end
disp(P)
plot(x0,err)
E3 = max(abs(err));
    