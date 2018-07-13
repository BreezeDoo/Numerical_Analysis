%用Newton迭代法求偏差点
function xx=Newton(pn,x)
xx = x;
while 1
    f1 = -50*xx/(1+25*xx^2)^2;
    f2 = (3750*xx^2-50)/(1+25*xx^2)^3;
    p1 = 8*pn(1)*xx^7+7*pn(2)*xx^6+6*pn(3)*xx^5+5*pn(4)*xx^4+4*pn(5)*xx^3+ ...
        3*pn(6)*xx^2+2*pn(7)*xx+pn(8);
    p2 = 56*pn(1)*xx^6+42*pn(2)*xx^5+30*pn(3)*xx^4+20*pn(4)*xx^3+12*pn(5)*xx^2+ ...
        6*pn(6)*xx+2*pn(7);
    temp = xx;
    xx = xx-(f1-p1)/(f2-p2);
    if abs(xx-temp)<1.0e-3
        break;
    end
end