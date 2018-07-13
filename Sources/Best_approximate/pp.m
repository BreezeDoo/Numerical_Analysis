%插值多项式pn
function y=pp(pn,x)
y = 0;
for i=1:9
    y = y+pn(i)*x^(9-i);
end