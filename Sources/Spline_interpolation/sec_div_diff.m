%求二阶差商
function f2=sec_div_diff(x,y)
deltx = diff(x);
delty = diff(y);
f1 = delty ./deltx;     %一阶
deltx2 = zeros(length(x)-2,1);
for i=1:length(x)-2
    deltx2(i) = x(i+2)-x(i);
end
delty2 = diff(f1);
f2 = delty2 ./deltx2;