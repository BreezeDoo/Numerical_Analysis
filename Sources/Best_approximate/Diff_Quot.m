%Çón½×²îÉÌ
function fn=Diff_Quot(x,y,n)
temp = y;
for i=1:n
    deltx = zeros(1,length(x)-i);
    delty = zeros(1,length(x)-i);
    for j=1:length(x)-i
        deltx(j) = x(j+i)-x(j);
    end
    delty = diff(temp);
    temp = delty ./deltx;
end
fn = temp;