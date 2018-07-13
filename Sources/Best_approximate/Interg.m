%复化梯形公式求积分
function T=Interg(P)
m = length(P);
[k,r] = deconv(P,[25,0,1]);
t0 = 0;
for i=1:2:m-2
    t0 = t0+2*k(i)/(m-1-i);
end
x = -1:0.01:1;
px = polyval(r,x);
f = 1 ./(1+25 .*x .^2) .*px;
t = 0;
for i=2:200
    t = t+f(i);
end
T = 0.005*(f(1)+f(21)+2*t)+t0;