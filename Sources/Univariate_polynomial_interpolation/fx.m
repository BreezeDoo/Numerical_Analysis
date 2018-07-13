%f(x)=sin(2x)函数定义
clear all;
x=-pi : 2*pi/29 : pi;
f=sin(2.*x);
plot(x,f)
grid on
title('f(x)=sin(2x)的图像')