//text.c
#include<stdio.h>
#include<math.h>
#include"remzf.h"

int main()
{
    int i;
    double a,b,eps,p[10],remzf(double);
    a=-1.0;
    b=1.0;
    eps=1.0e-10;  //精确度
    remz(a,b,p,9,eps,remzf);
    printf("\n");
    printf("各次系数为");
    for(i=0;i<=8;i++)
        printf("p(%2d)=%f\n",i,p[i]);

    printf("\n");
    printf("偏差为%e\n",p[9]);
    printf("\n");
    getch();
    return 0;
}

double remzf(double x)//需要逼近的函数
{
    double y;
    y=1/(1+25*x*x);
    return(y);
}
