#include"nhead.h"
#include<iostream>
#include"math.h"


int main(){

    double xl=0,xr=1,yd=0,yu=0.5,d=1.0/16.0,fac_d_h=0.9,k=2,Re=200;

    SPHfield f1(xl,xr,yd,yu,d,fac_d_h,k,Re);


    //流场信息赋值
    //int middle=f1.Npx_inner/2;

    int middle=f1.Npx_inner/2;

    for (int i = 0; i < f1.Npx_inner; i++)
    {
        for (int j = 0; j < f1.Npy_inner; j++)
        {
            f1.IP[i][j]->iniset(1,2*f1.IP[i][j]->x,2*f1.IP[i][j]->y,3*f1.IP[i][j]->x+2*f1.IP[i][j]->y);
        }
        
    }
    
    double a=f1.wallSPHestimate(f1.xl+d,f1.yd+d,0);


    double pause=1;


}