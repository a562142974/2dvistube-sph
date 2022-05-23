#include"nhead.h"


int main(){

    double xl=0,xr=1,yd=0,yu=0.5,d=1.0/16.0,fac_d_h=0.9,k=2,Re=200;

    SPHfield f1(xl,xr,yd,yu,d,fac_d_h,k,Re);


    //流场信息赋值
    int middle=f1.Npx_inner/2;

    for (int i = 0; i < f1.Npx_inner; i++)
    {
        for (int j = 0; j < f1.Npy_inner; j++)
        {
            f1.IP[i][j]->iniset(1,2*f1.IP[i][j]->x,2*f1.IP[i][j]->y,3*f1.IP[i][j]->x+2*f1.IP[i][j]->y);
        }
        
    }
    

    // //左半边
    // for (int i = 0; i < middle; i++)
    // {
    //     for (int j = 0; j < f1.Npy_inner; j++)
    //     {
    //         f1.IP[i][j]->iniset(2,1,1.5,1/gama/(gama-1));
    //     }
    // }

    // //右半边
    // for (int i = middle; i < f1.Npx_inner; i++)
    // {
    //     for (int j = 0; j < f1.Npy_inner; j++)
    //     {
    //         f1.IP[i][j]->iniset(3,1,1.5,1/gama/(gama-1));
    //     }
    // }




    
    f1.boundarycondition(1);    
    
    f1.updatedt();
    
    outputAVars(f1,"Vars.txt");
    outputparameters(f1,"para.txt");



    outputpositions(f1,2,3,"positions.txt");

}