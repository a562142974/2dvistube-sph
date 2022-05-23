#include"nhead.h"
#include<fstream>
#include<string>
#include<sstream>
#include<iostream>

const double xl=0;
const double xr=1;
const double yd=0;
const double yu=0.5;
const double d=1/20.0;
const double fac_d_h=1.1;
const double k=2;
const double Re=200;
const double eps=0.15;
const double h=d*fac_d_h;


std::string d2str(double t){
    std::stringstream ss;
    std::string stime;
    ss<<t;
    ss>>stime;
    return stime;
}

double caldt(SPHfield& f){
        double dt1=1,dt2=1,dt3=1,dt=1;
        for (int i = 0; i < f.Npx_inner; i++)
        {
            for (int j = 0; j < f.Npy_inner; j++)
            {
                dt1=0.25*h/sqrt(gama*f.IP[i][j]->p/f.IP[i][j]->rho)<dt1 ? 0.25*h/sqrt(gama*f.IP[i][j]->p/f.IP[i][j]->rho):dt1;
                dt2=0.25*sqrt(h/sqrt(pow(f.IP[i][j]->dtu,2)+pow(f.IP[i][j]->dtv,2)))<dt2 ? 0.25*sqrt(h/sqrt(pow(f.IP[i][j]->dtu,2)+pow(f.IP[i][j]->dtv,2))):dt2;
                dt3=eps*pow(h,2)*f.IP[i][j]->rho<dt3 ? eps*pow(h,2)*f.IP[i][j]->rho:dt3;
            }
        }

        dt=dt<dt1 ? dt : dt1;
        dt=dt<dt2 ? dt : dt2;
        dt=dt<dt3 ? dt : dt3;
        return dt;
};

void pre(SPHfield &f1,SPHfield &f2,double dt){
        for (int i = 0; i < f1.Npx_inner; i++)
        {
            for (int j = 0; j < f1.Npy_inner; j++)
            {
                f2.IP[i][j]->update(f1.IP[i][j]->rho+0.5*dt*f1.IP[i][j]->dtrho,\
                f1.IP[i][j]->u+0.5*dt*f1.IP[i][j]->dtu,\
                f1.IP[i][j]->v+0.5*dt*f1.IP[i][j]->dtv,\
                f1.IP[i][j]->e+0.5*dt*f1.IP[i][j]->dte);

                f2.IP[i][j]->x=f1.IP[i][j]->x+0.5*dt*f1.IP[i][j]->u;
                f2.IP[i][j]->y=f1.IP[i][j]->y+0.5*dt*f1.IP[i][j]->v;
            }
        }
        f2.get_particle_in_box();
        f2.boundarycondition(0);
        f2.updatedt();
};

void cor(SPHfield &f1,SPHfield &f2,double dt){
    for (int i = 0; i < f1.Npx_inner; i++)
        {
            for (int j = 0; j < f1.Npy_inner; j++)
            {
                f2.IP[i][j]->update(f1.IP[i][j]->rho+0.5*dt*(f1.IP[i][j]->dtrho+f2.IP[i][j]->dtrho),\
                f1.IP[i][j]->u+0.5*dt*(f1.IP[i][j]->dtu+f2.IP[i][j]->dtu),\
                f1.IP[i][j]->v+0.5*dt*(f1.IP[i][j]->dtv+f2.IP[i][j]->dtv),\
                f1.IP[i][j]->e+0.5*dt*(f1.IP[i][j]->dte+f2.IP[i][j]->dte));

                f2.IP[i][j]->x=f1.IP[i][j]->x+0.5*dt*(f1.IP[i][j]->u+f2.IP[i][j]->u);
                f2.IP[i][j]->y=f1.IP[i][j]->y+0.5*dt*(f1.IP[i][j]->v+f2.IP[i][j]->v);
            }
        }
        f2.get_particle_in_box();
        f2.boundarycondition(0);
};

void inter(SPHfield &f1,SPHfield &f2){
    double x,y;
        for (int i = 0; i < f1.Npx_inner; i++)
        {
            for (int j = 0; j < f1.Npy_inner; j++)
            {
                x=f1.IP[i][j]->x;
                y=f1.IP[i][j]->y;
                f1.IP[i][j]->update(f2.reSPHestimate(x,y,4),\
                f2.reSPHestimate(x,y,0),\
                f2.reSPHestimate(x,y,1),\
                f2.reSPHestimate(x,y,5));
            }
        }
        f1.boundarycondition(1);
};



int main(){
    //变量初始化
    SPHfield f1(xl,xr,yd,yu,d,fac_d_h,k,Re);
    SPHfield f2(xl,xr,yd,yu,d,fac_d_h,k,Re);


    outputparameters(f1,"para.txt");

    //初始流场赋值
    int middle=f1.Npx_inner/2;
    for (int i = 0; i < middle; i++)
    {
        for (int j = 0; j < f1.Npy_inner; j++)
        {
            f1.IP[i][j]->iniset(120/gama,0,0,1/gama/(gama-1));
            f2.IP[i][j]->iniset(120/gama,0,0,1/gama/(gama-1));
        }
    }
    for (int i = middle; i < f1.Npx_inner; i++)
    {
        for (int j = 0; j < f1.Npy_inner; j++)
        {
            f1.IP[i][j]->iniset(1.2/gama,0,0,1/gama/(gama-1));
            f2.IP[i][j]->iniset(1.2/gama,0,0,1/gama/(gama-1));
        }
    }

    f1.boundarycondition(1);

    double T=0.05;
    double t=0;
    double dt;
    std::vector<double> dts;
    std::vector<double> Markt{0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
    std::stringstream ss;


    


    int k=0,kt=0;
    while (t<1)
    {
        f1.updatedt();
        dt=caldt(f1);
        dts.push_back(dt);

        //预估步
        pre(f1,f2,dt);

        //矫正步
        cor(f1,f2,dt);

        //插值
        inter(f1,f2);

        //输出
        t+=dt;
        k+=1;
        std::cout<<"k="<<k<<std::endl;
        std::cout<<"t="<<t<<std::endl;
        std::cout<<std::endl;

        if (t>Markt[kt])
        {
            outputrho(f1,"rho"+d2str(Markt[kt])+".txt");
            kt+=1;
        }
    }
    
    outputrho(f1,"rho.txt");
    outputdt(dts,"dt.txt");
}