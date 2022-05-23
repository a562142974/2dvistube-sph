#include"math.h"
#include<vector>
#include<fstream>
#include<string>
#include<iomanip>


const double gama=1.4;
const double Pr=0.73;
const double PI=3.1415926;

//本文的无量纲化方式会导致e=T
class particle{
    public:
    double x,y,d;
    double m;
    double miu=1,k=miu*gama/Pr;
    double *Vars;//u,v,p,T,rho,e,txx,txy,tyy,qx,qy
    double u,v,p,T,rho,e,txx,txy,tyy,qx,qy;
    double ux,uy,vx,vy,Tx,Ty;
    double dtrho,dtu,dtv,dte;

    int isinner=0;

    particle* next;
    //默认构造函数
    particle(double nx,double ny,double nd);

    //初始化函数,仅用于内部粒子
    void iniset(double rho,double u,double v,double e);

    //更新函数，仅用于内部粒子
    void update(double rho,double u,double v,double e);

    //设置函数,仅用于边界
    void set(double u,double v,double p,double T);

    //根据导数更新粘性通量
    void setflux(double ux,double uy,double vx,double vy,double Tx,double Ty);

    //直接赋值导数
    void assders(double ux,double uy,double vx,double vy,double Tx,double Ty);
    
    //直接赋值粘性通量
    void assflux(double txx,double txy,double tyy,double qx,double qy);

    //时间导数赋值
    void assdt(double ndtrho,double ndtu,double ndtv,double ndte);
};

//默认构造函数
particle::particle(double nx,double ny,double nd):\
                    x(nx),y(ny),d(nd){
    Vars=new double[11];
};


//iniset
void particle::iniset(double nrho,double nu,double nv,double ne){
    rho=nrho;u=nu;v=nv;e=ne;
    m=rho*pow(d,2);
    T=e;
    p=(gama-1)*rho*e;

    Vars[0]=u;
    Vars[1]=v;
    Vars[2]=p;
    Vars[3]=T;
    Vars[4]=rho;
    Vars[5]=e;
};

//update
void particle::update(double nrho,double nu,double nv,double ne){
    rho=nrho;u=nu;v=nv;e=ne;
    T=e;
    p=(gama-1)*rho*e;

    Vars[0]=u;
    Vars[1]=v;
    Vars[2]=p;
    Vars[3]=T;
    Vars[4]=rho;
    Vars[5]=e;
};

//set,边界
void particle::set(double nu,double nv,double np,double nT){
    u=nu;v=nv;p=np;T=nT;
    e=T;
    rho=p/(gama-1)/e;
    m=rho*pow(d,2);
    
    Vars[0]=u;
    Vars[1]=v;
    Vars[2]=p;
    Vars[3]=T;
    Vars[4]=rho;
    Vars[5]=e;
};

//setflux,
void particle::setflux(double ux,double uy,double vx,double vy,double Tx,double Ty){
    txx=4.0/3.0*miu*ux-2.0/3.0*miu*vy;
    tyy=4.0/3.0*miu*vy-2.0/3.0*miu*ux;
    txy=miu*(uy+vx);
    qx=k*Tx;
    qy=k*Ty;

    Vars[6]=txx;
    Vars[7]=txy;
    Vars[8]=tyy;
    Vars[9]=qx;
    Vars[10]=qy;
};



//直接赋值导数
void particle::assders(double nux,double nuy,double nvx,double nvy,double nTx,double nTy){
    ux=nux;
    uy=nuy;
    vx=nvx;
    vy=nvy;
    Tx=nTx;
    Ty=nTy;
};



//直接赋值粘性通量
void particle::assflux(double ntxx,double ntxy,double ntyy,double nqx,double nqy){
    txx=ntxx;
    txy=ntxy;
    tyy=ntyy;
    qx=nqx;
    qy=nqy;

    Vars[6]=txx;
    Vars[7]=txy;
    Vars[8]=tyy;
    Vars[9]=qx;
    Vars[10]=qy;
};

//时间导数赋值
void particle::assdt(double ndtrho,double ndtu,double ndtv,double ndte){
    dtrho=ndtrho;
    dtu=ndtu;
    dtv=ndtv;
    dte=ndte;
};

//
class box{
    public:
    double l,r,u,d;
    particle* next=nullptr;
    box()=default;
};




//SPH场
class SPHfield{
    public:
    double xl,xr,yd,yu,d,fac_d_h,h,k;
    int N_particles;

    public:
    double inf_h;
    double Re;
    int N_layers_v_x,N_layers_v_y,Nbx,Nby,Npx,Npy,il,ir,jd,ju,Npx_inner,Npy_inner;
    double length_box,height_box,alphad;
    
    box** AB;
    public:
    particle*** AP;

    public:
    particle*** IP;

    //构造函数,k为inf_h与h之比
    SPHfield(double nxl,double nxr,double nyd,double nyu,double nd,double nfac_d_h,double nk,double nRe);

    //caldt
    void updatedt();


    //核函数
    double W(double dx,double dy);

    //核函数x一阶导数
    double dxW(double dx,double dy);

    //核函数y一阶导数
    double dyW(double dx,double dy);

    //get p in box 
    void get_particle_in_box();

    //边界处理
    void boundarycondition(int flag);

    //正常核函数0阶导估计
    double SPHestimate(double x,double y,int n);

    //正则核函数0阶导估计
    double reSPHestimate(double x,double y,int n);

    //墙壁附近截断支持域的核函数估计
    double wallSPHestimate(double x,double y,int n);

    //用于上境外虚粒子估计
    double innerSPHestimate(double x,double y,int n);

    //x一阶导估计
    double SPHestimatedx(particle* p,int n);

    //y一阶导估计
    double SPHestimatedy(particle* p,int n);

    //1/rho*dA/dx估计
    double SPHestimate_rhodx(particle* p,int n);

    //1/rho*dA/dy估计
    double SPHestimate_rhody(particle* p,int n);

    //rho*dA/dx估计
    double SPHestimaterhodx(particle* p,int n);

    //rho*dA/dy估计
    double SPHestimaterhody(particle* p,int n);

    //phi/rho*dA/dx估计
    double SPHestimatephirhodx(particle* p,int n,int k);

    //phi/rho*dA/dy估计
    double SPHestimatephirhody(particle* p,int n,int k);

    //找邻居
    std::vector<particle*> findneis(double x,double y);

    //找流体粒子邻居，忽略虚拟粒子和壁面粒子
    std::vector<particle*> findfluidneis(double x,double y);

    //找内部区域粒子
    std::vector<particle*> findinnerneis(double x,double y);

    //找盒
    int findi_b(double x);
    int findj_b(double y);

};


//核函数,采用B样条
double SPHfield::W(double dx,double dy){
    double q=sqrt(pow(dx,2)+pow(dy,2))/h;
    if(q>=0 && q<1){
        return alphad*(2.0/3.0-pow(q,2)+0.5*pow(q,3));
    }else if(q>=1 && q<2){
        return alphad*1.0/6.0*pow(2-q,3);
    }else{
        return 0;
    };
};

//核函数x一阶导数
double SPHfield::dxW(double dx,double dy){
    double q=sqrt(pow(dx,2)+pow(dy,2))/h;
    if(q>=0 && q<1){
        return alphad*dx/pow(h,2)*(-2+3.0/2.0*q);
    }else if(q>=1 && q<2){
        return alphad*dx/pow(h,2)*(-0.5/q)*pow(q-2,2);
    }else{
        return 0;
    };
};

//核函数y一阶导数
double SPHfield::dyW(double dx,double dy){
    double q=sqrt(pow(dx,2)+pow(dy,2))/h;
    if(q>=0 && q<1){
        return alphad*dy/pow(h,2)*(-2+3.0/2.0*q);
    }else if(q>=1 && q<2){
        return alphad*dy/pow(h,2)*(-0.5/q)*pow(q-2,2);
    }else{
        return 0;
    };
};

SPHfield::SPHfield(double nxl,double nxr,double nyd,double nyu,double nd,double nfac_h_d,double nk,double nRe):\
xl(nxl),xr(nxr),yd(nyd),yu(nyu),d(nd),fac_d_h(nfac_h_d),k(nk),Re(nRe){
    //常量计算
    h=d*fac_d_h;
    alphad=15.0/(7.0*PI*pow(h,2));
    inf_h=k*h;
    Nbx=floor((xr-xl)/inf_h)+2;
    Nby=floor((yu-yd)/inf_h)+2;
    length_box=(xr-xl)/(Nbx-2);
    height_box=(yu-yd)/(Nby-2);
    N_layers_v_x=floor(length_box/d);
    N_layers_v_y=floor(height_box/d);
    Npx_inner=(xr-xl)/d-1;
    Npy_inner=(yu-yd)/d-1;
    Npx=Npx_inner+2+2*N_layers_v_x;
    Npy=Npy_inner+2+2*N_layers_v_y;
    N_particles=Npx*Npy;
    il=N_layers_v_x;ir=Npx_inner+2+N_layers_v_x-1;
    jd=N_layers_v_y;ju=Npy_inner+2+N_layers_v_y-1;
    //分配AP
    AP= new particle**[Npx];
    for (int i = 0; i < Npx; i++)
    {
        AP[i]=new particle*[Npy];
        for (int j = 0; j < Npy; j++)
        {
            AP[i][j]=new particle(xl-N_layers_v_x*d+i*d,yd-N_layers_v_y*d+j*d,d);
        }
    }
    //设置IP
    IP=new particle**[Npx_inner];
    for (int i = 0; i < Npx_inner; i++)
    {
        IP[i]=new particle*[Npy_inner];
        for (int j = 0; j<Npy_inner; j++)
        {
            IP[i][j]=AP[i+N_layers_v_x+1][j+N_layers_v_y+1];
            IP[i][j]->isinner=1;
        }
    }
    //处理上中镜外粒子的inner
    for (int i = il+1; i <= ir-1; i++)
    {
        for (int j = ju; j <Npy; j++)
        {
            AP[i][j]->isinner=1;
        }
    }
    //初始化boxes
    AB=new box*[Nbx];
    for (int i = 0; i < Nbx; i++)
    {
        AB[i]=new box[Nby];
        for (int j = 0; j < Nby; j++)
        {
            AB[i][j].l=xl-length_box+i*length_box;
            AB[i][j].d=yd-height_box+j*height_box;
            AB[i][j].r=AB[i][j].l+length_box;
            AB[i][j].u=AB[i][j].d+height_box;
        }
    }
    //装盒
    get_particle_in_box();
};

//caldt
void SPHfield::updatedt(){
    double x,y;
    particle* p;
    for (int i = 0; i < Npx_inner; i++)
    {
        for (int j = 0; j < Npy_inner; j++)
        {
            p=IP[i][j];
            x=p->x;y=p->y;
            p->assdt(-SPHestimaterhodx(p,0)\
            -SPHestimaterhody(p,1),\

            -SPHestimate_rhodx(p,2)\
            +1/Re*SPHestimate_rhodx(p,6)\
            +1/Re*SPHestimate_rhody(p,7),\

            -SPHestimate_rhody(p,2)\
            +1/Re*SPHestimate_rhodx(p,7)\
            +1/Re*SPHestimate_rhody(p,8),\

            -SPHestimatephirhodx(p,0,2)\
            -SPHestimatephirhody(p,1,2)\
            +1/Re*SPHestimatephirhodx(p,0,6)\
            +1/Re*SPHestimatephirhody(p,0,7)\
            +1/Re*SPHestimatephirhodx(p,1,7)\
            +1/Re*SPHestimatephirhody(p,1,8)\
            +1/Re*SPHestimate_rhodx(p,9)\
            +1/Re*SPHestimate_rhody(p,10));
        }
    }
};

//装盒
void SPHfield::get_particle_in_box(){
    
    for (int i = 0; i <Nbx; i++)
    {
        for (int j = 0; j < Nby; j++)
        {   
            AB[i][j].next=nullptr;
        }
    }
    
    int m,n;
    for (int i = 0; i < Npx; i++)
    {
        for (int j = 0; j < Npy; j++)
        {
            m=findi_b(AP[i][j]->x);
            n=findj_b(AP[i][j]->y);
            AP[i][j]->next=AB[m][n].next;
            AB[m][n].next=AP[i][j];
        }
    }
};

//找流体邻居
std::vector<particle*> SPHfield::findfluidneis(double x,double y){
    int i_b,j_b,i_s,i_e,j_s,j_e;
    int this_Nbx=Nbx;
    int this_Nby=Nby;
    i_b=findi_b(x);
    j_b=findj_b(y);

    if (i_b==Nbx-1)
    {
        i_s=Nbx-2;i_e=Nbx-2;
    }else if (i_b==0)
    {
        i_s=1;i_e=1;
    }else{
        i_s=i_b-1;i_e=i_b+1;
    }
    
    if (j_b==Nby-1)
    {
        j_s=Nby-2;j_e=Nby-1;
    }else if (j_b==0)
    {
        j_s=1;j_e=1;
    }else{
        j_s=j_b-1;
        j_e=j_b+1;
    }

    std::vector<particle*> neis;
    particle* p;
    for (int i = i_s; i <= i_e; i++)
    {
        for (int j = j_s; j <= j_e; j++)
        {
            p=AB[i][j].next;
        while (p!=nullptr)
            {
                if(p->isinner==1 && sqrt(pow(p->x-x,2)+pow(p->y-y,2))<inf_h) neis.push_back(p);
                p=p->next;
            }
        }
    }
    return neis;
};

//找普通邻居,点有可能在边缘的盒子上，所以相邻的盒子就可能数目要少一些
std::vector<particle*> SPHfield::findneis(double x,double y){
    int i_b,j_b;
    i_b=findi_b(x);
    j_b=findj_b(y);
    std::vector<particle*> neis;
    particle* p;
    for (int i = i_b-(i_b!=0); i <= i_b+(i_b!=Nbx-1); i++)
    {
        for (int j = j_b-(j_b!=0); j <=j_b+(j_b!=Nby-1); j++)
        {
            p=AB[i][j].next;
            while(p!=nullptr){
                if(sqrt(pow(p->x-x,2)+pow(p->y-y,2))<inf_h) neis.push_back(p);
                p=p->next;
            }
        }
    }
    return neis;
};


//找内部区域粒子
std::vector<particle*> SPHfield::findinnerneis(double x,double y){
    int i_b,j_b;
    i_b=findi_b(x);
    j_b=findj_b(y);
    std::vector<particle*> neis;
    particle* p;
    for (int i = i_b-(i_b!=0); i <= i_b+(i_b!=Nbx-1); i++)
    {
        for (int j = j_b-(j_b!=0); j <=j_b+(j_b!=Nby-1); j++)
        {
            p=AB[i][j].next;
            while(p!=nullptr){
                if(sqrt(pow(p->x-x,2)+pow(p->y-y,2))<inf_h && (p->x-xl>0 && xr-p->x>0 && yu-p->y>0 && p->y-yd>0)) neis.push_back(p);
                p=p->next;
            }
        }
    }
    return neis;
};

//找盒
inline int SPHfield::findi_b(double x){
    return floor((x-(xl-length_box))/length_box)-(floor((x-(xl-length_box))/length_box)>Nbx-1);
}

inline int SPHfield::findj_b(double y){
    return floor((y-(yd-height_box))/height_box)-(floor((y-(yd-height_box))/height_box)>Nby-1);
}

//SPH核函数估计
double SPHfield::SPHestimate(double x,double y,int n){
    std::vector<particle*> neis=findneis(x,y);
    double re=0;double L=0;
    for (int j = 0; j < neis.size(); j++)
    {
        re+=neis[j]->m/neis[j]->rho*W(x-neis[j]->x,y-neis[j]->y)*neis[j]->Vars[n];
        //L+=neis[j]->m/neis[j]->rho*W(x-neis[j]->x,y-neis[j]->y);
    }
    return re;
}

//SPH正则核函数估计
double SPHfield::reSPHestimate(double x,double y,int n){
    std::vector<particle*> neis=findneis(x,y);

    double re=0;double L=0;
    for (int j = 0; j < neis.size(); j++)
    {
        re+=neis[j]->m/neis[j]->rho*W(x-neis[j]->x,y-neis[j]->y)*neis[j]->Vars[n];
        L+=neis[j]->m/neis[j]->rho*W(x-neis[j]->x,y-neis[j]->y);
    }
    return re/L;
};

//壁面附近核函数估计
double SPHfield::wallSPHestimate(double x,double y,int n){
    std::vector<particle*> neis=findfluidneis(x,y);
    double re=0;
    double L=0;
    for (int j = 0; j < neis.size(); j++)
    {
        re+=neis[j]->m/neis[j]->rho*W(x-neis[j]->x,y-neis[j]->y)*neis[j]->Vars[n];
        L+=neis[j]->m/neis[j]->rho*W(x-neis[j]->x,y-neis[j]->y);
    }
    return re/L;
};

double SPHfield::innerSPHestimate(double x,double y,int n){
    std::vector<particle*> neis=findinnerneis(x,y);
    double re=0;
    double L=0;
    for (int j = 0; j < neis.size(); j++)
    {
        re+=neis[j]->m/neis[j]->rho*W(x-neis[j]->x,y-neis[j]->y)*neis[j]->Vars[n];
        L+=neis[j]->m/neis[j]->rho*W(x-neis[j]->x,y-neis[j]->y);
    }
    return re/L;
}

//x普通一阶导估计
double SPHfield::SPHestimatedx(particle* p,int n){
    std::vector<particle*> neis=findneis(p->x,p->y);
    double re=0;
    for (int j = 0; j < neis.size(); j++)
    {
        re+=neis[j]->m/neis[j]->rho*(neis[j]->Vars[n]-p->Vars[n])*dxW(p->x-neis[j]->x,p->y-neis[j]->y);
    }
    return re;
};



//y一阶导估计
double SPHfield::SPHestimatedy(particle* p,int n){
    std::vector<particle*> neis=findneis(p->x,p->y);
    double re=0;
    for (int j = 0; j < neis.size(); j++)
    {
        re+=neis[j]->m/neis[j]->rho*(neis[j]->Vars[n]-p->Vars[n])*dyW(p->x-neis[j]->x,p->y-neis[j]->y);
    }
    return re;
};

//1/rho*dA/dx估计
double SPHfield::SPHestimate_rhodx(particle* p,int n){
    std::vector<particle*> neis=findneis(p->x,p->y);
    double re=0;
    for (int j = 0; j < neis.size(); j++)
    {
        re+=neis[j]->m*dxW(p->x-neis[j]->x,p->y-neis[j]->y)*(p->Vars[n]/pow(p->rho,2)+neis[j]->Vars[n]/pow(neis[j]->rho,2));
    }
    return re;
};

//1/rho*dA/dy估计
double SPHfield::SPHestimate_rhody(particle* p,int n){
    std::vector<particle*> neis=findneis(p->x,p->y);
    double re=0;
    for (int j = 0; j < neis.size(); j++)
    {
        re+=neis[j]->m*dyW(p->x-neis[j]->x,p->y-neis[j]->y)*(p->Vars[n]/pow(p->rho,2)+neis[j]->Vars[n]/pow(neis[j]->rho,2));
    }
    return re;
};

//rho*dA/dx估计
double SPHfield::SPHestimaterhodx(particle* p,int n){
    std::vector<particle*> neis=findneis(p->x,p->y);
    double re=0;
    for (int j = 0; j < neis.size(); j++)
    {
        re+=neis[j]->m*(neis[j]->Vars[n]-p->Vars[n])*dxW(p->x-neis[j]->x,p->y-neis[j]->y);
    }
    return re;
};


//rho*dA/dy估计
double SPHfield::SPHestimaterhody(particle* p,int n){
    std::vector<particle*> neis=findneis(p->x,p->y);
    double re=0;
    for (int j = 0; j < neis.size(); j++)
    {
        re+=neis[j]->m*(neis[j]->Vars[n]-p->Vars[n])*dyW(p->x-neis[j]->x,p->y-neis[j]->y);
    }
    return re;
};


//phi/rho*dA/dx估计
double SPHfield::SPHestimatephirhodx(particle* p,int n,int k){
    std::vector<particle*> neis=findneis(p->x,p->y);
    double re=0;
    for (int j = 0; j < neis.size(); j++)
    {
        re+=0.5*neis[j]->m*(neis[j]->Vars[n]-p->Vars[n])*dxW(p->x-neis[j]->x,p->y-neis[j]->y)*\
        (p->Vars[k]/pow(p->rho,2)+neis[j]->Vars[k]/pow(neis[j]->rho,2));
    }
    return re;
};

//phi/rho*dA/dy估计
double SPHfield::SPHestimatephirhody(particle* p,int n,int k){
    std::vector<particle*> neis=findneis(p->x,p->y);
    double re=0;
    for (int j = 0; j < neis.size(); j++)
    {
        re+=0.5*neis[j]->m*(neis[j]->Vars[n]-p->Vars[n])*dyW(p->x-neis[j]->x,p->y-neis[j]->y)*\
        (p->Vars[k]/pow(p->rho,2)+neis[j]->Vars[k]/pow(neis[j]->rho,2));
    }
    return re;
};

void SPHfield::boundarycondition(int flag){
    //flag用于确定当前处理的是规则粒子还是不规则粒子的边界，采用不同的方法处理边界外虚粒子

    double p,T;
    //下固壁
    for (int i = il+1; i <= ir-1; i++)
    {
        int j=jd;
        p=wallSPHestimate(AP[i][j]->x,AP[i][j]->y+d,2);
        T=wallSPHestimate(AP[i][j]->x,AP[i][j]->y+d,3);
        AP[i][j]->set(0,0,p,T);
    }


    double x,y;    double u,v;
    //中上壁外虚粒子，镜像粒子
    if (flag==1)
    {
        for (int i = il+1; i <= ir-1; i++)
        {
            for (int j = ju+1; j <= Npy-1; j++)
            {
                AP[i][j]->set(AP[i][ju-(j-ju)]->u,-AP[i][ju-(j-ju)]->v,\
                            AP[i][ju-(j-ju)]->p,AP[i][ju-(j-ju)]->T);
            }
            //上镜壁
            AP[i][ju]->set(0.5*(AP[i][ju+1]->u+AP[i][ju-1]->u),\
                        0.5*(AP[i][ju+1]->u+AP[i][ju-1]->v),\
                        0.5*(AP[i][ju+1]->u+AP[i][ju-1]->p),\
                        0.5*(AP[i][ju+1]->u+AP[i][ju-1]->T));
        }
    }else{
        //中上外虚粒子，镜像粒子,严格对称
        std::vector<particle*> neis;
        for (int i = il+1; i <=ir-1; i++)
        {
            for (int j = ju+1; j < Npy; j++)
            {
                x=AP[i][j]->x;
                y=yu-(AP[i][j]->y-yu);
                //处理中上壁外虚粒子时，上镜面还未处理，所以只搜索计算区域内部粒子
                //找邻居
                u=innerSPHestimate(x,y,0);
                v=innerSPHestimate(x,y,1);
                p=innerSPHestimate(x,y,2);
                T=innerSPHestimate(x,y,3);
                AP[i][j]->set(u,-v,p,T);
            }
            u=innerSPHestimate(AP[i][ju]->x,AP[i][ju]->y-d,0);
            v=innerSPHestimate(AP[i][ju]->x,AP[i][ju]->y-d,1);
            p=innerSPHestimate(AP[i][ju]->x,AP[i][ju]->y-d,2);
            T=innerSPHestimate(AP[i][ju]->x,AP[i][ju]->y-d,3);
            AP[i][ju]->set(0.5*(AP[i][ju+1]->u+u),\
                            0.5*(AP[i][ju+1]->v+v),\
                            0.5*(AP[i][ju+1]->p+p),\
                            0.5*(AP[i][ju+1]->T+T));
        }
    }
    
    
    //固壁
    
    //左固壁
    for (int j = jd; j <= Npy-1; j++)
    {
        int i=il;
        p=wallSPHestimate(AP[i][j]->x+d,AP[i][j]->y,2);
        T=wallSPHestimate(AP[i][j]->x+d,AP[i][j]->y,3);
        AP[i][j]->set(0,0,p,T);
    }

    //右固壁
    for (int j = jd; j <= Npy-1; j++)
    {
        int i=ir;
        p=wallSPHestimate(AP[i][j]->x-d,AP[i][j]->y,2);
        T=wallSPHestimate(AP[i][j]->x-d,AP[i][j]->y,3);
        AP[i][j]->set(0,0,p,T);
    }


    
    //角点，直接等于墙角点
    //左下
    for (int i = 0; i <= il; i++)
    {
        for (int j = 0; j <= jd; j++)
        {
            if(i!=il || j!=jd) AP[i][j]->set(AP[il][jd]->u,AP[il][jd]->v,AP[il][jd]->p,AP[il][jd]->T);
        }
    }

    //右下
    for (int i = ir; i <= Npx-1; i++)
    {
        for (int j = 0; j <= jd; j++)
        {
            if(i!=ir || j!=jd) AP[i][j]->set(AP[ir][jd]->u,AP[ir][jd]->v,AP[ir][jd]->p,AP[ir][jd]->T);
        }
    }

    //左上
    for (int i = 0; i <= il-1; i++)
    {
        for (int j = ju; j <= Npy-1; j++)
        {
            if(i!=il || j!=ju) AP[i][j]->set(AP[il][ju]->u,AP[il][ju]->v,AP[il][ju]->p,AP[il][ju]->T);
        }
    }

    //右上
    for (int i = ir+1; i <=Npx-1; i++)
    {
        for (int j = ju; j <= Npy-1; j++)
        {
            if(i!=ir || j!=ju) AP[i][j]->set(AP[ir][ju]->u,AP[ir][ju]->v,AP[ir][ju]->p,AP[ir][ju]->T);
        }
    }
    

    //壁面外虚粒子
    double uf,vf,pf,Tf;
    double dg,df;


    if (flag==1)
    {
            //左壁外虚粒子
        for (int j = jd+1; j <= ju-1; j++)
        {
            for (int i = 0; i <= il-1; i++)
            {
                dg=xl-(AP[i][j]->x);
                u=(1+dg/d)*AP[il][j]->u-dg/d*AP[il+1][j]->u;
                v=(1+dg/d)*AP[il][j]->v-dg/d*AP[il+1][j]->v;
                p=(1+dg/d)*AP[il][j]->p-dg/d*AP[il+1][j]->p;
                T=(1+dg/d)*AP[il][j]->T-dg/d*AP[il+1][j]->T;
                AP[i][j]->set(u,v,p,T);
            }   
        }
        
        //右壁外虚粒子
        for (int j = jd+1; j <= ju-1; j++)
        {
            for (int i = ir+1; i <= Npx-1; i++)
            {
                dg=(AP[i][j]->x)-xr;
                u=(1+dg/d)*AP[il][j]->u-dg/d*AP[ir-1][j]->u;
                v=(1+dg/d)*AP[il][j]->v-dg/d*AP[ir-1][j]->v;
                p=(1+dg/d)*AP[il][j]->p-dg/d*AP[ir-1][j]->p;
                T=(1+dg/d)*AP[il][j]->T-dg/d*AP[ir-1][j]->T;
                AP[i][j]->set(u,v,p,T);
            }   
        }

        //下壁外虚粒子
        for (int i = il+1; i <= ir-1; i++)
        {
            for (int j = 0; j <= jd-1; j++)
            {

            dg=yd-AP[i][j]->y;
                u=(1+dg/d)*AP[i][jd]->u-dg/d*AP[i][jd+1]->u;
                v=(1+dg/d)*AP[i][jd]->v-dg/d*AP[i][jd+1]->u;
                p=(1+dg/d)*AP[i][jd]->p-dg/d*AP[i][jd+1]->u;
                T=(1+dg/d)*AP[i][jd]->T-dg/d*AP[i][jd+1]->u;
                AP[i][j]->set(u,v,p,T);
            }
        }
    }else{
            //左壁外虚粒子
        for (int j = jd+1; j <= ju-1; j++)
        {
            y=AP[il][j]->y;
                    uf=wallSPHestimate(xl+d,y,0);
                    vf=wallSPHestimate(xl+d,y,1);
                    pf=wallSPHestimate(xl+d,y,2);
                    Tf=wallSPHestimate(xl+d,y,3);
            for (int i = 0; i <= il-1; i++)
            {
                dg=xl-(AP[i][j]->x);
                u=(1+dg/d)*AP[il][j]->u-dg/d*uf;
                v=(1+dg/d)*AP[il][j]->v-dg/d*vf;
                p=(1+dg/d)*AP[il][j]->p-dg/d*pf;
                T=(1+dg/d)*AP[il][j]->T-dg/d*Tf;
                AP[i][j]->set(u,v,p,T);
            }   
        }
        
        //右壁外虚粒子
        for (int j = jd+1; j <= ju-1; j++)
        {
            y=AP[ir][j]->y;
                    uf=wallSPHestimate(xr-d,y,0);
                    vf=wallSPHestimate(xr-d,y,1);
                    pf=wallSPHestimate(xr-d,y,2);
                    Tf=wallSPHestimate(xr-d,y,3);
            for (int i = ir+1; i <= Npx-1; i++)
            {
                dg=(AP[i][j]->x)-xr;
                u=(1+dg/d)*AP[il][j]->u-dg/d*uf;
                v=(1+dg/d)*AP[il][j]->v-dg/d*vf;
                p=(1+dg/d)*AP[il][j]->p-dg/d*pf;
                T=(1+dg/d)*AP[il][j]->T-dg/d*Tf;
                AP[i][j]->set(u,v,p,T);
            }   
        }

        //下壁外虚粒子
        for (int i = il+1; i <= ir-1; i++)
        {
            x=AP[i][jd]->x;
            uf=wallSPHestimate(x,yd+d,0);
            vf=wallSPHestimate(x,yd+d,1);
            pf=wallSPHestimate(x,yd+d,2);
            Tf=wallSPHestimate(x,yd+d,3);
            for (int j = 0; j <= jd-1; j++)
            {

            dg=yd-AP[i][j]->y;
                
                u=(1+dg/d)*AP[i][jd]->u-dg/d*uf;
                v=(1+dg/d)*AP[i][jd]->v-dg/d*vf;
                p=(1+dg/d)*AP[i][jd]->p-dg/d*pf;
                T=(1+dg/d)*AP[i][jd]->T-dg/d*Tf;
                AP[i][j]->set(u,v,p,T);
            }
        }
    }
    

    //计算内部粒子导数
    double ux,uy,vx,vy,Tx,Ty;
    for (int i = il; i <=ir; i++)
    {
        for (int j = jd; j <=ju; j++)
        {
            ux=SPHestimatedx(AP[i][j],0);
            uy=SPHestimatedy(AP[i][j],0);
            vx=SPHestimatedx(AP[i][j],1);
            vy=SPHestimatedy(AP[i][j],1);
            //这里强制T边界上梯度为0
            Tx=SPHestimatedx(AP[i][j],3)-(i==il || i==ir)*SPHestimatedx(AP[i][j],3);
            Ty=SPHestimatedy(AP[i][j],3)-(j==jd)*SPHestimatedy(AP[i][j],3);
            AP[i][j]->setflux(ux,uy,vx,vy,Tx,Ty);
            AP[i][j]->assders(ux,uy,vx,vy,Tx,Ty);
            //用于测试
        }
    }


    //更新边界上txx,txy,tyy,qx,qy
    //上镜面外粒子
    double txx,txy,tyy,qx,qy;
    if (flag==1)
    {
            for (int i = il+1; i <= ir-1; i++)
        {
            for (int j = ju+1; j < Npy; j++)
            {
                ux=SPHestimatedx(AP[i][ju-(j-ju)],0);
                uy=SPHestimatedy(AP[i][ju-(j-ju)],0);
                vx=SPHestimatedx(AP[i][ju-(j-ju)],1);
                vy=SPHestimatedy(AP[i][ju-(j-ju)],1);
                Tx=SPHestimatedx(AP[i][ju-(j-ju)],3);
                Ty=SPHestimatedy(AP[i][ju-(j-ju)],3);
                AP[i][j]->setflux(ux,-uy,-vx,vy,Tx,-Ty);
            }
        }
        
        //上镜面粒子
        for (int i = il+1; i <= ir-1; i++)
        {
            int j=ju;
            AP[i][j]->assflux(0.5*(AP[i][j+1]->Vars[6]+AP[i][j-1]->Vars[6]),\
                            0.5*(AP[i][j+1]->Vars[7]+AP[i][j-1]->Vars[7]),\
                            0.5*(AP[i][j+1]->Vars[8]+AP[i][j-1]->Vars[8]),\
                            0.5*(AP[i][j+1]->Vars[9]+AP[i][j-1]->Vars[9]),\
                            0.5*(AP[i][j+1]->Vars[10]+AP[i][j-1]->Vars[10]));
        }
    }else{
            //更新边界上txx,txy,tyy,qx,qy
        //上镜面外粒子
        
        for (int i = il+1; i < ir; i++)
        {
            x=AP[i][ju]->x;
            for (int j = ju+1; j < Npy; j++)
            {
                y=yu-(AP[i][j]->y-yu);
                txx=innerSPHestimate(x,y,6);
                txy=innerSPHestimate(x,y,7);
                tyy=innerSPHestimate(x,y,8);
                qx=innerSPHestimate(x,y,9);
                qy=innerSPHestimate(x,y,10);
                AP[i][j]->assflux(txx,-txy,tyy,qx,-qy);
            }
            txx=innerSPHestimate(x,yu-d,6);
            txy=innerSPHestimate(x,yu-d,7);
            tyy=innerSPHestimate(x,yu-d,8);
            qx=innerSPHestimate(x,yu-d,9);
            qy=innerSPHestimate(x,yu-d,10);
            AP[i][ju]->assflux(0.5*(txx+AP[i][ju+1]->txx),\
            0.5*(txy+AP[i][ju+1]->txy),\
            0.5*(tyy+AP[i][ju+1]->tyy),\
            0.5*(qx+AP[i][ju+1]->qx),\
            0.5*(qy+AP[i][ju+1]->qy));
        }
    }
    

    //固壁粒子

    //左右固壁
    
    for (int j = jd; j < Npy; j++)
    {   
        x=AP[il][jd]->x;
        //左固壁
        AP[il][j]->assflux(wallSPHestimate(x,AP[il][j]->y,6),\
        wallSPHestimate(x,AP[il][j]->y,7),\
        wallSPHestimate(x,AP[il][j]->y,8),\
        wallSPHestimate(x,AP[il][j]->y,9),\
        wallSPHestimate(x,AP[il][j]->y,10));

        x=AP[ir][jd]->x;
        //右固壁
        AP[ir][j]->assflux(wallSPHestimate(x,AP[ir][j]->y,6),\
        wallSPHestimate(x,AP[ir][j]->y,7),\
        wallSPHestimate(x,AP[ir][j]->y,8),\
        wallSPHestimate(x,AP[ir][j]->y,9),\
        wallSPHestimate(x,AP[ir][j]->y,10));
    }

    //下固壁
    for (int i = il+1; i <= ir-1; i++)
    {
        int j=jd;
        y=AP[il][jd]->y;
        AP[i][j]->assflux(\
        wallSPHestimate(AP[i][j]->x,y,6),\
        wallSPHestimate(AP[i][j]->x,y,7),\
        wallSPHestimate(AP[i][j]->x,y,8),\
        wallSPHestimate(AP[i][j]->x,y,9),\
        wallSPHestimate(AP[i][j]->x,y,10));
    }

    //角点
    //左上角
    for (int i = 0; i <= il-1; i++)
    {
        for (int j = ju; j <= Npy-1; j++)
        {
            AP[i][j]->assflux(AP[il][ju]->Vars[6],\
                            AP[il][ju]->Vars[7],\
                            AP[il][ju]->Vars[8],\
                            AP[il][ju]->Vars[9],\
                            AP[il][ju]->Vars[10]);
        }
    }

    //右上角
    for (int i = ir+1; i < Npx; i++)
    {
        for (int j = ju; j < Npy; j++)
        {
            AP[i][j]->assflux(AP[ir][ju]->Vars[6],\
                            AP[ir][ju]->Vars[7],\
                            AP[ir][ju]->Vars[8],\
                            AP[ir][ju]->Vars[9],\
                            AP[ir][ju]->Vars[10]);
        }
    }

    //左下角
    for (int i = 0; i <= il; i++)
    {
        for (int j = 0; j <= jd; j++)
        {
            if(i!=il || j!=jd){
                AP[i][j]->assflux(AP[il][jd]->Vars[6],\
                            AP[il][jd]->Vars[7],\
                            AP[il][jd]->Vars[8],\
                            AP[il][jd]->Vars[9],\
                            AP[il][jd]->Vars[10]);
            }
        }
    }

    //右下角
    for (int i = ir; i < Npx; i++)
    {
        for (int j = 0; j <= jd; j++)
        {
            AP[i][j]->assflux(AP[ir][jd]->Vars[6],\
                            AP[ir][jd]->Vars[7],\
                            AP[ir][jd]->Vars[8],\
                            AP[ir][jd]->Vars[9],\
                            AP[ir][jd]->Vars[10]);
        }
    }



    //壁外虚粒子
    if (flag==1)
    {
        //左壁外
        for (int j = jd+1; j <= ju-1; j++)
        {
            for (int i = 0; i <= il-1; i++)
            {
                dg=xl-(AP[i][j]->x);
                AP[i][j]->assflux((1+dg/d)*AP[il][j]->txx-dg/d*AP[il+1][j]->txx,\
                (1+dg/d)*AP[il][j]->txy-dg/d*AP[il+1][j]->txy,\
                (1+dg/d)*AP[il][j]->tyy-dg/d*AP[il+1][j]->tyy,\
                (1+dg/d)*AP[il][j]->qx-dg/d*AP[il+1][j]->qx,\
                (1+dg/d)*AP[il][j]->qy-dg/d*AP[il+1][j]->qy);
            }   
        }
        
        //右壁外虚粒子
        for (int j = jd+1; j <= ju-1; j++)
        {
            y=AP[ir][j]->y;
                    txx=wallSPHestimate(xr-d,y,6);
                    txy=wallSPHestimate(xr-d,y,7);
                    tyy=wallSPHestimate(xr-d,y,8);
                    qx=wallSPHestimate(xr-d,y,9);
                    qy=wallSPHestimate(xr-d,y,10);
            for (int i = ir+1; i <= Npx-1; i++)
            {
                dg=(AP[i][j]->x)-xr;
                AP[i][j]->assflux((1+dg/d)*AP[ir][j]->txx-dg/d*AP[ir-1][j]->txx,\
                (1+dg/d)*AP[ir][j]->txy-dg/d*AP[ir-1][j]->txy,\
                (1+dg/d)*AP[ir][j]->tyy-dg/d*AP[ir-1][j]->tyy,\
                (1+dg/d)*AP[ir][j]->qx-dg/d*AP[ir-1][j]->qx,\
                (1+dg/d)*AP[ir][j]->qy-dg/d*AP[ir-1][j]->qy);
            }   
        }

        //下壁外虚粒子
        for (int i = il+1; i <= ir-1; i++)
        {
            x=AP[i][jd]->x;
                    txx=wallSPHestimate(x,yd+d,6);
                    txy=wallSPHestimate(x,yd+d,7);
                    tyy=wallSPHestimate(x,yd+d,8);
                    qx=wallSPHestimate(x,yd+d,9);
                    qy=wallSPHestimate(x,yd+d,10);
            for (int j = 0; j <= jd-1; j++)
            {

            dg=yd-AP[i][j]->y;
            AP[i][j]->assflux((1+dg/d)*AP[i][jd]->txx-dg/d*AP[i][jd+1]->txx,\
                (1+dg/d)*AP[i][jd]->txy-dg/d*AP[i][jd+1]->txy,\
                (1+dg/d)*AP[i][jd]->tyy-dg/d*AP[i][jd+1]->tyy,\
                (1+dg/d)*AP[i][jd]->qx-dg/d*AP[i][jd+1]->qx,\
                (1+dg/d)*AP[i][jd]->qy-dg/d*AP[i][jd+1]->qy);
            }
        }
    }else{
            //壁外虚粒子
        //左壁外虚粒子
        for (int j = jd+1; j <= ju-1; j++)
        {
            y=AP[il][j]->y;
                    txx=wallSPHestimate(xl+d,y,6);
                    txy=wallSPHestimate(xl+d,y,7);
                    tyy=wallSPHestimate(xl+d,y,8);
                    qx=wallSPHestimate(xl+d,y,9);
                    qy=wallSPHestimate(xl+d,y,10);
            for (int i = 0; i <= il-1; i++)
            {
                dg=xl-(AP[i][j]->x);
                AP[i][j]->assflux((1+dg/d)*AP[il][j]->txx-dg/d*txx,\
                (1+dg/d)*AP[il][j]->txy-dg/d*txy,\
                (1+dg/d)*AP[il][j]->tyy-dg/d*tyy,\
                (1+dg/d)*AP[il][j]->qx-dg/d*qx,\
                (1+dg/d)*AP[il][j]->qy-dg/d*qy);
            }   
        }
        
        //右壁外虚粒子
        for (int j = jd+1; j <= ju-1; j++)
        {
            y=AP[ir][j]->y;
                    txx=wallSPHestimate(xr-d,y,6);
                    txy=wallSPHestimate(xr-d,y,7);
                    tyy=wallSPHestimate(xr-d,y,8);
                    qx=wallSPHestimate(xr-d,y,9);
                    qy=wallSPHestimate(xr-d,y,10);
            for (int i = ir+1; i <= Npx-1; i++)
            {
                dg=(AP[i][j]->x)-xr;
                AP[i][j]->assflux((1+dg/d)*AP[ir][j]->txx-dg/d*txx,\
                (1+dg/d)*AP[ir][j]->txy-dg/d*txy,\
                (1+dg/d)*AP[ir][j]->tyy-dg/d*tyy,\
                (1+dg/d)*AP[ir][j]->qx-dg/d*qx,\
                (1+dg/d)*AP[ir][j]->qy-dg/d*qy);
            }   
        }

        //下壁外虚粒子
        for (int i = il+1; i <= ir-1; i++)
        {
            x=AP[i][jd]->x;
                    txx=wallSPHestimate(x,yd+d,6);
                    txy=wallSPHestimate(x,yd+d,7);
                    tyy=wallSPHestimate(x,yd+d,8);
                    qx=wallSPHestimate(x,yd+d,9);
                    qy=wallSPHestimate(x,yd+d,10);
            for (int j = 0; j <= jd-1; j++)
            {

            dg=yd-AP[i][j]->y;
            AP[i][j]->assflux((1+dg/d)*AP[i][jd]->txx-dg/d*txx,\
                (1+dg/d)*AP[i][jd]->txy-dg/d*txy,\
                (1+dg/d)*AP[i][jd]->tyy-dg/d*tyy,\
                (1+dg/d)*AP[i][jd]->qx-dg/d*qx,\
                (1+dg/d)*AP[i][jd]->qy-dg/d*qy);
            }
        }
    }
};


void output(SPHfield& f){

    std::ofstream out("re.txt");
    for (int j = 0; j < f.Npy_inner; j++)
    {
        for (int i = 0; i < f.Npx_inner; i++)
        {
            out<<f.IP[i][j]->rho;
        }
        out<<std::endl;
    }
    
};

void outputIVars(SPHfield& f,std::string filename){

    std::ofstream out(filename);

    std::vector<std::string> names{"u","v","p","T","rho","e","txx","txy","tyy","qx","qy"};

    for (int k = 0; k < 11; k++)
    { 
        out<<names[k]<<std::endl;
        for (int j = 0; j < f.Npy_inner; j++)
        {
            for (int i = 0; i < f.Npx_inner; i++)
            {
                out<<f.IP[i][j]->Vars[k];
            }
            out<<std::endl;
        }
    }

    //
    out<<"dtrho"<<std::endl;
    for (int j = 0; j < f.Npy_inner; j++)
        {
            for (int i = 0; i < f.Npx_inner; i++)
            {
                out<<f.IP[i][j]->dtrho;
            }
            out<<std::endl;
        }

    
    out<<"dtu"<<std::endl;
    for (int j = 0; j < f.Npy_inner; j++)
        {
            for (int i = 0; i < f.Npx_inner; i++)
            {
                out<<f.IP[i][j]->dtu;
            }
            out<<std::endl;
        }


    out<<"dtv"<<std::endl;
    for (int j = 0; j < f.Npy_inner; j++)
        {
            for (int i = 0; i < f.Npx_inner; i++)
            {
                out<<f.IP[i][j]->dtv;
            }
            out<<std::endl;
        }


    out<<"dte"<<std::endl;
    for (int j = 0; j < f.Npy_inner; j++)
        {
            for (int i = 0; i < f.Npx_inner; i++)
            {
                out<<f.IP[i][j]->dte;
            }
            out<<std::endl;
        }
};

void outputAVars(SPHfield& f,std::string filename){

    std::ofstream out(filename);

    std::vector<std::string> names{"u","v","p","T","rho","e","txx","txy","tyy","qx","qy"};

    for (int k = 0; k < 11; k++)
    { 
        out<<names[k]<<std::endl;
        for (int j = 0; j < f.Npy; j++)
        {
            for (int i = 0; i < f.Npx; i++)
            {
                out<<std::setiosflags(std::ios::left)<<std::setw(15)<<f.AP[i][f.Npy-1-j]->Vars[k] <<" ";
            }
            out<<std::endl;
        }
    }


        out<<"ux"<<std::endl;
    for (int j = 0; j < f.Npy; j++)
        {
            for (int i = 0; i < f.Npx; i++)
            {
                out<<std::setiosflags(std::ios::left)<<std::setw(15)<<f.AP[i][f.Npy-1-j]->ux<<" ";
            }
            out<<std::endl;
        }

    
        out<<"uy"<<std::endl;
    for (int j = 0; j < f.Npy; j++)
        {
            for (int i = 0; i < f.Npx; i++)
            {
                out<<std::setiosflags(std::ios::left)<<std::setw(15)<<f.AP[i][f.Npy-1-j]->uy<<" ";
            }
            out<<std::endl;
        }

    
        out<<"vx"<<std::endl;
    for (int j = 0; j < f.Npy; j++)
        {
            for (int i = 0; i < f.Npx; i++)
            {
                out<<std::setiosflags(std::ios::left)<<std::setw(15)<<f.AP[i][f.Npy-1-j]->vx<<" ";
            }
            out<<std::endl;
        }

        out<<"vy"<<std::endl;
    for (int j = 0; j < f.Npy; j++)
        {
            for (int i = 0; i < f.Npx; i++)
            {
                out<<std::setiosflags(std::ios::left)<<std::setw(15)<<f.AP[i][f.Npy-1-j]->vy<<" ";
            }
            out<<std::endl;
        }


        out<<"Tx"<<std::endl;
    for (int j = 0; j < f.Npy; j++)
        {
            for (int i = 0; i < f.Npx; i++)
            {
                out<<std::setiosflags(std::ios::left)<<std::setw(15)<<f.AP[i][f.Npy-1-j]->Tx<<" ";
            }
            out<<std::endl;
        }


        out<<"Ty"<<std::endl;
    for (int j = 0; j < f.Npy; j++)
        {
            for (int i = 0; i < f.Npx; i++)
            {
                out<<std::setiosflags(std::ios::left)<<std::setw(15)<<f.AP[i][f.Npy-1-j]->Ty<<" ";
            }
            out<<std::endl;
        }


    //
    out<<"dtrho"<<std::endl;
    for (int j = 0; j < f.Npy; j++)
        {
            for (int i = 0; i < f.Npx; i++)
            {
                out<<std::setiosflags(std::ios::left)<<std::setw(15)<<f.AP[i][f.Npy-1-j]->dtrho<<" ";
            }
            out<<std::endl;
        }

    
    out<<"dtu"<<std::endl;
    for (int j = 0; j < f.Npy; j++)
        {
            for (int i = 0; i < f.Npx; i++)
            {
                out<<std::setiosflags(std::ios::left)<<std::setw(15)<<f.AP[i][f.Npy-1-j]->dtu<<" ";
            }
            out<<std::endl;
        }


    out<<"dtv"<<std::endl;
    for (int j = 0; j < f.Npy; j++)
        {
            for (int i = 0; i < f.Npx; i++)
            {
                out<<std::setiosflags(std::ios::left)<<std::setw(15)<<f.AP[i][f.Npy-1-j]->dtv<<" ";
            }
            out<<std::endl;
        }


    out<<"dte"<<std::endl;
    for (int j = 0; j < f.Npy; j++)
        {
            for (int i = 0; i < f.Npx; i++)
            {
                out<<std::setiosflags(std::ios::left)<<std::setw(15)<<f.AP[i][f.Npy-1-j]->dte<<" ";
            }
            out<<std::endl;
        }
};

//输出该SPHfield的所有点以及某个点的影响域以及背景网格信息的txt文件，交给python脚本处理
void outputpositions(SPHfield & f,int i,int j,std::string filename){
    //i,j 找邻居
    particle* ptofind=f.AP[i][j];
    std::vector<particle*> Neis=f.findneis(ptofind->x,ptofind->y);


    std::ofstream re(filename);
    re<<f.xl<<" "<<f.xr<<" "<<f.yd<<" "<<f.yu<<" "<<f.inf_h<<" "<<f.N_particles<<" "<<Neis.size()<<" "<<f.Nbx+1<<" "<<f.Nby+1<<std::endl;//常量


    for (int i = 0; i < f.Npx; i++)
    {
        for (int j = 0; j < f.Npy; j++)
        {
            re<<f.AP[i][j]->x<<" "<<f.AP[i][j]->y<<" "<<f.findi_b(f.AP[i][j]->x)<<" "<<f.findj_b(f.AP[i][j]->y)<<std::endl;
        }
        
    }
    

    for (int i = 0; i < Neis.size(); i++)
    {
        re<<Neis[i]->x<<" "<<Neis[i]->y<<" "<<f.findi_b(Neis[i]->x) <<" "<<f.findi_b(Neis[i]->y) <<std::endl;
    }
    re<<ptofind->x<<" "<<ptofind->y<<" "<<std::endl;
    for (int i = 0; i <f.Nbx; i++)
    {
        re<<f.AB[i][0].l<<std::endl;
    }
    re<<f.AB[f.Nbx-1][0].r<<std::endl;
    for (int j = 0; j < f.Nby; j++)
    {
        re<<f.AB[0][j].d<<std::endl;
    }
    re<<f.AB[0][f.Nby-1].u<<std::endl;
    re.close();
};

void outputparameters(SPHfield &f,std::string filename){
    std::ofstream re(filename);
    re<<"il="<<f.il<<std::endl;
    re<<"ir="<<f.ir<<std::endl;
    re<<"jd="<<f.jd<<std::endl;
    re<<"ju="<<f.ju<<std::endl;
    re<<"N_layers_v_x="<<f.N_layers_v_x<<std::endl;
    re<<"N_layers_v_y="<<f.N_layers_v_y<<std::endl;
    re<<"Nbx="<<f.Nbx<<std::endl;
    re<<"Nby="<<f.Nby<<std::endl;
    re<<"Npx="<<f.Npx<<std::endl;
    re<<"Npy="<<f.Npy<<std::endl;
    re<<"Npx_inner="<<f.Npx_inner<<std::endl;
    re<<"Npy_inner="<<f.Npy_inner<<std::endl;
    re<<"height_box="<<f.height_box<<std::endl;
    re<<"length_box="<<f.length_box<<std::endl;
    re<<"inf_h="<<f.inf_h<<std::endl;
}

void outputrho(SPHfield &f,std::string filename){
std::ofstream out(filename);
    out<<f.Npx_inner<<" "<<f.Npy_inner<<std::endl;;
    //位置
    for (int i = 0; i < f.Npx_inner; i++)
    {
        out<<f.IP[i][f.jd]->x<<" ";
    }
    out<<std::endl;
    for (int j = 0; j < f.Npy_inner; j++)
    {
        out<<f.IP[f.il][j]->y<<" ";
    }
    out<<std::endl;

    //密度值
    for (int j = 0; j < f.Npy_inner; j++)
    {
        for (int i = 0; i < f.Npx_inner; i++)
        {
            out<<std::setiosflags(std::ios::left)<<std::setw(15)<<f.IP[i][f.Npy_inner-1-j]->rho <<" ";
        }
        out<<std::endl;
    }
    out.close();
}


void outputdt(std::vector<double> dts,std::string filename){

    std::ofstream out(filename);

    for (int i = 0; i < dts.size(); i++)
    {
        out<<i<<" "<<dts[i]<<std::endl;
    }
    out.close();
    
}