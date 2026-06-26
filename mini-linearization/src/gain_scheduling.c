#include "gain_scheduling.h"
#include <string.h>
#include <assert.h>
#include <float.h>

GainSchedule *create_gain_schedule(int n_entries,int n_sched_vars){
    GainSchedule *gs=malloc(sizeof(GainSchedule)); memset(gs,0,sizeof(GainSchedule));
    gs->n_entries=n_entries;gs->n_scheduling_vars=n_sched_vars;
    gs->entries=calloc((size_t)n_entries,sizeof(GainScheduleEntry));
    gs->hysteresis_band=0.1; return gs;
}
void add_gain_schedule_entry(GainSchedule *gs,int index,const double *sv,
    const LinearizedSystem *sys,const double *K,const double *L,int ns,int ni,int no){
    GainScheduleEntry *e=&gs->entries[index]; e->n_scheduling_vars=gs->n_scheduling_vars;
    memcpy(e->scheduling_value,sv,(size_t)gs->n_scheduling_vars*sizeof(double));
    e->sys=(LinearizedSystem*)sys; e->n_states=ns;e->n_inputs=ni;e->n_outputs=no;
    e->K=malloc((size_t)(ni*ns)*sizeof(double)); memcpy(e->K,K,(size_t)(ni*ns)*sizeof(double));
    if(L){e->L=malloc((size_t)(ns*no)*sizeof(double));memcpy(e->L,L,(size_t)(ns*no)*sizeof(double));}
    if(gs->n_states==0){gs->n_states=ns;gs->n_inputs=ni;gs->n_outputs=no;}
}
void interpolate_gain_schedule(const GainSchedule *gs,const double *sv,
    double *K_out,double *L_out,GSInterpMethod method){
    int d=gs->n_scheduling_vars,ns=gs->n_states,ni=gs->n_inputs,no=gs->n_outputs;
    int ne=gs->n_entries; if(method==GS_NEAREST_NEIGHBOR||ne<=1){
        double bd=INFINITY;int best=0;for(int i=0;i<ne;i++){double dist=0;
            for(int j=0;j<d;j++){double df=sv[j]-gs->entries[i].scheduling_value[j];dist+=df*df;}
            if(dist<bd){bd=dist;best=i;}}
        memcpy(K_out,gs->entries[best].K,(size_t)(ni*ns)*sizeof(double));
        if(L_out&&gs->entries[best].L)memcpy(L_out,gs->entries[best].L,(size_t)(ns*no)*sizeof(double));
        return;
    }
    memset(K_out,0,(size_t)(ni*ns)*sizeof(double));
    if(L_out)memset(L_out,0,(size_t)(ns*no)*sizeof(double)); double ws=0;
    for(int i=0;i<ne;i++){double w=1.0;
        if(method==GS_BLENDING_FUNCTION){double d2=0;
            for(int j=0;j<d;j++){double df=sv[j]-gs->entries[i].scheduling_value[j];d2+=df*df;}
            w=exp(-d2/0.02);}
        else{for(int j=0;j<d;j++){double dv=fabs(sv[j]-gs->entries[i].scheduling_value[j]);
            w*=1.0/(1.0+dv*dv);}}
        for(int r=0;r<ni;r++)for(int c=0;c<ns;c++)K_out[r*ns+c]+=w*gs->entries[i].K[r*ns+c];
        if(L_out&&gs->entries[i].L)for(int r=0;r<ns;r++)for(int c=0;c<no;c++)L_out[r*no+c]+=w*gs->entries[i].L[r*no+c];
        ws+=w;
    }
    if(ws>1e-15){for(int i=0;i<ni*ns;i++)K_out[i]/=ws;
        if(L_out)for(int i=0;i<ns*no;i++)L_out[i]/=ws;}
}
GainScheduledController *create_gain_scheduled_controller(GainSchedule *gs,GSInterpMethod m,int ns){
    GainScheduledController *c=malloc(sizeof(GainScheduledController)); memset(c,0,sizeof(GainScheduledController));
    c->schedule=gs;c->method=m;c->current_K=malloc((size_t)(gs->n_inputs*gs->n_states)*sizeof(double));
    c->current_L=malloc((size_t)(gs->n_states*gs->n_outputs)*sizeof(double));
    c->observer_state=calloc((size_t)ns,sizeof(double)); return c;
}
double compute_gain_scheduled_control(GainScheduledController *ctrl,const double *sv,
    const double *y,double dt,const double *x_eq,const double *u_eq){
    GainSchedule *gs=ctrl->schedule; int ns=gs->n_states,ni=gs->n_inputs;
    interpolate_gain_schedule(gs,sv,ctrl->current_K,ctrl->current_L,ctrl->method);
    double u=0;
    if(y&&ctrl->current_L){
        for(int i=0;i<ni;i++){double s=0;
            for(int j=0;j<ns;j++)s+=ctrl->current_K[i*ns+j]*ctrl->observer_state[j]; u-=s;}
    }else if(y){
        for(int i=0;i<ni;i++){double s=0;
            for(int j=0;j<ns;j++)s+=ctrl->current_K[i*ns+j]*(y[j]-x_eq[j]); u-=s;}
    }
    if(gs->use_hysteresis)u=bumpless_transfer_filter(u,ctrl->prev_u,dt,0.1);
    ctrl->prev_u=u; if(u_eq)u+=u_eq[0]; return u;
}
double bumpless_transfer_filter(double un,double uo,double dt,double tau){
    if(dt<=0||tau<=0)return un; double a=1.0-exp(-dt/tau); return a*un+(1.0-a)*uo;
}
int hysteresis_scheduling(const double *sv,const double *pv,int ce,const GainSchedule *gs,double b){
    int d=gs->n_scheduling_vars; double dist=0;
    for(int j=0;j<d;j++){double df=sv[j]-pv[j];dist+=df*df;}
    if(sqrt(dist)<b)return ce; double bd=INFINITY;int best=ce;
    for(int i=0;i<gs->n_entries;i++){dist=0;
        for(int j=0;j<d;j++){double df=sv[j]-gs->entries[i].scheduling_value[j];dist+=df*df;}
        if(dist<bd){bd=dist;best=i;}}
    return best;
}
static void solve_are_iter(const double *A,const double *B,const double *Q,const double *R,
    int n,int m,double *P,double *K){
    for(int i=0;i<n*n;i++)P[i]=Q[i];
    double *Pn=malloc((size_t)(n*n)*sizeof(double)),*BK=malloc((size_t)(n*n)*sizeof(double));
    double Rinv=(R[0]>1e-15)?1.0/R[0]:1.0;
    for(int iter=0;iter<30;iter++){
        for(int i=0;i<m;i++)for(int j=0;j<n;j++){K[i*n+j]=0;
            for(int k=0;k<n;k++)K[i*n+j]+=Rinv*B[k*m+i]*P[k*n+j];}
        for(int i=0;i<n;i++)for(int j=0;j<n;j++){BK[i*n+j]=0;
            for(int k=0;k<m;k++)BK[i*n+j]+=B[i*m+k]*K[k*n+j];}
        double *Acl=malloc((size_t)(n*n)*sizeof(double));
        for(int i=0;i<n*n;i++)Acl[i]=A[i]-BK[i];
        double st=0.1;
        for(int i=0;i<n;i++)for(int j=0;j<n;j++){
            double res=0;for(int k=0;k<n;k++)res+=Acl[k*n+i]*P[k*n+j]+P[i*n+k]*Acl[k*n+j];
            res+=Q[i*n+j]+K[0]*R[0]*K[0]; Pn[i*n+j]=P[i*n+j]-st*res;
        }
        double diff=0;for(int i=0;i<n*n;i++){diff+=(Pn[i]-P[i])*(Pn[i]-P[i]);P[i]=Pn[i];}
        if(sqrt(diff)<1e-10)break; free(Acl);
    }
    free(Pn);free(BK);
}
void design_lqr_schedule(GainSchedule *gs,const double *Q,const double *R){
    int ns=gs->n_states,ni=gs->n_inputs;
    for(int i=0;i<gs->n_entries;i++){GainScheduleEntry *e=&gs->entries[i];
        if(!e->sys||!e->sys->A||!e->sys->B)continue;
        double *P=malloc((size_t)(ns*ns)*sizeof(double));
        double *K=malloc((size_t)(ni*ns)*sizeof(double));
        solve_are_iter(e->sys->A,e->sys->B,Q,R,ns,ni,P,K);
        if(e->K)free(e->K);e->K=K;free(P);}
}
void gain_schedule_to_lpv(const GainSchedule *gs,const double *sv,
    double *Ao,double *Bo,double *Co,double *Do){
    int ns=gs->n_states,ni=gs->n_inputs,no=gs->n_outputs,d=gs->n_scheduling_vars;
    memset(Ao,0,(size_t)(ns*ns)*sizeof(double));memset(Bo,0,(size_t)(ns*ni)*sizeof(double));
    if(Co)memset(Co,0,(size_t)(no*ns)*sizeof(double));
    if(Do)memset(Do,0,(size_t)(no*ni)*sizeof(double)); double ws=0;
    for(int i=0;i<gs->n_entries;i++){double w=1.0;
        for(int j=0;j<d;j++){double dv=fabs(sv[j]-gs->entries[i].scheduling_value[j]);
            w*=1.0/(1.0+dv*dv);}
        if(gs->entries[i].sys&&gs->entries[i].sys->A)
            for(int r=0;r<ns;r++)for(int c=0;c<ns;c++)Ao[r*ns+c]+=w*gs->entries[i].sys->A[r*ns+c];
        if(gs->entries[i].sys&&gs->entries[i].sys->B&&ni>0)
            for(int r=0;r<ns;r++)for(int c=0;c<ni;c++)Bo[r*ni+c]+=w*gs->entries[i].sys->B[r*ni+c];
        ws+=w;
    }
    if(ws>1e-15){for(int i=0;i<ns*ns;i++)Ao[i]/=ws;if(ni>0)for(int i=0;i<ns*ni;i++)Bo[i]/=ws;}
}
void design_pole_placement_schedule(GainSchedule *gs,const double *poles){
    int n=gs->n_states;
    for(int i=0;i<gs->n_entries;i++){GainScheduleEntry *e=&gs->entries[i];
        if(!e->sys||!e->sys->A||!e->sys->B||gs->n_inputs!=1)continue;
        double *Cmat=malloc((size_t)(n*n)*sizeof(double));
        for(int j=0;j<n;j++)Cmat[j]=e->sys->B[j];
        double *prev=malloc((size_t)n*sizeof(double));memcpy(prev,e->sys->B,(size_t)n*sizeof(double));
        for(int k=1;k<n;k++){double *AkB=malloc((size_t)n*sizeof(double));
            for(int r=0;r<n;r++){AkB[r]=0;for(int q=0;q<n;q++)AkB[r]+=e->sys->A[r*n+q]*prev[q];}
            for(int r=0;r<n;r++)Cmat[k*n+r]=AkB[r];
            memcpy(prev,AkB,(size_t)n*sizeof(double));free(AkB);
        }
        double *poly=malloc((size_t)(n+1)*sizeof(double));poly[0]=1;
        for(int k=0;k<n;k++){double p=-poles[i*n+k];
            for(int j=n;j>0;j--)poly[j]=poly[j]-p*poly[j-1];}
        for(int j=0;j<n;j++)e->K[j]=-poly[j+1];
        free(Cmat);free(prev);free(poly);
    }
}
void design_robust_gain_schedule(GainSchedule *gs,double alpha){
    double Qd[16]={1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1},Rd=1.0/alpha;
    design_lqr_schedule(gs,Qd,&Rd);
}
void free_gain_schedule(GainSchedule *gs){
    if(!gs)return;for(int i=0;i<gs->n_entries;i++){free(gs->entries[i].K);free(gs->entries[i].L);}
    free(gs->entries);free(gs);
}
void free_gain_scheduled_controller(GainScheduledController *c){
    if(c){free(c->current_K);free(c->current_L);free(c->observer_state);free(c);}
}
