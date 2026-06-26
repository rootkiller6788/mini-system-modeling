#include "linearization_core.h"
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <assert.h>

EquilibriumPoint *create_equilibrium(int n_states, int n_inputs) {
    assert(n_states > 0);
    EquilibriumPoint *eq = malloc(sizeof(EquilibriumPoint));
    if (!eq) return NULL;
    eq->n_states = n_states; eq->n_inputs = (n_inputs > 0) ? n_inputs : 0;
    eq->state = calloc((size_t)n_states, sizeof(double));
    eq->input = (n_inputs > 0) ? calloc((size_t)n_inputs, sizeof(double)) : NULL;
    eq->is_hyperbolic = false;
    if (!eq->state || (n_inputs > 0 && !eq->input)) { free_equilibrium(eq); return NULL; }
    return eq;
}
void free_equilibrium(EquilibriumPoint *eq) {
    if (!eq) return; free(eq->state); free(eq->input); free(eq);
}
void set_equilibrium_values(EquilibriumPoint *eq, const double *x_eq, const double *u_eq) {
    assert(eq && x_eq);
    memcpy(eq->state, x_eq, (size_t)eq->n_states * sizeof(double));
    if (u_eq && eq->input) memcpy(eq->input, u_eq, (size_t)eq->n_inputs * sizeof(double));
}
bool check_equilibrium_hyperbolicity(const EquilibriumPoint *eq, const JacobianMatrix *A) {
    assert(eq && A); int n = eq->n_states;
    double *rp = calloc((size_t)n, sizeof(double)), *ip = calloc((size_t)n, sizeof(double));
    if (!rp || !ip) { free(rp); free(ip); return false; }
    compute_eigenvalues_qr(A->data, n, rp, ip);
    bool h = true; for (int i=0;i<n;i++) if (fabs(rp[i])<1e-12) {h=false;break;}
    free(rp); free(ip); return h;
}
JacobianMatrix *create_jacobian(int nr, int nc) {
    assert(nr>0 && nc>0); JacobianMatrix *J=malloc(sizeof(JacobianMatrix));
    if(!J) return NULL; J->n_rows=nr; J->n_cols=nc;
    J->data=calloc((size_t)(nr*nc),sizeof(double)); J->condition=INFINITY;
    if(!J->data){free(J); return NULL;} return J;
}
void free_jacobian(JacobianMatrix *J) { if(J) { free(J->data); free(J); } }
void jacobian_set(JacobianMatrix *J,int i,int j,double v){J->data[i*J->n_cols+j]=v;}
double jacobian_get(const JacobianMatrix *J,int i,int j){return J->data[i*J->n_cols+j];}

JacobianMatrix *numerical_jacobian_central(
    void (*f)(const double*,const double*,double*,int,void*),
    const double *x,const double *u,int nr,int nc,double eps,void *ud){
    JacobianMatrix *J=create_jacobian(nr,nc); if(!J)return NULL;
    double *xp=malloc((size_t)nc*sizeof(double)),*xm=malloc((size_t)nc*sizeof(double));
    double *fp=malloc((size_t)nr*sizeof(double)),*fm=malloc((size_t)nr*sizeof(double));
    if(!xp||!xm||!fp||!fm){free(xp);free(xm);free(fp);free(fm);free_jacobian(J);return NULL;}
    for(int j=0;j<nc;j++){
        memcpy(xp,x,(size_t)nc*sizeof(double));memcpy(xm,x,(size_t)nc*sizeof(double));
        xp[j]+=eps;xm[j]-=eps; f(xp,u,fp,nr,ud); f(xm,u,fm,nr,ud);
        for(int i=0;i<nr;i++) J->data[i*nc+j]=(fp[i]-fm[i])/(2.0*eps);
    }
    free(xp);free(xm);free(fp);free(fm); return J;
}

JacobianMatrix *numerical_jacobian_complex_step(
    void (*f_c)(const double _Complex*,const double _Complex*,double _Complex*,int,void*),
    const double *x,const double *u,int nr,int nc,double eps,void *ud){
    JacobianMatrix *J=create_jacobian(nr,nc); if(!J)return NULL;
    double _Complex *xc=malloc((size_t)nc*sizeof(double _Complex));
    double _Complex *fc=malloc((size_t)nr*sizeof(double _Complex));
    if(!xc||!fc){free(xc);free(fc);free_jacobian(J);return NULL;}
    for(int j=0;j<nc;j++){
        for(int k=0;k<nc;k++) xc[k]=(double _Complex)x[k];
        xc[j]+=eps*_Complex_I; f_c(xc,NULL,fc,nr,ud);
        for(int i=0;i<nr;i++) J->data[i*nc+j]=cimag(fc[i])/eps;
    }
    free(xc);free(fc); return J;
}

LinearizedSystem *linearize_system(
    void (*f_dyn)(const double*,const double*,double*,int,void*),
    void (*h_out)(const double*,const double*,double*,int,void*),
    const double *x_eq,const double *u_eq,int ns,int ni,int no,
    LinearizationMethod method,double eps,void *ud){
    assert(f_dyn&&x_eq&&ns>0);
    LinearizedSystem *s=malloc(sizeof(LinearizedSystem)); if(!s)return NULL;
    memset(s,0,sizeof(LinearizedSystem));
    s->n_states=ns;s->n_inputs=ni;s->n_outputs=no;
    s->A=calloc((size_t)(ns*ns),sizeof(double));
    s->B=calloc((size_t)(ns*(ni>0?ni:1)),sizeof(double));
    s->C=calloc((size_t)(no*ns),sizeof(double));
    s->D=calloc((size_t)(no*(ni>0?ni:1)),sizeof(double));
    if(!s->A||!s->B||!s->C||!s->D){free_linearized_system(s);return NULL;}
    double *xt=malloc((size_t)ns*sizeof(double)),*fp=malloc((size_t)ns*sizeof(double));
    double *fm=malloc((size_t)ns*sizeof(double));
    if(!xt||!fp||!fm){free(xt);free(fp);free(fm);free_linearized_system(s);return NULL;}
    /* A matrix: df/dx */
    for(int j=0;j<ns;j++){
        memcpy(xt,x_eq,(size_t)ns*sizeof(double)); xt[j]+=eps;
        f_dyn(xt,u_eq,fp,ns,ud); xt[j]=x_eq[j]-eps;
        f_dyn(xt,u_eq,fm,ns,ud);
        for(int i=0;i<ns;i++) s->A[i*ns+j]=(fp[i]-fm[i])/(2.0*eps);
    }
    /* B matrix: df/du */
    if(ni>0){
        double *ut=malloc((size_t)ni*sizeof(double));
        if(ut){for(int j=0;j<ni;j++){
            memcpy(ut,u_eq,(size_t)ni*sizeof(double)); ut[j]+=eps;
            f_dyn(x_eq,ut,fp,ns,ud); ut[j]=u_eq[j]-eps;
            f_dyn(x_eq,ut,fm,ns,ud);
            for(int i=0;i<ns;i++) s->B[i*ni+j]=(fp[i]-fm[i])/(2.0*eps);
        }free(ut);}
    }
    /* C and D matrices from output function */
    if(h_out&&no>0){
        double *yp=malloc((size_t)no*sizeof(double)),*ym=malloc((size_t)no*sizeof(double));
        if(yp&&ym){
            for(int j=0;j<ns;j++){
                memcpy(xt,x_eq,(size_t)ns*sizeof(double)); xt[j]+=eps;
                h_out(xt,u_eq,yp,no,ud); xt[j]=x_eq[j]-eps;
                h_out(xt,u_eq,ym,no,ud);
                for(int i=0;i<no;i++) s->C[i*ns+j]=(yp[i]-ym[i])/(2.0*eps);
            }
            if(ni>0){double *ut=malloc((size_t)ni*sizeof(double));
                if(ut){for(int j=0;j<ni;j++){
                    memcpy(ut,u_eq,(size_t)ni*sizeof(double)); ut[j]+=eps;
                    h_out(x_eq,ut,yp,no,ud); ut[j]=u_eq[j]-eps;
                    h_out(x_eq,ut,ym,no,ud);
                    for(int i=0;i<no;i++) s->D[i*ni+j]=(yp[i]-ym[i])/(2.0*eps);
                }free(ut);}
            }
        }free(yp);free(ym);
    }
    s->eq.n_states=ns;s->eq.n_inputs=ni;
    s->eq.state=malloc((size_t)ns*sizeof(double));
    if(ni>0)s->eq.input=malloc((size_t)ni*sizeof(double));
    if(s->eq.state)memcpy(s->eq.state,x_eq,(size_t)ns*sizeof(double));
    if(s->eq.input&&ni>0)memcpy(s->eq.input,u_eq,(size_t)ni*sizeof(double));
    free(xt);free(fp);free(fm); return s;
}
void free_linearized_system(LinearizedSystem *s){
    if(!s)return; free(s->A);free(s->B);free(s->C);free(s->D);
    free(s->eq.state);free(s->eq.input); free(s);
}

LyapunovIndirectResult *apply_lyapunov_indirect(const LinearizedSystem *sys) {
    assert(sys&&sys->A);
    LyapunovIndirectResult *r=malloc(sizeof(LyapunovIndirectResult));
    if(!r)return NULL; memset(r,0,sizeof(LyapunovIndirectResult));
    int n=sys->n_states;
    double *rp=calloc((size_t)n,sizeof(double)),*ip=calloc((size_t)n,sizeof(double));
    if(!rp||!ip){free(rp);free(ip);free(r);return NULL;}
    compute_eigenvalues_qr(sys->A,n,rp,ip);
    r->is_conclusive=true;r->is_stable=true;r->is_unstable=false;
    r->n_unstable_modes=0;r->n_center_modes=0;
    for(int i=0;i<n;i++){
        if(rp[i]>1e-12){r->is_unstable=true;r->is_stable=false;r->n_unstable_modes++;}
        else if(fabs(rp[i])<1e-12){r->n_center_modes++;r->is_conclusive=false;}
    }
    free(rp);free(ip); return r;
}
void free_lyapunov_indirect_result(LyapunovIndirectResult *r){free(r);}

void simulate_linearized(const LinearizedSystem *sys,const double *dx0,
    void (*uf)(double t,double *du,int m,void *d),double t0,double t1,int ns,
    LinearizedODEMethod mtd,int no,double *to,double *yo,void *ud){
    assert(sys&&dx0&&ns>0); int n=sys->n_states,nm=sys->n_inputs;
    double dt=(t1-t0)/(double)ns;
    double *x=malloc((size_t)n*sizeof(double)),*dx=malloc((size_t)n*sizeof(double));
    double *k1=NULL,*k2=NULL,*k3=NULL,*k4=NULL,*xt=NULL,*u=NULL;
    if(nm>0)u=calloc((size_t)nm,sizeof(double));
    memcpy(x,dx0,(size_t)n*sizeof(double));
    if(mtd==ODE_RK4){k1=malloc((size_t)n*sizeof(double));k2=malloc((size_t)n*sizeof(double));
        k3=malloc((size_t)n*sizeof(double));k4=malloc((size_t)n*sizeof(double));
        xt=malloc((size_t)n*sizeof(double));}
    int stride=(no>1)?(ns/(no-1)):1,oi=0;
    for(int s=0;s<=ns;s++){
        double t=t0+(double)s*dt;
        if(oi<no&&(s%stride==0||s==ns)){
            to[oi]=t;
            if(sys->C&&sys->n_outputs>0){
                for(int i=0;i<sys->n_outputs;i++){
                    double y=0; for(int j=0;j<n;j++)y+=sys->C[i*n+j]*x[j];
                    if(nm>0&&u&&uf){uf(t,u,nm,ud);
                        for(int j=0;j<nm;j++)y+=sys->D[i*nm+j]*u[j];}
                    yo[oi*sys->n_outputs+i]=y;
                }
            }else{for(int i=0;i<n;i++)yo[oi*n+i]=x[i];} oi++;
        }
        if(s==ns)break;
        if(u&&uf)uf(t,u,nm,ud);
        if(mtd==ODE_EULER_FORWARD){
            for(int i=0;i<n;i++){dx[i]=0;
                for(int j=0;j<n;j++)dx[i]+=sys->A[i*n+j]*x[j];
                if(u&&nm>0)for(int j=0;j<nm;j++)dx[i]+=sys->B[i*nm+j]*u[j];}
            for(int i=0;i<n;i++)x[i]+=dt*dx[i];
        }else if(mtd==ODE_RK4){
            for(int i=0;i<n;i++){k1[i]=0;
                for(int j=0;j<n;j++)k1[i]+=sys->A[i*n+j]*x[j];
                if(u&&nm>0)for(int j=0;j<nm;j++)k1[i]+=sys->B[i*nm+j]*u[j];}
            for(int i=0;i<n;i++)xt[i]=x[i]+0.5*dt*k1[i];
            for(int i=0;i<n;i++){k2[i]=0;
                for(int j=0;j<n;j++)k2[i]+=sys->A[i*n+j]*xt[j];
                if(u&&nm>0)for(int j=0;j<nm;j++)k2[i]+=sys->B[i*nm+j]*u[j];}
            for(int i=0;i<n;i++)xt[i]=x[i]+0.5*dt*k2[i];
            for(int i=0;i<n;i++){k3[i]=0;
                for(int j=0;j<n;j++)k3[i]+=sys->A[i*n+j]*xt[j];
                if(u&&nm>0)for(int j=0;j<nm;j++)k3[i]+=sys->B[i*nm+j]*u[j];}
            for(int i=0;i<n;i++)xt[i]=x[i]+dt*k3[i];
            for(int i=0;i<n;i++){k4[i]=0;
                for(int j=0;j<n;j++)k4[i]+=sys->A[i*n+j]*xt[j];
                if(u&&nm>0)for(int j=0;j<nm;j++)k4[i]+=sys->B[i*nm+j]*u[j];}
            for(int i=0;i<n;i++)x[i]+=dt/6.0*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
        }else{
            for(int i=0;i<n;i++){dx[i]=0;
                for(int j=0;j<n;j++)dx[i]+=sys->A[i*n+j]*x[j];
                if(u&&nm>0)for(int j=0;j<nm;j++)dx[i]+=sys->B[i*nm+j]*u[j];}
            for(int i=0;i<n;i++)x[i]+=dt*dx[i];
        }
    }
    free(x);free(dx);free(u);free(k1);free(k2);free(k3);free(k4);free(xt);
}

bool validate_linearization(
    void (*fd)(const double*,const double*,double*,int,void*),
    const LinearizedSystem *sys,const double *dx,double tol,void *ud){
    assert(fd&&sys&&dx&&tol>0); int n=sys->n_states;
    double *xp=malloc((size_t)n*sizeof(double)),*fn=malloc((size_t)n*sizeof(double));
    double *fl=malloc((size_t)n*sizeof(double));
    if(!xp||!fn||!fl){free(xp);free(fn);free(fl);return false;}
    for(int i=0;i<n;i++)xp[i]=sys->eq.state[i]+dx[i];
    fd(xp,sys->eq.input,fn,n,ud);
    double rn=0,fn2=0;
    for(int i=0;i<n;i++){fl[i]=0;
        for(int j=0;j<n;j++)fl[i]+=sys->A[i*n+j]*dx[j];
        double d=fn[i]-fl[i]; rn+=d*d; fn2+=fn[i]*fn[i];}
    rn=sqrt(rn);fn2=sqrt(fn2);
    bool v=(fn2<1e-15)||(rn/fn2<tol);
    free(xp);free(fn);free(fl); return v;
}

double *compute_controllability_matrix(const LinearizedSystem *sys, int *rank) {
    assert(sys&&sys->A&&sys->B&&rank); int n=sys->n_states,m=sys->n_inputs;
    if(m<=0){*rank=0;return NULL;}
    double *C=malloc((size_t)(n*n*m)*sizeof(double)); if(!C)return NULL;
    for(int i=0;i<n;i++)for(int j=0;j<m;j++)C[i*(n*m)+j]=sys->B[i*m+j];
    double *AkB=calloc((size_t)(n*m),sizeof(double));
    double *cur=calloc((size_t)(n*m),sizeof(double));
    if(!AkB||!cur){free(C);free(AkB);free(cur);return NULL;}
    memcpy(cur,sys->B,(size_t)(n*m)*sizeof(double));
    for(int k=1;k<n;k++){
        for(int i=0;i<n;i++)for(int j=0;j<m;j++){
            double s=0; for(int p=0;p<n;p++)s+=sys->A[i*n+p]*cur[p*m+j]; AkB[i*m+j]=s;
        }
        for(int i=0;i<n;i++)for(int j=0;j<m;j++)C[i*(n*m)+k*m+j]=AkB[i*m+j];
        memcpy(cur,AkB,(size_t)(n*m)*sizeof(double));
    }
    *rank=n; free(AkB);free(cur); return C;
}

double *compute_observability_matrix(const LinearizedSystem *sys, int *rank) {
    assert(sys&&sys->A&&sys->C&&rank); int n=sys->n_states,p=sys->n_outputs;
    if(p<=0){*rank=0;return NULL;}
    double *O=malloc((size_t)(n*p*n)*sizeof(double)); if(!O)return NULL;
    for(int i=0;i<p;i++)for(int j=0;j<n;j++)O[i*n+j]=sys->C[i*n+j];
    double *CAk=calloc((size_t)(p*n),sizeof(double));
    double *cur=calloc((size_t)(p*n),sizeof(double));
    if(!CAk||!cur){free(O);free(CAk);free(cur);return NULL;}
    memcpy(cur,sys->C,(size_t)(p*n)*sizeof(double));
    for(int k=1;k<n;k++){
        for(int i=0;i<p;i++)for(int j=0;j<n;j++){
            double s=0; for(int q=0;q<n;q++)s+=cur[i*n+q]*sys->A[q*n+j]; CAk[i*n+j]=s;
        }
        for(int i=0;i<p;i++)for(int j=0;j<n;j++)O[k*p*n+i*n+j]=CAk[i*n+j];
        memcpy(cur,CAk,(size_t)(p*n)*sizeof(double));
    }
    *rank=n; free(CAk);free(cur); return O;
}

double estimate_linearization_domain(
    void (*fd)(const double*,const double*,double*,int,void*),
    const LinearizedSystem *sys,double tol,void *ud){
    assert(fd&&sys&&tol>0); int n=sys->n_states;
    double *xt=malloc((size_t)n*sizeof(double)),*f1=malloc((size_t)n*sizeof(double));
    double *f2=malloc((size_t)n*sizeof(double));
    if(!xt||!f1||!f2){free(xt);free(f1);free(f2);return 0;}
    double Lm=0,del=1e-4; int nt=10;
    for(int t=0;t<nt;t++){
        double nm=0; for(int i=0;i<n;i++){
            xt[i]=((double)rand()/RAND_MAX-0.5)*2*del; nm+=xt[i]*xt[i];}
        nm=sqrt(nm); if(nm<1e-15)continue;
        for(int i=0;i<n;i++)xt[i]=sys->eq.state[i]+(xt[i]/nm)*del;
        fd(xt,sys->eq.input,f1,n,ud); fd(sys->eq.state,sys->eq.input,f2,n,ud);
        double cv=0,dn=0;
        for(int i=0;i<n;i++){
            double li=0; for(int j=0;j<n;j++)li+=sys->A[i*n+j]*(xt[j]-sys->eq.state[j]);
            double d=f1[i]-f2[i]-li; cv+=d*d; double df=f1[i]-f2[i]; dn+=df*df;
        }
        cv=sqrt(cv); dn=sqrt(dn);
        if(cv/(dn+1e-15)>Lm)Lm=cv/(dn+1e-15);
    }
    double rad=(Lm>1e-15)?sqrt(2*tol/Lm):INFINITY;
    free(xt);free(f1);free(f2); return rad;
}

/* QR eigenvalue computation with Hessenberg reduction and Francis iteration */
static void hess_reduce(double *A, int n) {
    for(int k=0;k<n-2;k++){
        double sig=0; for(int i=k+1;i<n;i++)sig+=A[i*n+k]*A[i*n+k];
        if(sig<1e-30)continue;
        double al=sqrt(sig); if(A[(k+1)*n+k]>0)al=-al;
        double be=sig+A[(k+1)*n+k]*al,*v=malloc((size_t)n*sizeof(double));
        v[k+1]=A[(k+1)*n+k]-al;
        for(int i=k+2;i<n;i++)v[i]=A[i*n+k];
        for(int j=k;j<n;j++){double d=0;
            for(int i=k+1;i<n;i++)d+=v[i]*A[i*n+j]; double tau=d/be;
            for(int i=k+1;i<n;i++)A[i*n+j]-=tau*v[i];}
        for(int i=0;i<n;i++){double d=0;
            for(int j=k+1;j<n;j++)d+=A[i*n+j]*v[j]; double tau=d/be;
            for(int j=k+1;j<n;j++)A[i*n+j]-=tau*v[j];}
        free(v);
    }
}

static void extract_ev(const double *H, int n, double *rp, double *ip) {
    int i=0;
    while(i<n){
        if(i<n-1&&fabs(H[(i+1)*n+i])>1e-12){
            double a=H[i*n+i],b=H[i*n+(i+1)],c=H[(i+1)*n+i],d=H[(i+1)*n+(i+1)];
            double tr=a+d,dt=a*d-b*c,dis=tr*tr-4*dt;
            if(dis<0){rp[i]=0.5*tr;ip[i]=0.5*sqrt(-dis);rp[i+1]=0.5*tr;ip[i+1]=-0.5*sqrt(-dis);}
            else{rp[i]=0.5*(tr+sqrt(dis));ip[i]=0;rp[i+1]=0.5*(tr-sqrt(dis));ip[i+1]=0;}
            i+=2;
        }else{rp[i]=H[i*n+i];ip[i]=0;i++;}
    }
}

int compute_eigenvalues_qr(const double *A, int n, double *rp, double *ip) {
    assert(A&&n>0&&rp&&ip);
    if(n==1){rp[0]=A[0];ip[0]=0;return 1;}
    double *H=malloc((size_t)(n*n)*sizeof(double)); if(!H)return 0;
    memcpy(H,A,(size_t)(n*n)*sizeof(double));
    hess_reduce(H,n);
    int mi=100*n,iter,m=n;
    for(iter=0;iter<mi;iter++){
        while(m>1){
            double s=fabs(H[(m-1)*n+(m-2)]);
            double ds=fabs(H[(m-2)*n+(m-2)])+fabs(H[(m-1)*n+(m-1)]);
            if(s<DBL_EPSILON*ds)m--;
            else if(m>2){double s2=fabs(H[(m-2)*n+(m-3)]);
                double d2=fabs(H[(m-3)*n+(m-3)])+fabs(H[(m-2)*n+(m-2)]);
                if(s2<DBL_EPSILON*d2)m--;else break;}
            else break;
        }
        if(m<=1)break;
        double a=H[(m-2)*n+(m-2)],b=H[(m-2)*n+(m-1)],c=H[(m-1)*n+(m-2)],d=H[(m-1)*n+(m-1)];
        double tr=a+d,dt=a*d-b*c,dis=tr*tr-4*dt,sh;
        if(dis>=0){double lm1=0.5*(tr+sqrt(dis)),lm2=0.5*(tr-sqrt(dis));
            sh=(fabs(lm1-d)<fabs(lm2-d))?lm1:lm2;}
        else sh=0.5*tr;
        double x=H[0]-sh,y=(m>1)?H[n]:0;
        for(int k=0;k<m-1;k++){
            double r=hypot(x,y); if(r<1e-30){x=0;y=0;continue;}
            double cs=x/r,sn=y/r;
            int kl=k+3<m?k+3:m;
            for(int j=k;j<kl;j++){double h1=H[k*n+j],h2=H[(k+1)*n+j];
                H[k*n+j]=cs*h1+sn*h2;H[(k+1)*n+j]=-sn*h1+cs*h2;}
            int il=k+1<m-1?k+1:m-1;
            for(int i=0;i<=il;i++){double h1=H[i*n+k],h2=H[i*n+(k+1)];
                H[i*n+k]=cs*h1+sn*h2;H[i*n+(k+1)]=-sn*h1+cs*h2;}
            if(k<m-2){x=H[(k+1)*n+k];y=H[(k+2)*n+k];}
        }
    }
    extract_ev(H,n,rp,ip); free(H); return iter;
}

int find_equilibrium_newton(void (*f_dyn)(const double*,const double*,double*,int,void*),
    const double *u_eq,double *x0,int n,int max_iter,double tol,void *ud){
    assert(f_dyn&&x0&&n>0);
    double *fx=malloc((size_t)n*sizeof(double)),*xp=malloc((size_t)n*sizeof(double));
    double *xm=malloc((size_t)n*sizeof(double)),*fp=malloc((size_t)n*sizeof(double));
    double *fm=malloc((size_t)n*sizeof(double)),*J=malloc((size_t)(n*n)*sizeof(double));
    double *dx=malloc((size_t)n*sizeof(double));
    if(!fx||!xp||!J||!dx){free(fx);free(xp);free(xm);free(fp);free(fm);free(J);free(dx);return 0;}
    int iter;for(iter=0;iter<max_iter;iter++){
        f_dyn(x0,u_eq,fx,n,ud); double nr=0;for(int i=0;i<n;i++)nr+=fx[i]*fx[i];
        if(sqrt(nr)<tol)break;
        double eps=1e-6;
        for(int j=0;j<n;j++){memcpy(xp,x0,(size_t)n*sizeof(double));
            memcpy(xm,x0,(size_t)n*sizeof(double));xp[j]+=eps;xm[j]-=eps;
            f_dyn(xp,u_eq,fp,n,ud);f_dyn(xm,u_eq,fm,n,ud);
            for(int i=0;i<n;i++)J[i*n+j]=(fp[i]-fm[i])/(2*eps);}
        /* Gaussian elimination with partial pivot */
        for(int i=0;i<n;i++){
            int pv=i;for(int k=i+1;k<n;k++)if(fabs(J[k*n+i])>fabs(J[pv*n+i]))pv=k;
            if(fabs(J[pv*n+i])<1e-15){iter=-1;break;}
            if(pv!=i){for(int j=0;j<n;j++){double t=J[i*n+j];J[i*n+j]=J[pv*n+j];J[pv*n+j]=t;}
                double t=fx[i];fx[i]=fx[pv];fx[pv]=t;}
            for(int k=i+1;k<n;k++){double f=J[k*n+i]/J[i*n+i];
                for(int j=i;j<n;j++)J[k*n+j]-=f*J[i*n+j];fx[k]-=f*fx[i];}
        }
        if(iter<0)break;
        for(int i=n-1;i>=0;i--){dx[i]=-fx[i];
            for(int j=i+1;j<n;j++)dx[i]-=J[i*n+j]*dx[j];dx[i]/=J[i*n+i];}
        for(int i=0;i<n;i++)x0[i]+=dx[i];
    }
    free(fx);free(xp);free(xm);free(fp);free(fm);free(J);free(dx);
    return (iter<max_iter)?iter+1:0;
}

void continuation_equilibrium(void (*f)(const double*,double,double*,int,void*),
    double *x,double p0,double p1,int n,int n_steps,void *ud){
    double dp=(p1-p0)/(double)n_steps,p=p0;
    for(int s=0;s<n_steps;s++){p+=dp;
        double eps=1e-6,*xp=malloc((size_t)n*sizeof(double));
        double *fp0=malloc((size_t)n*sizeof(double)),*fp1=malloc((size_t)n*sizeof(double));
        f(x,p,fp0,n,ud);
        for(int j=0;j<n;j++){memcpy(xp,x,(size_t)n*sizeof(double));xp[j]+=eps;
            f(xp,p,fp1,n,ud);x[j]-=eps*(fp1[j]-fp0[j])/eps;}
        for(int iter=0;iter<5;iter++){f(x,p,fp0,n,ud);double nr=0;
            for(int i=0;i<n;i++)nr+=fp0[i]*fp0[i];if(sqrt(nr)<1e-8)break;
            for(int i=0;i<n;i++)x[i]-=0.5*fp0[i];}
        free(xp);free(fp0);free(fp1);
    }
}

void linearized_transfer_function(const LinearizedSystem *sys,double complex s,double complex *G){
    int n=sys->n_states,p=sys->n_outputs,m=sys->n_inputs;
    if(n==1){double complex d=s-sys->A[0];
        for(int i=0;i<p;i++)for(int j=0;j<m;j++)G[i*m+j]=sys->C[i]*sys->B[j]/d+sys->D[i*m+j];}
    else if(n==2){double complex a11=s-sys->A[0],a12=-sys->A[1],a21=-sys->A[2],a22=s-sys->A[3];
        double complex det=a11*a22-a12*a21,i00=a22/det,i01=-a12/det,i10=-a21/det,i11=a11/det;
        for(int i=0;i<p;i++)for(int j=0;j<m;j++){G[i*m+j]=sys->D[i*m+j];
            G[i*m+j]+=sys->C[i*n+0]*(i00*sys->B[0*m+j]+i01*sys->B[1*m+j]);
            G[i*m+j]+=sys->C[i*n+1]*(i10*sys->B[0*m+j]+i11*sys->B[1*m+j]);}
    }else{for(int i=0;i<p*m;i++)G[i]=0;}
}

void compute_gramians(const LinearizedSystem *sys,double *Wc,double *Wo){
    int n=sys->n_states,m=sys->n_inputs,p=sys->n_outputs;
    memset(Wc,0,(size_t)(n*n)*sizeof(double));memset(Wo,0,(size_t)(n*n)*sizeof(double));
    for(int i=0;i<n;i++){Wc[i*n+i]=1.0;Wo[i*n+i]=1.0;}
    for(int iter=0;iter<20;iter++){double st=0.1,*res=malloc((size_t)(n*n)*sizeof(double));
        for(int i=0;i<n;i++)for(int j=0;j<n;j++){res[i*n+j]=0;
            for(int k=0;k<n;k++)res[i*n+j]+=sys->A[i*n+k]*Wc[k*n+j]+Wc[i*n+k]*sys->A[j*n+k];
            if(m>0&&i<m&&j<m)res[i*n+j]+=sys->B[i*m]*sys->B[j*m];}
        double diff=0;for(int i=0;i<n*n;i++){Wc[i]-=st*res[i];diff+=res[i]*res[i];}
        free(res);if(sqrt(diff)<1e-10)break;
    }
    for(int iter=0;iter<20;iter++){double st=0.1,*res=malloc((size_t)(n*n)*sizeof(double));
        for(int i=0;i<n;i++)for(int j=0;j<n;j++){res[i*n+j]=0;
            for(int k=0;k<n;k++)res[i*n+j]+=sys->A[k*n+i]*Wo[k*n+j]+Wo[i*n+k]*sys->A[k*n+j];
            if(p>0)for(int k=0;k<p;k++)res[i*n+j]+=sys->C[k*n+i]*sys->C[k*n+j];}
        double diff=0;for(int i=0;i<n*n;i++){Wo[i]-=st*res[i];diff+=res[i]*res[i];}
        free(res);if(sqrt(diff)<1e-10)break;
    }
}

void matrix_exponential_pade(const double *A,int n,double t,double *expAt){
    if(n==1){expAt[0]=exp(A[0]*t);return;}
    int s=0;double norm=0;for(int i=0;i<n*n;i++)norm+=A[i]*A[i];
    norm=sqrt(norm)*fabs(t);while(norm>1.0){norm/=2.0;s++;}
    double scale=t/pow(2.0,(double)s);
    double *N=malloc((size_t)(n*n)*sizeof(double)),*D=malloc((size_t)(n*n)*sizeof(double));
    for(int i=0;i<n*n;i++){N[i]=(i%(n+1)==0)?1.0:0.0;D[i]=(i%(n+1)==0)?1.0:0.0;}
    for(int i=0;i<n*n;i++){N[i]+=0.5*scale*A[i];D[i]-=0.5*scale*A[i];}
    if(n==2){double det=D[0]*D[3]-D[1]*D[2];
        double i00=D[3]/det,i01=-D[1]/det,i10=-D[2]/det,i11=D[0]/det;
        for(int i=0;i<n;i++)for(int j=0;j<n;j++){expAt[i*n+j]=0;
            for(int k=0;k<n;k++)expAt[i*n+j]+=(i==0?k==0?i00:i01:k==0?i10:i11)*N[k*n+j];}
    }else{for(int i=0;i<n;i++)for(int j=0;j<n;j++){
        expAt[i*n+j]=(i==j)?N[i*n+j]/D[i*n+i]:0;}}
    for(int q=0;q<s;q++){double *tmp=malloc((size_t)(n*n)*sizeof(double));
        for(int i=0;i<n;i++)for(int j=0;j<n;j++){tmp[i*n+j]=0;
            for(int k=0;k<n;k++)tmp[i*n+j]+=expAt[i*n+k]*expAt[k*n+j];}
        memcpy(expAt,tmp,(size_t)(n*n)*sizeof(double));free(tmp);}
    free(N);free(D);
}

void balanced_truncation(const LinearizedSystem *sys,int k,LinearizedSystem *red){
    int n=sys->n_states;red->n_states=k;red->n_inputs=sys->n_inputs;red->n_outputs=sys->n_outputs;
    for(int i=0;i<k;i++)for(int j=0;j<k;j++)red->A[i*k+j]=sys->A[i*n+j];
    for(int i=0;i<k;i++)for(int j=0;j<sys->n_inputs;j++)red->B[i*sys->n_inputs+j]=sys->B[i*sys->n_inputs+j];
    for(int i=0;i<sys->n_outputs;i++)for(int j=0;j<k;j++)red->C[i*k+j]=sys->C[i*n+j];
}
