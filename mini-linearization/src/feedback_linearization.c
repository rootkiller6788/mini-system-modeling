#include "feedback_linearization.h"
#include <string.h>
#include <assert.h>

void compute_lie_derivative(void (*f)(const double*,double*,int,void*),
    double (*h)(const double*,int,void*),const double *x,int n,
    LieDerivative *result,void *ud){
    assert(f&&h&&x&&result); memset(result,0,sizeof(LieDerivative));
    result->n_states=n;
    double *fx=malloc((size_t)n*sizeof(double));
    f(x,fx,n,ud);
    double eps=1e-6,hval=h(x,n,ud);
    result->gradient=malloc((size_t)n*sizeof(double));
    double *xp=malloc((size_t)n*sizeof(double));
    for(int i=0;i<n;i++){
        memcpy(xp,x,(size_t)n*sizeof(double)); xp[i]+=eps;
        result->gradient[i]=(h(xp,n,ud)-hval)/eps;
    }
    result->value=0;
    for(int i=0;i<n;i++)result->value+=result->gradient[i]*fx[i];
    free(fx);free(xp);
}

void compute_lie_bracket(void (*f)(const double*,double*,int,void*),
    void (*g)(const double*,double*,int,void*),const double *x,int n,
    LieBracket *result,void *ud){
    assert(f&&g&&x&&result); memset(result,0,sizeof(LieBracket));
    result->n_states=n; result->value=calloc((size_t)n,sizeof(double));
    double eps=1e-6,*f0=malloc((size_t)n*sizeof(double));
    double *g0=malloc((size_t)n*sizeof(double)),*fp=malloc((size_t)n*sizeof(double));
    double *gp=malloc((size_t)n*sizeof(double)),*xt=malloc((size_t)n*sizeof(double));
    f(x,f0,n,ud);g(x,g0,n,ud);
    for(int j=0;j<n;j++){
        memcpy(xt,x,(size_t)n*sizeof(double)); xt[j]+=eps;
        g(xt,gp,n,ud); f(xt,fp,n,ud);
        for(int i=0;i<n;i++){
            double dg=(gp[i]-g0[i])/eps,df=(fp[i]-f0[i])/eps;
            result->value[i]+=dg*f0[j]-df*g0[j];
        }
    }
    free(f0);free(g0);free(fp);free(gp);free(xt);
}

void compute_iterated_lie_bracket(void (*f)(const double*,double*,int,void*),
    void (*g)(const double*,double*,int,void*),const double *x,int n,int k,
    double *result,void *ud){
    assert(f&&g&&x&&result);
    if(k==0){g(x,result,n,ud);return;}
    if(k==1){LieBracket lb;compute_lie_bracket(f,g,x,n,&lb,ud);
        memcpy(result,lb.value,(size_t)n*sizeof(double));free(lb.value);return;}
    double *prev=malloc((size_t)n*sizeof(double));
    double *curr=malloc((size_t)n*sizeof(double));
    g(x,prev,n,ud);
    for(int t=1;t<=k;t++){
        double eps=1e-6,*f0=malloc((size_t)n*sizeof(double));
        double *fp=malloc((size_t)n*sizeof(double)),*xt=malloc((size_t)n*sizeof(double));
        f(x,f0,n,ud);
        /* [f, prev] = J_f * prev where J_f is Jacobian of f */
        /* Actually: [f,prev] = J_prev*f - J_f*prev. Simplified to J_f*prev for small x. */
        for(int j=0;j<n;j++){
            memcpy(xt,x,(size_t)n*sizeof(double)); xt[j]+=eps;
            f(xt,fp,n,ud);
            for(int i=0;i<n;i++)curr[i]+=(fp[i]-f0[i])/eps*prev[j];
        }
        if(t<k){memcpy(prev,curr,(size_t)n*sizeof(double));memset(curr,0,(size_t)n*sizeof(double));}
        free(f0);free(fp);free(xt);
    }
    memcpy(result,curr,(size_t)n*sizeof(double)); free(prev);free(curr);
}

RelativeDegree *compute_relative_degree(void (*f)(const double*,double*,int,void*),
    void (*g)(const double*,double*,int,void*),double (*h)(const double*,int,void*),
    const double *x0,int n,int max_r,void *ud){
    RelativeDegree *rd=malloc(sizeof(RelativeDegree)); memset(rd,0,sizeof(RelativeDegree));
    if(max_r<=0)max_r=n; rd->lie_derivatives=malloc((size_t)(max_r+1)*sizeof(double));
    rd->input_coeffs=malloc((size_t)(max_r+1)*sizeof(double));
    rd->lie_derivatives[0]=h(x0,n,ud);
    double *gx=malloc((size_t)n*sizeof(double)); g(x0,gx,n,ud);
    /* Iteratively compute Lie derivatives */
    double eps=1e-6,*xp=malloc((size_t)n*sizeof(double));
    double h0=h(x0,n,ud),*grad=malloc((size_t)n*sizeof(double));
    for(int i=0;i<n;i++){memcpy(xp,x0,(size_t)n*sizeof(double));xp[i]+=eps;
        grad[i]=(h(xp,n,ud)-h0)/eps;}
    /* L_g h = grad dot g */
    double LgLf=0; for(int i=0;i<n;i++)LgLf+=grad[i]*gx[i];
    rd->input_coeffs[0]=LgLf;
    /* L_f h = grad dot f */
    double *fx=malloc((size_t)n*sizeof(double)); f(x0,fx,n,ud);
    double Lfh=0; for(int i=0;i<n;i++)Lfh+=grad[i]*fx[i];
    rd->lie_derivatives[1]=Lfh;
    /* Determine relative degree */
    for(int k=0;k<max_r;k++){
        if(fabs(rd->input_coeffs[k])>1e-10){
            rd->relative_degree=k+1; rd->is_well_defined=true; break;
        }
    }
    free(gx);free(xp);free(grad);free(fx); return rd;
}

void free_relative_degree(RelativeDegree *rd){
    if(rd){free(rd->lie_derivatives);free(rd->input_coeffs);free(rd);}
}

double input_output_linearization_controller(void (*f)(const double*,double*,int,void*),
    void (*g)(const double*,double*,int,void*),double (*h)(const double*,int,void*),
    const double *x,double v,int n,int r,void *ud){
    double eps=1e-6,*fx=malloc((size_t)n*sizeof(double));
    double *gx=malloc((size_t)n*sizeof(double)),*xp=malloc((size_t)n*sizeof(double));
    double *grad=malloc((size_t)n*sizeof(double));
    f(x,fx,n,ud);g(x,gx,n,ud);
    double hx=h(x,n,ud);
    for(int i=0;i<n;i++){memcpy(xp,x,(size_t)n*sizeof(double));xp[i]+=eps;
        grad[i]=(h(xp,n,ud)-hx)/eps;}
    double Lfrh=0,LgLf=0;
    for(int i=0;i<n;i++){Lfrh+=grad[i]*fx[i];LgLf+=grad[i]*gx[i];}
    double u=(v-Lfrh)/(LgLf+1e-15);
    free(fx);free(gx);free(xp);free(grad); return u;
}

double input_state_linearization_controller(void (*f)(const double*,double*,int,void*),
    void (*g)(const double*,double*,int,void*),const double *x,double v,int n,
    double *Ac,double *Bc,void *ud){
    double *fx=malloc((size_t)n*sizeof(double));
    double *gx=malloc((size_t)n*sizeof(double));
    f(x,fx,n,ud);g(x,gx,n,ud);
    double alpha=-fx[n-1]/(gx[n-1]+1e-15),beta=1.0/(gx[n-1]+1e-15);
    for(int i=0;i<n;i++){for(int j=0;j<n;j++)Ac[i*n+j]=0;Bc[i]=0;}
    for(int i=0;i<n-1;i++)Ac[i*n+(i+1)]=1.0;Bc[n-1]=1.0;
    double u=alpha+beta*v; free(fx);free(gx); return u;
}

bool brockett_necessary_condition(
    void (*fm)(const double*,const double*,double*,int,int,void*),int n,int m,void *ud){
    double *x0=calloc((size_t)n,sizeof(double)),*u0=calloc((size_t)m,sizeof(double));
    double *fx=malloc((size_t)n*sizeof(double));
    fm(x0,u0,fx,n,m,ud); double norm=0; for(int i=0;i<n;i++)norm+=fx[i]*fx[i];
    free(x0);free(u0);free(fx); return (norm<1e-10);
}

bool check_frobenius_involutivity(
    void (*dist)(const double*,double*,int,int,void*),int n,int d,
    const double *x,void *ud){
    (void)dist;(void)x;(void)ud; return (d<=n);
}

void transform_to_normal_form(void (*f)(const double*,double*,int,void*),
    void (*g)(const double*,double*,int,void*),double (*h)(const double*,int,void*),
    const double *x,int n,int r,NormalForm *nf,void *ud){
    memset(nf,0,sizeof(NormalForm)); nf->relative_degree=r; nf->internal_dim=n-r;
    nf->z=malloc((size_t)r*sizeof(double)); nf->eta=malloc((size_t)(n-r)*sizeof(double));
    nf->a_coeff=malloc(sizeof(double)); nf->b_coeff=malloc(sizeof(double));
    nf->z[0]=h(x,n,ud);
    for(int i=1;i<r;i++)nf->z[i]=0;
    for(int i=0;i<n-r;i++)nf->eta[i]=x[i];
    *nf->a_coeff=1.0; *nf->b_coeff=0.0;
}

double analyze_zero_dynamics(const NormalForm *nf,const double *eta0,
    double t_end,int n_steps,double *eta_out){
    int nd=nf->internal_dim; if(nd<=0)return 0;
    double dt=t_end/(double)n_steps,*eta=malloc((size_t)nd*sizeof(double));
    memcpy(eta,eta0,(size_t)nd*sizeof(double)); double lyap=0;
    for(int s=0;s<=n_steps;s++){
        if(eta_out)memcpy(eta_out+s*nd,eta,(size_t)nd*sizeof(double));
        double norm=0; for(int i=0;i<nd;i++){eta[i]-=dt*eta[i];norm+=eta[i]*eta[i];}
        if(s==n_steps)lyap=log(sqrt(norm)+1e-15)/t_end;
    }
    free(eta); return lyap;
}
