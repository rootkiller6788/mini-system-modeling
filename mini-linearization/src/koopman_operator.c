#include "koopman_operator.h"
#include <string.h>
#include <assert.h>
#include <float.h>

static int svd_decomp(double *A,int m,int n,double *U,double *S,double *Vt){
    memcpy(U,A,(size_t)(m*n)*sizeof(double));
    for(int j=0;j<n;j++){double nr=0;for(int i=0;i<m;i++)nr+=U[i*n+j]*U[i*n+j];S[j]=sqrt(nr);}
    for(int j=0;j<n;j++){if(S[j]>1e-15)for(int i=0;i<m;i++)U[i*n+j]/=S[j];}
    for(int i=0;i<n&&i<m;i++)Vt[i*n+i]=1.0;
    return (m<n?m:n);
}

void compute_dmd(const double *X,const double *Y,int n,int m,int r,double dt,
    DMDResult *result){
    assert(X&&Y&&result&&n>0&&m>1);
    memset(result,0,sizeof(DMDResult));
    int ms=m-1; if(r<=0||r>ms)r=(ms<n)?ms:n;
    result->n_state_dim=n;result->n_snapshots=m;result->r=r;
    double *U=malloc((size_t)(n*ms)*sizeof(double));
    double *Sd=malloc((size_t)ms*sizeof(double));
    double *Vt=malloc((size_t)(ms*ms)*sizeof(double));
    double *Xc=malloc((size_t)(n*ms)*sizeof(double));
    memcpy(Xc,X,(size_t)(n*ms)*sizeof(double));
    svd_decomp(Xc,n,ms,U,Sd,Vt);
    double *At=malloc((size_t)(r*r)*sizeof(double));
    memset(At,0,(size_t)(r*r)*sizeof(double));
    double *YV=malloc((size_t)(n*r)*sizeof(double));
    for(int i=0;i<n;i++)for(int j=0;j<r;j++){
        double s=0;for(int k=0;k<ms;k++)s+=Y[i*ms+k]*Vt[j*ms+k];YV[i*r+j]=s;
    }
    for(int i=0;i<r;i++)for(int j=0;j<r;j++){
        double s=0;for(int k=0;k<n;k++)s+=U[k*ms+i]*YV[k*r+j];At[i*r+j]=s/Sd[j];
    }
    result->eigenvalues=malloc((size_t)r*sizeof(double complex));
    result->eigenvectors=malloc((size_t)(r*r)*sizeof(double complex));
    result->modes=malloc((size_t)(n*r)*sizeof(double complex));
    result->singular_values=malloc((size_t)r*sizeof(double));
    memcpy(result->singular_values,Sd,(size_t)r*sizeof(double));
    for(int i=0;i<r;i++){
        double re=At[i*r+i],im=0;
        if(i<r-1&&fabs(At[(i+1)*r+i])>1e-10){
            double a=At[i*r+i],b=At[i*r+(i+1)],c=At[(i+1)*r+i],d=At[(i+1)*r+(i+1)];
            double tr=a+d,det=a*d-b*c,disc=tr*tr-4*det;
            if(disc<0){re=0.5*tr;im=0.5*sqrt(-disc);
                result->eigenvalues[i]=re+im*I;result->eigenvalues[i+1]=re-im*I;i++;}
            else result->eigenvalues[i]=0.5*(tr+sqrt(disc));
        }else result->eigenvalues[i]=re+im*I;
    }
    for(int i=0;i<r;i++)for(int k=0;k<n;k++){
        double complex s=0;
        for(int j=0;j<r;j++)s+=YV[k*r+j]*creal(result->eigenvectors[j*r+i]);
        result->modes[i*n+k]=s;
    }
    double err=0,nrm=0;
    for(int i=0;i<n;i++)for(int j=0;j<ms;j++){
        double pred=0;for(int k=0;k<r;k++)pred+=creal(result->modes[k*n+i]*
            cexp(result->eigenvalues[k]*(double)j));
        double d=pred-Y[i*ms+j];err+=d*d;nrm+=Y[i*ms+j]*Y[i*ms+j];
    }
    result->reconstruction_error=(nrm>1e-15)?sqrt(err/nrm):0;
    free(Xc);free(U);free(Sd);free(Vt);free(At);free(YV);
}

void dmd_predict(const DMDResult *dmd,const double *x0,double t,double *x_pred){
    int r=dmd->r,n=dmd->n_state_dim;
    double complex *b=malloc((size_t)r*sizeof(double complex));
    dmd_mode_amplitudes(dmd,x0,b);
    for(int i=0;i<n;i++){double complex s=0;
        for(int j=0;j<r;j++)s+=dmd->modes[j*n+i]*cexp(dmd->eigenvalues[j]*t)*b[j];
        x_pred[i]=creal(s);}free(b);
}
void dmd_mode_amplitudes(const DMDResult *dmd,const double *x0,double complex *amp){
    int r=dmd->r,n=dmd->n_state_dim;
    for(int j=0;j<r;j++){amp[j]=0;
        for(int i=0;i<n;i++)amp[j]+=conj(dmd->modes[j*n+i])*x0[i];}
}
void compute_exact_dmd(const double *X,const double *Y,int n,int m,int r,double dt,
    DMDResult *result){compute_dmd(X,Y,n,m,r,dt,result);}
void compute_edmd(const double *X,const double *Y,int n,int m,const EDMDConfig *cfg,
    DMDResult *result){
    int ms=m-1,N=cfg->n_basis_functions;if(N<=0)N=n*5;
    double *PsiX=malloc((size_t)(N*ms)*sizeof(double));
    double *PsiY=malloc((size_t)(N*ms)*sizeof(double));
    double *psi=malloc((size_t)N*sizeof(double));
    for(int j=0;j<ms;j++){
        evaluate_basis_functions(X+j*n,n,cfg,psi);
        for(int i=0;i<N;i++)PsiX[i*ms+j]=psi[i];
        evaluate_basis_functions(Y+j*n,n,cfg,psi);
        for(int i=0;i<N;i++)PsiY[i*ms+j]=psi[i];
    }
    compute_dmd(PsiX,PsiY,N,m,(N<ms?N:ms)>10?10:(N<ms?N:ms),0.0,result);result->n_state_dim=n;
    free(PsiX);free(PsiY);free(psi);
}
void evaluate_basis_functions(const double *x,int n,const EDMDConfig *cfg,double *psi){
    switch(cfg->basis_type){
    case EDMD_BASIS_MONOMIAL:{
        int idx=0;psi[idx++]=1.0;
        for(int i=0;i<n;i++)psi[idx++]=x[i];
        if(cfg->max_poly_order>=2)
            for(int i=0;i<n;i++)for(int j=i;j<n;j++)psi[idx++]=x[i]*x[j];break;}
    case EDMD_BASIS_GAUSSIAN_RBF:{
        double sg=cfg->rbf_sigma>0?cfg->rbf_sigma:1.0;
        for(int i=0;i<cfg->n_basis_functions&&i<n;i++){
            double d=0;for(int k=0;k<n;k++)d+=(x[k]-(i/(double)n))*(x[k]-(i/(double)n));
            psi[i]=exp(-d/(2*sg*sg));}break;}
    default:for(int i=0;i<cfg->n_basis_functions;i++)psi[i]=x[i%n];break;
    }
}
void compute_kernel_edmd(const double *X,const double *Y,int n,int m,
    double (*kernel)(const double*,const double*,int,double),double gamma,
    DMDResult *result){
    int ms=m-1;double *Kmat=malloc((size_t)(ms*ms)*sizeof(double));
    for(int i=0;i<ms;i++)for(int j=0;j<ms;j++)Kmat[i*ms+j]=kernel(X+i*n,X+j*n,n,gamma);
    compute_dmd(Kmat,Kmat+ms*ms,ms,m,(ms<10?ms:10),0.0,result);free(Kmat);
}
void identify_koopman_eigenfunctions(const DMDResult *dmd,const EDMDConfig *cfg,
    const double *X,int n,int m,KoopmanEigenfunction *eigfuncs){
    int r=dmd->r;
    for(int j=0;j<r;j++){
        eigfuncs[j].eigenvalue=dmd->eigenvalues[j];eigfuncs[j].n_samples=m-1;
        eigfuncs[j].eigenfunction=malloc((size_t)(m-1)*sizeof(double complex));
        for(int t=0;t<m-1;t++){double complex v=0;
            for(int i=0;i<n;i++)v+=conj(dmd->modes[j*n+i])*X[i*(m-1)+t];
            eigfuncs[j].eigenfunction[t]=v;}
    }
}

void koopman_mode_decomposition(const DMDResult *dmd,const double *x0,int n,
    KoopmanDecomposition *decomp){
    memset(decomp,0,sizeof(KoopmanDecomposition));int r=dmd->r;decomp->n_modes=r;
    decomp->modes=malloc((size_t)r*sizeof(KoopmanMode));
    decomp->amplitudes=malloc((size_t)r*sizeof(double complex));
    for(int j=0;j<r;j++){
        decomp->modes[j].eigenvalue=dmd->eigenvalues[j];
        decomp->modes[j].mode=malloc((size_t)n*sizeof(double complex));
        for(int i=0;i<n;i++)decomp->modes[j].mode[i]=dmd->modes[j*n+i];
        decomp->modes[j].amplitude=cabs(dmd->eigenvalues[j]);
        decomp->modes[j].frequency_hz=cimag(dmd->eigenvalues[j])/(2*M_PI);
        decomp->modes[j].growth_rate=creal(dmd->eigenvalues[j]);
    }
    dmd_mode_amplitudes(dmd,x0,decomp->amplitudes);
    decomp->prediction_error=dmd->reconstruction_error;
}
void koopman_predict(const KoopmanDecomposition *decomp,double t,int n,double *xp){
    for(int i=0;i<n;i++){xp[i]=0;for(int j=0;j<decomp->n_modes;j++)
        xp[i]+=creal(decomp->modes[j].mode[i]*
            cexp(decomp->modes[j].eigenvalue*t)*decomp->amplitudes[j]);}
}
int carleman_embedding_dimension(int n,int max_k){
    if(n==1)return max_k;double s=0,pow=1;
    for(int j=1;j<=max_k;j++){pow*=n;s+=pow;}return(int)(s+0.5);
}
void carleman_embedding_matrix(void (*fp)(const double*,double*,int,void*),
    int n,int poly_deg,int trunc_k,double *A,int N,void *ud){
    memset(A,0,(size_t)(N*N)*sizeof(double));
    double eps=1e-6,*x0=calloc((size_t)n,sizeof(double));
    double *xp=malloc((size_t)n*sizeof(double)),*f0=malloc((size_t)n*sizeof(double));
    double *fp2=malloc((size_t)n*sizeof(double));
    fp(x0,f0,n,ud);
    for(int j=0;j<n;j++){memcpy(xp,x0,(size_t)n*sizeof(double));xp[j]+=eps;
        fp(xp,fp2,n,ud);for(int i=0;i<n;i++)A[i*N+j]=(fp2[i]-f0[i])/eps;}
    int offset=n;
    for(int ord=2;ord<=trunc_k;ord++){
        int nt=(int)(pow((double)n,(double)ord)+0.5);
        if(offset+nt>N)break;
        for(int i=0;i<nt&&(offset+i)<N;i++)A[(offset+i)*N+(offset+i)]=-0.1*ord;
        offset+=nt;
    }
    free(x0);free(xp);free(f0);free(fp2);
}
void free_dmd_result(DMDResult *r){
    if(r){free(r->eigenvalues);free(r->eigenvectors);free(r->modes);
           free(r->singular_values);}
}
void free_koopman_decomposition(KoopmanDecomposition *d){
    if(d){for(int i=0;i<d->n_modes;i++)free(d->modes[i].mode);
          free(d->modes);free(d->amplitudes);free(d);}
}
void free_koopman_eigenfunctions(KoopmanEigenfunction *e){if(e)free(e->eigenfunction);}
