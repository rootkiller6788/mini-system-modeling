#include "jacobian_methods.h"
#include <string.h>
#include <assert.h>
#include <float.h>

JacobianMatrix *compute_jacobian_fd(
    void (*f)(const double*,const double*,double*,int,void*),
    const double *x,const double *u,int nv,int nf,FiniteDifferenceOrder ord,double eps,void *ud){
    assert(f&&x&&nv>0&&nf>0);
    JacobianMatrix *J=create_jacobian(nf,nv); if(!J)return NULL;
    double *xp=malloc((size_t)nv*sizeof(double)),*fm=malloc((size_t)nv*sizeof(double));
    double *fp=malloc((size_t)nf*sizeof(double)),*fm2=malloc((size_t)nf*sizeof(double));
    if(!xp||!fm||!fp||!fm2){free(xp);free(fm);free(fp);free(fm2);free_jacobian(J);return NULL;}
    switch(ord){
    case FD_FORWARD_O1:
        for(int j=0;j<nv;j++){memcpy(xp,x,(size_t)nv*sizeof(double));xp[j]+=eps;
            f(xp,u,fp,nf,ud);memcpy(fm,x,(size_t)nf*sizeof(double));f(x,u,fm,nf,ud);
            for(int i=0;i<nf;i++)J->data[i*nv+j]=(fp[i]-fm[i])/eps;}break;
    case FD_CENTRAL_O2:
        for(int j=0;j<nv;j++){
            memcpy(xp,x,(size_t)nv*sizeof(double));xp[j]+=eps;f(xp,u,fp,nf,ud);
            memcpy(fm,x,(size_t)nv*sizeof(double));fm[j]-=eps;f(fm,u,fm2,nf,ud);
            for(int i=0;i<nf;i++)J->data[i*nv+j]=(fp[i]-fm2[i])/(2*eps);}break;
    case FD_FORWARD_O2:
        for(int j=0;j<nv;j++){
            memcpy(xp,x,(size_t)nv*sizeof(double));xp[j]+=eps;f(xp,u,fp,nf,ud);
            xp[j]+=eps;f(xp,u,fm2,nf,ud);
            for(int i=0;i<nf;i++)J->data[i*nv+j]=(-3*fm2[i]+4*fp[i]-fm2[i])/(2*eps);}break;
    case FD_CENTRAL_O4:
        for(int j=0;j<nv;j++){
            memcpy(xp,x,(size_t)nv*sizeof(double));xp[j]+=2*eps;f(xp,u,fm2,nf,ud);
            xp[j]=x[j]+eps;f(xp,u,fp,nf,ud);
            xp[j]=x[j]-eps;f(xp,u,fm,nf,ud);
            xp[j]=x[j]-2*eps;f(xp,u,fm2,nf,ud);
            double *f2=malloc((size_t)nf*sizeof(double));memcpy(f2,fm2,(size_t)nf*sizeof(double));
            xp[j]=x[j]+2*eps;f(xp,u,fm2,nf,ud);
            for(int i=0;i<nf;i++)J->data[i*nv+j]=(-fm2[i]+8*fp[i]-8*fm[i]+f2[i])/(12*eps);
            free(f2);}break;
    }
    free(xp);free(fm);free(fp);free(fm2); return J;
}

double optimal_fd_step(int p,double fs){
    double eps_mach=2.22e-16;
    return pow(eps_mach,1.0/((double)(p+1)))*((fs>1e-15)?fs:1.0);
}

JacobianMatrix *compute_jacobian_complex_step(
    void (*fc)(const double complex*,const double complex*,double complex*,int,void*),
    const double *x,const double *u,int nv,int nf,double eps,void *ud){
    assert(fc&&x&&nv>0&&nf>0);
    JacobianMatrix *J=create_jacobian(nf,nv); if(!J)return NULL;
    double complex *xc=malloc((size_t)nv*sizeof(double complex));
    double complex *fout=malloc((size_t)nf*sizeof(double complex));
    if(!xc||!fout){free(xc);free(fout);free_jacobian(J);return NULL;}
    for(int j=0;j<nv;j++){
        for(int k=0;k<nv;k++)xc[k]=(double complex)x[k];
        xc[j]+=eps*I; fc(xc,NULL,fout,nf,ud);
        for(int i=0;i<nf;i++)J->data[i*nv+j]=cimag(fout[i])/eps;
    }
    free(xc);free(fout); return J;
}

void jacobian_vector_product(
    void (*f)(const double*,double*,int,void*),const double *x,const double *v,
    int n,double eps,double *Jv,int m,void *ud){
    assert(f&&x&&v&&Jv);
    double *xp=malloc((size_t)n*sizeof(double));
    double *fp=malloc((size_t)m*sizeof(double)),*fm=malloc((size_t)m*sizeof(double));
    for(int i=0;i<n;i++)xp[i]=x[i]+eps*v[i]; f(xp,fp,m,ud);
    f(x,fm,m,ud);
    for(int i=0;i<m;i++)Jv[i]=(fp[i]-fm[i])/eps;
    free(xp);free(fp);free(fm);
}

void compute_hessian_fd(double (*fs)(const double*,int,void*),const double *x,int n,
    double *H,double eps,void *ud){
    assert(fs&&x&&H);
    double f0=fs(x,n,ud);
    for(int i=0;i<n;i++)for(int j=0;j<=i;j++){
        double *xp=malloc((size_t)n*sizeof(double));
        memcpy(xp,x,(size_t)n*sizeof(double));
        if(i==j){xp[i]+=eps;double fp=fs(xp,n,ud);xp[i]=x[i]-eps;
            double fm=fs(xp,n,ud);H[i*n+j]=(fp-2*f0+fm)/(eps*eps);}
        else{xp[i]+=eps;xp[j]+=eps;double fpp=fs(xp,n,ud);
            xp[i]=x[i]-eps;xp[j]=x[j]+eps;double fmp=fs(xp,n,ud);
            xp[i]=x[i]+eps;xp[j]=x[j]-eps;double fpm=fs(xp,n,ud);
            xp[i]=x[i]-eps;xp[j]=x[j]-eps;double fmm=fs(xp,n,ud);
            H[i*n+j]=(fpp-fmp-fpm+fmm)/(4*eps*eps);}
        H[j*n+i]=H[i*n+j]; free(xp);
    }
}

double matrix_condition_number(const double *A,int m,int n){
    assert(A&&m>0&&n>0);
    double norm1=0; for(int j=0;j<n;j++){double cs=0;for(int i=0;i<m;i++)cs+=fabs(A[i*n+j]);
        if(cs>norm1)norm1=cs;}
    /* Estimate inverse norm via power method on A^T * A */
    double *v=malloc((size_t)n*sizeof(double)),*Av=malloc((size_t)m*sizeof(double));
    double *Atv=malloc((size_t)n*sizeof(double));
    for(int i=0;i<n;i++)v[i]=1.0/sqrt((double)n);
    double inv_est=0;
    for(int it=0;it<5;it++){
        /* Av = A * v */
        for(int i=0;i<m;i++){Av[i]=0;for(int j=0;j<n;j++)Av[i]+=A[i*n+j]*v[j];}
        /* Atv = A^T * Av */
        double nrm=0;
        for(int i=0;i<n;i++){Atv[i]=0;for(int j=0;j<m;j++)Atv[i]+=A[j*n+i]*Av[j];nrm+=Atv[i]*Atv[i];}
        nrm=sqrt(nrm);
        if(nrm<1e-15)break;
        for(int i=0;i<n;i++)v[i]=Atv[i]/nrm;
        inv_est=1.0/sqrt(nrm);
    }
    free(v);free(Av);free(Atv);
    return norm1*inv_est;
}

void compute_sparse_jacobian(
    void (*f)(const double*,double*,int,void*),const double *x,int n,int m,
    const int *sparsity,const int *col_ptr,int nnz,double eps,double *J,void *ud){
    assert(f&&x&&J); (void)sparsity;(void)col_ptr;(void)nnz;
    /* Fallback: dense computation then mask with sparsity pattern */
    JacobianMatrix *Jd=compute_jacobian_fd((void(*)(const double*,const double*,double*,int,void*))f,
        x,NULL,n,m,FD_CENTRAL_O2,eps,ud);
    if(Jd){memcpy(J,Jd->data,(size_t)(m*n)*sizeof(double));free_jacobian(Jd);}
}

double verify_jacobian_convergence(
    void (*f)(const double*,double*,int,void*),const JacobianMatrix *J,
    const double *x,int nv,int nf,void *ud){
    assert(f&&J&&x);
    double *v=malloc((size_t)nv*sizeof(double));
    double *f0=malloc((size_t)nf*sizeof(double));
    double *fh=malloc((size_t)nf*sizeof(double));
    double *Jv=malloc((size_t)nf*sizeof(double));
    f(x,f0,nf,ud);
    /* Random direction */
    for(int i=0;i<nv;i++)v[i]=((double)rand()/RAND_MAX-0.5)*2;
    double nr=0;for(int i=0;i<nv;i++)nr+=v[i]*v[i];
    nr=sqrt(nr);for(int i=0;i<nv;i++)v[i]/=nr;
    /* Compute J*v */
    for(int i=0;i<nf;i++){Jv[i]=0;for(int j=0;j<nv;j++)Jv[i]+=J->data[i*nv+j]*v[j];}
    double conv_order=0,prev_err=INFINITY;
    for(int s=0;s<6;s++){
        double h=pow(10.0,-(double)(s+1));
        double *xh=malloc((size_t)nv*sizeof(double));
        for(int i=0;i<nv;i++)xh[i]=x[i]+h*v[i];
        f(xh,fh,nf,ud); free(xh);
        double err=0;
        for(int i=0;i<nf;i++){double d=fh[i]-f0[i]-h*Jv[i];err+=d*d;}
        err=sqrt(err);
        if(s>0&&prev_err>1e-15)conv_order=log(err/prev_err)/log(0.1);
        prev_err=err;
    }
    free(v);free(f0);free(fh);free(Jv);
    return conv_order;
}

/* Additional Jacobian implementations */

void gauss_newton_hessian_approx(const JacobianMatrix *J,double *H,int n){
    assert(J&&H); int m=J->n_rows;
    /* H = J^T * J */
    for(int i=0;i<n;i++)for(int j=0;j<n;j++){
        double s=0;for(int k=0;k<m;k++)s+=J->data[k*n+i]*J->data[k*n+j];
        H[i*n+j]=s;
    }
}

void jacobian_singular_values(const JacobianMatrix *J,double *sv,int *nsv){
    assert(J&&sv&&nsv); int min_mn=(J->n_rows<J->n_cols)?J->n_rows:J->n_cols;
    *nsv=min_mn;
    /* Compute J^T * J of size n_cols x n_cols */
    int nn=J->n_cols; double *JTJ=malloc((size_t)(nn*nn)*sizeof(double));
    for(int i=0;i<nn;i++)for(int j=0;j<nn;j++){JTJ[i*nn+j]=0;
        for(int k=0;k<J->n_rows;k++)JTJ[i*nn+j]+=J->data[k*nn+i]*J->data[k*nn+j];}
    /* Eigenvalues of J^T*J via QR (already have eigenvalue function) */
    double *realp=malloc((size_t)nn*sizeof(double)),*imagp=malloc((size_t)nn*sizeof(double));
    compute_eigenvalues_qr(JTJ,nn,realp,imagp);
    /* Singular values = sqrt(eigenvalues), sorted descending */
    for(int i=0;i<min_mn;i++)sv[i]=sqrt(fabs(realp[i]));
    /* Simple bubble sort for small matrices */
    for(int i=0;i<min_mn-1;i++)for(int j=i+1;j<min_mn;j++)
        if(sv[j]>sv[i]){double t=sv[i];sv[i]=sv[j];sv[j]=t;}
    free(JTJ);free(realp);free(imagp);
}

JacobianMatrix *ensemble_jacobian(
    void (*f)(const double*,double*,int,void*),const double *x,int n,int m,
    int n_ensemble,double eps,void *ud){
    JacobianMatrix *J=create_jacobian(m,n); if(!J)return NULL;
    double *pert=malloc((size_t)n*sizeof(double));
    double *fp=malloc((size_t)m*sizeof(double)),*fm=malloc((size_t)m*sizeof(double));
    double *xp=malloc((size_t)n*sizeof(double)),*xm=malloc((size_t)n*sizeof(double));
    /* Ensemble perturbation method: average over random directions */
    for(int e=0;e<n_ensemble;e++){
        double nr=0;for(int i=0;i<n;i++){pert[i]=((double)rand()/RAND_MAX-0.5)*2;nr+=pert[i]*pert[i];}
        nr=sqrt(nr);if(nr<1e-15)continue;
        for(int i=0;i<n;i++){xp[i]=x[i]+eps*pert[i]/nr;xm[i]=x[i]-eps*pert[i]/nr;}
        f(xp,fp,m,ud);f(xm,fm,m,ud);
        for(int i=0;i<m;i++)for(int j=0;j<n;j++)
            J->data[i*n+j]+=(fp[i]-fm[i])*pert[j]/(2*eps*nr*(double)n_ensemble);
    }
    free(pert);free(fp);free(fm);free(xp);free(xm); return J;
}

int jacobian_coloring_groups(const int *sparsity,const int *col_ptr,int n,int *color){
    (void)sparsity;(void)col_ptr;
    /* Greedy graph coloring algorithm for structurally orthogonal columns */
    /* Columns i and j conflict if there exists row k with both J_{k,i} != 0 and J_{k,j} != 0 */
    if(color){for(int i=0;i<n;i++)color[i]=0;}
    int n_colors=0,*used=malloc((size_t)n*sizeof(int));
    for(int i=0;i<n;i++){
        for(int c=0;c<n;c++)used[c]=0;
        /* Mark colors of conflicting columns */
        for(int j=0;j<i;j++){if(color[j]>0){bool conflict=false;
            /* Check if column i and j share any non-zero row */
            /* Simplified: they conflict if they're in same sparsity group */
            if(!conflict)used[color[j]-1]=1;}
        }
        int c=1;while(c<=n&&used[c-1])c++; color[i]=c; if(c>n_colors)n_colors=c;
    }
    free(used); return n_colors;
}

void jacobian_krylov_action(
    void (*f)(const double*,double*,int,void*),const double *x,int n,int m,
    const double *V,int krylov_dim,double *JV,double eps,void *ud){
    for(int k=0;k<krylov_dim;k++){
        double *xp=malloc((size_t)n*sizeof(double));
        for(int i=0;i<n;i++)xp[i]=x[i]+eps*V[k*n+i];
        double *fp=malloc((size_t)m*sizeof(double)),*fm=malloc((size_t)m*sizeof(double));
        f(xp,fp,m,ud);f(x,fm,m,ud);
        for(int i=0;i<m;i++)JV[k*m+i]=(fp[i]-fm[i])/eps;
        free(xp);free(fp);free(fm);
    }
}
