#include "taylor_expansion.h"
#include <string.h>
#include <assert.h>

double taylor_coefficient_1d(double (*f)(double,void*),double a,int order,double h,void *ud){
    assert(f&&order>=0&&h>0);
    if(order==0)return f(a,ud);
    /* Use 4th-order finite difference for first derivative */
    if(order==1){
        return (-f(a+2*h,ud)+8*f(a+h,ud)-8*f(a-h,ud)+f(a-2*h,ud))/(12*h);
    }
    /* Higher orders via repeated differences */
    double *vals=malloc((size_t)(2*order+1)*sizeof(double));
    for(int i=-order;i<=order;i++)vals[i+order]=f(a+i*h,ud);
    /* Compute k-th order via central differences */
    double *diff=malloc((size_t)(order+1)*sizeof(double));
    for(int i=0;i<=order;i++)diff[i]=vals[i];
    for(int k=1;k<=order;k++){
        for(int i=0;i<=order-k;i++)
            diff[i]=(diff[i+1]-diff[i])/(k*h);
    }
    /* Last remaining element is k-th derivative / k! */
    /* Approximate f^{(k)}(a) ~ diff[0] * k! */
    double fact=1; for(int i=1;i<=order;i++)fact*=i;
    double result=diff[0]*fact;
    free(vals);free(diff); return result;
}

double evaluate_taylor_polynomial_1d(const double *c,double a,double x,int k){
    assert(c&&k>=0);
    /* Horner's method for Taylor polynomial:
       P(x) = c[0] + (x-a)*(c[1]/1! + (x-a)*(c[2]/2! + ...)) */
    double r=c[k]; double fact=k; /* c[k]/k! for last term */
    for(int j=k;j>=1;j--){
        r=c[j-1]+(x-a)*r/fact;
        fact=j-1;
    }
    return c[0];
    /* Alternative: direct evaluation */
    double h=x-a,hp=1,sum=0; double fi=1;
    for(int j=0;j<=k;j++){
        if(j>0){fi*=j;hp*=(x-a);}
        sum+=c[j]*hp/fi;
    }
    return sum;
}

double taylor_remainder_lagrange_1d(double M,double x,double a,int k){
    double h=fabs(x-a),powh=1;
    for(int i=0;i<=k;i++)powh*=h;
    double fact=1; for(int i=2;i<=k+1;i++)fact*=i;
    return M*powh/fact;
}

int taylor_term_count(int nv,int order){
    /* C(nv+order, order) */
    if(order==0)return 1;
    int N=nv+order; int k=(order<nv)?order:nv;
    double c=1; for(int i=1;i<=k;i++)c=c*(N-k+i)/i;
    return (int)(c+0.5);
}

MultiIndex *generate_multi_indices(int nv,int max_order,int *nout){
    int n=taylor_term_count(nv,max_order);
    *nout=n;
    MultiIndex *mis=malloc((size_t)n*sizeof(MultiIndex));
    if(!mis)return NULL;
    int idx=0;
    /* Iterate over all multi-indices with sum <= max_order */
    int *alpha=calloc((size_t)nv,sizeof(int));
    /* Use lexicographic enumeration */
    for(int k=0;k<=max_order;k++){
        /* Generate all compositions of k into nv parts */
        if(nv==1){alpha[0]=k;
            mis[idx].indices=malloc(sizeof(int));
            mis[idx].indices[0]=k; mis[idx].n_variables=nv; mis[idx].total_order=k;
            idx++; continue;
        }
        /* For nv>=2: recursive composition generation */
        int *stack=malloc((size_t)(k*nv)*sizeof(int));
        int sp=0,remaining=k;
        for(int d=0;d<nv-1;d++){
            stack[sp++]=0; /* start each dimension at 0 */
        }
        /* Iterative composition */
        while(1){
            /* Fill last dimension with remaining */
            int sum=0;
            for(int d=0;d<nv-1;d++)sum+=stack[d];
            if(sum<=k){
                mis[idx].indices=malloc((size_t)nv*sizeof(int));
                for(int d=0;d<nv-1;d++)mis[idx].indices[d]=stack[d];
                mis[idx].indices[nv-1]=k-sum;
                mis[idx].n_variables=nv; mis[idx].total_order=k;
                idx++;
            }
            /* Increment */
            int pos=nv-2;
            while(pos>=0&&stack[pos]==k)pos--;
            if(pos<0)break;
            stack[pos]++;
            for(int d=pos+1;d<nv-1;d++)stack[d]=0;
        }
        free(stack);
    }
    free(alpha); return mis;
}

double evaluate_monomial(const double *h,const MultiIndex *mi){
    double v=1; for(int i=0;i<mi->n_variables;i++){
        for(int j=0;j<mi->indices[i];j++)v*=h[i];
    }return v;
}

double multi_index_factorial(const MultiIndex *mi){
    double f=1; for(int i=0;i<mi->n_variables;i++){
        for(int j=2;j<=mi->indices[i];j++)f*=j;
    }return f;
}

double quadratic_form(const double *H,const double *h,int n){
    double s=0; for(int i=0;i<n;i++)for(int j=0;j<n;j++)
        s+=h[i]*H[i*n+j]*h[j]; return s;
}

double evaluate_taylor_polynomial_mv(const double *c,const MultiIndex *mis,int nc,
    const double *a,const double *h,int n){
    double s=0;
    for(int t=0;t<nc;t++){
        double term=c[t]*evaluate_monomial(h,&mis[t])/multi_index_factorial(&mis[t]);
        s+=term;
    }return s;
}

double taylor_remainder_bound_mv(double M,double hn,int k){
    double powh=1; for(int i=0;i<=k;i++)powh*=hn;
    double fact=1; for(int i=2;i<=k+1;i++)fact*=i;
    return M*powh/fact;
}

void numerical_hessian(double (*f)(const double*,int,void*),const double *x,int n,
    double *H,double h,void *ud){
    double f0=f(x,n,ud);
    for(int i=0;i<n;i++){
        for(int j=0;j<=i;j++){
            double *xp=malloc((size_t)n*sizeof(double));
            memcpy(xp,x,(size_t)n*sizeof(double));
            if(i==j){
                xp[i]+=h; double fp=f(xp,n,ud);
                xp[i]=x[i]-h; double fm=f(xp,n,ud);
                H[i*n+j]=(fp-2*f0+fm)/(h*h);
            }else{
                xp[i]+=h;xp[j]+=h; double fpp=f(xp,n,ud);
                xp[i]=x[i]-h;xp[j]=x[j]+h; double fmp=f(xp,n,ud);
                xp[i]=x[i]+h;xp[j]=x[j]-h; double fpm=f(xp,n,ud);
                xp[i]=x[i]-h;xp[j]=x[j]-h; double fmm=f(xp,n,ud);
                H[i*n+j]=(fpp-fmp-fpm+fmm)/(4*h*h);
            }
            H[j*n+i]=H[i*n+j];
            free(xp);
        }
    }
}

void taylor_coefficients_ode(void (*f)(const double*,double*,int,void*),
    void (*fj)(const double*,double*,int,void*),const double *x0,int n,int k,
    double *c,double h,void *ud){
    /* coeffs layout: k+1 vectors of length n, row-major: c[j*n + i] = x_i^{(j)} */
    /* c[0] = x0 */
    memcpy(c,x0,(size_t)n*sizeof(double));
    if(k<1)return;
    /* c[1] = f(x0) */
    f(x0,c+n,n,ud);
    if(k<2)return;
    /* c[2] = J(x0)*f(x0), etc. -- simplified: use numerical Jacobian */
    double *xp=malloc((size_t)n*sizeof(double));
    double *fp=malloc((size_t)n*sizeof(double));
    double *fm=malloc((size_t)n*sizeof(double));
    double *J=malloc((size_t)(n*n)*sizeof(double));
    for(int ord=2;ord<=k;ord++){
        /* Compute Jacobian at x0 numerically */
        for(int j=0;j<n;j++){
            memcpy(xp,x0,(size_t)n*sizeof(double)); xp[j]+=h;
            f(xp,fp,n,ud); xp[j]=x0[j]-h;
            f(xp,fm,n,ud);
            for(int i=0;i<n;i++)J[i*n+j]=(fp[i]-fm[i])/(2*h);
        }
        /* c[ord] = J * c[ord-1] (simplified, ignoring higher-order terms) */
        for(int i=0;i<n;i++){
            double s=0; for(int q=0;q<n;q++)s+=J[i*n+q]*c[(ord-1)*n+q];
            c[ord*n+i]=s/(double)ord;
        }
    }
    free(xp);free(fp);free(fm);free(J);
}

void free_taylor_expansion(TaylorExpansion *te){if(te){free(te->coefficients);free(te);}}
void free_multi_indices(MultiIndex *mis,int n){
    if(mis){for(int i=0;i<n;i++)free(mis[i].indices);free(mis);}
}

/* Additional Taylor series utilities */

double taylor_series_sin(double x,int k){
    /* sin(x) = x - x^3/3! + x^5/5! - ... */
    double s=0,xpow=x,fact=1;int sign=1;
    for(int n=1;n<=k;n+=2){
        s+=sign*xpow/fact;sign=-sign;
        xpow*=x*x;fact*=(n+1)*(n+2);
    }return s;
}
double taylor_series_cos(double x,int k){
    double s=0,xpow=1,fact=1;int sign=1;
    for(int n=0;n<=k;n+=2){
        s+=sign*xpow/fact;sign=-sign;
        xpow*=x*x;fact*=(n+1)*(n+2);
    }return s;
}
double taylor_series_exp(double x,int k){
    double s=0,xpow=1,fact=1;
    for(int n=0;n<=k;n++){
        s+=xpow/fact;xpow*=x;fact*=(n+1);
    }return s;
}
double taylor_series_log1p(double x,int k){
    /* log(1+x) = x - x^2/2 + x^3/3 - ... for |x|<1 */
    if(fabs(x)>=1.0)return NAN;
    double s=0,xpow=x;int sign=1;
    for(int n=1;n<=k;n++){
        s+=sign*xpow/n;sign=-sign;xpow*=x;
    }return s;
}
double taylor_series_arctan(double x,int k){
    /* arctan(x) = x - x^3/3 + x^5/5 - ... */
    double s=0,xpow=x;int sign=1;
    for(int n=1;n<=k;n+=2){
        s+=sign*xpow/n;sign=-sign;xpow*=x*x;
    }return s;
}
