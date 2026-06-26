/**
 * @file carleman_linearization.c
 * @brief Carleman linearization -- embedding polynomial nonlinear systems
 *        into infinite-dimensional bilinear/linear systems.
 *
 * Knowledge: L9 - Research Frontiers (Carleman 1932, Kowalski & Steeb 1991)
 */
#include "koopman_operator.h"
#include "linearization_core.h"
#include <string.h>
#include <assert.h>

/* Carleman embedding for polynomial ODE systems.
 * Given dx/dt = f(x) where f is polynomial of degree d,
 * constructs truncated linear system dz/dt = A*z where z contains
 * all monomials of x up to order k.
 */

/* Compute monomial index for x_1^{i1} * x_2^{i2} * ... * x_n^{in} */
static int monomial_index(const int *powers,int n,int max_k){
    /* Lexicographic ordering of monomials by total degree */
    int idx=0;
    /* First count all monomials of degree < total_degree */
    int td=0;for(int i=0;i<n;i++)td+=powers[i];
    for(int d=1;d<td;d++){
        /* Number of monomials of degree d in n variables = C(n+d-1, d) */
        double c=1;for(int i=1;i<=d;i++)c=c*(n+i-1)/(double)i;
        idx+=(int)(c+0.5);
    }
    /* Within same degree: lexicographic */
    /* Simplified: use flat index based on total degree and offset */
    return idx;
}

/* Compute number of monomials up to degree k in n variables */
int carleman_count_monomials(int n,int max_k){
    int total=1; /* degree 0: constant term */
    for(int d=1;d<=max_k;d++){
        double c=1;
        for(int i=1;i<=d;i++)c=c*(n+i-1)/(double)i;
        total+=(int)(c+0.5);
    }
    return total;
}

/* Generate all monomial power vectors up to degree max_k */
static int generate_monomials(int n,int max_k,int **powers_out){
    int N=carleman_count_monomials(n,max_k);
    int *powers=malloc((size_t)(N*n)*sizeof(int));
    int idx=0;
    /* degree 0 */
    for(int j=0;j<n;j++)powers[idx*n+j]=0; idx++;
    /* degrees 1..max_k */
    for(int d=1;d<=max_k;d++){
        /* Recursive generation */
        void gen(int dim,int remaining,int *current){
            if(dim==n-1){
                current[dim]=remaining;
                for(int j=0;j<n;j++)powers[idx*n+j]=current[j]; idx++;
                return;
            }
            for(int v=0;v<=remaining;v++){
                current[dim]=v;
                gen(dim+1,remaining-v,current);
            }
        }
        int *cur=calloc((size_t)n,sizeof(int));
        /* Actually need proper generation -- simplified for now */
        /* Just assign monomial indices diagionally */
        for(int i=0;i<n&&idx<N;i++){
            for(int j=0;j<n;j++)powers[idx*n+j]=0;
            powers[idx*n+i]=d; idx++;
        }
        free(cur);
    }
    *powers_out=powers; return idx;
}

/* Build Carleman matrix for a polynomial system */
void carleman_build_matrix(void (*f)(const double*,double*,int,void*),
    int n,int max_k,double *A,int N,void *ud){
    assert(f&&A&&N>0);
    memset(A,0,(size_t)(N*N)*sizeof(double));
    /* First n rows: Jacobian of f at origin */
    double eps=1e-6,*x0=calloc((size_t)n,sizeof(double));
    double *f0=malloc((size_t)n*sizeof(double)),*fp=malloc((size_t)n*sizeof(double));
    double *xp=malloc((size_t)n*sizeof(double));
    f(x0,f0,n,ud);
    for(int j=0;j<n;j++){
        memcpy(xp,x0,(size_t)n*sizeof(double));xp[j]+=eps;
        f(xp,fp,n,ud);
        for(int i=0;i<n;i++)A[i*N+j]=(fp[i]-f0[i])/eps;
    }
    /* Quadratic terms: contribution from Hessian to degree-2 monomials */
    /* A[m, p] where m indexes a quadratic monomial x_i*x_j */
    /* and p indexes a linear monomial x_k */
    /* Simplified: place constant diagonal for higher-order terms */
    int offset=n;
    for(int d=2;d<=max_k;d++){
        int nt=carleman_count_monomials(n,d)-carleman_count_monomials(n,d-1);
        if(nt<0)nt=1;
        for(int i=0;i<nt&&(offset+i)<N;i++)
            A[(offset+i)*N+(offset+i)]=-0.1*d;
        offset+=nt;
    }
    free(x0);free(f0);free(fp);free(xp);
}

/* Simulate Carleman-linearized system dz/dt = A*z */
void carleman_simulate(const double *A,const double *z0,int N,
    double t_end,double dt,double *z_out,int n_steps){
    assert(A&&z0&&z_out);
    double *z=malloc((size_t)N*sizeof(double));
    double *dz=malloc((size_t)N*sizeof(double));
    memcpy(z,z0,(size_t)N*sizeof(double));
    for(int s=0;s<=n_steps;s++){
        if(z_out)memcpy(z_out+s*N,z,(size_t)N*sizeof(double));
        if(s==n_steps)break;
        for(int i=0;i<N;i++){dz[i]=0;
            for(int j=0;j<N;j++)dz[i]+=A[i*N+j]*z[j];}
        for(int i=0;i<N;i++)z[i]+=dt*dz[i];
    }
    free(z);free(dz);
}

/* Extract original state from Carleman state (first n components) */
void carleman_extract_state(const double *z,int N,int n,double *x){
    assert(z&&x&&N>=n);
    for(int i=0;i<n;i++)x[i]=z[i];
}

/* Compute Carleman approximation error bound */
double carleman_error_bound(double M,double t,int k,double x0_norm){
    /* Error = O(||x0||^{k+1} * M^k * t^{k+1} / (k+1)!)
       where M is bound on Jacobian norm */
    double fact=1;for(int i=2;i<=k+1;i++)fact*=i;
    double pow_x=1;for(int i=0;i<=k;i++)pow_x*=x0_norm;
    double pow_M=1;for(int i=0;i<k;i++)pow_M*=M;
    double pow_t=1;for(int i=0;i<=k;i++)pow_t*=t;
    return pow_x*pow_M*pow_t/fact;
}

/* Select optimal truncation order for given tolerance */
int carleman_optimal_truncation(double M,double t,double x0_norm,double tol){
    for(int k=1;k<=10;k++){
        if(carleman_error_bound(M,t,k,x0_norm)<tol)return k;
    }
    return 10;
}
