#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "linearization_core.h"
#include "taylor_expansion.h"
#include "describing_function.h"
#include "jacobian_methods.h"
#include "koopman_operator.h"
#include "gain_scheduling.h"

static int tr=0,tp=0;
#define T(n) do{tr++;printf("  TEST %s... ",n);}while(0)
#define P() do{tp++;printf("PASS\n");}while(0)

static void dyn_lin1(const double*x,const double*u,double*dx,int n,void*d){
    (void)u;(void)d;dx[0]=-x[0];}
static void dyn_pend(const double*x,const double*u,double*dx,int n,void*d){
    (void)d;double g=9.81,l=1.0,b=0.1,m=1.0;
    dx[0]=x[1];dx[1]=(g/l)*sin(x[0])-b*x[1]+(u?u[0]:0)/(m*l*l);}
static void out_pend(const double*x,const double*u,double*y,int p,void*d){
    (void)u;(void)d;y[0]=x[0];}
static double sin_wrap(double x,void*d){(void)d;return sin(x);}

static void t1(void){T("create_equilibrium");
    EquilibriumPoint*e=create_equilibrium(2,1);assert(e);assert(e->n_states==2);
    free_equilibrium(e);P();}
static void t2(void){T("jacobian_create");
    JacobianMatrix*J=create_jacobian(3,2);assert(J);jacobian_set(J,0,0,1.5);
    assert(fabs(jacobian_get(J,0,0)-1.5)<1e-10);free_jacobian(J);P();}
static void t3(void){T("numerical_jacobian");
    double x[2]={0,0};JacobianMatrix*J=numerical_jacobian_central(dyn_pend,x,NULL,2,2,1e-6,NULL);
    assert(J);assert(fabs(jacobian_get(J,0,1)-1.0)<0.1);free_jacobian(J);P();}
static void t4(void){T("linearize_system");
    double xe[2]={0,0},ue[1]={0};
    LinearizedSystem*s=linearize_system(dyn_pend,out_pend,xe,ue,2,1,1,LINEARIZE_NUMERICAL_CENTRAL,1e-6,NULL);
    assert(s);assert(s->n_states==2);free_linearized_system(s);P();}
static void t5(void){T("lyapunov_indirect");
    double xe[2]={0,0},ue[1]={0};
    LinearizedSystem*s=linearize_system(dyn_pend,out_pend,xe,ue,2,1,1,LINEARIZE_NUMERICAL_CENTRAL,1e-6,NULL);
    assert(s);LyapunovIndirectResult*r=apply_lyapunov_indirect(s);
    assert(r);assert(r->is_unstable);free_lyapunov_indirect_result(r);free_linearized_system(s);P();}
static void t6(void){T("simulate_linearized");
    LinearizedSystem*s=malloc(sizeof(LinearizedSystem));
    s->n_states=1;s->n_inputs=0;s->n_outputs=1;
    s->A=malloc(sizeof(double));s->A[0]=-2;s->B=calloc(1,sizeof(double));
    s->C=malloc(sizeof(double));s->C[0]=1;s->D=calloc(1,sizeof(double));
    s->eq.n_states=1;s->eq.n_inputs=0;s->eq.state=calloc(1,sizeof(double));s->eq.input=NULL;
    double x0[1]={1},to[11],yo[11];
    simulate_linearized(s,x0,NULL,0,2,100,ODE_RK4,11,to,yo,NULL);
    assert(fabs(yo[10]-exp(-4))<0.01);free_linearized_system(s);P();}
static void t7(void){T("eigenvalues_qr");
    double A[4]={0,1,-2,-3},rp[2],ip[2];int it=compute_eigenvalues_qr(A,2,rp,ip);
    assert(it>0);assert(rp[0]<0&&rp[1]<0);P();}
static void t8(void){T("df_ideal_relay");
    DescribingFunction df=df_ideal_relay(1,2);assert(fabs(cabs(df.gain)-2/M_PI)<1e-10);P();}
static void t9(void){T("df_saturation");
    DescribingFunction df=df_saturation(1,0.5,1);double r=0.5,ex=(2/M_PI)*(asin(r)+r*sqrt(1-r*r));
    assert(fabs(cabs(df.gain)-ex)<1e-10);P();}
static void t10(void){T("taylor_coeff");
    double c1=taylor_coefficient_1d(sin_wrap,0,1,1e-4,NULL);
    assert(fabs(c1-1)<0.01);P();}
static void t11(void){T("controllability");
    double xe[1]={0},ue[1]={0};LinearizedSystem*s=linearize_system(dyn_lin1,NULL,xe,ue,1,1,1,LINEARIZE_NUMERICAL_CENTRAL,1e-6,NULL);
    assert(s);int rk=0;double*Cm=compute_controllability_matrix(s,&rk);
    assert(Cm);assert(rk==1);free(Cm);free_linearized_system(s);P();}
static void t12(void){T("taylor_series");
    double s=taylor_series_sin(0.5,5);assert(fabs(s-sin(0.5))<0.001);P();}
static void t13(void){T("carleman_dim");
    assert(carleman_embedding_dimension(2,2)==6);assert(carleman_embedding_dimension(3,2)==12);P();}
static void t14(void){T("DMD");
    double X[3]={1,0.5,0.25},Y[3]={0.5,0.25,0.125};DMDResult r;memset(&r,0,sizeof(r));
    compute_dmd(X,Y,1,4,1,1,&r);assert(r.r==1);assert(fabs(creal(r.eigenvalues[0])-0.5)<0.3);
    free_dmd_result(&r);P();}
static void t15(void){T("gain_schedule");
    GainSchedule*gs=create_gain_schedule(3,2);assert(gs);assert(gs->n_entries==3);
    free_gain_schedule(gs);P();}

int main(void){
    printf("=== mini-linearization Unit Tests ===\n\n");
    t1();t2();t3();t4();t5();t6();t7();t8();t9();t10();t11();t12();t13();t14();t15();
    printf("\n=== %d/%d tests passed ===\n",tp,tr);
    return (tp==tr)?0:1;
}
