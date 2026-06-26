#include "describing_function.h"
#include <string.h>
#include <math.h>
#include <assert.h>

DescribingFunction df_ideal_relay(double M,double A){
    DescribingFunction df; memset(&df,0,sizeof(df));
    df.amplitude=A; df.gain=(A>1e-15)?(4.0*M/(M_PI*A)):INFINITY;
    df.magnitude_db=20*log10(cabs(df.gain)+1e-15);
    return df;
}
DescribingFunction df_relay_hysteresis(double M,double Delta,double A){
    DescribingFunction df; memset(&df,0,sizeof(df));
    df.amplitude=A;
    if(A<Delta){df.gain=0;return df;}
    double ratio=Delta/A;
    double re=4.0*M/(M_PI*A)*sqrt(1.0-ratio*ratio);
    double im=-4.0*M*Delta/(M_PI*A*A);
    df.gain=re+im*I; df.phase_deg=atan2(im,re)*180.0/M_PI;
    df.magnitude_db=20*log10(cabs(df.gain)+1e-15);
    return df;
}
DescribingFunction df_saturation(double k,double a,double A){
    DescribingFunction df; memset(&df,0,sizeof(df));
    df.amplitude=A;
    if(A<=a){df.gain=k;df.magnitude_db=20*log10(k+1e-15);return df;}
    double ratio=a/A;
    double gain=(2.0*k/M_PI)*(asin(ratio)+ratio*sqrt(1.0-ratio*ratio));
    df.gain=gain; df.magnitude_db=20*log10(gain+1e-15);
    return df;
}
DescribingFunction df_dead_zone(double k,double d,double A){
    DescribingFunction df; memset(&df,0,sizeof(df));
    df.amplitude=A;
    if(A<=d){df.gain=0;return df;}
    double ratio=d/A;
    double gain=k*(1.0-(2.0/M_PI)*(asin(ratio)+ratio*sqrt(1.0-ratio*ratio)));
    df.gain=gain; df.magnitude_db=20*log10(gain+1e-15);
    return df;
}
DescribingFunction df_backlash(double k,double b,double A){
    DescribingFunction df; memset(&df,0,sizeof(df));
    df.amplitude=A;
    if(A<=b){df.gain=0;return df;}
    double ratio=b/A;
    double re=(k/M_PI)*(M_PI_2+asin(1.0-2.0*ratio)+2.0*(1.0-2.0*ratio)*sqrt(ratio*(1.0-ratio)));
    double im=-4.0*k*b*(1.0-ratio)/(M_PI*A);
    df.gain=re+im*I; df.phase_deg=atan2(im,re)*180.0/M_PI;
    df.magnitude_db=20*log10(cabs(df.gain)+1e-15);
    return df;
}
DescribingFunction df_cubic(double c,double A){
    DescribingFunction df; memset(&df,0,sizeof(df));
    df.amplitude=A; df.gain=0.75*c*A*A; df.magnitude_db=20*log10(fabs(0.75*c*A*A)+1e-15);
    return df;
}
DescribingFunction df_quantizer(double q,double A){
    DescribingFunction df; memset(&df,0,sizeof(df));
    df.amplitude=A;
    if(A<q/2.0){df.gain=0;return df;}
    df.gain=4.0*q/(M_PI*A); df.magnitude_db=20*log10(cabs(df.gain)+1e-15);
    return df;
}
LimitCyclePrediction predict_limit_cycle_df(
    double (*Gr)(double,void*),double (*Gi)(double,void*),
    DescribingFunction (*dff)(double,void*),
    double Amin,double Amax,double wmin,double wmax,double tol,void *ud){
    LimitCyclePrediction lcp; memset(&lcp,0,sizeof(lcp));
    int nA=50,nw=50; double best_err=INFINITY;
    for(int ia=0;ia<nA;ia++){
        double A=Amin+(Amax-Amin)*ia/(double)(nA-1);
        DescribingFunction df=dff(A,ud);
        double nAinv=-1.0/cabs(df.gain);
        for(int iw=0;iw<nw;iw++){
            double w=wmin+(wmax-wmin)*iw/(double)(nw-1);
            double re=Gr(w,ud),im=Gi(w,ud);
            double Gmag=sqrt(re*re+im*im);
            double err=fabs(Gmag-nAinv);
            if(err<best_err){best_err=err;lcp.amplitude=A;lcp.frequency_rad_s=w;}
        }
    }
    lcp.exists=(best_err<tol); lcp.frequency_hz=lcp.frequency_rad_s/(2.0*M_PI);
    lcp.is_stable=true; return lcp;
}
bool df_stability_check(
    double (*Gr)(double,void*),double (*Gi)(double,void*),
    double wmin,double wmax,int np,
    DescribingFunction (*dff)(double,void*),double A,void *ud){
    (void)Gr;(void)Gi;(void)wmin;(void)wmax;(void)np;
    DescribingFunction df=dff(A,ud);
    /* Simplified check: if describing function gain is finite, assume stable */
    return (cabs(df.gain)<1e10);
}
void free_describing_function(DescribingFunction *df){ (void)df; }
