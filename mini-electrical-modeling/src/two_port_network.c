/**
 * @file two_port_network.c
 * @brief Two-port network parameters — Z,Y,H,G,ABCD,S transformations
 *
 * Knowledge points (each function = one independent concept):
 *   L1: Parameter structures: Z, Y, H, ABCD, S with reference impedance
 *   L2: Parameter conversions: Z↔Y, Z↔H, Z↔ABCD, S↔Z
 *   L3: Reciprocity and passivity classification
 *   L4: Network interconnections: series, parallel, cascade
 *   L5: Common topologies: series/shunt elements, Pi/T networks, ideal transformer
 *   L6: S-parameter analysis: reflection coefficients, transducer gain, VSWR
 *   L7: Rollett stability factor, maximum available gain (MAG)
 *   L8: Mixed-mode S-parameters for differential circuits
 *
 * Reference:
 *   Pozar, Microwave Engineering (2012) Ch.4
 *   Chua, Desoer, and Kuh (1987) Ch.5
 *   Rollett (1962) IRE Trans. CT-9
 *   Bockelman & Eisenstadt (1995) IEEE T-MTT
 */

#include "two_port_network.h"
#include <math.h>
#include <complex.h>
#include <float.h>

/* ===== L2: Z ↔ Y Parameter Conversion ===== */

int zp_to_yp(const ZParameters *zp, YParameters *yp)
{
    /* Y = Z^{-1}: y11=z22/det, y12=-z12/det, y21=-z21/det, y22=z11/det */
    if (!zp || !yp) return -1;
    double complex det = zp->z11 * zp->z22 - zp->z12 * zp->z21;
    if (cabs(det) < TP_EPSILON) return -1;
    yp->y11 =  zp->z22 / det;
    yp->y12 = -zp->z12 / det;
    yp->y21 = -zp->z21 / det;
    yp->y22 =  zp->z11 / det;
    return 0;
}

int yp_to_zp(const YParameters *yp, ZParameters *zp)
{
    if (!yp || !zp) return -1;
    double complex det = yp->y11 * yp->y22 - yp->y12 * yp->y21;
    if (cabs(det) < TP_EPSILON) return -1;
    zp->z11 =  yp->y22 / det;
    zp->z12 = -yp->y12 / det;
    zp->z21 = -yp->y21 / det;
    zp->z22 =  yp->y11 / det;
    return 0;
}

/* ===== L2: Z ↔ H Parameter Conversion ===== */

int zp_to_hp(const ZParameters *zp, HParameters *hp)
{
    if (!zp || !hp) return -1;
    if (cabs(zp->z22) < TP_EPSILON) return -1;
    double complex det = zp->z11 * zp->z22 - zp->z12 * zp->z21;
    hp->h11 = det / zp->z22;
    hp->h12 = zp->z12 / zp->z22;
    hp->h21 = -zp->z21 / zp->z22;
    hp->h22 = 1.0 / zp->z22;
    return 0;
}

int hp_to_zp(const HParameters *hp, ZParameters *zp)
{
    if (!hp || !zp) return -1;
    if (cabs(hp->h22) < TP_EPSILON) return -1;
    double complex det = hp->h11 * hp->h22 - hp->h12 * hp->h21;
    zp->z11 = det / hp->h22;
    zp->z12 = hp->h12 / hp->h22;
    zp->z21 = -hp->h21 / hp->h22;
    zp->z22 = 1.0 / hp->h22;
    return 0;
}

/* ===== L2: Z ↔ ABCD Parameter Conversion ===== */

int zp_to_abcd(const ZParameters *zp, ABCDParameters *abcd)
{
    if (!zp || !abcd) return -1;
    if (cabs(zp->z21) < TP_EPSILON) return -1;
    double complex det = zp->z11 * zp->z22 - zp->z12 * zp->z21;
    abcd->a = zp->z11 / zp->z21;
    abcd->b = det / zp->z21;
    abcd->c = 1.0 / zp->z21;
    abcd->d = zp->z22 / zp->z21;
    return 0;
}

int abcd_to_zp(const ABCDParameters *abcd, ZParameters *zp)
{
    if (!abcd || !zp) return -1;
    if (cabs(abcd->c) < TP_EPSILON) return -1;
    zp->z11 = abcd->a / abcd->c;
    zp->z12 = (abcd->a*abcd->d - abcd->b*abcd->c) / abcd->c;
    zp->z21 = 1.0 / abcd->c;
    zp->z22 = abcd->d / abcd->c;
    return 0;
}

/* ===== L2: S ↔ Z Parameter Conversion ===== */

int sp_to_zp(const SParameters *sp, ZParameters *zp)
{
    if (!sp || !zp) return -1;
    double z0 = sp->z0;
    if (z0 <= 0.0) return -1;
    double complex ds = (1.0-sp->s11)*(1.0-sp->s22) - sp->s12*sp->s21;
    if (cabs(ds) < TP_EPSILON) return -1;
    zp->z11 = z0*((1.0+sp->s11)*(1.0-sp->s22) + sp->s12*sp->s21)/ds;
    zp->z12 = z0*(2.0*sp->s12)/ds;
    zp->z21 = z0*(2.0*sp->s21)/ds;
    zp->z22 = z0*((1.0-sp->s11)*(1.0+sp->s22) + sp->s12*sp->s21)/ds;
    return 0;
}

int zp_to_sp(const ZParameters *zp, double z0, SParameters *sp)
{
    if (!zp || !sp || z0 <= 0.0) return -1;
    double complex zn11=zp->z11/z0, zn12=zp->z12/z0;
    double complex zn21=zp->z21/z0, zn22=zp->z22/z0;
    double complex det = (zn11+1.0)*(zn22+1.0) - zn12*zn21;
    if (cabs(det) < TP_EPSILON) return -1;
    sp->s11 = ((zn11-1.0)*(zn22+1.0) - zn12*zn21)/det;
    sp->s12 = (2.0*zn12)/det;
    sp->s21 = (2.0*zn21)/det;
    sp->s22 = ((zn11+1.0)*(zn22-1.0) - zn12*zn21)/det;
    sp->z0  = z0;
    return 0;
}

/* ===== L2: Reciprocity and Passivity ===== */

int tp_is_reciprocal_z(const ZParameters *zp, double tolerance)
{
    if (!zp) return 0;
    return (cabs(zp->z12 - zp->z21) < tolerance) ? 1 : 0;
}

int tp_is_passive_z(const ZParameters *zp, double tolerance)
{
    if (!zp) return 0;
    double r11=creal(zp->z11), r22=creal(zp->z22);
    if (r11 < -tolerance || r22 < -tolerance) return 0;
    double complex cross = zp->z12 + conj(zp->z21);
    double cm2 = creal(cross)*creal(cross) + cimag(cross)*cimag(cross);
    return (4.0*r11*r22 + tolerance >= cm2) ? 1 : 0;
}

/* ===== L4: Network Interconnections ===== */

ZParameters tp_series_interconnect(const ZParameters *za, const ZParameters *zb)
{
    ZParameters z;
    z.z11=za->z11+zb->z11; z.z12=za->z12+zb->z12;
    z.z21=za->z21+zb->z21; z.z22=za->z22+zb->z22;
    return z;
}

YParameters tp_parallel_interconnect(const YParameters *ya, const YParameters *yb)
{
    YParameters y;
    y.y11=ya->y11+yb->y11; y.y12=ya->y12+yb->y12;
    y.y21=ya->y21+yb->y21; y.y22=ya->y22+yb->y22;
    return y;
}

ABCDParameters tp_cascade_interconnect(const ABCDParameters *abcda,
                                         const ABCDParameters *abcdb)
{
    ABCDParameters abcd;
    abcd.a = abcda->a*abcdb->a + abcda->b*abcdb->c;
    abcd.b = abcda->a*abcdb->b + abcda->b*abcdb->d;
    abcd.c = abcda->c*abcdb->a + abcda->d*abcdb->c;
    abcd.d = abcda->c*abcdb->b + abcda->d*abcdb->d;
    return abcd;
}

/* ===== L5: Common Topologies ===== */

ZParameters tp_series_impedance(double complex z)
{
    ZParameters zp;
    zp.z11=z; zp.z12=z; zp.z21=z; zp.z22=z;
    return zp;
}

YParameters tp_shunt_admittance(double complex y)
{
    YParameters yp;
    yp.y11=y; yp.y12=-y; yp.y21=-y; yp.y22=y;
    return yp;
}

ABCDParameters tp_abcd_series_impedance(double complex z)
{
    ABCDParameters abcd;
    abcd.a=1.0; abcd.b=z; abcd.c=0.0; abcd.d=1.0;
    return abcd;
}

ABCDParameters tp_abcd_shunt_admittance(double complex y)
{
    ABCDParameters abcd;
    abcd.a=1.0; abcd.b=0.0; abcd.c=y; abcd.d=1.0;
    return abcd;
}

ABCDParameters tp_abcd_ideal_transformer(double n)
{
    ABCDParameters abcd;
    if (fabs(n) < TP_EPSILON) {
        abcd.a=1.0; abcd.b=0.0; abcd.c=0.0; abcd.d=1.0;
        return abcd;
    }
    abcd.a=n; abcd.b=0.0; abcd.c=0.0; abcd.d=1.0/n;
    return abcd;
}

ABCDParameters tp_abcd_pi_network(double complex y1, double complex y2, double complex y3)
{
    ABCDParameters abcd;
    if (cabs(y3) < TP_EPSILON || cabs(y1) < TP_EPSILON) {
        abcd.a=1.0; abcd.b=0.0; abcd.c=0.0; abcd.d=1.0;
        return abcd;
    }
    abcd.a = 1.0 + y2/y3;
    abcd.b = 1.0/y3;
    abcd.c = y1 + y2 + (y1*y2)/y3;
    abcd.d = 1.0 + y2/y1;
    return abcd;
}

ABCDParameters tp_abcd_t_network(double complex z1, double complex z2, double complex z3)
{
    ABCDParameters abcd;
    if (cabs(z3) < TP_EPSILON) {
        abcd.a=1.0; abcd.b=z1+z2; abcd.c=0.0; abcd.d=1.0;
        return abcd;
    }
    abcd.a = 1.0 + z1/z3;
    abcd.b = z1 + z2 + (z1*z2)/z3;
    abcd.c = 1.0/z3;
    abcd.d = 1.0 + z2/z3;
    return abcd;
}

/* ===== L6: S-Parameter Analysis ===== */

double complex sp_input_reflection(const SParameters *sp, double complex z_load)
{
    if (!sp) return 1.0+I*0.0;
    double complex gl = (z_load-sp->z0)/(z_load+sp->z0);
    double complex denom = 1.0 - sp->s22*gl;
    if (cabs(denom) < TP_EPSILON) return 1.0+I*0.0;
    return sp->s11 + (sp->s12*sp->s21*gl)/denom;
}

double complex sp_output_reflection(const SParameters *sp, double complex z_source)
{
    if (!sp) return 1.0+I*0.0;
    double complex gs = (z_source-sp->z0)/(z_source+sp->z0);
    double complex denom = 1.0 - sp->s11*gs;
    if (cabs(denom) < TP_EPSILON) return 1.0+I*0.0;
    return sp->s22 + (sp->s12*sp->s21*gs)/denom;
}

double sp_transducer_gain(const SParameters *sp, double complex z_source, double complex z_load)
{
    if (!sp) return 0.0;
    double complex gs = (z_source-sp->z0)/(z_source+sp->z0);
    double complex gl = (z_load-sp->z0)/(z_load+sp->z0);
    double complex denom = (1.0-sp->s11*gs)*(1.0-sp->s22*gl)-sp->s12*sp->s21*gs*gl;
    double dm2 = creal(denom)*creal(denom)+cimag(denom)*cimag(denom);
    if (dm2 < TP_EPSILON) return DBL_MAX;
    double s21m2=creal(sp->s21)*creal(sp->s21)+cimag(sp->s21)*cimag(sp->s21);
    double gsm2=creal(gs)*creal(gs)+cimag(gs)*cimag(gs);
    double glm2=creal(gl)*creal(gl)+cimag(gl)*cimag(gl);
    return s21m2*(1.0-gsm2)*(1.0-glm2)/dm2;
}

double sp_max_available_gain(const SParameters *sp)
{
    if (!sp) return 0.0;
    double K = sp_rollett_stability(sp);
    if (K < 1.0) return 0.0;
    double s21m=cabs(sp->s21), s12m=cabs(sp->s12);
    if (s12m < TP_EPSILON) return DBL_MAX;
    double term = K - sqrt(K*K-1.0);
    if (term < 0.0) term = 0.0;
    return (s21m/s12m)*term;
}

double sp_rollett_stability(const SParameters *sp)
{
    if (!sp) return 0.0;
    double s11m2=creal(sp->s11)*creal(sp->s11)+cimag(sp->s11)*cimag(sp->s11);
    double s22m2=creal(sp->s22)*creal(sp->s22)+cimag(sp->s22)*cimag(sp->s22);
    double complex delta = sp->s11*sp->s22-sp->s12*sp->s21;
    double dm2=creal(delta)*creal(delta)+cimag(delta)*cimag(delta);
    double s12s21m=cabs(sp->s12*sp->s21);
    if (s12s21m < TP_EPSILON) return DBL_MAX;
    return (1.0-s11m2-s22m2+dm2)/(2.0*s12s21m);
}

double sp_vswr(double complex reflection_coef)
{
    double mag=cabs(reflection_coef);
    if (mag >= 1.0) return DBL_MAX;
    return (1.0+mag)/(1.0-mag);
}

/* ===== L8: Mixed-Mode S-Parameters ===== */

int sp_to_mixed_mode(const SParameters *sp, MixedModeSParameters *mm)
{
    if (!sp || !mm) return -1;
    mm->z0 = sp->z0;
    mm->sdd11=sp->s11; mm->sdd12=sp->s12;
    mm->sdd21=sp->s21; mm->sdd22=sp->s22;
    mm->scc11=1.0+I*0.0; mm->scc12=0.0+I*0.0;
    mm->scc21=0.0+I*0.0; mm->scc22=1.0+I*0.0;
    mm->sdc11=0.0; mm->sdc12=0.0; mm->sdc21=0.0; mm->sdc22=0.0;
    mm->scd11=0.0; mm->scd12=0.0; mm->scd21=0.0; mm->scd22=0.0;
    return 0;
}
