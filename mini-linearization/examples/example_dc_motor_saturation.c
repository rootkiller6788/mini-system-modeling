#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linearization_core.h"
#include "describing_function.h"

/* DC motor with saturation: di/dt = (V - R*i - Ke*w)/L, dw/dt = (Kt*Tsat(i) - B*w)/J */
static void motor_dyn(const double *x, const double *u, double *dx, int n, void *d) {
    (void)d; (void)n;
    double R = 0.5, L = 0.001, Ke = 0.01, Kt = 0.01, J = 0.001, B = 0.0001;
    double i_max = 10.0;
    double i = x[0], w = x[1];
    double i_sat = (i > i_max) ? i_max : ((i < -i_max) ? -i_max : i);
    dx[0] = (u[0] - R * i - Ke * w) / L;
    dx[1] = (Kt * i_sat - B * w) / J;
}
static void motor_out(const double *x, const double *u, double *y, int p, void *d) {
    (void)u; (void)d; (void)p;
    y[0] = x[1]; /* output = angular velocity */
}

int main(void) {
    printf("=== DC Motor Saturation Analysis ===\n\n");
    printf("Toyota Prius MG2 motor parameters:\n");
    printf("  R=0.5ohm, L=1mH, Ke=0.01Vs/rad, Kt=0.01Nm/A\n");
    printf("  J=0.001kg.m2, B=1e-4 Nms/rad, Imax=10A\n\n");

    double x_eq[2] = {0.0, 300.0}; /* 300 rad/s ~ 2865 RPM */
    double u_eq[1] = {4.5}; /* ~4.5V to maintain speed */

    printf("Linearizing at %.0f rad/s (%.0f RPM)...\n",
           x_eq[1], x_eq[1] * 60.0 / (2.0 * M_PI));

    LinearizedSystem *sys = linearize_system(
        motor_dyn, motor_out, x_eq, u_eq, 2, 1, 1,
        LINEARIZE_NUMERICAL_CENTRAL, 1e-6, NULL);

    printf("A = [[%.2f, %.2f], [%.2f, %.2f]]\n",
           sys->A[0], sys->A[1], sys->A[2], sys->A[3]);
    printf("B = [%.2f, %.2f]^T\n", sys->B[0], sys->B[1]);

    LyapunovIndirectResult *lyap = apply_lyapunov_indirect(sys);
    printf("System is %sstable\n", lyap->is_stable ? "" : "un");

    /* Describing function for saturation at large signal */
    printf("\nLarge-signal saturation analysis:\n");
    for (int k = 0; k < 4; k++) {
        double A = 2.0 + k * 4.0;
        DescribingFunction df = df_saturation(1.0, 5.0, A);
        printf("  Input amplitude=%.0fA: N(A)=%.3f (%.1fdB)\n",
               A, fabs(creal(df.gain)), df.magnitude_db);
    }

    free_lyapunov_indirect_result(lyap);
    free_linearized_system(sys);
    printf("\nDone.\n");
    return 0;
}
