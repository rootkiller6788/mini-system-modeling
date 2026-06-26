#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linearization_core.h"
#include "describing_function.h"

/* Magnetic levitation: dx1/dt = x2, dx2/dt = g - k*x3^2/(m*x1^2), dx3/dt = (u-R*x3)/L */
static void maglev_dyn(const double *x, const double *u, double *dx, int n, void *d) {
    (void)d; (void)n;
    double g = 9.81, k = 1e-4, m = 0.05, R = 1.0, L = 0.01;
    dx[0] = x[1];
    dx[1] = g - k * x[2] * x[2] / (m * x[0] * x[0]);
    dx[2] = (u[0] - R * x[2]) / L;
}
static void maglev_out(const double *x, const double *u, double *y, int p, void *d) {
    (void)u; (void)d; (void)p;
    y[0] = x[0]; /* output = ball position */
}

int main(void) {
    printf("=== Magnetic Levitation Linearization ===\n\n");

    /* Equilibrium: force balance => x3_eq = sqrt(m*g/k) * x1_eq */
    double x_eq[3] = {0.01, 0.0, 0.7};
    double u_eq[1] = {0.7}; /* Voltage to maintain equilibrium */

    printf("Linearizing at gap=%.1fmm, current=%.2fA\n",
           x_eq[0]*1000, x_eq[2]);
    LinearizedSystem *sys = linearize_system(
        maglev_dyn, maglev_out, x_eq, u_eq, 3, 1, 1,
        LINEARIZE_NUMERICAL_CENTRAL, 1e-6, NULL);
    printf("A matrix (3x3):\n");
    for (int i = 0; i < 3; i++) {
        printf("  ");
        for (int j = 0; j < 3; j++)
            printf("%8.2f ", sys->A[i*3+j]);
        printf("\n");
    }

    LyapunovIndirectResult *lyap = apply_lyapunov_indirect(sys);
    printf("Unstable modes: %d (maglev is inherently unstable)\n",
           lyap->n_unstable_modes);

    /* Stability analysis via describing function for relay control */
    printf("\nDescribing function analysis with relay (M=12V):\n");
    double A_test = 0.002; /* 2mm oscillation */
    DescribingFunction df = df_ideal_relay(12.0, A_test);
    printf("  N(A=%.3f) = %.3f + %.3fj\n",
           A_test, creal(df.gain), cimag(df.gain));
    printf("  Magnitude: %.1f dB\n", df.magnitude_db);
    printf("  Equivalent linear gain: %.4f\n", cabs(df.gain));

    free_lyapunov_indirect_result(lyap);
    free_linearized_system(sys);
    printf("\nDone.\n");
    return 0;
}
