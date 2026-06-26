#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linearization_core.h"

static void pend_dyn(const double *x, const double *u, double *dx, int n, void *d) {
    (void)d; (void)n;
    double g = 9.81, l = 1.0, b = 0.1, m = 1.0;
    dx[0] = x[1];
    dx[1] = (g/l) * sin(x[0]) - b * x[1] + (u ? u[0] : 0.0) / (m * l * l);
}
static void pend_out(const double *x, const double *u, double *y, int p, void *d) {
    (void)u; (void)d; (void)p;
    y[0] = x[0];
}

int main(void) {
    printf("=== Inverted Pendulum Linearization ===\n\n");

    double x_eq[2] = {0.0, 0.0};
    double u_eq[1] = {0.0};

    printf("Linearizing pendulum about upright equilibrium...\n");
    LinearizedSystem *sys = linearize_system(
        pend_dyn, pend_out, x_eq, u_eq, 2, 1, 1,
        LINEARIZE_NUMERICAL_CENTRAL, 1e-6, NULL);
    printf("A matrix:\n");
    printf("  [%.4f  %.4f]\n", sys->A[0], sys->A[1]);
    printf("  [%.4f  %.4f]\n", sys->A[2], sys->A[3]);
    printf("B matrix: [%.4f, %.4f]^T\n", sys->B[0], sys->B[1]);

    LyapunovIndirectResult *lyap = apply_lyapunov_indirect(sys);
    printf("\nLyapunov Indirect Method:\n");
    printf("  Unstable modes: %d\n", lyap->n_unstable_modes);
    printf("  Conclusive: %s\n", lyap->is_conclusive ? "yes" : "no");

    int rank_ctrb;
    double *Cmat = compute_controllability_matrix(sys, &rank_ctrb);
    printf("  Controllability rank: %d (system is %s)\n",
           rank_ctrb, rank_ctrb == 2 ? "controllable" : "NOT controllable");
    free(Cmat);

    double domain = estimate_linearization_domain(pend_dyn, sys, 0.05, NULL);
    printf("  Valid perturbation radius: ~%.4f rad (%.1f deg)\n",
           domain, domain * 180.0 / M_PI);

    free_lyapunov_indirect_result(lyap);
    free_linearized_system(sys);
    printf("\nDone.\n");
    return 0;
}
