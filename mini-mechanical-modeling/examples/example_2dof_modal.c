/**
 * @example example_2dof_modal.c
 * @brief 2-DOF system modal analysis --- natural frequencies, mode shapes, FRF
 *
 * Demonstrates: Generalized eigenvalue problem, mode orthogonality,
 * frequency response functions, modal superposition.
 * Reference: Ewins, Modal Testing (2000), Ch.1-2.
 */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "multi_dof_systems.h"
#include "mechanical_elements.h"

int main(void)
{
    printf("\n=== 2-DOF Modal Analysis ===\n\n");

    /* System: m1=10kg, m2=5kg, k1=500N/m, k2=300N/m, kc=100N/m */
    double m1 = 10.0, m2 = 5.0;
    double k1 = 500.0, k2 = 300.0, kc = 100.0;
    double b1 = 2.0, b2 = 1.0, bc = 0.5;

    MDOFSystem sys;
    mdof_build_2dof(&sys, m1, m2, k1, k2, kc, b1, b2, bc);

    printf("System Description:\n");
    printf("  Masses: m1=%.1f kg, m2=%.1f kg\n", m1, m2);
    printf("  Stiffness: k1=%.0f N/m (wall-m1), k2=%.0f N/m (m2-wall), kc=%.0f N/m (coupling)\n",
           k1, k2, kc);
    printf("  Damping: b1=%.1f, b2=%.1f, bc=%.1f N*s/m\n\n", b1, b2, bc);

    printf("Mass Matrix M:\n");
    printf("  [[%6.1f, %6.1f]\n", sys.M[0], sys.M[1]);
    printf("   [%6.1f, %6.1f]]\n\n", sys.M[2], sys.M[3]);

    printf("Stiffness Matrix K:\n");
    printf("  [[%7.0f, %7.0f]\n", sys.K[0], sys.K[1]);
    printf("   [%7.0f, %7.0f]]\n\n", sys.K[2], sys.K[3]);

    /* Eigenvalue solution */
    double evals[MDOF_MAX_DOF], evecs[MDOF_MAX_DOF*MDOF_MAX_DOF];
    int n_modes = mdof_eigen_solve(&sys, evals, evecs);

    printf("Eigenvalue Analysis:\n");
    printf("  Found %d modes\n\n", n_modes);

    for (int i = 0; i < n_modes; i++) {
        double wn_i = (evals[i] > 0) ? sqrt(evals[i]) : 0.0;
        double fn_i = wn_i / (2.0 * M_PI);
        printf("  Mode %d: omega_n = %.3f rad/s, f_n = %.3f Hz\n", i+1, wn_i, fn_i);
        printf("    Mode shape: [%.4f, %.4f]^T\n", evecs[i*2], evecs[i*2+1]);

        /* Modal parameters */
        ModalParams mp;
        mdof_modal_params(&sys, &evecs[i*2], i, &mp);
        printf("    Modal mass = %.4f, Modal stiffness = %.2f\n", mp.modal_mass, mp.modal_stiffness);
        printf("    Participation factor = %.4f, Effective mass = %.4f kg\n\n",
               mp.participation_factor, mp.effective_modal_mass);
    }

    /* Orthogonality check */
    double ortho = mdof_mass_ortho(&sys, evecs, &evecs[2]);
    printf("Mass Orthogonality phi_1^T*M*phi_2 = %.2e (should be ~0)\n\n", ortho);

    /* Modal Assurance Criterion */
    double mac = mdof_mac(evecs, &evecs[2], 2);
    printf("MAC(phi_1, phi_2) = %.6f (should be ~0 for distinct modes)\n\n", mac);

    double mac_self = mdof_mac(evecs, evecs, 2);
    printf("MAC(phi_1, phi_1) = %.6f (should be 1 for same mode)\n\n", mac_self);

    /* Frequency Response Function */
    ModalParams modes[2];
    mdof_modal_params(&sys, evecs, 0, &modes[0]);
    mdof_modal_params(&sys, &evecs[2], 1, &modes[1]);

    printf("Frequency Response H11 (displacement at m1 due to force at m1):\n");
    printf("  Freq [Hz]    |H11| [m/N]   Phase [deg]\n");
    double freqs[] = {0.5, 1.0, 1.5, 2.0, 3.0, 5.0};
    for (int i = 0; i < 6; i++) {
        double w = 2.0 * M_PI * freqs[i];
        double complex H = mdof_receptance_frf(&sys, modes, 2, 0, 0, w);
        double mag = cabs(H);
        double phase = atan2(cimag(H), creal(H)) * 180.0 / M_PI;
        printf("  %8.2f    %10.6f    %8.1f\n", freqs[i], mag, phase);
    }

    /* Rayleigh damping design */
    double alpha, beta;
    rayleigh_damping_coeffs(modes[0].natural_freq, modes[1].natural_freq,
                             0.02, 0.03, &alpha, &beta);
    printf("\nRayleigh Damping Coefficients:\n");
    printf("  alpha = %.6f, beta = %.6f\n", alpha, beta);
    printf("  zeta_1 = %.4f (target=0.02), zeta_2 = %.4f (target=0.03)\n",
           rayleigh_modal_damping(alpha, beta, modes[0].natural_freq),
           rayleigh_modal_damping(alpha, beta, modes[1].natural_freq));

    printf("\n=== Modal Analysis Complete ===\n");
    return 0;
}
