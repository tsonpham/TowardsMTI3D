//
//    SES3D_GRAD.C -- Frechet derivatives
//

#include <stdio.h>
#include "ses3d.h"

  void ses3d_grad(int it) {

#ifdef USE_GRAD_Q

    // These are the derivatives of the weights D with respect to alpha. 
    // They are hardcoded here, and may need to be changed.

      real dPdalpha[nrdiss];      
      dPdalpha[0] = -3.06;
      dPdalpha[1] = -1.54;
      dPdalpha[2] = 2.56;

	// Also, alpha and Q_kappa are hardcoded right now. 
	
	  real alpha = 0.3;
	  real Q_kappa = 57823.0;

#endif    // USE_GRAD_Q

    // Compute kernels by integrating on the fly.
      
      if (saving_vector[nt-(it+1)] > 0) {

        // read dynamic fields

        // read velocity fields

          ses3d_restore(vx_fw, it, tag_vx);
          ses3d_restore(vy_fw, it, tag_vy);
          ses3d_restore(vz_fw, it, tag_vz);

        // read strain rates

          ses3d_restore(exx_fw, it, tag_exx);
          ses3d_restore(eyy_fw, it, tag_eyy);
          ses3d_restore(ezz_fw, it, tag_ezz);

          ses3d_restore(exy_fw, it, tag_exy);
          ses3d_restore(exz_fw, it, tag_exz);
          ses3d_restore(eyz_fw, it, tag_eyz);

#ifdef DIAG_GRAD
          diag_restore_fw(it);
#endif

        // Preliminaries.

        // Adjoint strain tensor.

#ifdef USE_GPU
          compute_strain_tensor_GPU();
          
#else          
          FORALL_MD(IM, IE) {

              real exx_ad = dxux[IM][IE];
              real eyy_ad = dyuy[IM][IE];
              real ezz_ad = dzuz[IM][IE];

              real exy_ad = 0.5 * (dxuy[IM][IE] + dyux[IM][IE]);
              real exz_ad = 0.5 * (dxuz[IM][IE] + dzux[IM][IE]);
              real eyz_ad = 0.5 * (dyuz[IM][IE] + dzuy[IM][IE]);

            // Traces of forward and adjoint strain tensor.
		
              real tr_e_ad = exx_ad + eyy_ad + ezz_ad;
              real tr_e_fw = exx_fw[IM][IE] + eyy_fw[IM][IE] + ezz_fw[IM][IE];

            // Scalar product of forward and adjoint strain tensor.
		
              real e_ddot_e = 
                  exx_ad * exx_fw[IM][IE] + 
                  eyy_ad * eyy_fw[IM][IE] + 
                  ezz_ad * ezz_fw[IM][IE] + 
                  2.0 * exy_ad * exy_fw[IM][IE] + 
                  2.0 * exz_ad * exz_fw[IM][IE] + 
                  2.0 * eyz_ad * eyz_fw[IM][IE];

            // Frechet kernels for elastic parameters.
            // These are kernels for *absolute* perturbations.

              grad_cp[IM][IE] += 
                  samp_ad * dt * 
                      (2.0 * rho[IM][IE] * cp[IM][IE] * tr_e_fw * tr_e_ad) / Jac;

              grad_rho[IM][IE] += 
                  samp_ad * dt * (
                      vx_fw[IM][IE] * vx[IM][IE] + 
                      vy_fw[IM][IE] * vy[IM][IE] + 
                      vz_fw[IM][IE] * vz[IM][IE]) / Jac;
              grad_rho[IM][IE] += 
                  samp_ad * dt * (
                      (cp[IM][IE] * cp[IM][IE] - 2.0 * cs[IM][IE] * cs[IM][IE]) * 
                              tr_e_fw * tr_e_ad + 
                          2.0 * cs[IM][IE] * cs[IM][IE] * e_ddot_e) / Jac;

              grad_csh[IM][IE] += 
                  samp_ad * dt * (
                      4.0 * rho[IM][IE] * cs[IM][IE] * (
                          e_ddot_e - tr_e_ad * tr_e_fw - 
                          2.0 * exz_fw[IM][IE] * exz_ad - 
                          2.0 * eyz_fw[IM][IE] * eyz_ad)) / Jac;
              grad_csv[IM][IE] += 
                  samp_ad * dt * (
                      8.0 * rho[IM][IE] * cs[IM][IE] * (
                          exz_fw[IM][IE] * exz_ad + 
                          eyz_fw[IM][IE] * eyz_ad)) / Jac;

#ifdef USE_GRAD_Q

            // Frechet kernels for visco-elastic parameters.
            // These are kernels for *relative* perturbations.
            
              if (is_diss == 1) {

                  for (int k = 0; k < nrdiss; k++) {

                    // Shear Q

                      real Mtemp = 
                          Mxx[k][IM][IE] * exx_fw[IM][IE] + 
                          Myy[k][IM][IE] * eyy_fw[IM][IE] + 
                          Mzz[k][IM][IE] * ezz_fw[IM][IE] +
                          Mxy[k][IM][IE] * exy_fw[IM][IE] + 
                          Mxz[k][IM][IE] * exz_fw[IM][IE] + 
                          Myz[k][IM][IE] * eyz_fw[IM][IE] -
                          tr_e_fw * (
                              Mxx[k][IM][IE] +
                              Myy[k][IM][IE] +
                              Mzz[k][IM][IE]) / 3.0;

                      grad_Q_mu[IM][IE] += samp_ad * dt * Mtemp * tau_p[k] / Jac;
				      grad_alpha_mu[IM][IE] += 
				          samp_ad * dt * Mtemp * dPdalpha[k] * tau_p[k] / 
				              (D_p[k] * Jac);
				              
                    // Bulk Q

                      Mtemp  = 
                          (Mxx[k][IM][IE] + Myy[k][IM][IE] + Mzz[k][IM][IE]) * 
                              tr_e_fw;

                      grad_Q_kappa[IM][IE] += samp_ad * dt * Mtemp * tau_p[k] / Jac;
                      grad_alpha_kappa[IM][IE] +=
                          samp_ad * dt * Mtemp * dPdalpha[k] * tau_p[k] / 
                              (D_p[k] * Jac);

                      }     // for (int k ...)

                // Normalizations

                // ACHTUNG: The following code is located outside all
                //     loops and ifs in the original F90 version.
                //     Apparently, this was a bug.

                  if (it + 1 == nt) {

                      grad_Q_mu[IM][IE] = 
                          2.0 * mu[IM][IE] * grad_Q_mu[IM][IE] / QQ[IM][IE];
                      grad_alpha_mu[IM][IE] = -
                          2.0 * alpha * mu[IM][IE] * 
                              grad_alpha_mu[IM][IE] / QQ[IM][IE];
                      grad_Q_kappa[IM][IE] = 
                          kappa[IM][IE] * grad_Q_kappa[IM][IE] / Q_kappa;
#if 0    // ACHTUNG: Apparent bug in F90 version
                      grad_alpha_mu[IM][IE] = 
                          -alpha * kappa[IM][IE] * 
                              grad_alpha_kappa[IM][IE] / Q_kappa;
#else
                      grad_alpha_kappa[IM][IE] = 
                          -alpha * kappa[IM][IE] * 
                              grad_alpha_kappa[IM][IE] / Q_kappa;
#endif
          
                      }
                            
                  }     // if (is_diss ...)

#endif    //  USE_GRAD_Q
            
              }     // FORALL_MD (...)

#endif    // USE_GPU

          }     // if (saving_vector ...)
 
      }
