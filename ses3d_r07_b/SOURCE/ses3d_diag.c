//
//    SES3D_DIAG.C -- Diagnostics
//

#include <stdio.h>
#include <stdlib.h>
#include "ses3d.h"

  static void diag_fw(
          int it,
          char *tag,
          ARR_MD p_vx,
          ARR_MD p_vy,
          ARR_MD p_vz,
          ARR_MD p_exx,
          ARR_MD p_eyy,
          ARR_MD p_ezz,
          ARR_MD p_exy,
          ARR_MD p_exz,
          ARR_MD p_eyz) {

      real vx_min = MD_min_reduce(p_vx);
      real vx_max = MD_max_reduce(p_vx);
      real vy_min = MD_min_reduce(p_vy);
      real vy_max = MD_max_reduce(p_vy);
      real vz_min = MD_min_reduce(p_vz);
      real vz_max = MD_max_reduce(p_vz);
      real exx_min = MD_min_reduce(p_exx);
      real exx_max = MD_max_reduce(p_exx);
      real eyy_min = MD_min_reduce(p_eyy);
      real eyy_max = MD_max_reduce(p_eyy);
      real ezz_min = MD_min_reduce(p_ezz);
      real ezz_max = MD_max_reduce(p_ezz);
      real exy_min = MD_min_reduce(p_exy);
      real exy_max = MD_max_reduce(p_exy);
      real exz_min = MD_min_reduce(p_exz);
      real exz_max = MD_max_reduce(p_exz);
      real eyz_min = MD_min_reduce(p_eyz);
      real eyz_max = MD_max_reduce(p_eyz);

      SERIAL
          WRITELN("Iteration %d: %s", it+1, tag);
          WRITELN("    vx_fw:  min %g, max %g", vx_min, vx_max);
          WRITELN("    vy_fw:  min %g, max %g", vy_min, vy_max);   
          WRITELN("    vz_fw:  min %g, max %g", vz_min, vz_max);   
          WRITELN("    exx_fw: min %g, max %g", exx_min, exx_max);   
          WRITELN("    eyy_fw: min %g, max %g", eyy_min, eyy_max);
          WRITELN("    ezz_fw: min %g, max %g", ezz_min, ezz_max);   
          WRITELN("    exy_fw: min %g, max %g", exy_min, exy_max);   
          WRITELN("    exz_fw: min %g, max %g", exz_min, exz_max);   
          WRITELN("    eyz_fw: min %g, max %g", eyz_min, eyz_max);
      ENDSERIAL
      }

  void diag_store_fw(int it) {
      diag_fw(it, "store FW", 
          vx, vy, vz, 
          dxux, dyuy, dzuz, 
          tmp_md1, tmp_md2, tmp_md3);
      }
      
  void diag_restore_fw(int it) {
      diag_fw(it, "restore FW",
          vx_fw, vy_fw, vz_fw,
          exx_fw, eyy_fw, ezz_fw,
          exy_fw, exz_fw, eyz_fw);          
      }

  void diag_grad(int it) {

      real grad_rho_min = MD_min_reduce(grad_rho);
      real grad_rho_max = MD_max_reduce(grad_rho);
      real grad_cp_min = MD_min_reduce(grad_cp);
      real grad_cp_max = MD_max_reduce(grad_cp);
      real grad_csh_min = MD_min_reduce(grad_csh);
      real grad_csh_max = MD_max_reduce(grad_csh);
      real grad_csv_min = MD_min_reduce(grad_csv);
      real grad_csv_max = MD_max_reduce(grad_csv);
#ifdef USE_GRAD_Q
      real grad_Q_mu_min = MD_min_reduce(grad_Q_mu);
      real grad_Q_mu_max = MD_max_reduce(grad_Q_mu);
      real grad_alpha_mu_min = MD_min_reduce(grad_alpha_mu);
      real grad_alpha_mu_max = MD_max_reduce(grad_alpha_mu);
      real grad_Q_kappa_min = MD_min_reduce(grad_Q_kappa);
      real grad_Q_kappa_max = MD_max_reduce(grad_Q_kappa);
      real grad_alpha_kappa_min = MD_min_reduce(grad_alpha_kappa);
      real grad_alpha_kappa_max = MD_max_reduce(grad_alpha_kappa);
#endif

      SERIAL
          WRITELN("Iteration %d: gradients", it+1);
          WRITELN("    grad_rho:  min %g, max %g", 
              grad_rho_min, grad_rho_max);
          WRITELN("    grad_cp:   min %g, max %g", 
              grad_cp_min, grad_cp_max);
          WRITELN("    grad_csh:  min %g, max %g", 
              grad_csh_min, grad_csh_max);
          WRITELN("    grad_csv:  min %g, max %g", 
              grad_csv_min, grad_csv_max);
#ifdef USE_GRAD_Q
          WRITELN("    grad_Q_mu:         min %g, max %g", 
              grad_Q_mu_min, grad_Q_mu_max);
          WRITELN("    grad_alpha_mu:     min %g, max %g", 
              grad_alpha_mu_min, grad_alpha_mu_max);
          WRITELN("    grad_Q_kappa:      min %g, max %g", 
              grad_Q_kappa_min, grad_Q_kappa_max);
          WRITELN("    grad_alpha_kappa:  min %g, max %g", 
              grad_alpha_kappa_min, grad_alpha_kappa_max);
#endif
      ENDSERIAL
      }
