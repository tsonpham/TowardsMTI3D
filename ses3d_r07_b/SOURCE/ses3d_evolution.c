//
//    SES3D_EVOLUTION.C -- Main computations
//

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ses3d.h"

  static void make_single_force(int it) {

    // single point source, forward calculation

      if (so_idx_vector >= 0) {

          int IM = so_idx_vector;

          for (int IE = 0; IE < ED_MAX; IE++) {
              real delta = so_delta[IE];

            // add single force to stress divergence

              if (source_type == 1)
                  delta_sx[IM][IE] = -delta * so[it];
              else if (source_type == 2)
                  delta_sy[IM][IE] = -delta * so[it];
              else if (source_type == 3)
                  delta_sz[IM][IE] = -delta * so[it];
              }
          }
          
      else if (so_idx_tensor >= 0) {

          int IM = so_idx_tensor;

          for (int IE = 0; IE < ED_MAX; IE++) {
              real delta = so_delta[IE];

            // make source field

              src_xx[IM][IE] = -MOM_xx * delta * so[it];
              src_yy[IM][IE] = -MOM_yy * delta * so[it];
              src_zz[IM][IE] = -MOM_zz * delta * so[it];

              src_xy[IM][IE] = -MOM_xy * delta * so[it];
              src_yx[IM][IE] = -MOM_xy * delta * so[it];

              src_xz[IM][IE] = -MOM_xz * delta * so[it];
              src_zx[IM][IE] = -MOM_xz * delta * so[it];

              src_yz[IM][IE] = -MOM_yz * delta * so[it];
              src_zy[IM][IE] = -MOM_yz * delta * so[it];
              }
          }

    // adjoint point sources, adjoint calculation

      else if (adjoint_flag == 2)  {

        // loop through all the adjoint sources in this processor box

          for (int idx = 0; idx < nr_adsrc; idx++) {
              int IM = ad_stf_idx[idx];
              if (IM >= 0) {
                  for (int IE = 0; IE < ED_MAX; IE++) {
                      delta_sx[IM][IE] = 0.0;
                      delta_sy[IM][IE] = 0.0;
                      delta_sz[IM][IE] = 0.0;
                      }
                  }
              }
                  
          for (int idx = 0; idx < nr_adsrc; idx++) {
              int IM = ad_stf_idx[idx];

              if (IM >= 0) {

                  for (int IE = 0; IE < ED_MAX; IE++) {
                      real delta = ad_stf_delta[idx][IE];

                    // add single force to stress divergence

                      delta_sx[IM][IE] -= delta * ad_stf_x[idx][it];
                      delta_sy[IM][IE] -= delta * ad_stf_y[idx][it];
                      delta_sz[IM][IE] -= delta * ad_stf_z[idx][it];
                      }
                  }
              }
          }
  
      }

  static void communicate_global_acceleration(void) {    

      communicate_global(COMM_SX);
      communicate_global(COMM_SY);
      communicate_global(COMM_SZ);
      }
      
  void ses3d_evolution(int it) {

      make_single_force(it);

      FORALL_MD (IM, IE) {

          real exx = 0.0;    // (grad v)_(theta theta)
          real exy = 0.0;    // (grad v)_(theta phi)
          real exz = 0.0;    // (grad v)_(theta r)
          real eyy = 0.0;    // (grad v)_(phi phi)
          real eyx = 0.0;    // (grad v)_(phi theta)
          real eyz = 0.0;    // (grad v)_(phi r)
          real ezz = 0.0;    // (grad v)_(r r)
          real ezx = 0.0;    // (grad v)_(r theta)
          real ezy = 0.0;    // (grad v)_(r phi)

        // derivative proxies of velocity field
 
 #ifndef TEMP_OPT_EVOLUTION    // TODO: Revise this
          int i = ED_index1(IE);
          int j = ED_index2(IE);
          int k = ED_index3(IE);
      
          for (int q = 0; q <= lpd; q++) {
              real t1 = 2 * dl[q][i] / dx;
              real t2 = 2 * dl[q][j] / dy;
              real t3 = 2 * dl[q][k] / dz;
                
              int q1 = ED_offset(q, j, k);
              int q2 = ED_offset(i, q, k);
              int q3 = ED_offset(i, j, q); 

              exx += vx[IM][q1] * t1;
              exy += vy[IM][q1] * t1;
              exz += vz[IM][q1] * t1;

              eyy += vy[IM][q2] * t2;
              eyx += vx[IM][q2] * t2;
              eyz += vz[IM][q2] * t2;

              ezz += vz[IM][q3] * t3;
              ezx += vx[IM][q3] * t3;
              ezy += vy[IM][q3] * t3;
              }
#else

          for (int q = 0; q <= lpd; q++) {
              int q1 = offx[q][IE];
              int q2 = offy[q][IE];
              int q3 = offz[q][IE]; 

              real t1 = c1x[q][IE];
              real t2 = c1y[q][IE];
              real t3 = c1z[q][IE];
                
              exx += vx[IM][q1] * t1;
              exy += vy[IM][q1] * t1;
              exz += vz[IM][q1] * t1;

              eyy += vy[IM][q2] * t2;
              eyx += vx[IM][q2] * t2;
              eyz += vz[IM][q2] * t2;

              ezz += vz[IM][q3] * t3;
              ezx += vx[IM][q3] * t3;
              ezy += vy[IM][q3] * t3;
              }
#endif

        // complete velocity gradients in spherical coordinates

          real tr = 1.0 / r[IM][IE];
          real tcot = cot_theta[IM][IE];
          real tsin = 1.0 / sin_theta[IM][IE];

          exx = (vz[IM][IE] + exx) * tr;
          exy = exy * tr;
          exz = (exz - vx[IM][IE]) * tr;

          eyy = (eyy * tsin + vz[IM][IE] + vx[IM][IE] * tcot) * tr;
          eyx = (-vy[IM][IE] * tcot + eyx * tsin) * tr;
          eyz = (-vy[IM][IE] + eyz * tsin) * tr;

          ezz = -ezz;
          ezx = -ezx;
          ezy = -ezy;

        // integrate velocity gradients to displacement gradients

          dzuz[IM][IE] += dt * ezz;
          dyuy[IM][IE] += dt * eyy;
          dxux[IM][IE] += dt * exx;
          dzux[IM][IE] += dt * ezx;
          dxuz[IM][IE] += dt * exz;
          dzuy[IM][IE] += dt * ezy;
          dyuz[IM][IE] += dt * eyz;
          dxuy[IM][IE] += dt * exy;
          dyux[IM][IE] += dt * eyx;

        // build apml strain rates

          exx += (prof_z[IM][IE] + prof_y[IM][IE]) * dxux[IM][IE] * ispml;
          eyy += (prof_x[IM][IE] + prof_z[IM][IE]) * dyuy[IM][IE] * ispml;
          ezz += (prof_x[IM][IE] + prof_y[IM][IE]) * dzuz[IM][IE] * ispml;

          exy += (prof_y[IM][IE] + prof_z[IM][IE]) * dxuy[IM][IE] * ispml;
          eyx += (prof_x[IM][IE] + prof_z[IM][IE]) * dyux[IM][IE] * ispml;

          exz += (prof_y[IM][IE] + prof_z[IM][IE]) * dxuz[IM][IE] * ispml;
          ezx += (prof_x[IM][IE] + prof_y[IM][IE]) * dzux[IM][IE] * ispml;

          eyz += (prof_x[IM][IE] + prof_z[IM][IE]) * dyuz[IM][IE] * ispml;
          ezy += (prof_x[IM][IE] + prof_y[IM][IE]) * dzuy[IM][IE] * ispml;

          exy = (exy + eyx) / 2;
          exz = (exz + ezx) / 2;
          eyz = (eyz + ezy) / 2;

        // make strong form stress rates

          real sxx;
          real syy;
          real szz;
          real sxy;
          real sxz;
          real syz;

        // diagonal stress rates
          if (is_diss) {
            // dissipation on
              real Mdx = 0.0;
              real Mdy = 0.0;
              real Mdz = 0.0;

              for (int d = 0; d < nrdiss; d++) {
                  Mdx += Mxx[d][IM][IE];
                  Mdy += Myy[d][IM][IE];
                  Mdz += Mzz[d][IM][IE];
                  }

              real Mtemp = Mdx + Mdy + Mdz;

              real L = 
                  (kappa[IM][IE] - 2 * mu_tau[IM][IE] / 3) * (exx + eyy + ezz) - 
                      2 * mu[IM][IE] * tau[IM][IE] * Mtemp / 3;
              sxx = 
                  L + 2 * mu_tau[IM][IE] * exx + 
                      2 * mu[IM][IE] * tau[IM][IE] * Mdx + 
                      C[IM][IE] * ezz + 
                      A[IM][IE] * (eyy + exx);
              syy = 
                  L + 2 * mu_tau[IM][IE] * eyy + 
                      2 * mu[IM][IE] * tau[IM][IE] * Mdy + 
                      C[IM][IE] * ezz + 
                      A[IM][IE] * (eyy + exx);
              szz = 
                  L + 2 * mu_tau[IM][IE] * ezz + 
                      2 * mu[IM][IE] * tau[IM][IE] * Mdz + 
                      C[IM][IE] * (eyy + exx);
              }
          else {
            // dissipation off
              real L = lambda[IM][IE] * (exx + eyy + ezz);
              sxx = L + 2 * mu[IM][IE] * exx + C[IM][IE] * ezz + A[IM][IE] * (eyy + exx);
              syy = L + 2 * mu[IM][IE] * eyy + C[IM][IE] * ezz + A[IM][IE] * (eyy + exx);
              szz = L + 2 * mu[IM][IE] * ezz + C[IM][IE] * (eyy + exx);
              }
                
        // off-diagonal stress rates
          if (is_diss) {
            // dissipation on
              real Mdx = 0.0;
              real Mdy = 0.0;
              real Mdz = 0.0;

              for (int d = 0; d < nrdiss; d++) {
                  Mdx += Mxy[d][IM][IE];
                  Mdy += Mxz[d][IM][IE];
                  Mdz += Myz[d][IM][IE];
                  }

              sxy = 2 * mu_tau[IM][IE] * exy + 2 * mu[IM][IE] * tau[IM][IE] * Mdx;
              sxz = 2 * mu_tau[IM][IE] * exz + 2 * mu[IM][IE] * tau[IM][IE] * Mdy + 2 * B[IM][IE] * exz;
              syz = 2 * mu_tau[IM][IE] * eyz + 2 * mu[IM][IE] * tau[IM][IE] * Mdz + 2 * B[IM][IE] * eyz;
              }
          else {
            // dissipation off
              sxy = 2 * mu[IM][IE] * exy;
              sxz = 2 * mu[IM][IE] * exz + 2 * B[IM][IE] * exz;
              syz = 2 * mu[IM][IE] * eyz + 2 * B[IM][IE] * eyz;
              }

        // time extrapolation of memory variables

          if (is_diss) {
              for (int d = 0; d < nrdiss; d++) {
                  real L = dt * D_p[d] / tau_p[d];
                  Mxx[d][IM][IE] -= dt * Mxx[d][IM][IE] / tau_p[d] + L * exx;
                  Myy[d][IM][IE] -= dt * Myy[d][IM][IE] / tau_p[d] + L * eyy;
                  Mzz[d][IM][IE] -= dt * Mzz[d][IM][IE] / tau_p[d] + L * ezz;
                  Mxy[d][IM][IE] -= dt * Mxy[d][IM][IE] / tau_p[d] + L * exy;
                  Mxz[d][IM][IE] -= dt * Mxz[d][IM][IE] / tau_p[d] + L * exz;
                  Myz[d][IM][IE] -= dt * Myz[d][IM][IE] / tau_p[d] + L * eyz;
                  }
              }

        // march pml stresses (integrate stress rates)

          sxx_pml[IM][IE] += dt * (sxx - ispml * prof_x[IM][IE] * sxx_pml[IM][IE]);
          syy_pml[IM][IE] += dt * (syy - ispml * prof_y[IM][IE] * syy_pml[IM][IE]);
          szz_pml[IM][IE] += dt * (szz - ispml * prof_z[IM][IE] * szz_pml[IM][IE]);

          sxy_pml[IM][IE] += dt * (sxy - ispml * prof_x[IM][IE] * sxy_pml[IM][IE]);
          syx_pml[IM][IE] += dt * (sxy - ispml * prof_y[IM][IE] * syx_pml[IM][IE]);

          sxz_pml[IM][IE] += dt * (sxz - ispml * prof_x[IM][IE] * sxz_pml[IM][IE]);
          szx_pml[IM][IE] += dt * (sxz - ispml * prof_z[IM][IE] * szx_pml[IM][IE]);

          syz_pml[IM][IE] += dt * (syz - ispml * prof_y[IM][IE] * syz_pml[IM][IE]);
          szy_pml[IM][IE] += dt * (syz - ispml * prof_z[IM][IE] * szy_pml[IM][IE]);
          }

#if 0    // TODO: Revise this
      FORALL_MD (IM, IE) {

        // compute stress divergence (with moment tensor added to the stress tensor)
          real tsx = 0.0; 
          real tsy = 0.0; 
          real tsz = 0.0;

#ifndef TEMP_OPT_EVOLUTION    // TODO: Revise this
          int i = ED_index1(IE);
          int j = ED_index2(IE);
          int k = ED_index3(IE);

          real t1 = (2 / dx) * Jac * w[j] * w[k];
          real t2 = (2 / dy) * Jac * w[i] * w[k];
          real t3 = -(2 / dz) * Jac * w[i] * w[j];

          for (int q = 0; q <= lpd; q++) {
              int q1 = ED_offset(q, j, k);
              int q2 = ED_offset(i, q, k);
              int q3 = ED_offset(i, j, q); 

              real tx = w[q] * dl[i][q] * r[IM][q1] * sin_theta[IM][q1] * t1;
              real ty = w[q] * dl[j][q] * r[IM][q2] * t2;
              real tz = w[q] * dl[k][q] * r[IM][q3] * r[IM][q3] * sin_theta[IM][q3] * t3;

              tsx +=
                  (szx_pml[IM][q3] + src_zx[IM][q3]) * tz +
                  (syx_pml[IM][q2] + src_yx[IM][q2]) * ty +
                  (sxx_pml[IM][q1] + src_xx[IM][q1]) * tx;
              tsy +=
                  (szy_pml[IM][q3] + src_zy[IM][q3]) * tz +
                  (syy_pml[IM][q2] + src_yy[IM][q2]) * ty +
                  (sxy_pml[IM][q1] + src_xy[IM][q1]) * tx;
              tsz +=
                  (szz_pml[IM][q3] + src_zz[IM][q3]) * tz +
                  (syz_pml[IM][q2] + src_yz[IM][q2]) * ty +
                  (sxz_pml[IM][q1] + src_xz[IM][q1]) * tx;
              }
#else

          for (int q = 0; q <= lpd; q++) {
              int q1 = offx[q][IE];
              int q2 = offy[q][IE];
              int q3 = offz[q][IE]; 

              real tx = c2x[q][IE] * r[IM][q1] * sin_theta[IM][q1];
              real ty = c2y[q][IE] * r[IM][q2];
              real tz = c2z[q][IE] * r[IM][q3] * r[IM][q3] * sin_theta[IM][q3];

              tsx +=
                  (szx_pml[IM][q3] + src_zx[IM][q3]) * tz +
                  (syx_pml[IM][q2] + src_yx[IM][q2]) * ty +
                  (sxx_pml[IM][q1] + src_xx[IM][q1]) * tx;
              tsy +=
                  (szy_pml[IM][q3] + src_zy[IM][q3]) * tz +
                  (syy_pml[IM][q2] + src_yy[IM][q2]) * ty +
                  (sxy_pml[IM][q1] + src_xy[IM][q1]) * tx;
              tsz +=
                  (szz_pml[IM][q3] + src_zz[IM][q3]) * tz +
                  (syz_pml[IM][q2] + src_yz[IM][q2]) * ty +
                  (sxz_pml[IM][q1] + src_xz[IM][q1]) * tx;
              }
#endif

        // add external forces single force

          tsx += delta_sx[IM][IE];
          tsy += delta_sy[IM][IE];
          tsz += delta_sz[IM][IE];

          sx[IM][IE] = tsx;
          sy[IM][IE] = tsy;
          sz[IM][IE] = tsz;
          }

#else
      for (int IM = 0; IM < MD_length; IM++) {
          static real sh_txx[ED_MAX];
          static real sh_tyx[ED_MAX];
          static real sh_tzx[ED_MAX];
          static real sh_txy[ED_MAX];
          static real sh_tyy[ED_MAX];
          static real sh_tzy[ED_MAX];
          static real sh_txz[ED_MAX];
          static real sh_tyz[ED_MAX];
          static real sh_tzz[ED_MAX];
          
          for (int IE = 0; IE < ED_MAX; IE++) {
              real tx = r[IM][IE] * sin_theta[IM][IE];
              real ty = r[IM][IE];
              real tz = r[IM][IE] * r[IM][IE] * sin_theta[IM][IE];
          
              sh_txx[IE] = (sxx_pml[IM][IE] + src_xx[IM][IE]) * tx;
              sh_tyx[IE] = (syx_pml[IM][IE] + src_yx[IM][IE]) * ty;
              sh_tzx[IE] = (szx_pml[IM][IE] + src_zx[IM][IE]) * tz;
                  
              sh_txy[IE] = (sxy_pml[IM][IE] + src_xy[IM][IE]) * tx;
              sh_tyy[IE] = (syy_pml[IM][IE] + src_yy[IM][IE]) * ty;
              sh_tzy[IE] = (szy_pml[IM][IE] + src_zy[IM][IE]) * tz;

              sh_txz[IE] = (sxz_pml[IM][IE] + src_xz[IM][IE]) * tx;
              sh_tyz[IE] = (syz_pml[IM][IE] + src_yz[IM][IE]) * ty;
              sh_tzz[IE] = (szz_pml[IM][IE] + src_zz[IM][IE]) * tz;
              }

          for (int IE = 0; IE < ED_MAX; IE++) {
              real tsx = 0.0; 
              real tsy = 0.0; 
              real tsz = 0.0;

              int i = ED_index1(IE);
              int j = ED_index2(IE);
              int k = ED_index3(IE);

              for (int q = 0; q <= lpd; q++) {
                  int q1 = ED_offset(q, j, k);
                  int q2 = ED_offset(i, q, k);
                  int q3 = ED_offset(i, j, q); 

                  real tx = c2x[q][IE];
                  real ty = c2y[q][IE];
                  real tz = c2z[q][IE];

                  tsx += sh_txx[q1] * tx + sh_tyx[q2] * ty + sh_tzx[q3] * tz;
                  tsy += sh_txy[q1] * tx + sh_tyy[q2] * ty + sh_tzy[q3] * tz;
                  tsz += sh_txz[q1] * tx + sh_tyz[q2] * ty + sh_tzz[q3] * tz;
                  }

            // add external forces single force

              tsx += delta_sx[IM][IE];
              tsy += delta_sy[IM][IE];
              tsz += delta_sz[IM][IE];

              sx[IM][IE] = tsx;
              sy[IM][IE] = tsy;
              sz[IM][IE] = tsz;
              }
                    
          }

#endif    // ... TODO: Revise this

    // map local force vectors to global force vectors

      communicate_global_acceleration();

      FORALL_MD (IM, IE) {

        // make accelerations / divide weak stresses by mass matrix

          sx[IM][IE] /= MM[IM][IE];
          sy[IM][IE] /= MM[IM][IE];
          sz[IM][IE] /= MM[IM][IE];

        // time extrapolation of the displacement field

          vx[IM][IE] += dt * (sx[IM][IE] - ispml * prof[IM][IE] * vx[IM][IE]);
          vy[IM][IE] += dt * (sy[IM][IE] - ispml * prof[IM][IE] * vy[IM][IE]);
          vz[IM][IE] += dt * (sz[IM][IE] - ispml * prof[IM][IE] * vz[IM][IE]);

        // taper boundaries when pmls are turned off

          vx[IM][IE] *= taper[IM][IE];
          vy[IM][IE] *= taper[IM][IE];
          vz[IM][IE] *= taper[IM][IE];

          sxx_pml[IM][IE] *= taper[IM][IE]; 
          syy_pml[IM][IE] *= taper[IM][IE]; 
          szz_pml[IM][IE] *= taper[IM][IE];
          sxy_pml[IM][IE] *= taper[IM][IE]; 
          syx_pml[IM][IE] *= taper[IM][IE];
          szx_pml[IM][IE] *= taper[IM][IE]; 
          sxz_pml[IM][IE] *= taper[IM][IE];
          syz_pml[IM][IE] *= taper[IM][IE]; 
          szy_pml[IM][IE] *= taper[IM][IE];
          }

      }

  void record_seismograms(int it) {

      for (int idx = 0; idx < nr; idx++) {
          int IM = rec_idx[idx];
          if (IM >= 0) {
              real tx = 0.0;
              real ty = 0.0;
              real tz = 0.0;              

              for (int IE = 0; IE < ED_MAX; IE++) {
                  real d = rec_delta[idx][IE];
                  tx += vx[IM][IE] * d;
                  ty += vy[IM][IE] * d;
                  tz += vz[IM][IE] * d;
                  }

              seismogram_x[idx][it] = tx;
              seismogram_y[idx][it] = ty;
              seismogram_z[idx][it] = tz;              
              }  
          }
          
      }
