//
//    PROJECT.CPP -- Engine for projecting SES3D fields
//

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <assert.h>
#include <string>
#include "project.h"

//
//    Local macros
//

#define FOR2(I1, LO1, HI1, I2, LO2, HI2) \
    for (int I1 = (LO1); I1 <= (HI1); I1++) \
    for (int I2 = (LO2); I2 <= (HI2); I2++) 

#define FOR3(I1, LO1, HI1, I2, LO2, HI2, I3, LO3, HI3) \
    for (int I1 = (LO1); I1 <= (HI1); I1++) \
    for (int I2 = (LO2); I2 <= (HI2); I2++) \
    for (int I3 = (LO3); I3 <= (HI3); I3++) 

#define DEL_ARR(A) \
    if ((A) != NULL) { \
        delete[] (A); \
        (A) = NULL; \
        } 

#define DEL_RES(A) \
    DeleteResult(A); \
    (A) = NULL;

//
//    Local constants
//

  static const real PI = (real)3.1415926535898;

//
//    ProjectBase
//

// construction/destruction

  ProjectBase::ProjectBase() {
      rootPath = "../..";
      InitKnots();
      pp = 0;
      px = 0;
      py = 0;
      pz = 0;
      box = NULL;
      nsubvol = 0;
      subvol = NULL;
      x = NULL;
      y = NULL;
      z = NULL;
      bixMin = NULL;
      bixMax = NULL;
      biyMin = NULL;
      biyMax = NULL;
      bizMin = NULL;
      bizMax = NULL;
      bintx = NULL;
      binty = NULL;
      bintz = NULL;
      }
      
  ProjectBase::~ProjectBase() { 
      DeleteBox();
      DeleteBlocks();
      DeleteGrid();
      }
      
// interface

  void ProjectBase::SetRootPath(const char *path) {
      rootPath = path;
      }

// implementation: protected

  bool ProjectBase::Compute() {
      bool succ = TryCompute();
      DeleteBox();
      DeleteBlocks();
      DeleteResults();
      return succ;
      }

  bool ProjectBase::TryCompute() {
      if (!ReadBox())
          return false;
      if (!ReadBlocks())
          return false;
      InitResults();
      for (int ind = 0; ind < pp; ind++) {
          if (!ComputePe(ind))
              return false;
          }
      if (!WriteResults())
          return false;
      return true;
      }

  bool ProjectBase::ReadBox() {

      string fn = rootPath + "/MODELS/MODELS/boxfile";
      FILE *fp = fopen(fn.c_str(), "r");
      if (fp == NULL) {
          Log("Cannot open %s", fn.c_str());
          return false;
          }
      for (int i = 0; i < 14; i++)
          ReadLn(fp);
      pp = ReadInt(fp);      
      ReadLn(fp);
      px = ReadInt(fp);      
      ReadLn(fp);
      py = ReadInt(fp);      
      ReadLn(fp);
      pz = ReadInt(fp);      
      ReadLn(fp);
      ReadLn(fp);

      box = new Box[pp];

      for (int i = 0; i < pp; i++) {
          int rank = ReadInt(fp);
          ReadLn(fp);
          ReadLn(fp);           // ix/iy/iz_multi_loc
          int nxMin = ReadInt(fp);
          int nxMax = ReadInt(fp);
          ReadLn(fp);
          int nyMin = ReadInt(fp);
          int nyMax = ReadInt(fp);
          ReadLn(fp);
          int nzMin = ReadInt(fp);
          int nzMax = ReadInt(fp);
          ReadLn(fp);
          real xmin = ReadReal(fp);
          real xmax = ReadReal(fp);
          ReadLn(fp);
          real ymin = ReadReal(fp);
          real ymax = ReadReal(fp);
          ReadLn(fp);
          real zmin = ReadReal(fp);
          real zmax = ReadReal(fp);
          ReadLn(fp);
          ReadLn(fp);
          box[i].rank = rank;      
          box[i].nx = nxMax - nxMin;
          box[i].ny = nyMax - nyMin;
          box[i].nz = nzMax - nzMin;
          box[i].xmin = xmin;
          box[i].xmax = xmax;
          box[i].ymin = ymin;
          box[i].ymax = ymax;
          box[i].zmin = zmin;
          box[i].zmax = zmax;
          }

      fclose(fp);
      return true;
      }
      
  void ProjectBase::DeleteBox() {
      DEL_ARR(box);
      }

  bool ProjectBase::ReadBlocks() {

      string fn = rootPath + "/MODELS/MODELS_3D/block_x";
      FILE *fp = fopen(fn.c_str(), "r");
      if (fp == NULL) {
          Log("Cannot open %s", fn.c_str());
          return false;
          }
      nsubvol = ReadInt(fp);
      subvol = new Subvol[nsubvol];
      for (int i = 0; i < nsubvol; i++) {
          subvol[i].nbx = 0;
          subvol[i].nby = 0;
          subvol[i].nbz = 0;
          subvol[i].bxco = NULL;
          subvol[i].byco = NULL;
          subvol[i].bzco = NULL;
          }
      for (int i = 0; i < nsubvol; i++) {
          int nbx = ReadInt(fp);
          real *bxco = new real[nbx];
          subvol[i].nbx = nbx;
          subvol[i].bxco = bxco;
          for (int k = 0; k < nbx; k++)
              bxco[k] = ReadReal(fp);
          }
      fclose(fp);

      fn = rootPath + "/MODELS/MODELS_3D/block_y";
      fp = fopen(fn.c_str(), "r");
      if (fp == NULL) {
          Log("Cannot open %s", fn.c_str());
          return false;
          }
      int n = ReadInt(fp);
      if (n != nsubvol) {
          Log("Mismatched number of sub-volumes: %d <> %d", n, nsubvol);
          fclose(fp);
          return false;
          }
      for (int i = 0; i < nsubvol; i++) {
          int nby = ReadInt(fp);
          real *byco = new real[nby];
          subvol[i].nby = nby;
          subvol[i].byco = byco;
          for (int k = 0; k < nby; k++)
              byco[k] = ReadReal(fp);
          }
      fclose(fp);

      fn = rootPath + "/MODELS/MODELS_3D/block_z";
      fp = fopen(fn.c_str(), "r");
      if (fp == NULL) {
          Log("Cannot open %s", fn.c_str());
          return false;
          }
      n = ReadInt(fp);
      if (n != nsubvol) {
          Log("Mismatched number of sub-volumes: %d <> %d", n, nsubvol);
          fclose(fp);
          return false;
          }
      for (int i = 0; i < nsubvol; i++) {
          int nbz = ReadInt(fp);
          real *bzco = new real[nbz];
          subvol[i].nbz = nbz;
          subvol[i].bzco = bzco;
          for (int k = 0; k < nbz; k++)
              bzco[k] = ReadReal(fp);
          }
      fclose(fp);

      for (int i = 0; i < nsubvol; i++) {
          int nbx = subvol[i].nbx;
          real *bxco = subvol[i].bxco;
          subvol[i].dbx = bxco[1] - bxco[0];
          for (int k = 0; k < nbx; k++)
              bxco[k] *= PI / (real)180.0;
          int nby = subvol[i].nby;
          real *byco = subvol[i].byco;
          subvol[i].dby = byco[1] - byco[0];
          for (int k = 0; k < nby; k++)
              byco[k] *= PI / (real)180.0;
          int nbz = subvol[i].nbz;
          real *bzco = subvol[i].bzco;
          subvol[i].dbz = bzco[1] - bzco[0];
          for (int k = 0; k < nbz; k++)
              bzco[k] *= (real)1000.0;
          }

      return true;
      }

  void ProjectBase::DeleteBlocks() {
      if (subvol == NULL)
          return;
      for (int i = 0; i < nsubvol; i++) {
          if (subvol[i].bxco != NULL)
              delete[] subvol[i].bxco;
          if (subvol[i].byco != NULL)
              delete[] subvol[i].byco;
          if (subvol[i].bzco != NULL)
              delete[] subvol[i].bzco;      
          }
      delete[] subvol;
      subvol = NULL;
      }

  bool ProjectBase::ComputePe(int ind) {
      bool succ = TryComputePe(ind);
      DeleteFields();
      DeleteGrid();
      return succ;
      }

  bool ProjectBase::TryComputePe(int ind) {
      if (!ReadFields(ind))
          return false;
      for (int isubvol = 0; isubvol < nsubvol; isubvol++)
          Integrate(ind, isubvol);
      return true;
      }

  void ProjectBase::Integrate(int ind, int isubvol) {

      int nx = box[ind].nx;
      int ny = box[ind].ny;
      int nz = box[ind].nz;
      real xmin = box[ind].xmin;
      real xmax = box[ind].xmax;
      real ymin = box[ind].ymin;
      real ymax = box[ind].ymax;
      real zmin = box[ind].zmin;
      real zmax = box[ind].zmax;
      real dx = (xmax - xmin) / (nx + 1);
      real dy = (ymax - ymin) / (ny + 1);
      real dz = (zmax - zmin) / (nz + 1);

    // make coordinates

      x = new real[nx+1][LPD+1];
      y = new real[ny+1][LPD+1];
      z = new real[nz+1][LPD+1];
      
      FOR2 (i, 0, nx, n, 0, LPD)
          x[i][n] = xmin + i * dx + 
              (real)0.5 * ((real)1.0 + knots[n]) * dx;
      FOR2 (i, 0, ny, n, 0, LPD)
          y[i][n] = ymin + i * dy + 
              (real)0.5 * ((real)1.0 + knots[n]) * dy;
      FOR2 (i, 0, nz, n, 0, LPD)
          z[i][n] = zmax - i * dz - 
              (real)0.5 * ((real)1.0 + knots[n]) * dz;

    // find indices of integration limits

      int nbx = subvol[isubvol].nbx;
      int nby = subvol[isubvol].nby;
      int nbz = subvol[isubvol].nbz;
      real *bxco = subvol[isubvol].bxco;
      real *byco = subvol[isubvol].byco;
      real *bzco = subvol[isubvol].bzco;

      bixMin = new int[nbx-1];
      bixMax = new int[nbx-1];
      biyMin = new int[nby-1];
      biyMax = new int[nby-1];
      bizMin = new int[nbz-1];
      bizMax = new int[nbz-1];

      bintx = new real[(nbx-1)*(nx+1)][LPD+1];
      binty = new real[(nby-1)*(ny+1)][LPD+1];
      bintz = new real[(nbz-1)*(nz+1)][LPD+1];

    // x-direction

      for (int indx = 0; indx < nbx - 1; indx++) {
          if (bxco[indx] > xmax || bxco[indx+1] < xmin) {
              bixMin[indx] = -1;
              bixMax[indx] = -1;
              continue;
              }
          int ixMin = 0;
          for (int i = 0; i <= nx; i++) {
              if (x[i][0] <= bxco[indx])
                  ixMin = i;
              }
          int ixMax = nx;
          for (int i = nx; i >= 0; i--) {
              if (bxco[indx+1] <= x[i][LPD])
                  ixMax = i;
              }
        // integration limits
          for (int i = ixMin; i <= ixMax; i++) {
              real xlima;
              if (bxco[indx] > x[i][0])
                  xlima = (2 * bxco[indx] - x[i][LPD] - x[i][0]) / dx;
              else
                  xlima = -1.0;
              real xlimb;
              if (bxco[indx+1] < x[i][LPD])
                  xlimb = (2 * bxco[indx+1] - x[i][LPD] - x[i][0]) / dx;
              else
                  xlimb = 1.0;
              for (int l = 0; l <= LPD; l++)
                  bintx[indx*(nx+1)+i][l] = IntLag(l, xlima, xlimb);
              }
          bixMin[indx] = ixMin;
          bixMax[indx] = ixMax;
          }

    // y-direction
    
      for (int indy = 0; indy < nby - 1; indy++) {
          if (byco[indy] > ymax || byco[indy+1] < ymin) {
              biyMin[indy] = -1;
              biyMax[indy] = -1;
              continue;
              }
          int iyMin = 0;
          for (int j = 0; j <= ny; j++) {
              if (y[j][0] <= byco[indy])
                  iyMin = j;
              }
          int iyMax = ny;
          for (int j = ny; j >= 0; j--) {
              if (byco[indy+1] <= y[j][LPD])
                  iyMax = j;
              }
        // integration limits
          for (int j = iyMin; j <= iyMax; j++) {
              real ylima;
              if (byco[indy] > y[j][0])
                  ylima = (2 * byco[indy] - y[j][LPD] - y[j][0]) / dy;
              else
                  ylima = -1.0;
              real ylimb;
              if (byco[indy+1] < y[j][LPD])
                  ylimb = (2 * byco[indy+1] - y[j][LPD] - y[j][0]) / dy;
              else
                  ylimb = 1.0;
              for (int m = 0; m <= LPD; m++)
                  binty[indy*(ny+1)+j][m] = IntLag(m, ylima, ylimb);
              }
          biyMin[indy] = iyMin;
          biyMax[indy] = iyMax;
          }
     
    // z-direction

      for (int indz = 0; indz < nbz - 1; indz++) {
          if (bzco[indz] > zmax || bzco[indz+1] < zmin) {
              bizMin[indz] = -1;
              bizMax[indz] = -1;
              continue;
              }
          int izMin = 0;
          for (int k = 0; k <= nz; k++) {
              if (z[k][0] >= bzco[indz+1])
                  izMin = k;
              }
          int izMax = nz;
          for (int k = nz; k >= 0; k--) {
              if (bzco[indz] >= z[k][LPD])
                  izMax = k;
              }
        // integration limits
          for (int k = izMin; k <= izMax; k++) {
              real zlimb;
              if (bzco[indz] > z[k][LPD])
                  zlimb = (z[k][0] + z[k][LPD] - 2 * bzco[indz]) / dz;
              else
                  zlimb = 1.0;
              real zlima;
              if (bzco[indz+1] < z[k][0])
                  zlima = (z[k][0] + z[k][LPD] - 2 * bzco[indz+1]) / dz;
              else
                  zlima = -1.0;
              for (int n = 0; n <= LPD; n++)
                  bintz[indz*(nz+1)+k][n] = IntLag(n, zlima, zlimb);
              }
          bizMin[indz] = izMin;
          bizMax[indz] = izMax;
          }

    // perform integration

      real dbx = subvol[isubvol].dbx;
      real dby = subvol[isubvol].dby;
      real dbz = subvol[isubvol].dbz;
      real dbn = (real)1.0 / (dbx * dby * dbz); 

      FOR3 (indx, 0, nbx - 2, indy, 0, nby - 2, indz, 0, nbz - 2) {
          int ixMin = bixMin[indx];
          int iyMin = biyMin[indy];
          int izMin = bizMin[indz];
          if (ixMin < 0 || iyMin < 0 || izMin < 0)
              continue;
          int ixMax = bixMax[indx];
          int iyMax = biyMax[indy];
          int izMax = bizMax[indz];
          FOR3 (i, ixMin, ixMax, j, iyMin, iyMax, k, izMin, izMax)
              IntegrateCell(isubvol, indx, indy, indz, i, j, k, 
                  nbx, nby, nbz, nx, ny, nz, dbn);
          }
      }

  void ProjectBase::DeleteGrid() {
      DEL_ARR(x);
      DEL_ARR(y);
      DEL_ARR(z);
      DEL_ARR(bixMin);
      DEL_ARR(bixMax);
      DEL_ARR(biyMin);
      DEL_ARR(biyMax);
      DEL_ARR(bizMin);
      DEL_ARR(bizMax);
      DEL_ARR(bintx);
      DEL_ARR(binty);
      DEL_ARR(bintz);
      }

  real ProjectBase::IntLag(int idx, real limLeft, real limRight) {
      real z0[LPD];
      real a[LPD+1];

    // make vector of negative roots
      int j = 0;
      for (int i = 0; i <= LPD; i++) {
          if (i != idx) {
              z0[j] = -knots[i];
              j++;
              }
          }

    // compute polynomial coefficients
      if (LPD == 1) {
          a[1] = (real)1.0;
          a[0] = z0[0];
          }
      else if (LPD == 2) {
          a[2] = (real)1.0;
          a[1] = z0[0] + z0[1];
          a[0] = z0[0] * z0[1];
          }
      else if (LPD > 2) {
        // initialisation of the iteration at degree 2
          a[2] = (real)1.0;
          a[1] = z0[0] + z0[1];
          a[0] = z0[0] * z0[1];
        // successively compute coefficients for higher polynomials
          for (int i = 3; i <= LPD; i++) {
            // compute the remaining coefficients
              for (int k = i - 2; k >= 1; k--)
                  a[k] = a[k] * z0[i-1] + a[k-1];
            // set the trivial coefficients for degree i >= 3
              a[i] = (real)1.0;
              a[i-1] = (real)0.0;
              a[0] = (real)1.0;
              for (int k = 0; k < i; k++) {
                  a[i-1] += z0[k];
                  a[0] *= z0[k];
                  }
              }
          }

    // compute constant factor
      real norm = (real)1.0;
      for (int i = 0; i <= LPD; i++) {
          if (i != idx)
              norm *= knots[idx] - knots[i];
          }

    // evaluate the integral
      real intLag = (real)0.0;
      real cr = limRight;
      real cl = limLeft;
      for (int i = 0; i <= LPD; i++) {
          intLag += (cr - cl) * a[i] / (i + 1);
          cr *= limRight;
          cl *= limLeft;
          }

      intLag /= norm;
      return intLag;
      }
      
  real **ProjectBase::NewResult() {
      real **data = new real *[nsubvol];
      for (int i = 0; i < nsubvol; i++) {
          int n = 
              (subvol[i].nbx - 1) *
              (subvol[i].nby - 1) *
              (subvol[i].nbz - 1);
          real *p = new real[n];
          for (int k = 0; k < n; k++)
              p[k] = (real)0.0;
          data[i] = p;
          }
      return data;
      }
      
  void ProjectBase::DeleteResult(real **data) {
      if (data == NULL)
          return;
      for (int i = 0; i < nsubvol; i++)
          delete[] data[i];
      delete[] data;
      }

  bool ProjectBase::WriteResult(const char *fn, real **data) {
      FILE *fp = fopen(fn, "w");
      if (fp == NULL) {
          Log("Cannot open %s", fn);
          return false;
          }
      fprintf(fp, "%d\n", nsubvol);
      for (int i = 0; i < nsubvol; i++) {
          int n = 
              (subvol[i].nbx - 1) *
              (subvol[i].nby - 1) *
              (subvol[i].nbz - 1);
          fprintf(fp, "%d\n", n);
          real *p = data[i];  
          for (int k = 0; k < n; k++)
              fprintf(fp, REAL_FMT"\n", p[k]);
          }
      fclose(fp);
      return true;
      }

  bool ProjectBase::ReadField(
          const char *fn, int ind, real (*field)[ED]) {
      int size =
          (box[ind].nx + 1) *
          (box[ind].ny + 1) *
          (box[ind].nz + 1) * ED;
      real *tmpF90 = new real[size];
      FILE *fp = fopen(fn, "rb");
      if (fp == NULL) {
          Log("Cannot open %s", fn);
          delete[] tmpF90;
          return false;
          }
    // unformatted Fortran file assumed
      int n0;
      if (fread(&n0, sizeof(int), 1, fp) != 1) {
          Log("Missing data header: file %s", fn);
          fclose(fp);
          delete[] tmpF90;
          return false;
          }
      int n = fread(tmpF90, sizeof(real), size, fp);
      if (n != size) {
          Log("Invalid data size: file %s: want %d, got %d", fn, size, n);
          fclose(fp);
          delete tmpF90;
          return false;
          }
      int n1;
      if (fread(&n1, sizeof(int), 1, fp) != 1) {
          Log("Missing data trailer: file %s", fn);
          fclose(fp);
          delete tmpF90;
          return false;
          }
      if (size * sizeof(real) != n0 || size * sizeof(real) != n1) {
          Log("Invalid data format: file %s: want %d, got (%d, %d)", 
              fn, size, n0, n1);
          fclose(fp);
          delete[] tmpF90;
          return false;
          }
      fclose(fp);
      F90ToField(ind, tmpF90, field);
      delete[] tmpF90;
      return true;
      }


  void ProjectBase::F90ToField(int ind, real *f90, real (*field)[ED]) {
      int nx = box[ind].nx;
      int ny = box[ind].ny;
      int nz = box[ind].nz;
      int im = 0;
      FOR3 (i, 0, nx, j, 0, ny, k, 0, nz) {
          int ie = 0;
          FOR3 (l, 0, LPD, m, 0, LPD, n, 0, LPD) {
              int ib = n * (LPD + 1) + m;
              ib = ib * (LPD + 1) + l;
              ib = ib * (nz + 1) + k;
              ib = ib * (ny + 1) + j;
              ib = ib * (nx + 1) + i;
              field[im][ie] = f90[ib];
              ie++;
              }
          im++;
          }
      }

  int ProjectBase::ReadInt(FILE *fp) {
      int n;
      fscanf(fp, "%d", &n);
      return n;
      }
      
  real ProjectBase::ReadReal(FILE *fp) {
      double v;
      fscanf(fp, "%lg", &v);
      return (real)v;
      }

  void ProjectBase::ReadLn(FILE *fp) {
      for ( ; ; ) {
          int ch = fgetc(fp);
          if (ch == EOF || ch == '\n')
              break;
          }
      }

  void ProjectBase::Log(const char *fmt, ...) {
      va_list arg;
      va_start(arg, fmt);
      vfprintf(stderr, fmt, arg);
      va_end(arg);
      fputc('\n', stderr);
      }

// implementation: private

  void ProjectBase::InitKnots() {
      if (LPD == 2) {
          knots[0] = (real)-1.0;
          knots[1] = (real)0.0;
          knots[2] = (real)1.0;
          }
      else if (LPD == 3) {
          knots[0] = (real)-1.0;
          knots[1] = (real)-0.4472135954999579;
          knots[2] = (real)0.4472135954999579;
          knots[3] = (real)1.0;
          }
      else if (LPD == 4) {
          knots[0] = (real)-1.0;
          knots[1] = (real)-0.6546536707079772;
          knots[2] = (real)0.0;
          knots[3] = (real)0.6546536707079772;
          knots[4] = (real)1.0;
          }
      else if (LPD == 5) {
          knots[0] = (real)-1.0;
          knots[1] = (real)-0.7650553239294647;
          knots[2] = (real)-0.2852315164806451;
          knots[3] = (real)0.2852315164806451;
          knots[4] = (real)0.7650553239294647;
          knots[5] = (real)1.0;
          }
      else if (LPD == 6) {
          knots[0] = (real)-1.0;
          knots[1] = (real)-0.8302238962785670;
          knots[2] = (real)-0.4688487934707142;
          knots[3] = (real)0.0;
          knots[4] = (real)0.4688487934707142;
          knots[5] = (real)0.8302238962785670;
          knots[6] = (real)1.0;
          }
      else if (LPD == 7) {
          knots[0] = (real)-1.0;
          knots[1] = (real)-0.8717401485096066;
          knots[2] = (real)-0.5917001814331423;
          knots[3] = (real)-0.2092992179024789;
          knots[4] = (real)0.2092992179024789;
          knots[5] = (real)0.5917001814331423;
          knots[6] = (real)0.8717401485096066;
          knots[7] = (real)1.0;
          }
      else
          assert(false);
      }

//
//    ProjectKernel
//

// construction/destruction

  ProjectKernel::ProjectKernel() { 
      isDiss = false;
      uCsv = NULL;
      uCsh = NULL;
      uCp = NULL;
      uRho = NULL;
      uQMu = NULL;
      uQKappa = NULL;
      uAlphaMu = NULL;
      uAlphaKappa = NULL;
      gradientCsv = NULL;
      gradientCsh = NULL;
      gradientCp = NULL;
      gradientRho = NULL;
      gradientQMu = NULL;
      gradientQKappa = NULL;
      gradientAlphaMu = NULL;
      gradientAlphaKappa = NULL;
      }
  
  ProjectKernel::~ProjectKernel() { 
      DeleteResults();
      DeleteFields();
      }

// interface

  bool ProjectKernel::Run(
          const char *fnGrad, const char *fnOutput, bool isDiss) {
      this->fnGrad = fnGrad;
      this->fnOutput = fnOutput;
      this->isDiss = isDiss;
      if (!Compute())
          return false;
      return true;
      }

// overrides

  void ProjectKernel::InitResults() {
      gradientCsv = NewResult();
      gradientCsh = NewResult();
      gradientCp = NewResult();
      gradientRho = NewResult();
      if (isDiss) {
          gradientQMu = NewResult();
          gradientQKappa = NewResult();
          gradientAlphaMu = NewResult();
          gradientAlphaKappa = NewResult();
          }
      }
      
  void ProjectKernel::DeleteResults() {
      DEL_RES(gradientCsv);
      DEL_RES(gradientCsh);
      DEL_RES(gradientCp);
      DEL_RES(gradientRho);
      DEL_RES(gradientQMu);
      DEL_RES(gradientQKappa);
      DEL_RES(gradientAlphaMu);
      DEL_RES(gradientAlphaKappa);
      }

  bool ProjectKernel::WriteResults() {
      string fn = fnOutput + "/gradient_csh";
      if (!WriteResult(fn.c_str(), gradientCsh))
          return false;
      fn = fnOutput + "/gradient_csv";
      if (!WriteResult(fn.c_str(), gradientCsv))
          return false;
      fn = fnOutput + "/gradient_cp";
      if (!WriteResult(fn.c_str(), gradientCp))
          return false;
      fn = fnOutput + "/gradient_rho";
      if (!WriteResult(fn.c_str(), gradientRho))
          return false;
      if (isDiss) {
          fn = fnOutput + "/gradient_Q_mu";
          if (!WriteResult(fn.c_str(), gradientQMu))
              return false;
          fn = fnOutput + "/gradient_Q_kappa";
          if (!WriteResult(fn.c_str(), gradientQKappa))
              return false;
          fn = fnOutput + "/gradient_alpha_mu";
          if (!WriteResult(fn.c_str(), gradientAlphaMu))
              return false;
          fn = fnOutput + "/gradient_alpha_kappa";
          if (!WriteResult(fn.c_str(), gradientAlphaKappa))
              return false;
          }
      return true;
      }

  bool ProjectKernel::ReadFields(int ind) {
      int n = 
          (box[ind].nx + 1) * 
          (box[ind].ny + 1) * 
          (box[ind].nz + 1);
      uCsv = new real[n][ED];
      uCsh = new real[n][ED];
      uCp = new real[n][ED];
      uRho = new real[n][ED];
      char sfx[64+1];
      sprintf(sfx, "%d", box[ind].rank-1);
      string fn = fnGrad + "/grad_csv_" + sfx;
      if (!ReadField(fn.c_str(), ind, uCsv))
          return false;
      fn = fnGrad + "/grad_csh_" + sfx;
      if (!ReadField(fn.c_str(), ind, uCsh))
          return false;
      fn = fnGrad + "/grad_cp_" + sfx;
      if (!ReadField(fn.c_str(), ind, uCp))
          return false;
      fn = fnGrad + "/grad_rho_" + sfx;
      if (!ReadField(fn.c_str(), ind, uRho))
          return false;
      if (isDiss) {
          uQMu = new real[n][ED];
          uQKappa = new real[n][ED];
          uAlphaMu = new real[n][ED];
          uAlphaKappa = new real[n][ED];
          fn = fnGrad + "/grad_Q_mu_" + sfx;
          if (!ReadField(fn.c_str(), ind, uQMu))
              return false;
          fn = fnGrad + "/grad_Q_kappa_" + sfx;
          if (!ReadField(fn.c_str(), ind, uQKappa))
              return false;
          fn = fnGrad + "/grad_alpha_mu_" + sfx;
          if (!ReadField(fn.c_str(), ind, uAlphaMu))
              return false;
          fn = fnGrad + "/grad_alphaKappa_" + sfx;
          if (!ReadField(fn.c_str(), ind, uAlphaKappa))
              return false;
          }
      return true;
      }

  void ProjectKernel::DeleteFields() {
      DEL_ARR(uCsv);
      DEL_ARR(uCsh);
      DEL_ARR(uCp);
      DEL_ARR(uRho);
      DEL_ARR(uQMu);
      DEL_ARR(uQKappa);
      DEL_ARR(uAlphaMu);
      DEL_ARR(uAlphaKappa);
      }

  void ProjectKernel::IntegrateCell(
          int isubvol,
          int indx, 
          int indy, 
          int indz, 
          int i, 
          int j, 
          int k,
          int nbx,
          int nby,
          int nbz,
          int nx,
          int ny,
          int nz,
          real dbn) {

      real *px = bintx[indx*(nx+1)+i];
      real *py = binty[indy*(ny+1)+j];
      real *pz = bintz[indz*(nz+1)+k];

      int im = (i * (ny + 1) + j) * (nz + 1) + k; 
      real *pCsv = uCsv[im];
      real *pCsh = uCsh[im];
      real *pCp = uCp[im];
      real *pRho = uRho[im];

      real csv = (real)0.0;
      real csh = (real)0.0;
      real cp = (real)0.0;
      real rho = (real)0.0;

      int ie = 0;
      FOR3 (l, 0, LPD, m, 0, LPD, n, 0, LPD) {
          real coeff = px[l] * py[m] * pz[n] * dbn;
          csv += coeff * pCsv[ie];
          csh += coeff * pCsh[ie];
          cp += coeff * pCp[ie];
          rho += coeff * pRho[ie];
          ie++;      
          }

      int ib = (indx * (nby - 1) + indy) * (nbz - 1) + indz;
      gradientCsv[isubvol][ib] += csv;
      gradientCsh[isubvol][ib] += csh;
      gradientCp[isubvol][ib] += cp;
      gradientRho[isubvol][ib] += rho;

      if (isDiss) {
          real *pQMu = uQMu[im];
          real *pQKappa = uQKappa[im];
          real *pAlphaMu = uAlphaMu[im];
          real *pAlphaKappa = uAlphaKappa[im];
      
          real qMu = (real)0.0;
          real qKappa = (real)0.0;
          real alphaMu = (real)0.0;
          real alphaKappa = (real)0.0;

          ie = 0;
          FOR3 (l, 0, LPD, m, 0, LPD, n, 0, LPD) {
              real coeff = px[i] * py[j] * pz[k] * dbn;
              qMu += coeff * pQMu[ie];
              qKappa += coeff * pQKappa[ie];
              alphaMu += coeff * pAlphaMu[ie];
              alphaKappa += coeff * pAlphaKappa[ie];
              ie++;      
              }

          gradientQMu[isubvol][ib] += qMu;
          gradientQKappa[isubvol][ib] += qKappa;
          gradientAlphaMu[isubvol][ib] += alphaMu;
          gradientAlphaKappa[isubvol][ib] += alphaKappa;
          }
      }

//
//    ProjectModel
//

// construction/destruction

  ProjectModel::ProjectModel() {
      lambda = NULL;
      mu = NULL;
      B = NULL;
      rhoinv = NULL;
      normalisation = NULL;
      pCsv = NULL;
      pCsh = NULL;
      pCp = NULL;
      pRho = NULL;
      pLambda = NULL;
      pMu = NULL;
      pB = NULL;
      pRhoinv = NULL;
      pNormalisation = NULL;
      }
      
  ProjectModel::~ProjectModel() {
      DeleteResults();
      DeleteFields();
      }
      
// interface

  bool ProjectModel::Run(const char *fnModel, const char *fnOutput) {
      this->fnModel = fnModel;
      this->fnOutput = fnOutput;
      if (!Compute())
          return false;
      return true;
      }

// overrides

  void ProjectModel::InitResults() {
      pCsv = NewResult();
      pCsh = NewResult();
      pCp = NewResult();
      pRho = NewResult();
      pLambda = NewResult();
      pMu = NewResult();
      pB = NewResult();
      pRhoinv = NewResult();
      pNormalisation = NewResult();
      }
      
  void ProjectModel::DeleteResults() {
      DEL_RES(pCsv);
      DEL_RES(pCsh);
      DEL_RES(pCp);
      DEL_RES(pRho);
      DEL_RES(pLambda);
      DEL_RES(pMu);
      DEL_RES(pB);
      DEL_RES(pRhoinv);
      DEL_RES(pNormalisation);
      }
      
  bool ProjectModel::WriteResults() {

      for (int i = 0; i < nsubvol; i++) {
          int n = 
              (subvol[i].nbx - 1) * 
              (subvol[i].nby - 1) * 
              (subvol[i].nbz - 1);
          real *sNormalisation = pNormalisation[i];
          real *sMu = pMu[i];
          real *sLambda = pLambda[i];
          real *sB = pB[i];
          real *sRhoinv = pRhoinv[i];
          real *sCsv = pCsv[i];
          real *sCsh = pCsh[i];
          real *sCp = pCp[i];
          real *sRho = pRho[i];
          for (int k = 0; k < n; k++) {
              real vNormalisation = sNormalisation[k];
              real vMu = sMu[k] / vNormalisation;
              real vLambda = sLambda[k] / vNormalisation;
              real vB = sB[k] / vNormalisation;
              real vRhoinv = sRhoinv[k] / vNormalisation;
              sCsv[k] = (real)0.001 * sqrt((vMu+vB)/vRhoinv);
              sCsh[k] = (real)0.001 * sqrt(vMu*vRhoinv);
              sCp[k] = (real)0.001 * sqrt((vLambda+2*vMu)*vRhoinv);
              sRho[k] = (real)1.0 / ((real)1000.0 * vRhoinv);
              }
          }
          
      string fn = fnOutput + "/vsh";
      if (!WriteResult(fn.c_str(), pCsh))
          return false;
      fn = fnOutput + "/vsv";
      if (!WriteResult(fn.c_str(), pCsv))
          return false;
      fn = fnOutput + "/vp";
      if (!WriteResult(fn.c_str(), pCp))
          return false;
      fn = fnOutput + "/rho";
      if (!WriteResult(fn.c_str(), pRho))
          return false;

      return true;
      }
      
  bool ProjectModel::ReadFields(int ind) {
      int n = 
          (box[ind].nx + 1) * 
          (box[ind].ny + 1) * 
          (box[ind].nz + 1);
      lambda = new real[n][ED];
      mu = new real[n][ED];
      B = new real[n][ED];
      rhoinv = new real[n][ED];
      normalisation = new real[n][ED];
      char sfx[64+1];
      sprintf(sfx, "%d", box[ind].rank-1);
      string fn = fnModel + "/lambda" + sfx;
      if (!ReadField(fn.c_str(), ind, lambda))
          return false;
      fn = fnModel + "/mu" + sfx;
      if (!ReadField(fn.c_str(), ind, mu))
          return false;
      fn = fnModel + "/B" + sfx;
      if (!ReadField(fn.c_str(), ind, B))
          return false;
      fn = fnModel + "/rhoinv" + sfx;
      if (!ReadField(fn.c_str(), ind, rhoinv))
          return false;
      for (int im = 0; im < n; im++) {
          for (int ie = 0; ie < ED; ie++)
              normalisation[im][ie] = (real)1.0;
          }
      return true;
      }
      
  void ProjectModel::DeleteFields() {
      DEL_ARR(lambda);
      DEL_ARR(mu);
      DEL_ARR(B);
      DEL_ARR(rhoinv);
      DEL_ARR(normalisation);
      }
      
  void ProjectModel::IntegrateCell(
          int isubvol,
          int indx, 
          int indy, 
          int indz, 
          int i, 
          int j, 
          int k,
          int nbx,
          int nby,
          int nbz,
          int nx,
          int ny,
          int nz,
          real dbn) {

      real *px = bintx[indx*(nx+1)+i];
      real *py = binty[indy*(ny+1)+j];
      real *pz = bintz[indz*(nz+1)+k];

      int im = (i * (ny + 1) + j) * (nz + 1) + k; 
      real *sLambda = lambda[im];
      real *sMu = mu[im];
      real *sB = B[im];
      real *sRhoinv = rhoinv[im];
      real *sNormalisation = normalisation[im];

      real vLambda = (real)0.0;
      real vMu = (real)0.0;
      real vB = (real)0.0;
      real vRhoinv = (real)0.0;
      real vNormalisation = (real)0.0;

      int ie = 0;
      FOR3 (l, 0, LPD, m, 0, LPD, n, 0, LPD) {
          real coeff = px[l] * py[m] * pz[n] * dbn;
          vLambda += coeff * sLambda[ie];
          vMu += coeff * sMu[ie];
          vB += coeff * sB[ie];
          vRhoinv += coeff * sRhoinv[ie];
          vNormalisation += coeff * sNormalisation[ie];
          ie++;      
          }

      int ib = (indx * (nby - 1) + indy) * (nbz - 1) + indz;
      pLambda[isubvol][ib] += vLambda;
      pMu[isubvol][ib] += vMu;
      pB[isubvol][ib] += vB;
      pRhoinv[isubvol][ib] += vRhoinv;
      pNormalisation[isubvol][ib] += vNormalisation;
      }

//
//    ProjectVelocity
//

// construction/destruction

  ProjectVelocity::ProjectVelocity() {
      it = 0;
      vx = NULL;
      vy = NULL;
      vz = NULL;
      pVx = NULL;
      pVy = NULL;
      pVz = NULL;
      }
      
  ProjectVelocity::~ProjectVelocity() {
      DeleteResults();
      DeleteFields();
      }

// interface

  bool ProjectVelocity::Run(
          const char *fnGrad, const char *fnOutput, int it) {
      this->fnGrad = fnGrad;
      this->fnOutput = fnOutput;
      this->it = it;
      if (!Compute())
          return false;
      return true;
      }

// overrides

  void ProjectVelocity::InitResults() {
      pVx = NewResult();
      pVy = NewResult();
      pVz = NewResult();
      }
      
  void ProjectVelocity::DeleteResults() {
      DEL_RES(pVx);
      DEL_RES(pVy);
      DEL_RES(pVz);
      }

  bool ProjectVelocity::WriteResults() {
      char sfx[64+1];
      sprintf(sfx, "%d", it);
      string fn = fnOutput + "/vx_" + sfx;
      if (!WriteResult(fn.c_str(), pVx))
          return false;
      fn = fnOutput + "/vy_" + sfx;
      if (!WriteResult(fn.c_str(), pVy))
          return false;
      fn = fnOutput + "/vz_" + sfx;
      if (!WriteResult(fn.c_str(), pVz))
          return false;
      return true;
      }
      
  bool ProjectVelocity::ReadFields(int ind) {
      int n = 
          (box[ind].nx + 1) * 
          (box[ind].ny + 1) * 
          (box[ind].nz + 1);
      vx = new real[n][ED];
      vy = new real[n][ED];
      vz = new real[n][ED];
      char sfx[64+1];
      sprintf(sfx, "%d_%d", box[ind].rank-1, it);
      string fn = fnGrad + "/vx_" + sfx;
      if (!ReadField(fn.c_str(), ind, vx))
          return false;
      fn = fnGrad + "/vy_" + sfx;
      if (!ReadField(fn.c_str(), ind, vy))
          return false;
      fn = fnGrad + "/vz_" + sfx;
      if (!ReadField(fn.c_str(), ind, vz))
          return false;
      return true;
      }
      
  void ProjectVelocity::DeleteFields() {
      DEL_ARR(vx);
      DEL_ARR(vy);
      DEL_ARR(vz);
      }
      
  void ProjectVelocity::IntegrateCell(
          int isubvol,
          int indx, 
          int indy, 
          int indz, 
          int i, 
          int j, 
          int k,
          int nbx,
          int nby,
          int nbz,
          int nx,
          int ny,
          int nz,
          real dbn) {

      real *px = bintx[indx*(nx+1)+i];
      real *py = binty[indy*(ny+1)+j];
      real *pz = bintz[indz*(nz+1)+k];

      int im = (i * (ny + 1) + j) * (nz + 1) + k; 
      real *sVx = vx[im];
      real *sVy = vy[im];
      real *sVz = vz[im];

      real vVx = (real)0.0;
      real vVy = (real)0.0;
      real vVz = (real)0.0;

      int ie = 0;
      FOR3 (l, 0, LPD, m, 0, LPD, n, 0, LPD) {
          real coeff = px[l] * py[m] * pz[n] * dbn;
          vVx += coeff * sVx[ie];
          vVy += coeff * sVy[ie];
          vVz += coeff * sVz[ie];
          ie++;      
          }

      int ib = (indx * (nby - 1) + indy) * (nbz - 1) + indz;
      pVx[isubvol][ib] += vVx;
      pVy[isubvol][ib] += vVy;
      pVz[isubvol][ib] += vVz;
      }

