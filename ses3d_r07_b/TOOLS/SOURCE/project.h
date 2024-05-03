//
//    PROJECT.H -- Engine for projecting SES3D fields
//

  using std::string;

  typedef float real;

#define REAL_FMT "%g"

  class ProjectBase {
  public:
      ProjectBase();
      virtual ~ProjectBase();
  public:
      void SetRootPath(const char *path);
  protected:
      virtual void InitResults() = 0;
      virtual void DeleteResults() = 0;
      virtual bool WriteResults() = 0;
      virtual bool ReadFields(int ind) = 0;
      virtual void DeleteFields() = 0;
      virtual void IntegrateCell(
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
          real dbn) = 0;
  protected:
    // configurable parameters
      static const int LPD = 4;
      static const int ED = (LPD + 1) * (LPD + 1) * (LPD + 1);
  protected:
      bool Compute();
      bool TryCompute();
      bool ReadBox();
      void DeleteBox();
      bool ReadBlocks();
      void DeleteBlocks();
      bool ComputePe(int ind);
      bool TryComputePe(int ind);
      void Integrate(int ind, int isubvol);
      void DeleteGrid();
      real IntLag(int idx, real limLeft, real limRight);
      real **NewResult();
      void DeleteResult(real **data);
      bool WriteResult(const char *fn, real **data);
      bool ReadField(const char *fn, int ind, real (*field)[ED]);
      void F90ToField(int ind, real *f90, real (*field)[ED]);
      int ReadInt(FILE *fp);
      real ReadReal(FILE *fp);
      void ReadLn(FILE *fp);
      void Log(const char *fmt, ...);
  private:
      void InitKnots();
  protected:
      struct Box {
          int rank;
          int nx;
          int ny;
          int nz;
          real xmin;
          real xmax;
          real ymin;
          real ymax;
          real zmin;
          real zmax;
          };
      struct Subvol {
          int nbx;
          int nby;
          int nbz;
          real *bxco;
          real *byco;
          real *bzco;
          real dbx;
          real dby;
          real dbz;
          };
  protected:
      string rootPath;
      real knots[7+1];
      int pp;
      int px;
      int py;
      int pz;
      Box *box;
      int nsubvol;
      Subvol *subvol;
      real (*x)[LPD+1];
      real (*y)[LPD+1];
      real (*z)[LPD+1];
      int *bixMin;
      int *bixMax;
      int *biyMin;
      int *biyMax;
      int *bizMin;
      int *bizMax;
      real (*bintx)[LPD+1];
      real (*binty)[LPD+1];
      real (*bintz)[LPD+1];
      };
      
  class ProjectKernel: public ProjectBase {
  public:
      ProjectKernel();
      ~ProjectKernel();
  public:
      bool Run(const char *fnGrad, const char *fnOutput, bool isDiss);
  protected:
      void InitResults();
      void DeleteResults();
      bool WriteResults();
      bool ReadFields(int ind);
      void DeleteFields();
      void IntegrateCell(
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
          real dbn);
  private:
      string fnGrad;
      string fnOutput;
      bool isDiss;
      real (*uCsv)[ED];
      real (*uCsh)[ED];
      real (*uCp)[ED];
      real (*uRho)[ED];
      real (*uQMu)[ED];
      real (*uQKappa)[ED];
      real (*uAlphaMu)[ED];
      real (*uAlphaKappa)[ED];
      real **gradientCsv;
      real **gradientCsh;
      real **gradientCp;
      real **gradientRho;
      real **gradientQMu;
      real **gradientQKappa;
      real **gradientAlphaMu;
      real **gradientAlphaKappa;
      };

  class ProjectModel: public ProjectBase {
  public:
      ProjectModel();
      ~ProjectModel();
  public:
      bool Run(const char *fnModel, const char *fnOutput);
  protected:
      void InitResults();
      void DeleteResults();
      bool WriteResults();
      bool ReadFields(int ind);
      void DeleteFields();
      void IntegrateCell(
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
          real dbn);
  private:
      string fnModel;
      string fnOutput;
      real (*lambda)[ED];
      real (*mu)[ED];
      real (*B)[ED];
      real (*rhoinv)[ED];
      real (*normalisation)[ED];
      real **pCsv;
      real **pCsh;
      real **pCp;
      real **pRho;
      real **pLambda;
      real **pMu;
      real **pB;
      real **pRhoinv;
      real **pNormalisation;
      };

  class ProjectVelocity: public ProjectBase {
  public:
      ProjectVelocity();
      ~ProjectVelocity();
  public:
      bool Run(const char *fnGrad, const char *fnOutput, int it);
  protected:
      void InitResults();
      void DeleteResults();
      bool WriteResults();
      bool ReadFields(int ind);
      void DeleteFields();
      void IntegrateCell(
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
          real dbn);
  private:
      string fnGrad;
      string fnOutput;
      int it;
      real (*vx)[ED];
      real (*vy)[ED];
      real (*vz)[ED];
      real **pVx;
      real **pVy;
      real **pVz;
      };

