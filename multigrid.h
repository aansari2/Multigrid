int iterM;


void BC(int M, int N, double phi[M+2][N+2]);

/* secondary smoother*/
void smooth2(int M, int N, double phi[M+2][N+2], double d2phi[M+2][N+2], double rho[M+2][N+2], double dx);

/* Primary smoother*/
void smooth(int M, int N, double phi[M+2][N+2], double d2phi[M+2][N+2], double rho[M+2][N+2], double dx);


/*primary calcResid*/
void residual(int M, int N, double r[M+2][N+2], double phi[M+2][N+2], double d2phi[M+2][N+2], double rho[M+2][N+2], double dx);

/* restricting in Multigrid*/
void Restrict(int M, int N, double rhs[M/2+2][N/2+2], double r[M+2][N+2]);

/* prolonging in Multigrid*/
void prolong(int M, int N, double eps[M+2][N+2], double epsc[2*M+2][2*N+2], double alpha);

int log2int(int val);

double maxabs(int M, int N, double r[M+2][N+2]);

double minabs(int M, int N, double r[M+2][N+2]);

void printArray(int M, int N, double a[M][N]);

void zerofy(int M, int N, double a[M][N]);

void meanzero(int M, int N, double phi[M+2][N+2]);

void VCycle(int M, int N, double phi[M+2][N+2], double div[M+2][N+2], double rho[M+2][N+2], double dx);

void WCycle(int M, int N, double phi[M+2][N+2], double div[M+2][N+2], double rho[M+2][N+2], double dx);

void FCycle(int M, int N, double phi[M+2][N+2], double div[M+2][N+2], double rho[M+2][N+2], double dx);