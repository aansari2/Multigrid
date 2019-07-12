#include "mex.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "stdlib.h"
#include "multigrid.h"

int bc = 1, iter = 1, iterT = 1;
double beta;

void BC(int M, int N, double phi[M+2][N+2]){
    int i,j;
    for(i = 0; i < M+2; i++){
        phi[i][N+1] = bc * phi[i][N];
        phi[i][0] = bc * phi[i][1];
    }
    for(j = 0; j < N+2; j++){
        phi[M+1][j] = bc * phi[M][j];
        phi[0][j] = bc * phi[1][j];
    }
}

/* secondary smoother*/
void smooth2(int M, int N, double phi[M+2][N+2], double d2phi[M+2][N+2], double rho[M+2][N+2], double dx){
    int i,j,k;
    double (*p)[N+2]; p = malloc((M+2) * sizeof *p);
    for (k = 0; k < iterM; k++){
        for(i = 1; i < M+1; i++){
            for(j = 1; j < N+1; j++){
                double rho_ij = rho[i][j];
                double rho_ipj = rho[i+1][j] + rho_ij;
                double rho_imj = rho[i-1][j] + rho_ij;
                double rho_ijp = rho[i][j+1] + rho_ij;
                double rho_ijm = rho[i][j-1] + rho_ij;
                double im = 1.0/(dx*dx)*(
                        1.0/(rho_ipj)+1.0/(rho_imj)+1.0/(rho_ijp)+1.0/(rho_ijm));
                double pg = 1/(dx*dx)*(phi[i+1][j]/rho_ipj+phi[i-1][j]/rho_imj+
                        phi[i][j+1]/rho_ijp+phi[i][j-1]/rho_ijm)-d2phi[i][j];
                //double beta = 1.0;
                p[i][j] = 1.0/im*pg;//beta/im*pg+(1-beta)*phi[i][j];
            }
        }
        BC(M,N,p);
        memcpy(phi, p, sizeof (double) * (M+2) * (N+2));

        //if ((im>(0.1/(dx*dx))) && (im<(2.0/(dx*dx))))
    }
    free(p);
}

/* Primary smoother*/
void smooth(int M, int N, double phi[M+2][N+2], double d2phi[M+2][N+2], double rho[M+2][N+2], double dx){
    int i,j,k;
    
    for (k = 0; k < iterM; k++){
        for(i = 1; i < M+1; i++){
            for(j = 1; j < N+1; j++){
                double rho_ij = rho[i][j];
                double rho_ipj = rho[i+1][j] + rho_ij;
                double rho_imj = rho[i-1][j] + rho_ij;
                double rho_ijp = rho[i][j+1] + rho_ij;
                double rho_ijm = rho[i][j-1] + rho_ij;
                double im = 1.0/(dx*dx)*(
                        1.0/(rho_ipj)+1.0/(rho_imj)+1.0/(rho_ijp)+1.0/(rho_ijm));
                double pg = 1/(dx*dx)*(phi[i+1][j]/rho_ipj+phi[i-1][j]/rho_imj+
                        phi[i][j+1]/rho_ijp+phi[i][j-1]/rho_ijm)-d2phi[i][j];
                //double beta = 1.0;
                phi[i][j] = 1.0/im*pg;//beta/im*pg+(1-beta)*phi[i][j];
            }
        }
        BC(M,N,phi);
        //if ((im>(0.1/(dx*dx))) && (im<(2.0/(dx*dx))))
    }
}


/*primary calcResid*/
void residual(int M, int N, double r[M+2][N+2], double phi[M+2][N+2], double d2phi[M+2][N+2], double rho[M+2][N+2], double dx){
    int i,j,  q = 0,  m = M,  n = N;
    //#pragma omp parallel for private(j)
    for(i = 1; i < m+1; i++){
        for(j = 1; j < n+1; j++){
            double rho_ij = rho[i][j];
            double phi_ij = phi[i][j];
            double pg = 1.0/(dx*dx)*(
                    (phi[i+1][j]-phi_ij)/(rho[i+1][j]+rho_ij)-
                    (phi_ij-phi[i-1][j])/(rho_ij+rho[i-1][j])+
                    (phi[i][j+1]-phi_ij)/(rho[i][j+1]+rho_ij)-
                    (phi_ij-phi[i][j-1])/(rho[i][j-1]+rho_ij));
            r[i][j] = d2phi[i][j] - pg;
        }
    }
    BC(M,N,r);
}

/* restricting in Multigrid*/
void Restrict(int M, int N, double rhs[M/2+2][N/2+2], double r[M+2][N+2]){
    int i,j, m = M/2,  n = N/2;
    //#pragma omp parallel for private(j)
    for(i = 1; i < m+1; i++){
        for(j = 1; j < n+1; j++){
            rhs[i][j] = 0.25*(r[2*i][2*j]+
                    r[2*i-1][2*j]+
                    r[2*i][2*j-1]+
                    r[2*i-1][2*j-1]);
        }
    }
    BC(m,n,rhs);
}


/* prolonging in Multigrid*/
void prolong(int M, int N, double eps[M+2][N+2], double epsc[2*M+2][2*N+2], double alpha){
    int i,j, mn = M, nn = N; 
    //double cart = ((double) (log2int(M)+1))/8.0; alpha = 0.05;
    //double cart = 0.03;
    for(i = 1; i < mn+1; i++){
        for(j = 1; j < nn+1; j++){
            double val = eps[i][j];
            epsc[2*i-1][2*j-1] += alpha*(0.5625*val+0.1875*(eps[i-1][j]+eps[i][j-1])+0.0625*eps[i-1][j-1]);
            epsc[2*i-1][2*j] += alpha*(0.5625*val+0.1875*(eps[i-1][j]+eps[i][j+1])+0.0625*eps[i-1][j+1]);
            epsc[2*i][2*j-1] += alpha*(0.5625*val+0.1875*(eps[i+1][j]+eps[i][j-1])+0.0625*eps[i+1][j-1]);
            epsc[2*i][2*j] += alpha*(0.5625*val+0.1875*(eps[i+1][j]+eps[i][j+1])+0.0625*eps[i+1][j+1]);
        }
    }
    int m = mn * 2;
    int n = nn * 2;
    BC(m,n,epsc);
}
double maxabs(int M, int N, double r[M+2][N+2]){
    double maximum = -1e300;
    int i,j;
    for(i = 0; i < M+2; i++){
        for(j = 0; j < N+2; j++){
            if (fabs(r[i][j])>maximum){
                maximum = fabs(r[i][j]);
            }
        }
    }
    return maximum;
}

double minabs(int M, int N, double r[M+2][N+2]){
    double minimum = 1e300;
    int i,j;
    for(i = 0; i < M+2; i++){
        for(j = 0; j < N+2; j++){
            if ((r[i][j])<minimum){
                minimum = (r[i][j]);
            }
        }
    }
    return minimum;
}

int log2int(int val){
    if (val == 1) return 0;
    int ret = 0;
    while (val > 1){
        val >>= 1;
        ret++;
    }
    return ret;
}

void printArray(int M, int N, double a[M][N]){
    int j,k;
    for (k=0;k<M;k++){
        for (j=0;j<N;j++){
            printf("%.2e\t",a[k][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void zerofy(int M, int N, double a[M][N]){
    int i,j;
    for (i = 0; i < M; i++){
        for (j = 0; j < N; j++)
            a[i][j] = 0;
    }
}

void meanzero(int M, int N, double phi[M+2][N+2]){
    double average = 0;
    int i,j;
    for(i = 0; i < M+2; i++){
        for(j = 0; j < N+2; j++){
            average = average + phi[i][j];
        }
    }
    average = average/((M+2)*(N+2));
    for(i = 0; i < M+2; i++){
        for(j = 0; j < N+2; j++){
            phi[i][j] = phi[i][j] - average;
        }
    }
}

void VCycle(int M, int N, double phi[M+2][N+2], double div[M+2][N+2], double rho[M+2][N+2], double dx){
    smooth(M, N, phi, div, rho, dx);
    
    double (*r)[N+2]; r = malloc((M+2) * sizeof *r);
    residual(M, N, r, phi, div, rho, dx);
    
    
    double (*rhs)[N/2+2]; rhs = malloc((M/2+2) * sizeof *rhs);
    Restrict(M, N, rhs, r);
    double (*rho_2h)[N/2+2]; rho_2h = malloc((M/2+2) * sizeof *rho_2h);
    Restrict(M, N, rho_2h, rho);
    
    
    double (*eps)[N/2+2]; eps = malloc((M/2+2) * sizeof *eps); zerofy(M/2+2,N/2+2,eps);
    //printArray(M/2+2,N/2+2,rhs);
    
    if (((M/2)%2) || ((N/2)%2)){
        smooth(M/2, N/2, eps, rhs, rho_2h, 2*dx);}
    else{
        VCycle(M/2, N/2, eps, rhs, rho_2h, 2*dx);
    }
    
    prolong(M/2,N/2,eps,phi,1.0);
    smooth(M, N, phi, div, rho, dx);
}

void WCycle(int M, int N, double phi[M+2][N+2], double div[M+2][N+2], double rho[M+2][N+2], double dx){
    smooth(M, N, phi, div, rho, dx);
    
    double (*r)[N+2]; r = malloc((M+2) * sizeof *r);
    residual(M, N, r, phi, div, rho, dx);
    
    double (*rhs)[N/2+2]; rhs = malloc((M/2+2) * sizeof *rhs);
    Restrict(M, N, rhs, r);
    
    double (*rho_2h)[N/2+2]; rho_2h = malloc((M/2+2) * sizeof *rho_2h);
    Restrict(M, N, rho_2h, rho);
    double alpha = 1;//pow(minabs(M/2,N/2,rho_2h)/maxabs(M/2,N/2,rho_2h),0.5);

    double (*eps)[N/2+2]; eps = malloc((M/2+2) * sizeof *eps); zerofy(M/2+2,N/2+2,eps);    
    if (((M/2)%2) || ((N/2)%2)){
        iterT = iterM; iterM = (M/2)*(N/2);
        smooth(M/2, N/2, eps, rhs, rho_2h, 2*dx);
        iterM = iterT;
    }
    else{
        WCycle(M/2, N/2, eps, rhs, rho_2h, 2*dx);
    }
    prolong(M/2,N/2,eps,phi,alpha);
    smooth(M, N, phi, div, rho, dx);
    residual(M, N, r, phi, div, rho, dx);
    Restrict(M, N, rhs, r);
    zerofy(M/2+2,N/2+2,eps);    
    if (((M/2)%2) || ((N/2)%2)){
        iterT = iterM; iterM = (M/2)*(N/2);
        smooth(M/2, N/2, eps, rhs, rho_2h, 2*dx);
        iterM = iterT;}
    else
        WCycle(M/2, N/2, eps, rhs, rho_2h, 2*dx);
    
    prolong(M/2,N/2,eps,phi,alpha);
    smooth(M, N, phi, div, rho, dx);
    free(r); free(rhs); free(rho_2h); free(eps);
}


void FCycle(int M, int N, double phi[M+2][N+2], double div[M+2][N+2], double rho[M+2][N+2], double dx){
    smooth(M, N, phi, div, rho, dx);
    
    double (*r)[N+2]; r = malloc((M+2) * sizeof *r);
    residual(M, N, r, phi, div, rho, dx);
    
    double (*rhs)[N/2+2]; rhs = malloc((M/2+2) * sizeof *rhs);
    Restrict(M, N, rhs, r);
    
    double (*rho_2h)[N/2+2]; rho_2h = malloc((M/2+2) * sizeof *rho_2h);
    Restrict(M, N, rho_2h, rho);
    
    double (*eps)[N/2+2]; eps = malloc((M/2+2) * sizeof *eps); zerofy(M/2+2,N/2+2,eps);
    double alpha = pow(minabs(M/2,N/2,rho_2h)/maxabs(M/2,N/2,rho_2h),2.0);

    if (((M/2)%2) || ((N/2)%2)){
        smooth(M/2, N/2, eps, rhs, rho_2h, 2*dx);}
    else{
        FCycle(M/2, N/2, eps, rhs, rho_2h, 2*dx);
    }
    prolong(M/2,N/2,eps,phi,1.0);
    smooth(M, N, phi, div, rho, dx);
    residual(M, N, r, phi, div, rho, dx);
    Restrict(M, N, rhs, r);
    zerofy(M/2+2,N/2+2,eps);
    if (((M/2)%2) || ((N/2)%2))
        smooth(M/2, N/2, eps, rhs, rho_2h, 2*dx);
    else
        VCycle(M/2, N/2, eps, rhs, rho_2h, 2*dx);
    prolong(M/2,N/2,eps,phi,1.0);
    smooth(M, N, phi, div, rho, dx);

    free(r); free(rhs); free(rho_2h); free(eps);
}