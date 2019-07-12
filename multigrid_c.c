#include "mex.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "stdlib.h"
#include "multigrid.h"


typedef void (*f)(int, int N, double (*)[N+2], double(*)[N+2], double(*)[N+2], double);  //declare typdef

f func[5] = {&VCycle, &FCycle, &WCycle, &smooth, &smooth2};      //make array func of type f,
//the pointer to a function

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %
// % MULTIGRID: Calculates the variable coefficient poisson equation using V,W 
// % or F cycle in 2D with nuemnann boundary condition of d(phi(x,y))/d(n) 
// % equal to zero for M by M mesh points for equal dx and dy mesh size                                                                                                              
// % PDE being solved:
// % 
// %                                      yMMNNNNNNNNMd`     o                                     
// %                                .o:    sMm.     :y`   - .oo:    :o-                      `/+so`
// % -oooooooooooooooooo-         `yMd:     +MN-   /s    ym ms:N+   -dMh`                   sMmyo+`
// % `hMMMNmmmmmmmmmmNMh`        `mMo        :NM/`o+     M+ M+ sN     /Mm.                 sMy`    
// %   yMMm.        `hs          dMo          -mMh/      yh/Ny/ds      +Mm`  `........ `oyyMMhyyo  
// %    oMMN:      .do    `oo`  /MN            .o-        -/Nh/-        mM+  oddddddds `++mMy++/:  
// %     /NMN+    -d/     .hh.  yMy                         h:          sMd               NM.      
// %      :NMMo  /d:            hMs      oooooooooooooooooooooooooo     oMd  .-------.   :Md       
// %       .mMMyod.             oMd                .:/:`                yMy  +hhhhhhh/   yM+       
// %        `hMMh`              -MM.              sms+ym-              `NM:             .MN`       
// %         `//`                oMm.            :M/   Mo             .hMs              -y/        
// %                              +NNo.          +My::sm.           .oNNo                          
// %                               `/y+          oM/oo/`            /y+`                           
// %                                             yd                                                
// %
// %   • phi: 2D array of cell centered solution variable including ghost cells, containing 
// %   upon call the initial guess of the solution. Size = (M+2) x (N+2)
// %
// %   • f: 2D array of cell centered PDE right hand side including ghost cells (even 
// %   though these are not necessarily defined). Size = (M+2) x (N+2)
// %
// %   • rho: 2D array of cell centered PDE variable coefficents including ghost cells 
// %   (even though these are not necessarily defined). Size = (M+2) x (N+2)
// %
// %   • h: mesh spacing
// %
// %   • nIterMax: integer number of maximum V-cycle iterations to be performed
// %
// %   • alpha: convergence threshold for residual
// %
// %   • smoothCounts: Number of smoothing operations in pre and post smoothing
// %
// %   • cycleType: 1 for V-Cycle, 2 for F-Cycle and 3 for W-Cycle, 4 for Gauss Seidel,
// %   5 for Jacobi
// %
// % OutPuts:
// %   • phi: 2D array of cell centered solution variable including ghost cells, containing 
// %   the solution
// %
// %   • resid: array of residuals for each iteration
// %
// % Note:
// %   • Ensure prime factors of M and N is mostly comprised of 2. Any non-two 
// %   factors will slow down convergence. 
// %
// % Author: Adil Ansari, MS
// % 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (nrhs==0) {mexPrintf("Test\n"); return;}
    /* Macros for the ouput and input arguments */
    double *B, *C, *pin, *d2pin, *rhoin, dx, thresh;
    int i, j, k, maxit, choice;
    int M, N;
    pin = mxGetPr(prhs[0]);
    d2pin = mxGetPr(prhs[1]);
    rhoin =  mxGetPr(prhs[2]);
    dx =  mxGetScalar(prhs[3]);
    M = mxGetM(prhs[0])-2; /* Get the dimensions of A */
    N = mxGetN(prhs[0])-2;
    maxit = mxGetScalar(prhs[4]);
    thresh = mxGetScalar(prhs[5]);
    iterM =  mxGetScalar(prhs[6]);
    choice =  mxGetScalar(prhs[7])-1;

    plhs[0] = mxCreateDoubleMatrix(M+2, N+2, mxREAL); /* Create the output matrix */
    plhs[1] = mxCreateDoubleMatrix(1, maxit+1, mxREAL); /* Create the output matrix */
    B = mxGetPr(plhs[0]);
    C = mxGetPr(plhs[1]); /* Get the pointer to the data of B */
    if (choice<0 || choice>4){
        choice = 3;
    }
    
    double (*phi)[N+2];  phi = malloc((M+2) * sizeof *phi);
    double (*d2phi)[N+2];d2phi = malloc((M+2) * sizeof *d2phi);
    double (*rhom)[N+2]; rhom = malloc((M+2) * sizeof *rhom);
    double (*r)[N+2]; r = malloc((M+2) * sizeof *r);

    for(j = 0; j < N+2; j++){
        for(i = 0; i < M+2; i++){
            phi  [i][j] = pin  [i + (M+2)*j];
            d2phi[i][j] = d2pin[i + (M+2)*j];
            rhom [i][j] = rhoin[i + (M+2)*j];
        }
    }
    
    /*********************************************************
     *                  Begin Computations
     **********************************************************/
    residual(M, N, r, phi, d2phi, rhom, dx);
    C[0] = maxabs(M, N, r);
    for (k = 1; k <= maxit; k++){
        func[choice](M, N, phi, d2phi, rhom, dx);
        if (!(k%5)) meanzero(M,N,phi);
        residual(M, N, r, phi, d2phi, rhom, dx);
        C[k] = maxabs(M, N, r)/C[0];
        if (C[k] < thresh) break;
    }
    C[0] = 1;
    /***********************************************************
     *                     End Computation                     *
     ***********************************************************/
    for(i = 0; i < M+2; i++){
        for(j = 0; j < N+2; j++)
            B[i + (M+2)*j] = phi[i][j];
    }
    free(phi); free(r); free(rhom); free(d2phi);
    return;
}