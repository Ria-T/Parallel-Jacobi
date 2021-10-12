/************************************************************
 * Program to solve a finite difference
 * discretization of the screened Poisson equation:
 * (d2/dx2)u + (d2/dy2)u - alpha u = f
 * with zero Dirichlet boundary condition using the iterative
 * Jacobi method with overrelaxation.
 *
 * RHS (source) function
 *   f(x,y) = -alpha*(1-x^2)(1-y^2)-2*[(1-x^2)+(1-y^2)]
 *
 * Analytical solution to the PDE
 *   u(x,y) = (1-x^2)(1-y^2)
 *
 * Current Version: Christian Iwainsky, RWTH Aachen University
 * MPI C Version: Christian Terboven, RWTH Aachen University, 2006
 * MPI Fortran Version: Dieter an Mey, RWTH Aachen University, 1999 - 2005
 * Modified: Sanjiv Shah,        Kuck and Associates, Inc. (KAI), 1998
 * Author:   Joseph Robicheaux,  Kuck and Associates, Inc. (KAI), 1998
 *
 * Unless READ_INPUT is defined, a meaningful input dataset is used (CT).
 *
 * Input : n     - grid dimension in x direction
 *         m     - grid dimension in y direction
 *         alpha - constant (always greater than 0.0)
 *         tol   - error tolerance for the iterative solver
 *         relax - Successice Overrelaxation parameter
 *         mits  - maximum iterations for the iterative solver
 *
 * On output
 *       : u(n,m)       - Dependent variable (solution)
 *       : f(n,m,alpha) - Right hand side function
 *
 *************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define BLOCK_SIZE 256

// fy, fx values are standard, but the program recalculates them at each iteration.
// Having an array with all the values pre-calculated once saves a lot of computing time in the iteration
int calculate_fX_fY_arrays(double xStart, double yStart, int n, int m,  double deltaX, double deltaY, double **fX, double **fY){
	*fX = (double*)calloc(n, sizeof(double));
	*fY = (double*)calloc(m, sizeof(double));
	if( *fX == NULL || *fY == NULL) return -1;
	for (int x = 0; x < n; x++){
		(*fX)[x] = xStart + x*deltaX;
	}
	for (int y = 0; y < m; y++){
		(*fY)[y] = yStart + y*deltaY;
	}
	return 0;

}

/*************************************************************
 * Performs one iteration of the Jacobi method and computes
 * the residual value.
 *
 * NOTE: u(0,*), u(maxXCount-1,*), u(*,0) and u(*,maxYCount-1)
 * are BOUNDARIES and therefore not part of the solution.
 *************************************************************/

__global__ void one_jacobi_iteration(double xStart, double yStart,
                            int maxXCount, int maxYCount,
                            double *src, double *dst,
                            double deltaX, double deltaY,
                            double alpha, double omega,
                            double const * const cx, double const * const cy, double const * const cc,
                            double *fX, double *fY,
                            double *iterationError)
{
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    double error = 0.0;
    double updateVal;
    double f;
    // Coefficients
    /*double cx = 1.0/(deltaX*deltaX);
    double cy = 1.0/(deltaY*deltaY);
    double cc = -2.0*cx-2.0*cy-alpha;*/
    if ((0<x && x<maxXCount-1) && (0<y && y<maxYCount-1)) /* if x and y are within the boundaries */
    {
        //fX = xStart + (x-1)*deltaX; Pre calculated to save time
        f = -alpha*(1.0-fX[x-1]*fX[x-1])*(1.0-fY[y-1]*fY[y-1]) - 2.0*(1.0-fX[x-1]*fX[x-1]) - 2.0*(1.0-fY[y-1]*fY[y-1]);
        updateVal = ((SRC(x-1,y) + SRC(x+1,y))*(*cx) +
                	(SRC(x,y-1) + SRC(x,y+1))*(*cy) +
                	SRC(x,y)*(*cc) - f
					)/(*cc);
        DST(x,y) = SRC(x,y) - omega*updateVal;
        error += updateVal*updateVal;
    }

    *iterationError = sqrt(error)/((maxXCount-2)*(maxYCount-2));
}



/**********************************************************
 * Checks the error between numerical and exact solutions
 **********************************************************/
double checkSolution(double xStart, double yStart,
                     int maxXCount, int maxYCount,
                     double *u,
                     double deltaX, double deltaY,
                     double alpha)
{
#define U(XX,YY) u[(YY)*maxXCount+(XX)]
    int x, y;
    double fX, fY;
    double localError, error = 0.0;

    for (y = 1; y < (maxYCount-1); y++)
    {
        fY = yStart + (y-1)*deltaY;
        for (x = 1; x < (maxXCount-1); x++)
        {
            fX = xStart + (x-1)*deltaX;
            localError = U(x,y) - (1.0-fX*fX)*(1.0-fY*fY);
            error += localError*localError;
        }
    }
    return sqrt(error)/((maxXCount-2)*(maxYCount-2));
}


int main(int argc, char **argv)
{
    int n, m, mits;
    double alpha, tol, relax;
    double maxAcceptableError;
    double error, *d_error;
    double *u, *u_old, *tmp;
    double *d_u, *d_u_old;
    int allocCount;
    int iterationCount, maxIterationCount;

//    printf("Input n,m - grid dimension in x,y direction:\n");
    scanf("%d,%d", &n, &m);
//    printf("Input alpha - Helmholtz constant:\n");
    scanf("%lf", &alpha);
//    printf("Input relax - successive over-relaxation parameter:\n");
    scanf("%lf", &relax);
//    printf("Input tol - error tolerance for the iterative solver:\n");
    scanf("%lf", &tol);
//    printf("Input mits - maximum solver iterations:\n");
    scanf("%d", &mits);


    printf("-> %d, %d, %g, %g, %g, %d\n", n, m, alpha, relax, tol, mits);

    allocCount = (n+2)*(m+2);
    // Those two calls also zero the boundary elements
    u = 	(double*)calloc(allocCount, sizeof(double)); //reverse order
    u_old = (double*)calloc(allocCount, sizeof(double));

//    printf("allocCount=%d u=%p u_old=%p\n", allocCount, u, u_old);

    if (u == NULL || u_old == NULL)
    {
        printf("Not enough memory for two %ix%i matrices\n", n+2, m+2);
        exit(1);
    }
    maxIterationCount = mits;
    maxAcceptableError = tol;

    // Solve in [-1, 1] x [-1, 1]
    double xLeft = -1.0, xRight = 1.0;
    double yBottom = -1.0, yUp = 1.0;

    double deltaX = (xRight-xLeft)/(n-1);
    double deltaY = (yUp-yBottom)/(m-1);

    iterationCount = 0;
    error = HUGE_VAL;
    clock_t start = clock(), diff;


    // Jacobi iteration Coefficient variables (values are calculated once -> faster execution of inline function)
    double JIV_cx = 1.0/(deltaX*deltaX);
    double JIV_cy = 1.0/(deltaY*deltaY);
    double JIV_cc = -2.0*JIV_cx-2.0*JIV_cy-alpha;
	double *fX = NULL,*fY = NULL;

	if( calculate_fX_fY_arrays(xLeft,yBottom,n,m,deltaX,deltaY,&fX,&fY) != 0)
		{ printf("Can not calculate the table due to memory limitations\n" ); exit(1); }

    size_t arraySizeInBytes = allocCount * sizeof(double);
    cudaMalloc(&d_u, arraySizeInBytes); 
    cudaMalloc(&d_u_old, arraySizeInBytes); 
    cudaMalloc(&d_error, sizeof(double));

    cudaMemcpy(d_u, u, arraySizeInBytes, cudaMemcpyHostToDevice); 
    cudaMemcpy(d_u_old, u_old, arraySizeInBytes, cudaMemcpyHostToDevice); 

    /* Iterate as long as it takes to meet the convergence criterion */
    while (iterationCount < maxIterationCount && error > maxAcceptableError)
    {
        //printf("Iteration %i", iterationCount);

        one_jacobi_iteration<<<(n*m+BLOCK_SIZE) / BLOCK_SIZE, BLOCK_SIZE>>>(xLeft, yBottom,
                                                                                    n+2, m+2,
                                                                                    u_old, u,
                                                                                    deltaX, deltaY,
                                                                                    alpha, relax,
                                                                                    &JIV_cx, &JIV_cy, &JIV_cc,
                                                                                    fX,fY,
                                                                                    d_error);
        cudaDeviceSynchronize();

        //printf("\tError %g\n", error);
        iterationCount++;
        // Swap the buffers
        tmp = u_old;
        u_old = u;
        u = tmp;
    }

    // Copy results back to host.
    // d_u_old holds the solution after the most recent buffers swap.
    cudaMemcpy(u_old, d_u_old, arraySizeInBytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(&error, d_error, sizeof(double), cudaMemcpyDeviceToHost);

    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("Residual %g\n",error);

    free(fX);
    free(fY);

    // u_old holds the solution
    double absoluteError = checkSolution(xLeft, yBottom,
                                         n+2, m+2,
                                         u_old,
                                         deltaX, deltaY,
                                         alpha);
    printf("The error of the iterative solution is %g\n", absoluteError);

    cudaFree(d_u); 
    cudaFree(d_u_old); 
    cudaFree(d_error);
    free(u);
    free(u_old);

    return 0;
}
