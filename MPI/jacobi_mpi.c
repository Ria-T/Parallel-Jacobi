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
#include <mpi.h>

#include "partitioning.h"

/*************************************************************
 * Performs one iteration of the Jacobi method and computes
 * the residual value.
 *
 * NOTE: u(0,*), u(maxXCount-1,*), u(*,0) and u(*,maxYCount-1)
 * are BOUNDARIES and therefore not part of the solution.
 *************************************************************/
inline double one_jacobi_iteration(double xStart, double yStart,
                            int maxXCount, int maxYCount,
                            double *src, double *dst,
                            double deltaX, double deltaY,
                            double alpha, double omega,
                            double const * const cx, double const * const cy, double const * const cc)
{
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]
    int x, y;
    double fX, fY;
    double error = 0.0;
    double updateVal;
    double f;
    // Coefficients
    /*double cx = 1.0/(deltaX*deltaX);
    double cy = 1.0/(deltaY*deltaY);
    double cc = -2.0*cx-2.0*cy-alpha;*/

    for (y = 1; y < (maxYCount-1); y++)
    {
        fY = yStart + (y-1)*deltaY;
        for (x = 1; x < (maxXCount-1); x++)
        {
            fX = xStart + (x-1)*deltaX;
            f = -alpha*(1.0-fX*fX)*(1.0-fY*fY) - 2.0*(1.0-fX*fX) - 2.0*(1.0-fY*fY);
            updateVal = (	(SRC(x-1,y) + SRC(x+1,y))*(*cx) +
                			(SRC(x,y-1) + SRC(x,y+1))*(*cy) +
                			SRC(x,y)*(*cc) - f
						)/(*cc);
            DST(x,y) = SRC(x,y) - omega*updateVal;
            error += updateVal*updateVal;
        }
    }
    return sqrt(error)/((maxXCount-2)*(maxYCount-2));
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
    double error;
    double *u, *u_old, *tmp;
    int allocCount;
    int iterationCount, maxIterationCount;
    double t1, t2;

    int world_size, rank;

    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
//    printf("Input n,m - grid dimension in x,y direction:\n");
    scanf("%d,%d", &n, &m);
//    printf("Input alpha - Helmholtz constant:\n");
    scanf("%lf", &alpha);
//    printf("Input relax - successive over-relaxation parameter:\n");
    scanf("%lf", &relax);
//    printf("Input tol - error tolerance for the iterrative solver:\n");
    scanf("%lf", &tol);
//    printf("Input mits - maximum solver iterations:\n");
    scanf("%d", &mits);
    }

    // Proccess 0 reads the input, the rest get the input after proccess 0 bcasts it
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&relax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mits, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int sizee[2], coords[2];
    sizee[0] = n;
    sizee[1] = m;
    get_local_table(sizee, coords, &rank, world_size);
    return 0;

    printf("%d out of %d\n-> %d, %d, %g, %g, %g, %d\n gonna take %d out of %d, and %d out of %d\n"
    ,rank, world_size, n, m, alpha, relax, tol, mits, n/world_size, n, m/world_size, m);

    // Now each proccess will get it's "working" area of the table
    MPI_Barrier(MPI_COMM_WORLD);
    int startRow,endRow;
    if(rank == 0){
        startRow = 0;
        endRow = round(m/world_size);//m/world_size-1;
        MPI_Send(&endRow, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }else if(rank != world_size-1){
        MPI_Recv(&startRow, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        startRow++; // Starting where the previous ended
        endRow = startRow + round(m/world_size);//m/world_size-1;
        MPI_Send(&endRow, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }else{
        MPI_Recv(&startRow, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        startRow++;
        endRow = m;
    }

    int bold[2][2], size[2];
    size[0] = n;
    size[1] = m;
    bold[0][0] = 0; bold[0][1] = startRow;
    bold[1][0] = n-1; bold[1][1] = endRow;
    draw_table(bold, size , rank);

    int original_n = n, original_m = m;
    n = n;
    m = endRow - startRow;
    printf("(%d) small range: %d,%d\n",rank,n,m);

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

    //MPI_Init(NULL,NULL); should it get up?
    t1 = MPI_Wtime();

    /* Iterate as long as it takes to meet the convergence criterion */
    while (iterationCount < maxIterationCount && error > maxAcceptableError)
    {    	
        //printf("Iteration %i", iterationCount);

        error = one_jacobi_iteration(xLeft, yBottom,
                                     n+2, m+2,
                                     u_old, u,
                                     deltaX, deltaY,
                                     alpha, relax,
                                     &JIV_cx, &JIV_cy, &JIV_cc);

        //printf("\tError %g\n", error);
        iterationCount++;
        // Swap the buffers
        tmp = u_old;
        u_old = u;
        u = tmp;
    }

    t2 = MPI_Wtime();
    printf( "Iterations=%3d Elapsed MPI Wall time is %f\n", iterationCount, t2 - t1 ); 
    MPI_Finalize();
    
    
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("Residual %g\n",error);

    // u_old holds the solution after the most recent buffers swap
    double absoluteError = checkSolution(xLeft, yBottom,
                                         n+2, m+2,
                                         u_old,
                                         deltaX, deltaY,
                                         alpha);
    printf("The error of the iterative solution is %g\n", absoluteError);

    return 0;
}
