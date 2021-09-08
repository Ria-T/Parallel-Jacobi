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
#include "helpers.h"

/*************************************************************
 * Performs one iteration of the Jacobi method and computes
 * the residual value.
 *
 * NOTE: u(0,*), u(maxXCount-1,*), u(*,0) and u(*,maxYCount-1)
 * are BOUNDARIES and therefore not part of the solution.
 *************************************************************/
/*inline*/ double one_jacobi_iteration(double xStart, double yStart,
                            int maxXCount, int maxYCount,
                            double *src, double *dst,
                            double deltaX, double deltaY,
                            double alpha, double omega,
                            double const * const cx, double const * const cy, double const * const cc,
                            int *neighbors)
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


    MPI_Comm cart_comm;
    int size[2], *neighbors, **coordinates;
    get_local_table(&n, &m, &coordinates, &rank, world_size,&cart_comm);
    n = coordinates[1][0]-coordinates[0][0];
    m = coordinates[1][1]-coordinates[0][1];
    size[0] = n;
    size[1] = m;
    printf("%d:%dx%d\n",rank,n,m);
    neighbors = get_neighbors(&cart_comm);

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

    MPI_Request SRequests[4],RRequests[4];
    int requestsCount = 0;
    MPI_Datatype table_row, table_column;
    MPI_Type_vector(n, 1, 1, MPI_DOUBLE, &table_row);
    MPI_Type_commit(&table_row);
    MPI_Type_vector(m, 1, n+2, MPI_DOUBLE, &table_column);
    MPI_Type_commit(&table_column);

    //MPI_Init(NULL,NULL); should it get up?
    t1 = MPI_Wtime();

    //init_debug();

    /* Iterate as long as it takes to meet the convergence criterion */
    while (iterationCount < maxIterationCount && error > maxAcceptableError)
    {
        //printf("Iteration %i\n", iterationCount);

        if(rank > -1){
            int l;
	
			// second row
            for(l=0; l<n+2; l++){
                u_old[1*(n+2)+l]=0.11+rank;
            }
			
			// pre last row
            for(l=0; l<n+2; l++){
                u_old[(m)*(n+2)+l]=0.44+rank;
            }
			
			// left column + 1
			for(l=0; l<m+2; l++){
                u_old[l*(n+2)+1]=0.22+rank;
            }
			
			// right column - 1
			for(l=0; l<m+2; l++){
                u_old[(l)*(n+2)+n]=0.33+rank;
            }

        }
        requestsCount = 0;

        write_table(u_old,size,rank,neighbors);

        if(neighbors[UP] >= 0){
            MPI_Irecv(u_old, 1, table_row, neighbors[UP], 0, cart_comm, &RRequests[requestsCount]);
            MPI_Isend(&u_old[n+2], 1, table_row, neighbors[UP], 0, cart_comm, &SRequests[requestsCount]);
            printf("%d got/sent line from/to %d\n",rank,neighbors[UP]);
            requestsCount++;
        }

        if(neighbors[DOWN] >= 0){
            MPI_Irecv(&u_old[(m+1)*(n+2)], 1, table_row, neighbors[DOWN], 0, cart_comm, &RRequests[requestsCount]);
            MPI_Isend(&u_old[n+2], 1, table_row, neighbors[DOWN], 0, cart_comm, &SRequests[requestsCount]);
            requestsCount++;
            printf("%d got/sent line from/to %d\n",rank,neighbors[DOWN]);
        }

        if(neighbors[LEFT] >= 0){
            MPI_Irecv(u_old, 1, table_column, neighbors[LEFT], 0, cart_comm,&RRequests[requestsCount]);
            MPI_Isend(&u_old[n+2+1], 1, table_column, neighbors[LEFT], 0, cart_comm,&SRequests[requestsCount]);
            requestsCount++;
            //printf("%d got/sent line from/to %d\n",rank,neighbors[LEFT]);
        }

        if(neighbors[RIGHT] >= 0){
            MPI_Irecv(&u_old[n+1], 1, table_column, neighbors[RIGHT], 0, cart_comm,&RRequests[requestsCount]);
            MPI_Isend(&u_old[n], 1, table_column, neighbors[RIGHT], 0, cart_comm,&SRequests[requestsCount]);
            requestsCount++;
            //printf("%d got/sent line from/to %d\n",rank,neighbors[RIGHT]);
        }

//        MPI_Wait(RRequests, MPI_STATUS_IGNORE); 
        MPI_Waitall(requestsCount, RRequests, MPI_STATUS_IGNORE);
        /*if(neighbors[UP] >= 0 ){
            MPI_Wait(RRequests + UP, MPI_STATUS_IGNORE); 
        }
        if(neighbors[DOWN] >= 0 ){
            MPI_Wait(RRequests + DOWN, MPI_STATUS_IGNORE);
        }*/

        //printf("1\n");

        /*error = one_jacobi_iteration(xLeft, yBottom,
                                     n+2, m+2,
                                     u_old, u,
                                     deltaX, deltaY,
                                     alpha, relax,
                                     &JIV_cx, &JIV_cy, &JIV_cc, neighbors);*/

        //printf("\tError %g\n", error);
        iterationCount++;
        // Swap the buffers
        tmp = u_old;
        u_old = u;
        u = tmp;

        //printf("2\n");
//        MPI_Wait(SRequests, MPI_STATUS_IGNORE); 
        //printf("rank %d wait for %d\n",rank,requestsCount);
        MPI_Waitall(requestsCount, SRequests, MPI_STATUS_IGNORE);
        /*if(neighbors[UP] >= 0 ){
           MPI_Wait(SRequests + UP, MPI_STATUS_IGNORE); 
        }*/
        
        //printf("3\n");
        /*if(neighbors[DOWN] >= 0 ){
            MPI_Wait(SRequests + DOWN, MPI_STATUS_IGNORE);
        }*/
        
       //printf("4\n");
       write_table(u,size,rank,neighbors);
       break;
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

    free(neighbors);
    free(coordinates);
    return 0;
}
