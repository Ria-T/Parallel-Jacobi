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

    for (y = 2; y < (maxYCount-2); y++)
    {
        fY = yStart + (y-1)*deltaY;
        for (x = 2; x < (maxXCount-2); x++)
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
    return error;//sqrt(error)/((maxXCount-2)*(maxYCount-2));
}

/*inline*/ double one_halo_jacobi_iteration(double xStart, double yStart,
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

    for(y=1; y<maxYCount-1; y+=maxYCount-3){
        fY = yStart + (y-1)*deltaY;
        for(x=1; x<maxXCount-1; x++){
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

    for(y=2; y<maxYCount-2; y++){
        fY = yStart + (y-1)*deltaY;
        for(x=1; x<maxXCount-1; x+=maxXCount-3){
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

    return error;//sqrt(error)/((maxXCount-2)*(maxYCount-2));
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
    int n, m, actual_n, actual_m, mits, *card_cords;
    double alpha, tol, relax;
    double maxAcceptableError;
    double error, globalError;
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
    printf("-> %d, %d, %g, %g, %g, %d\n", n, m, alpha, relax, tol, mits);
    }

    // Proccess 0 reads the input, the rest get the input after proccess 0 bcasts it
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&relax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mits, 1, MPI_INT, 0, MPI_COMM_WORLD);

    actual_n = n;
    actual_m = m;

    MPI_Comm cart_comm;
    int size[2], *neighbors, **coordinates, *topology_dims;
    get_local_table(&n, &m, &topology_dims, &coordinates, &rank, world_size,&cart_comm);
    n = coordinates[1][0]-coordinates[0][0]+1;
    m = coordinates[1][1]-coordinates[0][1]+1;
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

    card_cords = malloc(sizeof(int)*2);
    MPI_Cart_coords(cart_comm, rank, 2, card_cords);

    calculate_range(&xLeft, &xRight, card_cords[1], topology_dims[0]);
    calculate_range(&yBottom, &yUp, card_cords[0], topology_dims[1]);
    if(rank==0){
        xLeft = -1.0;
        xRight = 1.0;
        yBottom = 0.0;
        yUp = 1.0;
    }else{
        xLeft = -1.0;
        xRight = 1.0;
        yBottom = -1.0;
        yUp = 0.0;
    }
    printf("(%d:%d,%d) out of [-1,1] I get [%f,%f] & [%f,%f]\n",rank,card_cords[0],card_cords[1],xLeft,xRight,yBottom,yUp);

    free(topology_dims);
    free(card_cords);

    double deltaX = (xRight-xLeft)/(actual_n-1);//(n-1);
    double deltaY = (yUp-yBottom)/(actual_m-1);//(m-1);

    iterationCount = 0;
    error = HUGE_VAL;
    clock_t start = clock(), diff;
    

    // Jacobi iteration Coefficient variables (values are calculated once -> faster execution of inline function)
    double JIV_cx = 1.0/(deltaX*deltaX);
    double JIV_cy = 1.0/(deltaY*deltaY);
    double JIV_cc = -2.0*JIV_cx-2.0*JIV_cy-alpha;

    MPI_Request SRequests[4],RRequests[4];
    int RrequestsCount = 0,SrequestsCount = 0;
    MPI_Datatype table_row, table_column;
    MPI_Type_vector(n, 1, 1, MPI_DOUBLE, &table_row);
    MPI_Type_commit(&table_row);
    MPI_Type_vector(m, 1, n+2, MPI_DOUBLE, &table_column);
    MPI_Type_commit(&table_column);

    //MPI_Init(NULL,NULL); should it get up?
    t1 = MPI_Wtime();

    //init_debug();

    /* Iterate as long as it takes to meet the convergence criterion */
    int kkk=2;
    while (iterationCount < maxIterationCount && error > maxAcceptableError && kkk > 0)
    {
        kkk--;
        //printf("Iteration %i\n", iterationCount);

        #ifdef DEBUG

        if(rank > -1){
            int l,kk;
	
			// second row
            for(l=1; l<n+1; l++){
                u_old[1*(n+2)+l]=0.1+rank+((float)( (l-1)%10 )-1.0)/100.0;
            }
			
			// pre last row
            for(l=2; l<n+1; l++){
                u_old[(m)*(n+2)+l]=0.2+rank+((float)( (l-2)%10 )-2.0)/100.0;
            }
			
			// left column + 1
			for(l=1; l<m+1; l++){
                u_old[l*(n+2)+1]=0.3+rank+((float)( (l-1)%10 )-1.0)/100.0;
            }
			
			// right column - 1
			for(l=1; l<m+1; l++){
                u_old[(l)*(n+2)+n]=0.4+rank+((float)( (l-1)%10 )-1.0)/100.0;
            }


            //center
            for(l=2; l<n; l++){
                for(kk=2; kk<m; kk++){
                    u_old[(kk)*(n+2)+l]=0.91;
                }
            }

        }
        write_table(u_old,size,rank,neighbors);
        #endif

        RrequestsCount = 0;
        SrequestsCount = 0;

        if(neighbors[UP] >= 0){
            MPI_Irecv(&u_old[1], 1, table_row, neighbors[UP], 0, cart_comm, &RRequests[RrequestsCount]);
            RrequestsCount++;
        }

        if(neighbors[DOWN] >= 0){
            MPI_Irecv(&u_old[(m+1)*(n+2)+1], 1, table_row, neighbors[DOWN], 0, cart_comm, &RRequests[RrequestsCount]);
            RrequestsCount++;
        }

        if(neighbors[LEFT] >= 0){
            MPI_Irecv(&u_old[(n+2)], 1, table_column, neighbors[LEFT], 0, cart_comm,&RRequests[RrequestsCount]);
            RrequestsCount++;
        }

        if(neighbors[RIGHT] >= 0){
            MPI_Irecv(&u_old[2*n+3], 1, table_column, neighbors[RIGHT], 0, cart_comm,&RRequests[RrequestsCount]);
            RrequestsCount++;
        }

        if(neighbors[UP] >= 0){
            MPI_Isend(&u_old[(n+2)+1], 1, table_row, neighbors[UP], 0, cart_comm, &SRequests[SrequestsCount]);
            SrequestsCount++;
        }

        if(neighbors[DOWN] >= 0){
            MPI_Isend(&u_old[(m)*(n+2)+1], 1, table_row, neighbors[DOWN], 0, cart_comm, &SRequests[SrequestsCount]);
            SrequestsCount++;
        }

        if(neighbors[LEFT] >= 0){
            MPI_Isend(&u_old[(n+2)+1], 1, table_column, neighbors[LEFT], 0, cart_comm,&SRequests[SrequestsCount]);
            SrequestsCount++;
        }

        if(neighbors[RIGHT] >= 0){
            MPI_Isend(&u_old[(n+2)+n], 1, table_column, neighbors[RIGHT], 0, cart_comm,&SRequests[SrequestsCount]);
            SrequestsCount++;
        }

        //printf("1\n");

        #ifdef DEBUG

        #else
        error = one_jacobi_iteration(xLeft, yBottom,
                                     n+2, m+2,
                                     u_old, u,
                                     deltaX, deltaY,
                                     alpha, relax,
                                     &JIV_cx, &JIV_cy, &JIV_cc);
        #endif

//        MPI_Wait(RRequests, MPI_STATUS_IGNORE); 
        MPI_Waitall(RrequestsCount, RRequests, MPI_STATUS_IGNORE);
        /*if(neighbors[UP] >= 0 ){
            MPI_Wait(RRequests + UP, MPI_STATUS_IGNORE); 
        }
        if(neighbors[DOWN] >= 0 ){
            MPI_Wait(RRequests + DOWN, MPI_STATUS_IGNORE);
        }*/

        #ifdef DEBUG

        #else
        error += one_halo_jacobi_iteration(xLeft, yBottom,
                                    n+2, m+2,
                                    u_old, u,
                                    deltaX, deltaY,
                                    alpha, relax,
                                    &JIV_cx, &JIV_cy, &JIV_cc);
        #endif

        //printf("\tError %g\n", error);
        iterationCount++;
        // Swap the buffers
        tmp = u_old;
        u_old = u;
        u = tmp;

        printf("%d: reducing %f\n",rank,error);
        MPI_Reduce(&error, &globalError, 1, MPI_DOUBLE, MPI_SUM, 0, cart_comm);

        if(rank == 0){
            error = sqrt(globalError)/(actual_n*actual_m);
            printf("\tsqrt(%f)/(%d*%d)=%g\n",globalError,actual_n,actual_m,error);
            //printf("\tError %g\n", error);
        }
        MPI_Bcast(&error, 1, MPI_DOUBLE, 0, cart_comm);

        //printf("2\n");
//        MPI_Wait(SRequests, MPI_STATUS_IGNORE); 
        //printf("rank %d wait for %d\n",rank,requestsCount);
        MPI_Waitall(SrequestsCount, SRequests, MPI_STATUS_IGNORE);
        /*if(neighbors[UP] >= 0 ){
           MPI_Wait(SRequests + UP, MPI_STATUS_IGNORE); 
        }*/
        
        //printf("3\n");
        /*if(neighbors[DOWN] >= 0 ){
            MPI_Wait(SRequests + DOWN, MPI_STATUS_IGNORE);
        }*/
        
       //printf("4\n");
       #ifdef DEBUG
       write_table(u,size,rank,neighbors);
       break;
       #endif
       write_table(u,size,rank,neighbors);
    }

    t2 = MPI_Wtime();
    printf( "Iterations=%3d Elapsed MPI Wall time is %f\n", iterationCount, t2 - t1 ); 
    MPI_Finalize();

    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;

    #ifdef DEBUG

    #else
    write_table(u_old,size,rank,neighbors);
    #endif

    if(rank == 0){
        printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
        printf("Residual %g\n",error);

        // u_old holds the solution after the most recent buffers swap
        double absoluteError = checkSolution(xLeft, yBottom,
                                            n+2, m+2,
                                            u_old,
                                            deltaX, deltaY,
                                            alpha);
        printf("The error of the iterative solution is %g\n", absoluteError);
    }

    free(neighbors);
    free(coordinates);
    return 0;
}
