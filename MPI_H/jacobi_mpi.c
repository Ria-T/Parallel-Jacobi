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
#include <omp.h>

#include "partitioning.h"

// fy, fx values are standard for all itterations. Though they are recallculated at each itteration/jacobi_itteration calling.
// Having an arrey with all the values pre-calculated once saves a lot of computing time while itterating/executing jacobi_itteration.
// The following function does just the above, allocating two arreys and saving in them all the posible needed fX, fY values.
int calculate_fX_fY_arreys(double xStart, double yStart, int n, int m, double deltaX, double deltaY, double **fX, double **fY){
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
 * 
 * NOTE FOR EXAMINATOR: In the mpi program there are two jacobi_iteration functions, 
 * the first caclulates the "white" part of the table, and the second (one_halo_jacobi_iteration) the halos.
 *************************************************************/
inline double one_jacobi_iteration(int maxXCount, int maxYCount,
                            double *src, double *dst,
                            double alpha, double omega,
                            double const * const cx, double const * const cy, double const * const cc,
                            double *fX, double *fY)
{
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]
    int x, y;
    //double fX, fY; As mentioned above we now use the tables, to avoid calculating the same values multiple times
    double error = 0.0;
    double updateVal;
    double f;
    // Coefficients, are passed by refernece and are calculated once outside of the calculation loop
    /*double cx = 1.0/(deltaX*deltaX);
    double cy = 1.0/(deltaY*deltaY);
    double cc = -2.0*cx-2.0*cy-alpha;*/

    #pragma omp parallel default(none) \
    shared(maxXCount, maxYCount, src, dst, omega, alpha, cx, cy, cc, fX, fY) \
    private(x, y, f, updateVal) \
    reduction (+: error)
    {
        #pragma omp for collapse(2) schedule(static)
        for (y = 2; y < (maxYCount-2); y++){
            for (x = 2; x < (maxXCount-2); x++){
                f = -alpha*(1.0-fX[x-1]*fX[x-1])*(1.0-fY[y-1]*fY[y-1]) - 2.0*(1.0-fX[x-1]*fX[x-1]) - 2.0*(1.0-fY[y-1]*fY[y-1]);

                //#pragma omp critical
                updateVal = (	(SRC(x-1,y) + SRC(x+1,y))*(*cx) +
                                (SRC(x,y-1) + SRC(x,y+1))*(*cy) +
                                SRC(x,y)*(*cc) - f
                            )/(*cc);

                //#pragma omp critical
                DST(x,y) = SRC(x,y) - omega*updateVal;

                //#pragma omp atomic
                error += updateVal*updateVal;
            }
        }
    }

    // The error is simply returned since the exprasion to find this itteration error is done in the main loop by process 0 for all processes
    return error;
    //sqrt(error)/((maxXCount-2)*(maxYCount-2));
}

// This function is the same as the above but, as mentioned above, it calculates the halo values of the local table
inline double one_halo_jacobi_iteration(int maxXCount, int maxYCount,
                            double *src, double *dst,
                            double alpha, double omega,
                            double const * const cx, double const * const cy, double const * const cc,
                            double *fX, double *fY)
{
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]
    int x, y;
    //double fX, fY;
    double error = 0.0;
    double updateVal;
    double f;
    // Coefficients
    /*double cx = 1.0/(deltaX*deltaX);
    double cy = 1.0/(deltaY*deltaY);
    double cc = -2.0*cx-2.0*cy-alpha;*/

    // There are two for loops, one for rows, and one for columns

    #pragma omp parallel default(none) \
    shared(maxXCount, maxYCount, src, dst, omega, alpha, cx, cy, cc, fX, fY) \
    private(x, y, f, updateVal) \
    reduction (+: error)
    {

        #pragma omp for collapse(2) schedule(static)
        for(y=1; y<maxYCount-1; y+=maxYCount-3){
            for(x=1; x<maxXCount-1; x++){
                f = -alpha*(1.0-fX[x-1]*fX[x-1])*(1.0-fY[y-1]*fY[y-1]) - 2.0*(1.0-fX[x-1]*fX[x-1]) - 2.0*(1.0-fY[y-1]*fY[y-1]);

                //#pragma omp critical
                updateVal = (	(SRC(x-1,y) + SRC(x+1,y))*(*cx) +
                                (SRC(x,y-1) + SRC(x,y+1))*(*cy) +
                                SRC(x,y)*(*cc) - f
                            )/(*cc);

                //#pragma omp critical
                DST(x,y) = SRC(x,y) - omega*updateVal;

                error += updateVal*updateVal;
            }
        }

        #pragma omp for collapse(2) schedule(static)
        for(y=2; y<maxYCount-2; y++){
            for(x=1; x<maxXCount-1; x+=maxXCount-3){
                f = -alpha*(1.0-fX[x-1]*fX[x-1])*(1.0-fY[y-1]*fY[y-1]) - 2.0*(1.0-fX[x-1]*fX[x-1]) - 2.0*(1.0-fY[y-1]*fY[y-1]);

                //#pragma omp critical
                updateVal = (	(SRC(x-1,y) + SRC(x+1,y))*(*cx) +
                                (SRC(x,y-1) + SRC(x,y+1))*(*cy) +
                                SRC(x,y)*(*cc) - f
                            )/(*cc);

                //#pragma omp critical
                DST(x,y) = SRC(x,y) - omega*updateVal;

                error += updateVal*updateVal;
            }
        }

    }

    return error;//sqrt(error)/((maxXCount-2)*(maxYCount-2));
}



/**********************************************************
 * Checks the error between numerical and exact solutions
 * 
 * NOTE FOR EXAMINATOR: The only change is that (like in the once_itteration_jacobi), the error value is returned,
 * and the expresion is calculated once for all processes by process 0
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
    return error;//sqrt(error)/((maxXCount-2)*(maxYCount-2));
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

    // Proccess 0 reads the input, other processes get it from proccess 0 after it reads it and bcasts it
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&relax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mits, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // The size of the whole table are needed and saved here, n,m will become the size of the local table in the following lines
    actual_n = n;
    actual_m = m;

    MPI_Comm cart_comm;
    int *neighbors, **coordinates, *topology_dims;
    get_local_table(&n, &m, &topology_dims, &coordinates, &rank, world_size,&cart_comm);
    n = coordinates[1][0]-coordinates[0][0]+1;
    m = coordinates[1][1]-coordinates[0][1]+1;
    //printf("%d:%dx%d\n",rank,n,m);
    neighbors = get_neighbors(&cart_comm);

    allocCount = (n+2)*(m+2);
    // Those two calls also zero the boundary elements
    u = 	(double*)calloc(allocCount, sizeof(double)); //reverse order
    u_old = (double*)calloc(allocCount, sizeof(double));
    
    //printf("allocCount=%d u=%p u_old=%p\n", allocCount, u, u_old);
    
    if (u == NULL || u_old == NULL)
    {
        printf("Not enough memory for two %ix%i matrices\n", n+2, m+2);
        exit(1);
    }
    maxIterationCount = mits;
    maxAcceptableError = tol;

    // Solve in [-1, 1] x [-1, 1]
    double actual_xLeft = -1.0, actual_xRight = 1.0;
    double actual_yBottom = -1.0, actual_yUp = 1.0;

    double xLeft = -1.0, xRight = 1.0;
    double yBottom = -1.0, yUp = 1.0;

    card_cords = malloc(sizeof(int)*2);
    MPI_Cart_coords(cart_comm, rank, 2, card_cords);

    double deltaX = (actual_xRight-actual_xLeft)/(actual_n-1);//(n-1);
    double deltaY = (actual_yUp-actual_yBottom)/(actual_m-1);//(m-1);

    // Calculating the local xLeft,xRight,yBottom,yUp values for the size of the local table
    calculate_xy_range(actual_xLeft,actual_xRight, &xLeft, &xRight, deltaX, coordinates[0][0], coordinates[1][0]);
    calculate_xy_range(actual_yBottom,actual_yUp, &yBottom, &yUp, deltaY, coordinates[0][1], coordinates[0][0]);

    //printf("(%d:%d,%d) out of [-1,1] I get [%f,%f] & [%f,%f]\n",rank,card_cords[0],card_cords[1],xLeft,xRight,yBottom,yUp);

    free(card_cords);

    iterationCount = 0;
    error = HUGE_VAL;
    clock_t start = clock(), diff;
    

    // Jacobi iteration Coefficient variables (values are calculated once -> faster execution of inline function)...
    double JIV_cx = 1.0/(deltaX*deltaX);
    double JIV_cy = 1.0/(deltaY*deltaY);
    double JIV_cc = -2.0*JIV_cx-2.0*JIV_cy-alpha;
    double *JIV_fX = NULL, *JIV_fY = NULL;

    //...The same goes for the jacobi iteration fX, fY variables
    if( calculate_fX_fY_arreys(xLeft, yBottom, n, m, deltaX, deltaY, &JIV_fX, &JIV_fY) != 0 )
        printf("Not enough memory to calculate fY & fX values!\n");

    MPI_Request SRequests[4],RRequests[4];
    int RrequestsCount = 0,SrequestsCount = 0;
    MPI_Datatype table_row, table_column;

    // Creating two MPI types, one for the row and the column of tha small/"local" table.
    MPI_Type_vector(n, 1, 1, MPI_DOUBLE, &table_row);
    MPI_Type_commit(&table_row);
    MPI_Type_vector(m, 1, n+2, MPI_DOUBLE, &table_column);
    MPI_Type_commit(&table_column);

    t1 = MPI_Wtime();

    /* Iterate as long as it takes to meet the convergence criterion */
    while (iterationCount < maxIterationCount && error > maxAcceptableError)
    {
        //printf("Iteration %i\n", iterationCount);

        RrequestsCount = 0;
        SrequestsCount = 0;

        // First recieving table data from all neighbors (or simply those that exist)
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


        // And then sending table data to all neighbors (or simply those that exist)
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

        // First calling one_jacobi_iteration, to caluculate only the white part of the table
        error = one_jacobi_iteration(n+2, m+2,
                                     u_old, u,
                                     alpha, relax,
                                     &JIV_cx, &JIV_cy, &JIV_cc, 
                                     JIV_fX, JIV_fY);

        // Befora we calculate the halo values, we have to wait for all sends to finish
        MPI_Waitall(RrequestsCount, RRequests, MPI_STATUS_IGNORE);

        // And then calculating the values for the halos with the according function
        error += one_halo_jacobi_iteration(n+2, m+2,
                                    u_old, u,
                                    alpha, relax,
                                    &JIV_cx, &JIV_cy, &JIV_cc, 
                                    JIV_fX, JIV_fY);

        //printf("\tError %g\n", error);
        iterationCount++;
        // Swap the buffers
        tmp = u_old;
        u_old = u;
        u = tmp;

        // In order to find the error for this itteration we first have to collect all the error values
        MPI_Reduce(&error, &globalError, 1, MPI_DOUBLE, MPI_SUM, 0, cart_comm);

        // Procces 0 calulates the error for all procceses
        if(rank == 0){
            error = sqrt(globalError)/(actual_n*actual_m);
            //printf("\tsqrt(%f)/(%d*%d)=%g\n",globalError,actual_n,actual_m,error);
            //printf("\tError %g\n", error);
        }
        // And then broadcasts the resuled error to all the other processes
        MPI_Bcast(&error, 1, MPI_DOUBLE, 0, cart_comm);

        // Finaly we wait for all the sending of data to complete
        MPI_Waitall(SrequestsCount, SRequests, MPI_STATUS_IGNORE);
    }

    t2 = MPI_Wtime();

    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;

    free(u);
    free(JIV_fX); free(JIV_fY);

    int *iterationCounts = NULL;
    if( rank == 0 ){
        iterationCounts = calloc( world_size, sizeof(int) );
        // Procces 0 prints the total number of itterations
        printf( "Iterations=%3d Elapsed MPI Wall time is %f\n", iterationCount, t2 - t1 ); 
    }
    // And gathers how many itterations the other procceses have done in order to print any anomalies
    MPI_Gather(&iterationCount, 1, MPI_INT, iterationCounts, 1, MPI_INT, 0, cart_comm);

    // Each procces caluculates the absoluteError of it's local table
    double absoluteError = checkSolution(xLeft, yBottom,
                                        n+2, m+2,
                                        u_old,
                                        deltaX, deltaY,
                                        alpha);

    double total_absoluteError;
    // And procces 0 collects all the error values, to calculate the total error and print it
    MPI_Reduce(&absoluteError, &total_absoluteError, 1, MPI_DOUBLE, MPI_SUM, 0, cart_comm);

    if(rank == 0){
        // Printing if any procceses did a different of itterations
        for(int i = 0; i<world_size; i++)
            if(iterationCount != iterationCounts[i])
                printf("Process %d specificaly did %d itterations\n",i,iterationCounts[i]);
        
        free(iterationCounts);

        printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
        printf("Residual %g\n",error);
        // u_old holds the solution after the most recent buffers swap
        total_absoluteError = sqrt(total_absoluteError)/(actual_n*actual_m);
        printf("The error of the iterative solution is %g\n", total_absoluteError);
    }

    free(u_old);
    free(neighbors);
    free(coordinates);
    free(topology_dims);
    MPI_Finalize();
    return 0;
}