#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <signal.h>
#include <execinfo.h>
#include <mpi.h>

//Signal handler, useful for debugging
#define UNUSED(x) (void)(x)
#define BACKTRACE_DEPTH 30
void signal_handler(int signal) {
    void* array[BACKTRACE_DEPTH];
    size_t size = backtrace(array, BACKTRACE_DEPTH);

    fprintf(stderr, "%s\n", strsignal(signal));
    backtrace_symbols_fd(array, size, 2);
    exit(1);
}

//Number of NODES / INTERVALS / SUBDIVISIONS in the solver
#define NNODE 101
#define NSUB (NNODE - 1)


/* Sets up the geometry of the problem depending on the Initial Boundary Conditions (IBC)
 * The right / left trim represent values at either end of the tridiagonal matrix that -
 * are not solved for since they are given by the IBC. Only the threads with ranks -
 * that are at either end of the solver need these values set. All other threads see 0.
 *
 * Mode 1, left = normal      right = derivative
 * Mode 2, left = derivative  right = normal
 * Mode 3, left = normal      right = normal
 * Mode 4, left = derivative  right = derivative
 */

void geometry(int ibc, int comm_rank, int comm_size, int* left_trim, int *right_trim) {
    switch (ibc) {
        case 1:
            *left_trim = 1;
            *right_trim = 0;
            break;
        case 2:
            *left_trim = 0;
            *right_trim = 1;
            break;
        case 3:
            *left_trim = 1;
            *right_trim = 1;
            break;
        case 4:
            *left_trim = 0;
            *right_trim = 0;
            break;
    }

    //If it's not the first process
    if (comm_rank != 0) *left_trim = 0;

    //If it's not the last process
    if (comm_rank != comm_size - 1) *right_trim = 0;
}


//These can be arbitrary functions of position for modelling the ODE
//In the example code they're just static values
static inline double ff(double x) { UNUSED(x); return 0.0;}
static inline double pp(double x) { UNUSED(x); return 1.0;}
static inline double qq(double x) { UNUSED(x); return 0.0;}

//Static lookup table needed in the following function to avoid uneeded branching
static const int sign[] = {1, -1, -1, 1};


/* This is the core computation of the assemble function. Computes the integral
 * over the interval for the linear basis functions. Factored out to simplify
 * the calculation at the workload boundary between threads.
*/
static inline void integrateBasis(
    double* aleft,
    double* adiag,
    double* aright,
    double* f,
    double ul,
    double ur,
    int i,
    int j,
    int ibc,
    double quadPosition,
    double interval_size,
    int working_size,
    int comm_rank,
    int comm_size
) {
    int node = i + j;

    //Phi is always 0.5 if the spacing between nodes is regular. Factored out.
    f[node] += interval_size * ff(quadPosition) * 0.5;

    //Since each thread's (node) is relative, detecting the overall boundary requires checking the rank too
    if (node == 0 && comm_rank == 0) {
        f[node] = f[node] - pp(0.0) * ul;
    } else if (node == working_size - 1 && comm_rank == comm_size - 1) {
        f[node] = f[node] + pp(1.0) * ur;
    }

    //For each basis function
    for (int m = 0; m < 2; ++m) {
        int node2  = i + m;
        double aij = interval_size * (qq(quadPosition) * 0.25 + sign[m + 2*j] * pp(quadPosition) / interval_size / interval_size);

        //Once again the thread relative variables (node) require checking the thread rank for overall boundary cases
        //Makes the conditions more complicated but simplifies the f[node], adiag[node], etc.
        if (
            node2 == 0
            && comm_rank == 0
            && (ibc == 1 || ibc == 3)
        ) {
                f[node] -= aij * ul;
        } else if (
            node2 == working_size - 1
            && comm_rank == comm_size - 1
            && (ibc == 2 || ibc == 3)
        ) {
                f[node] -= aij * ur;
        } else {
            if (m == j) {
                adiag[node] += aij;
            } else if (m < j) {
                aleft[node] += aij;
            } else {
                aright[node] += aij;
            }
        }
    }
}

/* The assemble function puts together the tridiagonal matrix for solving the
 * system. We factored out the inner most loop of the assemble function. Since each
 * interval affects the values of its neighbouring nodes, when splitting the
 * workload between seperate proccesses a slight overcalculation is needed at the
 * boundaries to compensate for the calculations carried out on other threads
 * this is required for each thread to be able to solve its own piece.
 */

void assemble(
    double* aleft,
    double* adiag,
    double* aright,
    double* f,
    double ul,
    double ur,
    double xl,
    double xr,
    int ibc,
    int thread_start,
    int working_size,
    int comm_rank,
    int comm_size
) {
    double interval_size = (xr - xl) / NSUB;

    //Iterate over the intervals (one less than #nodes)
    for (int i = 0; i < working_size - 1; ++i) {

        //While the iteration variables are relative, the quad position is global.
        double quadPosition = (i + thread_start + 0.5) * interval_size;

        //For the two nodes either side of the interval i.
        for (int j = 0; j < 2; ++j) {
            integrateBasis(aleft, adiag, aright, f, ul, ur, i, j, ibc, quadPosition, interval_size, working_size, comm_rank, comm_size);
        }
    }

    //If you're not the lowest rank, calculate the lower threads upper node
    if (comm_rank != 0) {
        integrateBasis(aleft, adiag, aright, f, ul, ur, -1, 1, ibc, (thread_start - 0.5) * interval_size, interval_size, working_size, comm_rank, comm_size);
    }

    //If you're not the highest rank, calculate the upper threads lower node
    if (comm_rank != comm_size - 1) {
        integrateBasis(aleft, adiag, aright, f, ul, ur, working_size - 1, 0, ibc, (working_size + thread_start - 0.5)*interval_size, interval_size, working_size, comm_rank, comm_size);
    }
}

/* A standard serial tridiagonal solver, needed to solve the reduced system on
 * each thread. The existing inplace solver is destructive to the inputs. This
 * one uses more memory but preserves the inputs.
 */
double* serial_solve(const double* const a, const double* const b, const double* const c, const double* const r, int n) {
    double* cTmp = calloc(n, sizeof(double));
    double* rTmp = calloc(n, sizeof(double));
    double* sol  = calloc(n, sizeof(double));

    cTmp[0] = c[0] / b[0];

    for (int i = 1; i < n; ++i) {
        cTmp[i] = c[i] / (b[i] - a[i]*cTmp[i-1]);
    }

    rTmp[0] = r[0] / b[0];

    for (int i = 1; i < n; ++i) {
        rTmp[i] = (r[i] - a[i] * rTmp[i-1]) / (b[i] - a[i] * cTmp[i-1]);
    }

    sol[n-1] = rTmp[n-1];

    for (int i = n - 2; i >= 0; --i) {
        sol[i] = rTmp[i] - cTmp[i] * sol[i+1];
    }

    free(cTmp);
    free(rTmp);

    return sol;
}

/* Solves a tridiagonal matrix in parallel. Each thread solves a subsection of
 * the problem, all the threads communicate a reduced matrix based on the
 * LU decomp values at the boundaries between sub problems, solving this reduced
 * system that only scales with the number of threads produces a correction to
 * the inital solution of the sub problem.
 * Based on this paper: http://www.mcs.anl.gov/~zippy/publications/partrid/partrid.html
*/
void parallel_solve(
    const double* const a,
    const double* const b,
    const double* const c,
    const double* const r,
    double* s,
    int size,
    int comm_rank,
    int comm_size
) {
    //TODO try to figure out a more inplace version of this algorithm.
    //TODO investigate why example in the paper doesn't work.
    double* w   = calloc(size, sizeof(double));
    double* y   = calloc(size, sizeof(double));
    double* xR  = calloc(size, sizeof(double));
    double* xLH = calloc(size, sizeof(double));
    double* wUH = calloc(size, sizeof(double));
    double* xUH = calloc(size, sizeof(double));

    w[0] = c[0] / b[0];
    for (int i = 1; i < size; ++i) w[i] = c[i] / (b[i] - a[i] * w[i - 1]);

    y[0] = r[0] / b[0];
    for (int i = 1; i < size; ++i) y[i] = (r[i] - a[i] * y[i - 1]) / (b[i] - a[i] * w[i - 1]);

    xR[size-1] = y[size-1];
    for (int i = size - 2; i >= 0; --i) xR[i] = y[i] - w[i] * xR[i+1];

    xLH[size-1] = -w[size-1];
    for (int i = size - 2; i >= 0; --i) xLH[i] = -w[i] * xLH[i+1];

    wUH[size - 1] = a[size-1] / b[size-1];
    for (int i = size - 2; i >= 0; --i) wUH[i] = a[i] / (b[i] - c[i] * wUH[i+1]);

    xUH[0] = -wUH[0];
    for (int i = 1; i < size; ++i) xUH[i] = -wUH[i] * xUH[i-1];

    //Setup the reduced global system
    //Should really by 2 * comm_size - 2, but the MPI_send would be more complicated.
    int reduced_size = 2 * comm_size;

    double* reducedA = calloc(reduced_size, sizeof(double));
    double* reducedB = calloc(reduced_size, sizeof(double));
    double* reducedC = calloc(reduced_size, sizeof(double));
    double* reducedR = calloc(reduced_size, sizeof(double));

    double tempA[] = {-1,     xUH[size - 1]};
    double tempB[] = {xUH[0], xLH[size - 1]};
    double tempC[] = {xLH[0], -1};
    double tempR[] = {-xR[0], -xR[size - 1]};

    //Each thread builds a picture of the reduced system.
    MPI_Allgather(tempA, 2, MPI_DOUBLE, reducedA, 2, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(tempB, 2, MPI_DOUBLE, reducedB, 2, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(tempC, 2, MPI_DOUBLE, reducedC, 2, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(tempR, 2, MPI_DOUBLE, reducedR, 2, MPI_DOUBLE, MPI_COMM_WORLD);

    //Solve the reduced system ignoring the boundary values at each end since they're not coupled with further threads.
    double* reduced_solution = serial_solve(reducedA + 1, reducedB + 1, reducedC + 1, reducedR + 1, 2 * comm_size - 2);

    //Not the highest rank
    double coLH = (comm_rank != comm_size - 1) ? reduced_solution[comm_rank * 2] : 0;

    //Not the lowest rank
    double coUH = (comm_rank != 0) ? reduced_solution[comm_rank * 2 - 1] : 0;

    //Correct initial solution with values from the reduced global system.
    for (int i = 0; i < size; i++) {
        s[i] = xR[i] + coLH*xLH[i] + coUH*xUH[i];
    }
}

FILE* fp = NULL;

int main(int argc, char** argv) {

    //Setup signal handlers
    signal(SIGSEGV, signal_handler);
    signal(SIGABRT, signal_handler);

    MPI_Init(&argc, &argv);

    //Get thread dimensions.
    int comm_rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    //Boundary values
    double ul = 0.0;
    double ur = 1.0;

    //Position of left and right boundaries
    double xl = 0.0;
    double xr = 1.0;

    //Boundary problem type
    int ibc = 1;

    //Bits trimmed off either side of the matrix depending on the IBC.
    int left_trim;
    int right_trim;
    geometry(ibc, comm_rank, comm_size, &left_trim, &right_trim);


    //Work out this threads workload based on the thread dimensions
    //TODO change workload breakup to be vectorizable.
    int thread_start = (comm_rank + 0) * NNODE / comm_size;
    int thread_end   = (comm_rank + 1) * NNODE / comm_size;
    int working_size = thread_end - thread_start;

    //TODO align this memory for vector instructions.
    double* aleft  = calloc(working_size, sizeof(double));
    double* adiag  = calloc(working_size, sizeof(double));
    double* aright = calloc(working_size, sizeof(double));
    double* f      = calloc(working_size, sizeof(double));
    double* s      = calloc(working_size, sizeof(double));

    //Setup timing for the compute functions.
    struct timeval start, end;
    gettimeofday(&start, NULL);

    //Assemble the system to be solved
    assemble(aleft, adiag, aright, f, ul, ur, xl, xr, ibc, thread_start, working_size, comm_rank, comm_size);

    //Solve the sytem
    parallel_solve(
        aleft  + left_trim,
        adiag  + left_trim,
        aright + left_trim,
        f      + left_trim,
        s      + left_trim,
        working_size - left_trim - right_trim,
        comm_rank,
        comm_size
    );

    gettimeofday(&end, NULL);
    double delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
             end.tv_usec - start.tv_usec) / 1.e6;
    printf("%d, %d, %f\n", comm_size, NSUB, delta);

    for (int i = 0; i < working_size; i++) {
        printf("%3d: %f\n", thread_start + i, s[i]);
    }

    MPI_Finalize();

    return 0;
}
