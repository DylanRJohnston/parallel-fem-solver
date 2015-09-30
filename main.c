#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/time.h>

#define MAT_SIZE 100000000
#define NUM_THREADS 4
#define WORK_SIZE (MAT_SIZE / NUM_THREADS)
#define REDUCED_SIZE (2 * NUM_THREADS - 2)

double* tridiagonal_solve(const double* const a, const double* const b, const double* const c, const double* const r, int n) {
    double* cTmp = calloc(n,  sizeof(double));
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

void solve(
    const double* const a,
    const double* const b,
    const double* const c,
    const double* const r,
    double* s,
    int size,
    int pid,
    double* reducedA,
    double* reducedB,
    double* reducedC,
    double* reducedR
) {

    double* w = calloc(size, sizeof(double));
    w[0] = c[0] / b[0];
    for (int i = 1; i < size; ++i) w[i] = c[i] / (b[i] - a[i] * w[i - 1]);

    double* y = calloc(size, sizeof(double));
    y[0] = r[0] / b[0];
    for (int i = 1; i < size; ++i) y[i] = (r[i] - a[i] * y[i - 1]) / (b[i] - a[i] * w[i - 1]);

    double* xR  = calloc(size, sizeof(double));
    xR[size-1] = y[size-1];
    for (int i = size - 2; i >= 0; --i) xR[i] = y[i] - w[i] * xR[i+1];

    double* xLH = calloc(size, sizeof(double));
    xLH[size-1] = -w[size-1];
    for (int i = size - 2; i >= 0; --i) xLH[i] = -w[i] * xLH[i+1];

    double* wUH = calloc(size, sizeof(double));
    wUH[size - 1] = a[size-1] / b[size-1];
    for (int i = size - 2; i >= 0; --i) wUH[i] = a[i] / (b[i] - c[i] * wUH[i+1]);

    double* xUH = calloc(size, sizeof(double));
    xUH[0] = -wUH[0];
    for (int i = 1; i < size; ++i) xUH[i] = -wUH[i] * xUH[i-1];

    int lowerIndex = pid * 2-1;
    if (lowerIndex >= 0) {
        reducedA[lowerIndex] = -1;
        reducedB[lowerIndex] = xUH[0];
        reducedC[lowerIndex] = xLH[0];
        reducedR[lowerIndex] = -xR[0];
    }

    int upperIndex = pid * 2;
    if (upperIndex < REDUCED_SIZE) {
        reducedA[upperIndex] = xUH[size - 1];
        reducedB[upperIndex] = xLH[size - 1];
        reducedC[upperIndex] = -1;
        reducedR[upperIndex] = -xR[size - 1];
    }

    #pragma omp barrier
    double* reduced_solution = tridiagonal_solve(reducedA, reducedB, reducedC, reducedR, REDUCED_SIZE);

    double coLH = pid != NUM_THREADS - 1 ? reduced_solution[pid*2]     : 0;
    double coUH = pid != 0         ? reduced_solution[pid*2 - 1] : 0;

    for (int i = 0; i < size; i++) {
        s[i + size*pid] = xR[i] + coLH*xLH[i] + coUH*xUH[i];
    }
}

int main(void) {
    // const double a[] = {0.0,9.56926,9.96821,9.22169,2.86898,6.25015,1.40662,4.10846,8.39917,5.04963};
    // const double b[] = {9.91092,8.11347,1.99982,9.07977,3.11674,5.30049,2.55666,3.01966,5.73383,2.05164};
    // const double c[] = {4.17284,5.81033,5.52629,6.55847,2.22006,6.39581,1.12284,4.85695,3.41811,0.0};
    // const double r[] = {8.43678, 4.59904, 7.20084, 1.53088, 7.56609, 2.8395, 5.57885, 4.17973, 7.71119, 7.35723};
    double* a = calloc(MAT_SIZE, sizeof(double));
    double* b = calloc(MAT_SIZE, sizeof(double));
    double* c = calloc(MAT_SIZE, sizeof(double));
    double* r = calloc(MAT_SIZE, sizeof(double));
    double* s = calloc(MAT_SIZE, sizeof(double));

    double* reducedA = malloc(REDUCED_SIZE * sizeof(double));
    double* reducedB = malloc(REDUCED_SIZE * sizeof(double));
    double* reducedC = malloc(REDUCED_SIZE * sizeof(double));
    double* reducedR = malloc(REDUCED_SIZE * sizeof(double));

    srand(time(NULL));
    for (int i = 0; i < MAT_SIZE; ++i) {
        a[i] = 1.0 + (double) rand() / (RAND_MAX) * 9.0;
        b[i] = 1.0 + (double) rand() / (RAND_MAX) * 9.0;
        c[i] = 1.0 + (double) rand() / (RAND_MAX) * 9.0;
        r[i] = 1.0 + (double) rand() / (RAND_MAX) * 9.0;
    }

    a[0] = 0.0;
    c[MAT_SIZE-1] = 0.0;

    struct timeval start, diff;
    gettimeofday(&start, NULL);
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        int thread_num = omp_get_thread_num();
        int offset = thread_num * WORK_SIZE;
        solve(
            a+offset,
            b+offset,
            c+offset,
            r+offset,
            s,
            WORK_SIZE,
            thread_num,
            reducedA,
            reducedB,
            reducedC,
            reducedR
        );
    }
    gettimeofday(&diff, NULL);
    printf("Parallel takes %f\n", ((diff.tv_sec - start.tv_sec) * 1e6 + (diff.tv_usec - start.tv_usec))/1e6);

    gettimeofday(&start, NULL);
    double* other_solution = tridiagonal_solve(a, b, c, r, MAT_SIZE);
    gettimeofday(&diff, NULL);
    printf("Linear takes %f\n", ((diff.tv_sec - start.tv_sec) * 1e6 + (diff.tv_usec - start.tv_usec))/1e6);

    other_solution++;

    /*printf("Solution is \n");
    for (int i = 0; i < MAT_SIZE; ++i) {
        printf("%d: %f %f\n", i, s[i], other_solution[i]);
    }*/
}
