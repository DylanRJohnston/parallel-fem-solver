#include <stdio.h>
#include <stdlib.h>


/*
  
*/


const int a[] = {0, 1, 2, 3,  4,  5,  6,  7,  8, 9};
const int b[] = {1, 2, 3, 4,  5,  6,  7,  8,  9, 10};
const int c[] = {3, 5, 7, 9, 11, 13, 15, 17, 19, 0};

int main(int argc, char** argv) {

}

void forward_elimination() {
    xuh[1] = c[1] / b[1];
    xlh[1] = r[1] / b[1];

    for (int i = 2; i <= M; ++i) {
        denom = b[i] - c[i] * xuh[i-1];

        xuh[i] = c[i] / denom;
        xlh[i] = (r[i] - a[i] * xlh[i - 1]) / denom;
    }
}

void back_substitution() {
    xR[M] = xlh[M];
    xlh[M] = -xuh[M];
    xuh[M] = a[M]/b[M];

    for (int i = M-1; i >= 1; --i) {
        xR[i]  = xlh[i] - xuh[i] * xR[i+1];
        xlh[i] = -xuh[i] * xlh[i+1];
        denom = b[i] - c[i] * xuh[i+1];
        xuh[i] = -a[i] / denom;
    }
}

void forward_substitution() {
    xuh[1] = -xuh[1];

    for (int i = 2; i <= M; ++i) {
        xuh[i] = -xuh[i]*xuh[i-1];
    }
}
