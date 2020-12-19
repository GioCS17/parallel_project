#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <complex>

using namespace std;

void BitReversal(int n, double xReal[], double xImg[], double log2_of_n);
void FillArray(int n, double xReal[], double xImg[], double yReal[], double yImg[]);
void CalcFFT(int log2_of_n, double xReal[], double xImg[], double wReal[], double wImg[], int N);
void GetTwiddleFactors(int n, int N, double wReal[], double wImg[]);
long long cont = 0;

int main(int argc, char** argv) {
    if (argc != 3) {
        cerr << "usage: FFT <size of array> <number of threads>" << endl;
        exit(1);
    }

    int n = atoi(argv[1]); //Assign size of array
    if (!(n > 0) || !((n & (n - 1)) == 0)) {
        cerr << "usage: <size of array> must be power of 2" << endl;
        exit(1);
    }

    int p = atoi(argv[2]); //Assign number of threads to use
    double log2_of_n = log10(n) / log10(2); //Set for use later
    double *xReal = new double[n];
    double *xImg = new double[n];
    double *yReal = new double[n];
    double *yImg = new double[n];
    double *wReal = new double[(int) log2_of_n];
    double *wImg = new double[(int) log2_of_n];
    double *W_OracleR = new double[n / 2];
    double *W_OracleI = new double[n / 2];
    omp_set_dynamic(0);
    omp_set_num_threads(p); //set number of threads

    FillArray(n, xReal, xImg, yReal, yImg);
    timespec sTime;
    double t1 = omp_get_wtime();
    BitReversal(n, xReal, xImg, log2_of_n);
    GetTwiddleFactors(log2_of_n, n, wReal, wImg);
    CalcFFT(log2_of_n, xReal, xImg, wReal, wImg, n);
    timespec fTime;
    double t2 = omp_get_wtime();
    cout<<"Num threads::"<<p<<std::endl;
    cout<<"Size::"<<n<<" - Time::"<<t2-t1<<std::endl;
    cout<<"FLOPS::"<<cont<<std::endl;
    cout<<"\n";

    return 0;
}

void BitReversal(int n, double xReal[], double xImg[], double log2_of_n) {
    int i;
    double temp;
    #pragma omp parallel for default(none) private(i, temp) shared(n, xReal, xImg, log2_of_n, cont)
    for (i = 1; i < n; i++) {
        int k, rev = 0;
        int inp = i;
        for (k = 0; k < log2_of_n; k++) {
            rev = (rev << 1) | (inp & 1);
            inp >>= 1;
	    cont = cont  + 4;
        }

        if (rev <= i) continue; //Skip if already done
        temp = xReal[i]; //Store into temp values and swap
        xReal[i] = xReal[rev];
        xReal[rev] = temp;
        temp = xImg[i];
        xImg[i] = xImg[rev];
        xImg[rev] = temp;
    }
}

void CalcFFT(int log2_of_n, double xReal[], double xImg[], double wReal[], double wImg[], int N) {
    int n, d, i, k;
    double temp_w, temp_x;
    for (n = 1; n < log2_of_n + 1; n++) {   
	//Iterate through time slice
        d = pow(2, n); // log(log(N)) flops 
        #pragma omp parallel for default(none) private(k, i, temp_w, temp_x) shared(n, N, d, xReal, xImg, wReal, wImg, log2_of_n, cont)
        for (k = 0; k < (d / 2); k++) {   
	//Iterate through even and odd elements
            for (i = k; i < N; i += d) { //Butterfly operation
                temp_w = wReal[n - 1] * xReal[i + (d / 2)]; //Multiply by twiddle factor
                temp_x = xReal[i];  //Store in temp's and restore with correct values
                xReal[i] = temp_w + temp_x;
                xReal[i + (d / 2)] = temp_x - temp_w;
                temp_w = wImg[n - 1] * xImg[i + (d / 2)]; //Multiply by twiddle factor
                temp_x = xImg[i]; //Store in temp's and restore with correct values
                xImg[i] = temp_w + temp_x;
                xImg[i + (d / 2)] = temp_x - temp_w;
		cont = cont + 10;
            }
        }
    }
}

void FillArray(int n, double xReal[], double xImg[], double yReal[], double yImg[]) {
    int i;
    struct drand48_data buffer;
    srand48_r(time(NULL) ^ omp_get_thread_num(), &buffer); //Seed for each thread
    //#pragma omp parallel for default(none) private(i) shared(x,n)
    for (i = 0; i < n; i++) {
        drand48_r(&buffer, &xReal[i]); //Store random double into xReal
        yReal[i] = xReal[i]; //Copy
        drand48_r(&buffer, &xImg[i]); //Store random double into xImg
        yImg[i] = xImg[i]; //Copy
    }
}

void GetTwiddleFactors(int log2_of_n, int n, double wReal[], double wImg[]) {
    int i;
#pragma omp parallel for default(none) private(i) shared(n, wReal, wImg, log2_of_n, cont)
    for (i = 0; i < log2_of_n; i++) {
        wReal[i] = cos(((double) i * 2.0 * M_PI) / ((double) n));
        wImg[i] = sin(((double) i * 2.0 * M_PI) / ((double) n));
	cont = cont + 8;
    }
}
