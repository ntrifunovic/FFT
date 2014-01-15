#include "Maxfiles.h"
#include "MaxSLiCInterface.h"

#include <cstdio>
#include <cstdlib>
#include<cmath>
#include<cstdio>
#include<ctime>
#include<algorithm>

using namespace std;

#define PI 3.14159265359

// Okretanje bitova u broju
unsigned bit_reverse(register unsigned x) {
  x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
  x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
  x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
  x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));

  return((x >> 16) | (x << 16));
}

// Cooley-Tukey radix-2 inplace FFT algoritam
void fft(float real[], float imag[], unsigned n) {
  ///////////////////////////////////////////////
  //  Preuredjivanje niza u bit reverse order  //
  ///////////////////////////////////////////////
  unsigned log_n = __builtin_ctz(n);

  for(unsigned i = 0; i < n; i++) {
    unsigned j = bit_reverse(i) >> ( sizeof(i)*8 - log_n);

    if(i > j) {
      swap(real[i], real[j]);
      swap(imag[i], imag[j]);
    }
  }

  for(unsigned i = 1; i < n; i <<= 1) {
    for(unsigned k = 0; k < i; k++) {
	  float w_real = cos(-2*PI*k/(i<<1));
      float w_imag = sin(-2*PI*k/(i<<1));

      for(unsigned j = 0; j < n; j += i<<1) {
    	float temp_real = real[j+k+i]*w_real - imag[j+k+i]*w_imag;
    	float temp_imag = real[j+k+i]*w_imag + imag[j+k+i]*w_real;

        real[j+k+i] = real[j+k] - temp_real;
        imag[j+k+i] = imag[j+k] - temp_imag;

        real[j+k] += temp_real;
        imag[j+k] += temp_imag;
      }
    }
  }
}


//#define PRINT_OUTPUT

const unsigned long long max_broj_primera = 10000000;

int main(int argcnt, char* args[]) {

	unsigned long long broj_primera = argcnt == 2 ? atoll(args[1]) : 1;

	if(broj_primera > max_broj_primera) {
		printf("Maksimalni broj primera je %llu:\n", max_broj_primera);

		broj_primera = max_broj_primera;
	}

	printf("Broj primera: %llu\nDuzina sekvence: %u\n", broj_primera, FFT_n);

	static float real[max_broj_primera][FFT_n];
	static float imag[max_broj_primera][FFT_n];

	for(unsigned long long i = 0; i < broj_primera ; i++)
		for(unsigned long long j = 0; j < FFT_n; j++) {
			real[i][j] = rand()%10;
			imag[i][j] = rand()%10;
		}

	printf("Zavrseno genersianje random sekvenci\n");

#ifdef PRINT_OUTPUT
	for(unsigned long long k = 0; k < broj_primera; k++)
		for(unsigned long long i = 0; i < FFT_n; i++)
			printf("%.3f + %.3fj\n", real[k][i], imag[k][i]);
#endif

	static float o_real[max_broj_primera][FFT_n];
	static float o_imag[max_broj_primera][FFT_n];

	clock_t startTime;

	printf("DFE:\n");

	startTime = clock();


	FFT2_nonblock(broj_primera, (float *)o_imag, (float *)o_real);
	FFT(broj_primera, (float *)imag, (float *)real);

	printf("Vreme racunanja FFT-a: %f Normalizovano vreme: %e\n", ((double)(clock()-startTime))/CLOCKS_PER_SEC, ((double)(clock()-startTime))/CLOCKS_PER_SEC/broj_primera);

#ifdef PRINT_OUTPUT
	for(unsigned long long k = 0; k < broj_primera; k++)
		for(unsigned long long i = 0; i < FFT_n; i++)
			printf("%.3f + %.3fj\n", o_real[k][i], o_imag[k][i]);
#endif

	printf("CPU:\n");

	startTime = clock();

	for(unsigned long long k = 0; k < broj_primera; k++){
		fft(real[k], imag[k], FFT_n);
	}

	printf("Vreme racunanja FFT-a: %f Normalizovano vreme: %e\n", ((double)(clock()-startTime))/CLOCKS_PER_SEC, ((double)(clock()-startTime))/CLOCKS_PER_SEC/broj_primera);

#ifdef PRINT_OUTPUT
	for(unsigned long long k = 0; k < broj_primera; k++)
		for(unsigned long long i = 0; i < FFT_n; i++)
			printf("%.3f + %.3fj\n", real[k][i], imag[k][i]);
#endif

	for(unsigned long long k = 0; k < broj_primera; k++)
		for(unsigned long long i = 0; i < FFT_n; i++)
			if(real[k][i] - o_real[k][i] > 0.1e-10 || imag[k][i] - o_imag[k][i] > 0.1e-10) {
				return 1;
			}

	printf("Kraj");

	return 0;
}
