#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include "Excercise1.h"
#include "Excercise2.h"
#include "Excercise5.h"

double excercise3_serial(int dim) {
	double x, pi, sum = 0.0;
	double step = 1.0 / (double) dim;
	double begin_time = omp_get_wtime();
	for (int i = 0; i < dim; i++) {
		x = (i + 0.5) * step;
		sum += 4.0 / (1.0 + x * x);
	}
	pi = step * sum;
	double end_time = omp_get_wtime() - begin_time;
	printf("*SERIAL EXECUTION TERMINATED IN: %.8g\tVALUE OF PI= %.16g*\n",
	      end_time, pi);
	return end_time;
}
double excercise3_critical(int dim, int numThreads) {
	double pi = 0.0, begin_time, end_time;
	double step = 1.0 / (double) dim;
	omp_set_num_threads(numThreads);
	begin_time = omp_get_wtime();
#pragma omp parallel if(dim > 10000)
	{
		double x, sum;
#pragma omp for
		for (int i = 0; i < dim; i++) {
			x = (i + 0.5) * step;
			sum += 4.0 / (1.0 + x * x);
		}

#pragma omp critical
		pi += step * sum;
	}
	end_time = omp_get_wtime() - begin_time;
	printf("*CRITICAL EXECUTION TERMINATED IN: %.8g\tVALUE OF PI= %.16g*\n",
	      end_time, pi);
	return end_time;
}
double excercise3_padding(int dim, int numThreads) {
	double pi = 0.0, begin_time, end_time, x;
	double pi_padding[numThreads][8];
	double step = 1.0 / (double) dim;
	omp_set_num_threads(numThreads);
	begin_time = omp_get_wtime();
#pragma omp parallel if(dim > 10000)
	{
		int id = omp_get_thread_num();
#pragma omp for
		for (int i = 0; i < dim; i++) {
			x = (i + 0.5) * step;
			pi_padding[id][0] += 4.0 / (1.0 + x * x);
		}
	}
	for (int i = 0; i < numThreads; i++)
		pi += pi_padding[i][0] * step;

	end_time = omp_get_wtime() - begin_time;
	printf("*PADDING EXECUTION TERMINATED IN: %.8g\tVALUE OF PI= %.16g*\n",
	      end_time, pi);
	return end_time;
}
double excercise3_reduction(int dim, int numThreads) {
	int i;
	double x, pi, sum = 0.0;
	double step = 1.0 / (double) dim;
	omp_set_num_threads(numThreads);
	double begin_time = omp_get_wtime();
#pragma omp parallel for private(i) reduction (+:sum) if(dim > 10000)
	for (i = 0; i < dim; i++) {
		x = (i + 0.5) * step;
		sum += 4.0 / (1.0 + x * x);
	}
	pi = step * sum;
	double end_time = omp_get_wtime() - begin_time;
	printf("*REDUCTION EXECUTION TERMINATED IN: %.8g\tVALUE OF PI= %.16g*\n",
	      end_time, pi);
	return end_time;
}
double excercise3_serial_monte_carlo(int dim, int numThreads) {
	int i, nCirc = 0;
	double r = 1.0, pi, x, y;
	omp_set_num_threads(numThreads);
	unsigned int seed = (int) r * 2;
	double begin_time = omp_get_wtime();
	for (i = 0; i < dim; i++) {
		x = (double) rand_r(&seed);
		y = (double) rand_r(&seed) / RAND_MAX;
		if (pow(x, 2) + pow(y, 2) <= pow(r, 2))
			nCirc++;
	}
	pi = 4.0 * ((double) nCirc / (double) dim);
	double end_time = omp_get_wtime() - begin_time;
	printf(
	      "*SERIAL MONTE CARLO EXECUTION TERMINATED IN: %.8g\tVALUE OF PI= %.16g*\n",
	      end_time, pi);
	return end_time;
}
double excercise3_parallel_monte_carlo(int dim, int numThreads) {
	int i, nCirc = 0;
	double r = 1.0, pi, x, y;
	unsigned int seed = 0;
	double begin_time = omp_get_wtime();
#pragma omp parallel for private(x,y,i) reduction (+:nCirc) if(dim > 10000000)
	for (i = 0; i < dim; i++) {
		x = (double) rand_r(&seed) / RAND_MAX;
		y = (double) rand_r(&seed) / RAND_MAX;
		if (pow(x, 2) + pow(y, 2) <= pow(r, 2))
			nCirc++;
	}
	pi = 4.0 * ((double) nCirc / (double) dim);
	double end_time = omp_get_wtime() - begin_time;
	printf(
	      "*PARALLEL MONTE CARLO EXECUTION TERMINATED IN: %.8g\tVALUE OF PI= %.16g*\n",
	      end_time, pi);
	return end_time;
}
void excercise3(int dim, int numThreads) {
	printf("*********************\nBEGIN EXCERCISE 3\n");
	excercise3_serial(dim);
	excercise3_critical(dim, numThreads);
	excercise3_padding(dim, numThreads);
	excercise3_reduction(dim, numThreads);
	excercise3_serial_monte_carlo(dim, numThreads);
	excercise3_parallel_monte_carlo(dim, numThreads);
	printf("END EXCERCISE 3\n*********************\n");
}

void excercise4Reduction(int *vector, int sizeVector, int toSearch) {
	int parallel_find = 0, a;
	double time_begin = omp_get_wtime();
#pragma omp parallel for shared(vector) private(a) reduction (+:parallel_find)
	for (a = 0; a < sizeVector; a++)
		if (vector[a] == toSearch)
			parallel_find++;
	double time_end = omp_get_wtime() - time_begin;
	printf(
	      "***PARALLEL (reduction) -> Finished in: %.8g\tFOUNDED: %d times (SIZE OF VECTOR:%d)\n",
	      time_end, parallel_find, sizeVector);
}

void excercise4Atomic(int *vector, int sizeVector, int toSearch) {
	int parallel_find = 0;
	double time_begin = omp_get_wtime();
#pragma omp parallel
	{
#pragma omp  for
		for (int i = 0; i < sizeVector; i++) {
			if (vector[i] == toSearch)
#pragma omp atomic
				parallel_find++;
		}

	}
	double time_end = omp_get_wtime() - time_begin;
	printf(
	      "***PARALLEL (atomic) -> Finished in: %.8g\tFOUNDED: %d times (SIZE OF VECTOR:%d)\n",
	      time_end, parallel_find, sizeVector);
}

void excercise4(int dim, int numThreads) {
	printf("*********************\nBEGIN EXCERCISE 4\n");
	int a;
	double time_begin, time_end;
	srand(time(NULL));
	printf("START GENERATION VECTOR\n");
	int sizeVector = dim;
	int * vector = (int*) malloc(sizeVector * sizeof(int));
	for (a = 0; a < sizeVector; a++)
		vector[a] = rand() % 100;
	printf("END GENERATION VECTOR\n");
	int toSearch = rand() % 100;
	time_begin = omp_get_wtime();
	int find = 0;
	for (a = 0; a < sizeVector; a++)
		if (vector[a] == toSearch)
			find++;
	time_end = omp_get_wtime() - time_begin;
	printf(
	      "***SERIAL -> Finished in: %.8g\tFOUNDED: %d times (SIZE OF VECTOR:%d)\n",
	      time_end, find, sizeVector);

	excercise4Reduction(vector, sizeVector, toSearch);
	excercise4Atomic(vector, sizeVector, toSearch);
	free(vector);
	printf("END EXCERCISE 4\n*********************\n");
}

int main() {
	printf(
	      "Insert number of bonus:\n 0 (for sum of Vector)\n 1 (for 2DMatrixGeneration) \n 2 (Exercise 2)\n 3 (Calculation of PI)\n");
	printf("4 (for find an element in a vector)\n 5 (for Game of Life)\n ");
	printf(
	      "followed by number of thread and dimension of container in form: excercise-numThreads-dimension \n EXAMPLE: 3-6-10000\n");
	int toLaunch, num_thread = 0, dim;
	const int numberExcercise = 5;
	const int max_num_threads = omp_get_max_threads();
//	scanf("%d-%d-%d", &toLaunch, &num_thread, &dim);
//	while(toLaunch>numberExcercise || toLaunch<1 || max_num_threads < num_thread || num_thread==0)
//	{
//		printf("ERROR INPUT!!!\nEXAMPLE: 3-6-10000\n");
//		scanf("%d-%d-%d", &toLaunch, &num_thread, &dim);
//	}
	num_thread = max_num_threads - 3;
	dim = 10;
	toLaunch = 5;
	printf("THREADS USED: %d\n", num_thread);
	switch (toLaunch) {
		case 1: {
			printf("***Exercise 1*** \n");
			Excercise1 ex1(dim, num_thread);
			ex1.execute();
		}
			break;
		case 2: {
			printf("***Exercise 2***\n");
			Excercise2 ex2(dim, num_thread);
			ex2.execute();
		}
			break;
		case 3: {
			printf("***Exercise 3*** \n");
			excercise3(dim, num_thread);
		}
			break;
		case 4: {
			printf("***Exercise 4*** \n");
			excercise4(dim, num_thread);
		}
			break;
		case 5: {
			printf("***Exercise 5*** \n");
			Excercise5 ex5(dim, 20, num_thread);
			ex5.execute();
		}
			break;
		default:
			printf("NUMBER NOT KNOW");
	}
	printf("IRON MAN MUORE\n");
	return 0;
}
