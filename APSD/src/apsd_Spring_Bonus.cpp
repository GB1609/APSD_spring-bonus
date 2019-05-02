#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include "Excercise1.h"
#include "Excercise2.h"
#include "Excercise3.h"
#include "Excercise5.h"

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
	dim = 100000000;
	toLaunch = 3;
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
			Excercise3 ex3(dim,num_thread);
			ex3.execute();
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
