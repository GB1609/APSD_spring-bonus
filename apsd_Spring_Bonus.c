#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <stdlib.h>

void sumOfVector() {
	srand(time(NULL));
	printf("START GENERATION VECTOR\n");
	int sizeVector = INT_MAX;
	int * vettore;
	int sommaSeriale = 0;
	vettore = (int*) malloc(sizeVector * sizeof(int));
	for (int a = 0; a < sizeVector; a++)
		vettore[a] = rand() % 10;
	printf("END GENERATION VECTOR\n");
	printf("BEGIN SERIAL\n");
	double time_begin = omp_get_wtime();
	for (int a = 0; a < sizeVector; a++)
		sommaSeriale += vettore[a];
	double time_serial = omp_get_wtime() - time_begin;
	printf("SERIAL TIME= %.8g \n SUM: %d\n", time_serial, sommaSeriale);
	int numThreadMax = omp_get_max_threads();
	int numthread = numThreadMax;
	printf("BEGIN PARALLEL\n");
	printf("Number of threads %d , USED: %d \n", numThreadMax, numthread);
	omp_set_num_threads(numthread);
	time_begin = omp_get_wtime();
	int final_sum = 0;
	int index;
#pragma omp parallel for shared(vettore) private(index) reduction(+:final_sum)
	for (index = 0; index < sizeVector; index++)
		final_sum += vettore[index];
	double time_parallel = omp_get_wtime() - time_begin;
	printf("TEMPO PARALLELO= %.8g \n SOMMA= %d \n", time_parallel, final_sum);
	printf("END SUM OF VECTOR \n");
}

void matrixGeneration2D() {

}

int main() {
	printf(
			"Insert number of bonus:\n 0 (for sum of Vector)\n 1 (for 2DMatrixGeneration) \n 2 (Exercise 2)\n 3 (Calculation of PI)\n 4 (for find an element in a vector)\n 5 (for Game of Life)\n ");
	int toLaunch;
	scanf("%d",&toLaunch);

	switch (toLaunch) {
	case 0:
		sumOfVector();
		break;
	case 1:
		matrixGeneration2D();
		break;
	default:
		printf("NUMBER NOT KNOW");
	}

	return 0;
}
