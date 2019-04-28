#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

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
void printMatrix(double **matrix, int dim) {
	printf("%d", dim);
	for (int a = 0; a < dim; a++) {
		for (int b = 0; b < dim; b++) {
			printf("%.8g \t", matrix[a][b]);
		}
		printf("\n");
	}
}
int checkMatrix(double **matrix1, double **matrix2, int dim, int numThreads) {

	int finalSum = 0, a, b;
#pragma omp parallel for shared(matrix1,matrix2) private(a,b) reduction(+:finalSum)
	for (a = 0; a < dim; a++)
		for (b = 0; b < dim; b++)
			if (matrix1[a][b] != matrix2[a][b])
				finalSum++;

	return finalSum;

}
void matrixGeneration2D(int dim, int numThreads) {
	double **matrixA = (double **) malloc(dim * sizeof(double));
	double **matrixB = (double **) malloc(dim * sizeof(double));
	for (int i = 0; i < dim; i++) {
		matrixA[i] = (double *) malloc(dim * sizeof(double));
		matrixB[i] = (double *) malloc(dim * sizeof(double));
	}
	int a, b;
	double time_serial, time_begin, time_serial_check;
	printf("***********************\nBEGIN SEQUENTIAL EXERCISE 1\n");
	time_begin = omp_get_wtime();
	for (a = 0; a < dim; a++)
		for (b = 0; b < dim; b++) {
			matrixA[a][b] = 5 * pow(a, 3) + 5 * M_PI * pow(b, 6);
			matrixB[a][b] = (10 / 3) * (matrixA[a][b]);
		}
//	printMatrix(matrixA, dim);
	time_serial = omp_get_wtime() - time_begin;
	printf(
			"END SEQUENTIAL EXERCISE 1\nTIME EXECUTION SEQ: %8g\n***********************",
			time_serial);
	printf("***********************\nBEGIN PARALLEL EXERCISE 1\n");
	time_begin = omp_get_wtime();
	double **matrixAParallel = (double **) malloc(dim * sizeof(double));
	double **matrixBParallel = (double **) malloc(dim * sizeof(double));
	for (int i = 0; i < dim; i++) {
		matrixAParallel[i] = (double *) malloc(dim * sizeof(double));
		matrixBParallel[i] = (double *) malloc(dim * sizeof(double));
	}
	omp_set_num_threads(numThreads);

#pragma omp parallel private(a, b)
	{
#pragma omp for schedule(static)
		for (a = 0; a < dim; a++) {
			for (b = 0; b < dim; b++) {
				matrixAParallel[a][b] = 5 * pow(a, 3) + 5 * M_PI * pow(b, 6);
				matrixBParallel[a][b] = (10 / 3) * (matrixA[a][b]);
			}
		}
	}
//	printMatrix(matrixA,dim);

	time_serial = omp_get_wtime() - time_begin;
	time_begin = omp_get_wtime();
	if (checkMatrix(matrixA, matrixAParallel, dim, numThreads) == 0
			&& checkMatrix(matrixB, matrixBParallel, dim, numThreads) == 0) {
		time_serial_check = omp_get_wtime() - time_begin;
		printf(
				"END PARALLEL EXERCISE 1\nTIME EXECUTION PARALLEL: %8g\n***********************",
				time_serial);
		printf("TIME CHECKING: %8g\n***********************",
				time_serial_check);
	} else
		printf("END PARALLEL EXERCISE 1 WITH ERROR\n\n***********************");
}

int main() {
	printf(
			"Insert number of bonus:\n 0 (for sum of Vector)\n 1 (for 2DMatrixGeneration) \n 2 (Exercise 2)\n 3 (Calculation of PI)\n 4 (for find an element in a vector)\n 5 (for Game of Life)\n ");
//	int toLaunch;
//	scanf("%d", &toLaunch);
	int numThreadMax = omp_get_max_threads();
	int numthread = numThreadMax;
//	switch (toLaunch) {
//	case 0:
//		printf("SUM OF VECTOR\n");
//		sumOfVector();
//		break;
//	case 1:
//		printf("Exercise one \n");

	matrixGeneration2D(10000, numthread);
//		break;
//	default:
//		printf("NUMBER NOT KNOW");
//	}
	return 0;
}
