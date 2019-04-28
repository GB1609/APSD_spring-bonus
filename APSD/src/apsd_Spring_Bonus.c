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
void printMatrixEx1(double **matrix, int dim) {
	printf("%d", dim);
	for (int a = 0; a < dim; a++) {
		for (int b = 0; b < dim; b++) {
			printf("%.8g \t", matrix[a][b]);
		}
		printf("\n");
	}
}
int checkMatrixEx1(double **matrix1, double **matrix2, int dim, int numThreads) {

	int finalSum = 0, a, b;
#pragma omp parallel for shared(matrix1,matrix2) private(a,b) reduction(+:finalSum)
	for (a = 0; a < dim; a++)
		for (b = 0; b < dim; b++)
			if (matrix1[a][b] != matrix2[a][b])
				finalSum++;

	return finalSum;

}

int checkVectorsEx2(int *vetor1, int *vector2, int dim, int numThreads) {

	int finalSum = 0, a;
#pragma omp parallel for shared(vetor1,vector2) private(a) reduction(+:finalSum)
	for (a = 0; a < dim; a++)
		if (vetor1[a] != vector2[a])
			finalSum++;
	return finalSum;

}

void excerise1(int dim, int numThreads) {
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
	if (checkMatrixEx1(matrixA, matrixAParallel, dim, numThreads) == 0
			&& checkMatrixEx1(matrixB, matrixBParallel, dim, numThreads) == 0) {
		time_serial_check = omp_get_wtime() - time_begin;
		printf(
				"END PARALLEL EXERCISE 1\nTIME EXECUTION PARALLEL: %8g\n***********************",
				time_serial);
		printf("TIME CHECKING: %8g\n***********************",
				time_serial_check);
	} else
		printf("END PARALLEL EXERCISE 1 WITH ERROR\n\n***********************");
	free(matrixA);
	free(matrixB);
	free(matrixBParallel);
	free(matrixAParallel);
	printf("END EXCERCISE 1\n*********************\n");
}

void excercise2(int dim, int numThreads) {

	printf("*********************BEGIN EXCERCISE 2\n");
	int a;
	double time_begin, time_end;
	srand(time(NULL));
	printf("START GENERATION VECTOR\n");
	int sizeVector = dim;
	int * vectorA;
	int* vectorC;
	int* vectorCParallel;
	int* vectorB;
	vectorA = (int*) malloc(sizeVector * sizeof(int));
	vectorB = (int*) malloc(sizeVector * sizeof(int));
	vectorC = (int*) malloc(sizeVector * sizeof(int));
	vectorCParallel = (int*) malloc(sizeVector * sizeof(int));
	for (a = 0; a < sizeVector; a++) {
		vectorA[a] = rand() % 10;
		vectorB[a] = rand() % 10;
	}
	printf("END GENERATION VECTORS A & B\n");

	time_begin = omp_get_wtime();

	for (a = 0; a < sizeVector; a++) {
		vectorC[a] = vectorA[a] + vectorB[a];
	}
	int serial_sum = 0;
	for (a = 0; a < sizeVector; a++) {
		serial_sum += vectorC[a] * vectorB[a];
	}
	time_end = omp_get_wtime() - time_begin;
	printf(
			"*****TIME SERIAL EXECUTION= %.8g\tWITH SCALAR PRODUCT SERIAL= %d*****\n",
			time_end, serial_sum);
	omp_set_num_threads(numThreads);
	time_begin = omp_get_wtime();
	int parallel_sum = 0, partial_sum = 0, thread_n;
#pragma omp parallel for shared(vectorA,vectorB,vectorCParallel) private(a)
	for (a = 0; a < sizeVector; a++) {
		vectorCParallel[a] = vectorA[a] + vectorB[a];
	}
#pragma omp parallel private(a,thread_n,partial_sum)
	{
		partial_sum = 0;
		thread_n = omp_get_thread_num();

#pragma omp for reduction(+:parallel_sum)
		for (a = 0; a < sizeVector; a++) {
			parallel_sum += (vectorCParallel[a] * vectorB[a]);
			partial_sum = parallel_sum;
		}
	}
	time_end = omp_get_wtime() - time_begin;
	printf(
			"*****TIME PARALLEL EXECUTION= %.8g\tWITH SCALAR PRODUCT PARALLEL= %d*****\n",
			time_end, parallel_sum);

	printf("DIFFERENCE BEETWEN VECTORC AND VECTORCPARALLEL= %d\n",
			checkVectorsEx2(vectorC, vectorCParallel, sizeVector, numThreads));
	printf("DIFFERENCE BEETWEN SUM AND PARALLEL SUM= %d\n",
			parallel_sum - serial_sum);
	free(vectorA);
	free(vectorB);
	free(vectorC);
	free(vectorCParallel);
	printf("END EXCERCISE 2\n*********************\n");
}

void excercise3(int dim, int numThreads) {
	printf("*********************\nBEGIN EXCERCISE 3\n");
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
	int parallel_find = 0, a;
	double time_begin = omp_get_wtime();
#pragma omp parallel
	{
		int id = omp_get_thread_num();
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
void excercise5(int dim, int numThreads) {
	printf("*********************\nBEGIN EXCERCISE 3\n");
	printf("END EXCERCISE 3\n*********************\n");
}

int main() {
	printf(
			"Insert number of bonus:\n 0 (for sum of Vector)\n 1 (for 2DMatrixGeneration) \n 2 (Exercise 2)\n 3 (Calculation of PI)\n 4 (for find an element in a vector)\n 5 (for Game of Life)\n ");
	int toLaunch = 4;
//	scanf("%d", &toLaunch);
	int max_num_threads = omp_get_max_threads();
	int num_thread = max_num_threads;
	int dim = 10000;
	switch (toLaunch) {
	case 0:
		printf("SUM OF VECTOR\n");
		sumOfVector();
		break;
	case 1:
		printf("Exercise 1 \n");
		excerise1(dim, num_thread);
		break;
	case 2:
		printf("Exercise 2 \n");
		excercise2(dim, num_thread);
		break;
	case 3:
		printf("Exercise 3 \n");
		excercise3(dim, num_thread);
		break;
	case 4:
		printf("Exercise 4 \n");
		excercise4(dim, num_thread);
		break;
	case 5:
		printf("Exercise 5 \n");
		excercise5(dim, num_thread);
		break;
	default:
		printf("NUMBER NOT KNOW");
	}
	return 0;
}
