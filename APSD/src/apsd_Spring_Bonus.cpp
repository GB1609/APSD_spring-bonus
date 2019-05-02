#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include "Excercise1.h"
#include "Excercise2.h"

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

void printMatrixEx5(int **matrix, int dim) {
	printf("DIMENSION OF MATRIX ====> %d\n", dim);
	for (int a = 0; a < dim; a++) {
		for (int b = 0; b < dim; b++) {
			printf("%d \t", matrix[a][b]);
		}
		printf("\n");
	}
}

int checkMatrixEx5(int **matrix1, int **matrix2, int dim, int numThreads) {
	int finalSum = 0, a, b;
	for (a = 0; a < dim; a++)
		for (b = 0; b < dim; b++)
			if (matrix1[a][b] != matrix2[a][b])
				finalSum++;
	return finalSum;
}

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

void generate_matrix_exercise_5(int ** matrix, int ** s_matrix, int ** p_matrix,
        int dim) {
	int cont = 0;
	for (int a = 0; a < dim; a++) {
		for (int b = 0; b < dim; b++) {
			matrix[a][b] = rand() % 2;
			s_matrix[a][b] = matrix[a][b];
			p_matrix[a][b] = matrix[a][b];
		}
		cont++;
		for (int s = 0; s < cont; s++)
			printf("*");
		printf("\n");
	}
}

int module(int x, int dim) {
	return (x % dim + dim) % dim;
}

int count_live_cell_ex5(int ** matrix, int row, int col, int dim,
        int toroidale) {
	int cont = 0;

	int up = row - 1, down = row + 1, left = col - 1, right = col + 1,
	        up_toroidale = module((row - 1), dim), left_toroidale = module(
	                (col - 1), dim), down_toroidale = module((row + 1), dim),
	        right_toroidale = module((col + 1), dim);
	//UP
	if (up > 0 && matrix[up][col])
		cont++;
	else if (toroidale && matrix[up_toroidale][col])
		cont++;
	//LEFT
	if (left > 0 && matrix[row][left])
		cont++;
	else if (toroidale && matrix[row][left_toroidale])
		cont++;
	//RIGHT
	if (right < dim && matrix[row][right])
		cont++;
	else if (toroidale && matrix[row][right_toroidale])
		cont++;
	//DOWN
	if (down < dim && matrix[down][col])
		cont++;
	else if (toroidale && matrix[down_toroidale][col])
		cont++;

	//DOWN_LEFT
	if (down < dim && left > 0 && matrix[down][left])
		cont++;
	else if (toroidale && matrix[down_toroidale][left_toroidale])
		cont++;
	//up_right
	if (up > 0 && right < dim && matrix[up][right])
		cont++;
	else if (toroidale && matrix[up_toroidale][right_toroidale])
		cont++;
	//UP_LEFT
	if (up > 0 && left > 0 && matrix[up][left])
		cont++;
	else if (toroidale && matrix[up_toroidale][left_toroidale])
		cont++;
	//up_right
	if (up > 0 && right < dim && matrix[up][right])
		cont++;
	else if (toroidale && matrix[up_toroidale][right_toroidale])
		cont++;
	return cont;
}

int** generateTemp(int dim) {
	int** temp = (int**) malloc(dim * sizeof(int*));
	for (int t = 0; t < dim; t++)
		temp[t] = (int*) malloc(dim * sizeof(int));
	return temp;
}

void update_temp_matrix_Ex5(int dim, int** matrix, int **temp) {
	for (int c1 = 0; c1 < dim; c1++)
		for (int c2 = 0; c2 < dim; c2++)
			temp[c1][c2] = matrix[c1][c2];
}

void update_matrix_ex5(int x, int j, int live, int** matrix, int** temp) {
	if (!temp[x][j] && live == 3)
		matrix[x][j] = 1;
	else if (temp[x][j] && (live < 2 || live > 3))
		matrix[x][j] = 0;
}

int serial_ex5(int numIteration, int ** matrix, int dim) {
	double begin_time = omp_get_wtime();
	int** temp = generateTemp(dim);
	for (int a = 0; a < numIteration; a++) {
		update_temp_matrix_Ex5(dim, matrix, temp);
		for (int x = 0; x < dim; x++)
			for (int j = 0; j < dim; j++) {
				int live = count_live_cell_ex5(temp, x, j, dim, 1); //ultimo valore indica se toroidale o no
//				printf("sono x=%d y =%d e intorno ho %d vive \n", x, j, live);
				update_matrix_ex5(x, j, live, matrix, temp);
			}
//		printf("*****TURN %d*****", a);
//		printMatrixEx5(matrix, dim);
	}
	free(temp);
	double end_time = omp_get_wtime() - begin_time;
	printf("*SERIAL GOL EXECUTION TERMINATED IN: %.8g\n", end_time);
	return end_time;
}
double parallel_ex5(int numIteration, int ** matrix, int dim, int nThreads) {
	omp_set_num_threads(nThreads);
	double begin_time = omp_get_wtime();
	int x, j;
	int** temp = generateTemp(dim);
	for (int a = 0; a < numIteration; a++) {
		update_temp_matrix_Ex5(dim, matrix, temp);
#pragma omp parallel for private (x,j) schedule(dynamic,3)
		for (x = 0; x < dim; x++) {
			for (j = 0; j < dim; j++) {
				int live = count_live_cell_ex5(temp, x, j, dim, 1); //ultimo valore indica se toroidale o no
//					printf("sono x=%d y =%d e intorno ho %d vive \n", x, j,
//							live);
				update_matrix_ex5(x, j, live, matrix, temp);
			}
		}

//		printf("*****TURN %d*****", a);
//		printMatrixEx5(matrix, dim);
	}
	free(temp);
	double end_time = omp_get_wtime() - begin_time;
	printf("*PARALLEL EXECUTION TERMINATED IN: %.8g\n", end_time);
	return end_time;
}

void excercise5(int numIteration, int numThreads) {
	int dim = 20;
//	numIteration = 1;
	printf("*********************\nBEGIN EXCERCISE 5\n");
	printf("*BEGIN GENERATION MATRIX*\n");
	int **matrix = (int **) malloc(dim * sizeof(int*));
	int **serial_matrix = (int **) malloc(dim * sizeof(int*));
	int **parallel_matrix = (int **) malloc(dim * sizeof(int*));
	for (int i = 0; i < dim; i++) {
		matrix[i] = (int *) malloc(dim * sizeof(int));
		serial_matrix[i] = (int *) malloc(dim * sizeof(int));
		parallel_matrix[i] = (int *) malloc(dim * sizeof(int));
	}
	generate_matrix_exercise_5(matrix, serial_matrix, parallel_matrix, dim);
	printf("*END GENERATION MATRIX*\nORIGINAL MATRIX\n");
//	printMatrixEx5(matrix, dim);
	serial_ex5(numIteration, serial_matrix, dim);
	parallel_ex5(numIteration, parallel_matrix, dim, 2);
	printf("DiFF============>%d\n",
	        checkMatrixEx5(serial_matrix, parallel_matrix, dim, numThreads));
	printf("END EXCERCISE 5\n*********************\n");
	free(matrix);
	free(serial_matrix);
	free(parallel_matrix);
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
	num_thread = max_num_threads-3;
	dim = 100000000;
	toLaunch = 2;
	printf("THREADS USED: %d",num_thread);
	switch (toLaunch) {
		case 0:
			printf("SUM OF VECTOR\n");
			sumOfVector();
			break;
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
			excercise5(dim, num_thread);
		}
			break;
		default:
			printf("NUMBER NOT KNOW");
	}
	printf("IRON MAN MUORE\n");
	return 0;
}
