/*
 * Excercise1.h
 *
 *  Created on: 2 mag 2019
 *      Author: gb1609
 */

#ifndef EXCERCISE1_H_
#define EXCERCISE1_H_

#include <omp.h>

class Excercise1 {
private:
	int dim;
	int numThread;
	double **matrixA;
	double **matrixB;
	double **matrixAParallel = (double **) malloc(dim * sizeof(double));
	double **matrixBParallel = (double **) malloc(dim * sizeof(double));
public:
	Excercise1(int d, int nt): dim(d),numThread(nt) {
		matrixA = (double **) malloc(dim * sizeof(double*));
		matrixB = (double **) malloc(dim * sizeof(double*));
		matrixAParallel = (double **) malloc(dim * sizeof(double));
		matrixBParallel = (double **) malloc(dim * sizeof(double));
		for (int i = 0; i < dim; i++) {
			matrixA[i] = (double *) malloc(dim * sizeof(double));
			matrixB[i] = (double *) malloc(dim * sizeof(double));
			matrixAParallel[i] = (double *) malloc(dim * sizeof(double));
			matrixBParallel[i] = (double *) malloc(dim * sizeof(double));
		}
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

	void execute()
	{
		double serial=serial_execute();
		double parallel=parallel_execute();
		check_equality();
		printf("SPEED_UP:%.4g\n",serial/parallel);
		free(matrixA);
		free(matrixB);
		free(matrixAParallel);
		free(matrixBParallel);
	}

	double serial_execute() {
		int a, b;
		double time_serial, time_begin, time_serial_check;
		time_begin = omp_get_wtime();
		for (a = 0; a < dim; a++)
		for (b = 0; b < dim; b++) {
			matrixA[a][b] = 5 * pow(a, 3) + 5 * M_PI * pow(b, 6);
			matrixB[a][b] = (10 / 3) * (matrixA[a][b]);
		}
		//	printMatrix(matrixA, dim);
		time_serial = omp_get_wtime() - time_begin;
		printf(
				"*TIME SERIAL EXECUTION: %8g\n*",
				time_serial);
		return time_serial;
	}

	double parallel_execute() {
		int a,b;
		double end_time,begin_time;
		omp_set_num_threads(numThread);
		begin_time=omp_get_wtime();
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
		end_time = omp_get_wtime() - begin_time;
		printf(
				"*TIME PARALLEL EXECUTION: %8g\n*",
				end_time);
		return end_time;
	}

	void check_equality()
	{
		double end_time,begin_time;
		begin_time=omp_get_wtime();
		int finalSum = 0, a, b;
#pragma omp parallel for shared(matrixA,matrixAParallel) private(a,b) reduction(+:finalSum)
			for (a = 0; a < dim; a++)
				for (b = 0; b < dim; b++)
					if (matrixA[a][b] != matrixAParallel[a][b])
						finalSum++;
		int sum=finalSum;
		finalSum=0;
#pragma omp parallel for shared(matrixB,matrixBParallel) private(a,b) reduction(+:finalSum)
			for (a = 0; a < dim; a++)
			for (b = 0; b < dim; b++)
			if (matrixB[a][b] != matrixBParallel[a][b])
			finalSum++;

		finalSum+=sum;
		if(finalSum>0)
		printf("*ERROR");
		else
		printf("*CORRECT, SAME VALUE");
		end_time=omp_get_wtime()-begin_time;
		printf("TIME CHECKING:%.8g\n",end_time);
	}

	void printMatrixEx1(double **matrix) {
		printf("%d", dim);
		for (int a = 0; a < dim; a++) {
			for (int b = 0; b < dim; b++) {
				printf("%.8g \t", matrix[a][b]);
			}
			printf("\n");
		}
	}

};

#endif /* EXCERCISE1_H_ */


