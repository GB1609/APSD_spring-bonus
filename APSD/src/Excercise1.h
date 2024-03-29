/*
 * Excercise1.h
 *
 *  Created on: 2 mag 2019
 *      Author: gb1609
 */

#ifndef EXCERCISE1_H_
#define EXCERCISE1_H_

class Excercise1
{
  private:
  int dim;
  int num_threads;
  double serial;
  double **matrixA;
  double **matrixB;
  double **parallel_matrixA = (double **) malloc(dim * sizeof(double));
  double **parallel_matrixB = (double **) malloc(dim * sizeof(double));
  FileWriter fw;
  public:
  Excercise1(int d, FileWriter fw): dim(d), fw(fw)
  {
    serial=0.0;
    num_threads=0;
    matrixA = (double **) malloc(dim * sizeof(double*));
    matrixB = (double **) malloc(dim * sizeof(double*));
    parallel_matrixA = (double **) malloc(dim * sizeof(double));
    parallel_matrixB = (double **) malloc(dim * sizeof(double));
    for (int i = 0; i < dim; i++)
    {
      matrixA[i] = (double *) malloc(dim * sizeof(double));
      matrixB[i] = (double *) malloc(dim * sizeof(double));
      parallel_matrixA[i] = (double *) malloc(dim * sizeof(double));
      parallel_matrixB[i] = (double *) malloc(dim * sizeof(double));
    }
  }

  void printMatrix(double **matrix, int dim)
  {
    printf("%d", dim);
    for (int a = 0; a < dim; a++)
    {
      for (int b = 0; b < dim; b++)
      {
	printf("%.8g \t", matrix[a][b]);
      }
      printf("\n");
    }
  }

  void execute(int num_thread)
  {
    num_threads=num_thread;
    double parallel=parallel_execute();
    check_equality();
    double speed_up = serial / parallel;
    printf("NT: %d,SERIAL: %.8g,PARALLEL: %.8g,SPEED_UP:%.4g\n",num_thread,serial,parallel,speed_up);
    fw.update_string(num_threads,serial,parallel,speed_up,"ex1 static");
  }

  void write()
  {
    fw.write();
  }

  void clean()
  {
    free(matrixA);
    free(matrixB);
    free(parallel_matrixA);
    free(parallel_matrixB);
  }

  void serial_execute()
  {
    int a, b;
    double time_serial, time_begin, time_serial_check;
    time_begin = omp_get_wtime();
    for (a = 0; a < dim; a++)
    for (b = 0; b < dim; b++)
    {
      matrixA[a][b] = 5 * pow(a, 3) + 5 * M_PI * pow(b, 6);
      matrixB[a][b] = (10 / 3) * (matrixA[a][b]);
    }
    //	printMatrix(matrixA, dim);
    time_serial = omp_get_wtime() - time_begin;
    printf(
	"*TIME SERIAL EXECUTION: %8g\n*",
	time_serial);
    serial=time_serial;
  }

  double parallel_execute()
  {
    int a,b;
    double end_time,begin_time;
    omp_set_num_threads(num_threads);
    begin_time=omp_get_wtime();
#pragma omp parallel private(a, b)
    {
#pragma omp for
      for (a = 0; a < dim; a++)
      {
	for (b = 0; b < dim; b++)
	{
	  parallel_matrixA[a][b] = 5 * pow(a, 3) + 5 * M_PI * pow(b, 6);
	  parallel_matrixB[a][b] = (10 / 3) * (matrixA[a][b]);
	}
      }
    }
    end_time = omp_get_wtime() - begin_time;
    return end_time;
  }

  void check_equality()
  {
    double end_time,begin_time;
    begin_time=omp_get_wtime();
    int finalSum = 0, a, b;
#pragma omp parallel for shared(matrixA,parallel_matrixA) private(a,b) reduction(+:finalSum)
    for (a = 0; a < dim; a++)
    for (b = 0; b < dim; b++)
    if (matrixA[a][b] != parallel_matrixA[a][b])
    finalSum++;
    int sum=finalSum;
    finalSum=0;
#pragma omp parallel for shared(matrixB,parallel_matrixB) private(a,b) reduction(+:finalSum)
    for (a = 0; a < dim; a++)
    for (b = 0; b < dim; b++)
    if (matrixB[a][b] != parallel_matrixB[a][b])
    finalSum++;

    finalSum+=sum;
    if(finalSum>0)
    printf("*ERROR\n");
    else
    printf("*CORRECT, SAME VALUES\t");
    end_time=omp_get_wtime()-begin_time;
    printf("TIME CHECKING: %.8g*\n",end_time);
  }

};

#endif /* EXCERCISE1_H_ */

