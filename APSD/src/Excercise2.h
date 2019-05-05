/*
 * Excercise2.h
 *
 *  Created on: 2 mag 2019
 *      Author: gb1609
 */

#ifndef EXCERCISE2_H_
#define EXCERCISE2_H_

#include <omp.h>

class Excercise2
{
  private:
  FileWriter fw;
  int dim;
  int num_threads;
  double serial;
  double *vectorA, *vectorB, *vectorC,*vectorCParallel;
  public:
  Excercise2(int d,FileWriter fw): dim(d), fw(fw)
  {
    serial=0.0;num_threads=0;
    vectorA = (double*) malloc(dim * sizeof(double));
    vectorB = (double*) malloc(dim * sizeof(double));
    vectorC = (double*) malloc(dim * sizeof(double));
    vectorCParallel = (double*) malloc(dim * sizeof(double));
    for (int a = 0; a < dim; a++)
    {
      vectorA[a] = (double) (pow((rand() % 1000)/INT_MAX,6)+1.0);
      vectorB[a] = (double)(pow((rand() % 1000)/INT_MAX,5)+1.0);
    }

  }

  void write()
  {
    fw.write();
  }

  void clean()
  {
    free(vectorA);
    free(vectorB);
    free(vectorC);
    free(vectorCParallel);
  }

  void serial_execution()
  {
    double begin=omp_get_wtime();
    int a;
    for (a = 0; a < dim; a++)
    {
      vectorC[a] = pow(pow(vectorA[a] + vectorB[a],3)*3.14,3); // CAMBIATO PERCHE ALTRIMENTI CI METTEVA TROPPO POCO
    }
    int serial_sum = 0;
    for (a = 0; a < dim; a++)
    {
      serial_sum += vectorC[a] * vectorB[a];
    }
    double end=omp_get_wtime()-begin;
    serial= end;
  }

  double parallel_execution()
  {
    omp_set_num_threads(num_threads);
    int a;
    double begin=omp_get_wtime();
    int parallel_sum = 0, partial_sum, thread_n;
#pragma omp parallel for shared(vectorA,vectorB,vectorCParallel) private(a)
    for (a = 0; a < dim; a++)
    {
      vectorCParallel[a] = pow(pow(vectorA[a] + vectorB[a],3)*3.14,3);
    }
#pragma omp parallel private(a,thread_n,partial_sum)
    {
      partial_sum = 0;
      thread_n = omp_get_thread_num();

#pragma omp for private(a) reduction(+:parallel_sum)
      for (a = 0; a < dim; a++)
      {
	parallel_sum += (vectorCParallel[a] * vectorB[a]);
	partial_sum = parallel_sum;
      }
    }
    double end=omp_get_wtime()-begin;
    return end;
  }

  void execute(int nt)
  {
    num_threads=nt;
    double parallel=parallel_execution();
    double speed_up=serial/parallel;
    if(checkVectors()>0)
    printf("ERROR IN C\n");
    else
    printf("CORRECT \n");
    printf("NT: %d,SERIAL: %.8g,PARALLEL: %.8g,SPEED_UP:%.4g\n",num_threads,serial,parallel,speed_up);
    fw.update_string(num_threads,serial,parallel,speed_up,"ex2 parallel");
  }

  int checkVectors()
  {
    int finalSum = 0, a;
#pragma omp parallel for shared(vectorC,vectorCParallel) private(a) reduction(+:finalSum)
    for (a = 0; a < dim; a++)
    if (vectorC[a] != vectorCParallel[a])
    finalSum++;
    return finalSum;

  }
};

#endif /* EXCERCISE2_H_ */
