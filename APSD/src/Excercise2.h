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
  int *vectorA, *vectorB, *vectorC,*vectorCParallel;
  public:
  Excercise2(int d,int nt,FileWriter fw): dim(d),num_threads(nt), fw(fw)
  {
    vectorA = (int*) malloc(dim * sizeof(int));
    vectorB = (int*) malloc(dim * sizeof(int));
    vectorC = (int*) malloc(dim * sizeof(int));
    vectorCParallel = (int*) malloc(dim * sizeof(int));
    for (int a = 0; a < dim; a++)
    {
      vectorA[a] = rand() % 10;
      vectorB[a] = rand() % 10;
    }

  }

  double serial_execution()
  {
    double begin=omp_get_wtime();
    int a;
    for (a = 0; a < dim; a++)
    {
      vectorC[a] = vectorA[a] + vectorB[a];
    }
    int serial_sum = 0;
    for (a = 0; a < dim; a++)
    {
      serial_sum += vectorC[a] * vectorB[a];
    }
    double end=omp_get_wtime()-begin;
    printf("*SERIAL EXECUTION TIME:%.8g*\tSCALAR PRODUCT: %d*\n",end,serial_sum);
    return end;
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
      vectorCParallel[a] = vectorA[a] + vectorB[a];
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
    printf("*PARALLEL EXECUTION TIME:%.8g\tSCALAR PRODUCT: %d*\n",end,parallel_sum);
    return end;
  }

  void execute()
  {
    double serial=serial_execution();
    double parallel=parallel_execution();
    if(checkVectors()>0)
    printf("ERROR IN C");
    else
    printf("SPEED_UP:%.4g\n",serial/parallel);
    free(vectorA);
    free(vectorB);
    free(vectorC);
    free(vectorCParallel);
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
