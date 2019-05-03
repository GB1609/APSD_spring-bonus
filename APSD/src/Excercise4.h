/*
 * Excercise4.h
 *
 *  Created on: 2 mag 2019
 *      Author: gb1609
 */

#ifndef EXCERCISE4_H_
#define EXCERCISE4_H_

class Excercise4
{
  private:
  int dim;
  int num_threads;
  int *vector;
  int toSearch;
  FileWriter fw;
  public:
  Excercise4(int d, int nt,FileWriter fw):dim(d),num_threads(nt), fw(fw)
  {
    vector = (int*) malloc(dim * sizeof(int));
    for (int a = 0; a < dim; a++)
    vector[a] = rand() % 1000;
    toSearch = rand() % 1000;
  }

  double serial_find_count_execute()
  {
    double time_begin = omp_get_wtime();
    int find = 0;
    for (int a = 0; a < dim; a++)
    if (vector[a] == toSearch)
    find++;
    double time_end = omp_get_wtime() - time_begin;
    printf(
	"*FIND AND COUNT SERIAL EXECUTION TIME: %.8g\tFOUNDED: %d times (SIZE OF VECTOR:%d)*\n",
	time_end, find, dim);
    return time_end;
  }

  double parallel_reduction_find_count_execute()
  {
    int parallel_find = 0, a;
    double time_begin = omp_get_wtime();
#pragma omp parallel for shared(vector) private(a) reduction (+:parallel_find)
    for (a = 0; a < dim; a++)
    if (vector[a] == toSearch)
    parallel_find++;
    double time_end = omp_get_wtime() - time_begin;
    printf(
	"*FIND AND COUNT PARALLEL REDUCITON EXECUTION TIME: %.8g\tFOUNDED: %d times (SIZE OF VECTOR:%d)*\n",
	time_end, parallel_find, dim);
    return time_end;
  }

  double parallel_atomic_find_count_execute()
  {
    int parallel_find = 0;
    double time_begin = omp_get_wtime();
#pragma omp parallel
    {
#pragma omp  for
      for (int i = 0; i < dim; i++)
      {
	if (vector[i] == toSearch)
#pragma omp atomic
	parallel_find++;
      }

    }
    double time_end = omp_get_wtime() - time_begin;
    printf(
	"*FIND AND COUNT PARALLEL ATOMIC EXECUTION TIME: %.8g\tFOUNDED: %d times (SIZE OF VECTOR:%d)*\n",
	time_end, parallel_find, dim);
    return time_end;
  }

  bool serial_find_only()
  {
    double begin_time=omp_get_wtime();
    for (int a=0;a<dim;a++)
    if(vector[a]==toSearch)
    {
      return true;
    }
    return false;
  }

  int parallel_find_only_execute()
  {
    int a;
    int parallel_find =0;
#pragma omp parallel shared(vector,parallel_find) private(a)
    {
      for (a = 0; a < dim && parallel_find==0; a++)
      if (vector[a] == toSearch)
      {
#pragma omp critical
	parallel_find= 1;
      }

    }
    return parallel_find;
  }
  void execute()
  {
    double serial_fc=serial_find_count_execute();
    double parallel_fc_reduction=parallel_reduction_find_count_execute();
    double parallel_fc_atomic=parallel_atomic_find_count_execute();
    printf("SPEED UP REDUCTION: %.8g\nSPEED UP ATOMIC: %.8g\n",serial_fc/parallel_fc_atomic,serial_fc/parallel_fc_atomic);

    double begin_time,serial_time,parallel_time;
    begin_time=omp_get_wtime();
    if(serial_find_only())
    printf("*SERIAL ONLY FIND ===>>>> FOUND IN: ");
    else
    printf("*SERIAL ONLY FIND ===>>>> NOT FOUND IN: ");
    serial_time=omp_get_wtime()-begin_time;
    printf("%.8g*\n",serial_time);
    begin_time=omp_get_wtime();
    if(parallel_find_only_execute()>0)
    printf("*PARALLEL ONLY FIND ===>>>> FOUND IN: ");
    else
    printf("*PARALLEL ONLY FIND ===>>>> NOT FOUND IN: ");
    parallel_time=omp_get_wtime()-begin_time;
    printf("%.8g*\n",parallel_time);
    printf("*SPEED UP ONLY FIND %.8g*\n",serial_time/parallel_time);

  }
};

#endif /* EXCERCISE4_H_ */
