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
  double serial_count;
  double serial_only;
  FileWriter fw;
  public:
  Excercise4(int d, FileWriter fw):dim(d), fw(fw)
  {
    serial_count=0.0;
    serial_only=0.0;
    num_threads=0;
    vector = (int*) malloc(dim * sizeof(int));
    for (int a = 0; a < dim; a++)
    vector[a] = rand() % 10000;
    toSearch = rand() % 10000;
  }

  void clean()
  {
    free(vector);
  }

  void write()
  {
    fw.write();
  }

  void serial_find_count_execute()
  {
    double time_begin = omp_get_wtime();
    int find = 0;
    for (int a = 0; a < dim; a++)
    if (vector[a] == toSearch)
    find++;
    double time_end = omp_get_wtime() - time_begin;
    serial_count=time_end;
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
    return time_end;
  }

  void serial_find_only()
  {
    double begin_time=omp_get_wtime();
    for (int a=0;a<dim;a++)
    if(vector[a]==toSearch)
    {
      serial_only=omp_get_wtime()-begin_time;
      return;
    }
    serial_only=omp_get_wtime()-begin_time;
    return;
  }

  double parallel_find_only_execute()
  {
    double begin=omp_get_wtime(),end;
    int a;
    bool find=false;
    double parallel_find =0.0;
    int slice=dim/num_threads;
#pragma omp parallel shared(find,begin) private(a)
    {
#pragma omp for
      for (a = 0; a < dim; a++)
      {
	if (!find && vector[a] == toSearch)
	{
#pragma omp critical
	  { find=true;
	    end=omp_get_wtime()-begin;
	  }
	  if(find)
	  a=dim;
	}
      }
    }
    if(!find)
    end=omp_get_wtime()-begin;
    return end;
  }
  void execute(int nt)
  {
    num_threads=nt;
    double parallel_fc_reduction=parallel_reduction_find_count_execute();
    double red_speed_up=serial_count/parallel_fc_reduction;
    double parallel_fc_atomic=parallel_atomic_find_count_execute();
    double atomic_speed_up=serial_count/parallel_fc_atomic;
    printf("REDUCTION - NT: %d,SERIAL: %.8g, PARALLEL: %.8g,SPEED UP: %.8g\n",num_threads,serial_count,parallel_fc_reduction,red_speed_up);
    printf("ATOMIC - NT: %d,SERIAL: %.8g, PARALLEL: %.8g,SPEED UP: %.8g\n",num_threads,serial_count,parallel_fc_atomic,atomic_speed_up);
    double pfo= parallel_find_only_execute();
    double pfo_speed_up=serial_only/pfo;
    printf("ONLY FIND - NT: %d,SERIAL: %.8g, PARALLEL: %.8g,SPEED UP: %.8g\n",num_threads,serial_only,pfo,pfo_speed_up);
    fw.update_string(num_threads,serial_count,parallel_fc_atomic,atomic_speed_up,"ex4 COUNT ATOMIC");
    fw.update_string(num_threads,serial_count,parallel_fc_reduction,red_speed_up,"ex4 COUNT REDUCTION");
    fw.update_string(num_threads,serial_only,pfo,pfo_speed_up,"ex4 FIND ONLY");
  }
};

#endif /* EXCERCISE4_H_ */
