/*
 * Excercise3.h
 *
 *  Created on: 2 mag 2019
 *      Author: gb1609
 */

#ifndef EXCERCISE3_H_
#define EXCERCISE3_H_

class Excercise3
{
  private:
  int dim;
  int num_threads;
  double serial;
  double serial_montecarlo;
  FileWriter fw;
  public:
  Excercise3(int d,FileWriter fw):dim(d), fw(fw)
  {
    num_threads=0;
    serial=0.0;
    serial_montecarlo=0.0;
  };

  void write()
  { fw.write();};

  void serial_execute()
  {
    double x, pi, sum = 0.0;
    double step = 1.0 / (double) dim;
    double begin_time = omp_get_wtime();
    for (int i = 0; i < dim; i++)
    {
      x = (i + 0.5) * step;
      sum += 4.0 / (1.0 + x * x);
    }
    pi = step * sum;
    double end_time = omp_get_wtime() - begin_time;
    serial= end_time;
  }
  double critical_execute()
  {
    double pi = 0.0, begin_time, end_time;
    double step = 1.0 / (double) dim;
    omp_set_num_threads(num_threads);
    begin_time = omp_get_wtime();
#pragma omp parallel if(dim > 10000)
    {
      double x, sum;
#pragma omp for
      for (int i = 0; i < dim; i++)
      {
	x = (i + 0.5) * step;
	sum += 4.0 / (1.0 + x * x);
      }

#pragma omp critical
      pi += step * sum;
    }
    end_time = omp_get_wtime() - begin_time;
    return end_time;
  }
  double padding_execute()
  {
    double pi = 0.0, begin_time, end_time, x;
    double pi_padding[num_threads][8];
    double step = 1.0 / (double) dim;
    omp_set_num_threads(num_threads);
    begin_time = omp_get_wtime();
#pragma omp parallel if(dim > 10000)
    {
      int id = omp_get_thread_num();
#pragma omp for
      for (int i = 0; i < dim; i++)
      {
	x = (i + 0.5) * step;
	pi_padding[id][0] += 4.0 / (1.0 + x * x);
      }
    }
    for (int i = 0; i < num_threads; i++)
    pi += pi_padding[i][0] * step;
    end_time = omp_get_wtime() - begin_time;
    return end_time;
  }
  double parallel_reduction()
  {
    int i;
    double x, pi, sum = 0.0;
    double step = 1.0 / (double) dim;
    omp_set_num_threads(num_threads);
    double begin_time = omp_get_wtime();
#pragma omp parallel for private(i) reduction (+:sum) if(dim > 10000)
    for (i = 0; i < dim; i++)
    {
      x = (i + 0.5) * step;
      sum += 4.0 / (1.0 + x * x);
    }
    pi = step * sum;
    double end_time = omp_get_wtime() - begin_time;
    return end_time;
  }
  void serial_monte_carlo_execute()
  {
    int i, nCirc = 0;
    double r = 1.0, pi, x, y;
    unsigned int seed = 0;
    double begin_time = omp_get_wtime();
    for (i = 0; i < dim; i++)
    {
      x = (double) rand_r(&seed) / RAND_MAX;
      y = (double) rand_r(&seed) / RAND_MAX;
      if (x*x + y*y <= r*r)
      nCirc+=1;
    }
    pi = 4.0 * ((double) nCirc / (double) dim);
    double end_time = omp_get_wtime() - begin_time;
    serial_montecarlo= end_time;
  }
  double parallel_reduction_monte_carlo()
  {
    int i, nCirc = 0;
    double r = 1.0, pi, x, y;
    unsigned int seed = 0;
    double begin_time = omp_get_wtime();
    omp_set_num_threads(num_threads);
#pragma omp parallel private(x, y)
    {
#pragma omp for reduction (+:nCirc)
      for (i = 0; i < dim; i++)
      {
	x = (double) rand_r(&seed) / RAND_MAX;
	y = (double) rand_r(&seed) / RAND_MAX;
	if (x*x + y*y <= r*r)
	nCirc+=1;
      }
    }
    pi = 4.0 * ((double) nCirc / (double) dim);
    return omp_get_wtime() - begin_time;
  }
  double parallel_critical_monte_carlo()
  {
    int i, nCirc = 0;
    double r = 1.0, pi, x, y;
    unsigned int seed = 0;
    double begin_time = omp_get_wtime();
    omp_set_num_threads(num_threads);
#pragma omp parallel for private(x,y)
    for (i = 0; i < dim; i++)
    {
      x = (double) rand_r(&seed) / RAND_MAX;
      y = (double) rand_r(&seed) / RAND_MAX;
      if (pow(x, 2) + pow(y, 2) <= pow(r, 2))
#pragma omp critical
      nCirc+=1;
    }
    pi = 4.0 * ((double) nCirc / (double) dim);
    return omp_get_wtime() - begin_time;
  }
  void execute(int nt)
  {
    num_threads=nt;
    double critical=critical_execute();
    double padding=padding_execute();
    double reduction=parallel_reduction();
    double critila_sp_up = serial / critical;
    double padding_speed_up = serial / padding;
    double reduction_speed_up = serial / reduction;
    double par_red_mc=parallel_reduction_monte_carlo();
    double reduction_mc_speed_up = serial_montecarlo / par_red_mc;
    printf("CRITICAL NT: %d,SERIAL: %.8g, PARALLEL: %.8g,SPEED UP MC: %.8g\n",num_threads,serial,critical,critila_sp_up);
    printf("PADDING NT: %d,SERIAL: %.8g, PARALLEL: %.8g,SPEED UP MC: %.8g\n",num_threads,serial,padding,padding_speed_up);
    printf("REDUCTION NT: %d,SERIAL: %.8g, PARALLEL: %.8g,SPEED UP MC: %.8g\n",num_threads,serial,reduction,reduction_speed_up);
    printf("REDUCTION NT: %d,SERIAL: %.8g, PARALLEL: %.8g,SPEED UP MC: %.8g\n",num_threads,serial_montecarlo,par_red_mc,reduction_mc_speed_up);
    fw.update_string(num_threads,serial_montecarlo,par_red_mc,reduction_mc_speed_up,"ex3 reduction MONTECARLO");
    fw.update_string(num_threads,serial,critical,critila_sp_up,"ex3 CIRTICAL");
    fw.update_string(num_threads,serial,padding,padding_speed_up,"ex3 PADDING");
    fw.update_string(num_threads,serial,reduction,reduction_speed_up,"ex3 REDUCTION");

  }
};

#endif /* EXCERCISE3_H_ */
