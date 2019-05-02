/*
 * Excercise3.h
 *
 *  Created on: 2 mag 2019
 *      Author: gb1609
 */

#ifndef EXCERCISE3_H_
#define EXCERCISE3_H_

class Excercise3 {
	private:
	int dim;
	int num_threads;
	public:
	Excercise3(int d, int nt):dim(d),num_threads(nt) {}

	double excercise3_serial() {
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

	double excercise3_critical() {
		double pi = 0.0, begin_time, end_time;
		double step = 1.0 / (double) dim;
		omp_set_num_threads(num_threads);
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
	double excercise3_padding() {
		double pi = 0.0, begin_time, end_time, x;
		double pi_padding[num_threads][8];
		double step = 1.0 / (double) dim;
		omp_set_num_threads(num_threads);
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
		for (int i = 0; i < num_threads; i++)
		pi += pi_padding[i][0] * step;

		end_time = omp_get_wtime() - begin_time;
		printf("*PADDING EXECUTION TERMINATED IN: %.8g\tVALUE OF PI= %.16g*\n",
				end_time, pi);
		return end_time;
	}
	double excercise3_reduction() {
		int i;
		double x, pi, sum = 0.0;
		double step = 1.0 / (double) dim;
		omp_set_num_threads(num_threads);
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
	double excercise3_serial_monte_carlo() {
		int i, nCirc = 0;
		double r = 1.0, pi, x, y;
		omp_set_num_threads(num_threads);
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
	double excercise3_parallel_monte_carlo() {
		int i, nCirc = 0;
		double r = 1.0, pi, x, y;
		unsigned int seed = 0;
		double begin_time = omp_get_wtime();
		omp_set_num_threads(num_threads);
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

	void execute()
	{
		double serial=excercise3_serial();
		double critical=excercise3_critical();
		double padding=excercise3_padding();
		double reduction=excercise3_reduction();
		printf("SPEED UP CRITICAL: %.8g\nSPEED UP PADDING: %.8g\nSPEED UP REDUCTION: %.8g\n",critical,padding,reduction);

		double serial_mc=excercise3_serial_monte_carlo();
		double parallel_mc=excercise3_parallel_monte_carlo();
		printf("SPEED UP MC: %.8g\n",parallel_mc);
	}
};

#endif /* EXCERCISE3_H_ */
