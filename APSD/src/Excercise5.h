/*
 * Excercise5.h
 *
 *  Created on: 2 mag 2019
 *      Author: gb1609
 */

#ifndef EXCERCISE5_H_
#define EXCERCISE5_H_

#include <omp.h>

class Excercise5
{
  private:
  int dim,num_steps,num_threads;
  int **serial_matrix,**parallel_matrix;
  FileWriter fw;

  public:
  Excercise5(int d,int ns,int nt,FileWriter fw): dim(d),num_steps(ns),num_threads(nt),fw(fw)
  {
    printf("*BEGIN GENERATION MATRIX*\n");
    serial_matrix = (int **) malloc(dim * sizeof(int*));
    parallel_matrix = (int **) malloc(dim * sizeof(int*));
    for (int i = 0; i < dim; i++)
    {
      serial_matrix[i] = (int *) malloc(dim * sizeof(int));
      parallel_matrix[i] = (int *) malloc(dim * sizeof(int));
    }
    generate_matrix();
    printf("*END GENERATION MATRIX*\nORIGINAL MATRIX\n");
  }

  void execute()
  {
    double serial=serial_execution();
    double parallel=parallel_execution();
    if(checkMatrix()>0)
    printf("*ERROR FINAL VALUES*\n");
    else
    printf("*CORRECT VALUES\t");
    printf("SPEED_UP:%.4g*\n",serial/parallel);
    free(serial_matrix);free(parallel_matrix);
  }
  double serial_execution()
  {
    double begin_time = omp_get_wtime();
    int** temp = generate_temp_matrix();
    for (int a = 0; a < num_steps; a++)
    {
      update_temp_matrix(serial_matrix, temp);
      for (int x = 0; x < dim; x++)
      {
	for (int j = 0; j < dim; j++)
	{
	  int live = count_live_cell(temp, x, j, dim, 1); //ultimo valore indica se toroidale o no
	  //				printf("sono x=%d y =%d e intorno ho %d vive \n", x, j, live);
	  update_values_matrix(x, j, live, serial_matrix, temp);
	}}
      //		printf("*****TURN %d*****", a);
      //		printMatrixEx5(matrix, dim);
    }
    free(temp);
    double end_time = omp_get_wtime() - begin_time;
    printf("*SERIAL GOL EXECUTION TERMINATED IN: %.8g\n", end_time);
    return end_time;
  }

  double parallel_execution()
  {
    omp_set_num_threads(num_threads);
    double begin_time = omp_get_wtime();
    int x, j;
    int** temp = generate_temp_matrix();
    for (int a = 0; a < num_steps; a++)
    {
      update_temp_matrix(parallel_matrix, temp);
#pragma omp parallel for private (x,j) schedule(dynamic,3)
      for (x = 0; x < dim; x++)
      {
	for (j = 0; j < dim; j++)
	{
	  int live = count_live_cell(temp, x, j, dim, 1); //ultimo valore indica se toroidale o no
	  //					printf("sono x=%d y =%d e intorno ho %d vive \n", x, j,
	  //							live);
	  update_values_matrix(x, j, live, parallel_matrix, temp);
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

  int module(int x)
  {
    return (x % dim + dim) % dim;
  }

  void generate_matrix()
  {
    int cont = 0;
    for (int a = 0; a < dim; a++)
    {
      for (int b = 0; b < dim; b++)
      {
	serial_matrix[a][b] = rand() % 2;
	parallel_matrix[a][b] = serial_matrix[a][b];
      }
      cont++;
      for (int s = 0; s < cont; s++)
      printf("*");
      printf("\n");
    }
  }

  void printMatrix(int** matrix)
  {
    printf("DIMENSION OF MATRIX ====> %d\n", dim);
    for (int a = 0; a < dim; a++)
    {
      for (int b = 0; b < dim; b++)
      {
	printf("%d \t", matrix[a][b]);
      }
      printf("\n");
    }
  }

  int checkMatrix()
  {
    int finalSum = 0, a, b;
    for (a = 0; a < dim; a++)
    for (b = 0; b < dim; b++)
    if (serial_matrix[a][b] != parallel_matrix[a][b])
    finalSum++;
    return finalSum;
  }

  int** generate_temp_matrix()
  {
    int **temp = (int**) malloc(dim * sizeof(int*));
    for (int t = 0; t < dim; t++)
    temp[t] = (int*) malloc(dim * sizeof(int));
    return temp;
  }

  void update_temp_matrix(int** matrix, int **temp)
  {
    for (int c1 = 0; c1 < dim; c1++)
    for (int c2 = 0; c2 < dim; c2++)
    temp[c1][c2] = matrix[c1][c2];
  }

  void update_values_matrix(int x, int j, int live, int** matrix, int** temp)
  {
    if (!temp[x][j] && live == 3)
    matrix[x][j] = 1;
    else if (temp[x][j] && (live < 2 || live > 3))
    matrix[x][j] = 0;
  }

  int count_live_cell(int ** matrix, int row, int col, int dim,
      int toroidale)
  {
    int cont = 0;

    int up = row - 1, down = row + 1, left = col - 1, right = col + 1,
    up_toroidale = module(row - 1), left_toroidale = module(
	col - 1), down_toroidale = module(row + 1),
    right_toroidale = module(col + 1);
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

};

#endif /* EXCERCISE5_H_ */
