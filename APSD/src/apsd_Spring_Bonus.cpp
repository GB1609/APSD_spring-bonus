#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
#include "FileWriter.h"
#include "Excercise1.h"
#include "Excercise2.h"
#include "Excercise3.h"
#include "Excercise4.h"
#include "Excercise5.h"

int
main ()
{
  printf (
      "Insert number of bonus:\n 0 (for sum of Vector)\n 1 (for 2DMatrixGeneration) \n 2 (Exercise 2)\n 3 (Calculation of PI)\n");
  printf ("4 (for find an element in a vector)\n 5 (for Game of Life)\n ");
  printf (
      "followed by number of thread and dimension of container in form: excercise-numThreads-dimension \n EXAMPLE: 3-6-10000\n");
  int toLaunch, num_thread = 0, dim;
  const int numberExcercise = 5;
  const int max_num_threads = omp_get_max_threads ();
//	scanf("%d-%d-%d", &toLaunch, &num_thread, &dim);
//	while(toLaunch>numberExcercise || toLaunch<1 || max_num_threads < num_thread || num_thread==0)
//	{
//		printf("ERROR INPUT!!!\nEXAMPLE: 3-6-10000\n");
//		scanf("%d-%d-%d", &toLaunch, &num_thread, &dim);
//	}
  num_thread = max_num_threads - 3;
  dim = 20000;
  toLaunch = 5;
  printf ("THREADS USED: %d\n", num_thread);
  FileWriter fw (toLaunch);
  switch (toLaunch)
  {
    case 1:
      printf ("***Exercise 1*** \n");
      Excercise1 (dim, num_thread, fw).execute ();
      break;
    case 2:
      printf ("***Exercise 2***\n");
      Excercise2 (dim, num_thread, fw).execute ();
      break;
    case 3:
      printf ("***Exercise 3*** \n");
      Excercise3 (dim, num_thread, fw).execute ();
      break;
    case 4:
      printf ("***Exercise 4*** \n");
      Excercise4 (dim, num_thread, fw).execute ();
      break;
    case 5:
      printf ("***Exercise 5*** \n");
      Excercise5 (dim, 1, num_thread, fw).execute ();
      break;
    default:
      printf ("NUMBER NOT KNOW");
  }
  return 0;
}
