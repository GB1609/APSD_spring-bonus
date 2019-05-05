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
  int toLaunch, num_thread, dim;
  const int max_num_threads = omp_get_max_threads ();
  dim = 10000000;
  toLaunch = 3;
  FileWriter fw (toLaunch);
  switch (toLaunch)
  {
    case 1:
    {
      printf ("***Exercise 1*** \n");
      Excercise1 ex1 = Excercise1 (dim, fw);
      ex1.serial_execute ();
      for (num_thread = 2; num_thread < max_num_threads + 1; num_thread++)
	ex1.execute (num_thread);
      ex1.write ();
      ex1.clean ();
    }
      break;
    case 2:
    {
      printf ("***Exercise 2***\n");
      Excercise2 ex2 = Excercise2 (dim, fw);
      ex2.serial_execution ();
      for (num_thread = 2; num_thread < max_num_threads + 1; num_thread++)
	ex2.execute (num_thread);
      ex2.write ();
      ex2.clean ();
    }
      break;
    case 3:
    {
      printf ("***Exercise 3*** \n");
      Excercise3 excercise3 = Excercise3 (dim, fw);
      excercise3.serial_execute ();
      excercise3.serial_monte_carlo_execute ();
      for (num_thread = 2; num_thread < max_num_threads + 1; num_thread++)
	excercise3.execute (num_thread);
      excercise3.write ();
    }
      break;
    case 4:
    {
      printf ("***Exercise 4*** \n");
      Excercise4 ex4 (dim, fw);
      ex4.serial_find_count_execute ();
      ex4.serial_find_only ();
      for (num_thread = 2; num_thread < max_num_threads + 1; num_thread++)
	ex4.execute (num_thread);
      ex4.write ();
      ex4.clean ();
    }
      break;
    case 5:
    {
      printf ("***Exercise 5*** \n");
      Excercise5 ex5 (dim, 100, fw);
      ex5.serial_execute ();
      for (num_thread = 2; num_thread < max_num_threads + 1; num_thread++)
      {
	ex5.execute (num_thread);
	ex5.reset_parallel_matrix ();
      }
      ex5.write ();
      ex5.clean ();
    }
      break;
    default:
      printf ("NUMBER NOT KNOW");
  }
  return 0;
}
