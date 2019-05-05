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
  dim = 400000000;
  toLaunch = 2;
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
      printf ("***Exercise 3*** \n");
      Excercise3 (dim, num_thread, fw).execute ();
      break;
    case 4:
    {
      printf ("***Exercise 4*** \n");
      Excercise4 ex4 (dim, num_thread, fw);
      ex4.execute ();
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
