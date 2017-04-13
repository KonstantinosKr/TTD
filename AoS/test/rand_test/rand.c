#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "loba.h"
#include "error.h"
#include "tribal_ispc.h"

/* migrate triangles "in-place" to new ranks */
static void migrate_triangles (int size, int nt, REAL *t[3][3], REAL *v[3], int *rank)
{
  /* TODO */
}

int main (int argc, char **argv)
{
  REAL *t[3][3]; /* triangles */
  REAL *v[3]; /* velocities */
  unsigned int *tid; /* triangle identifiers */
  REAL lo[3] = {0, 0, 0}; /* lower corner */
  REAL hi[3] = {1, 1, 1}; /* upper corner */
  int nt; /* number of triangles */
  int *rank; /* migration ranks */
  int size; /* buffers size */
  int i, j;
  int myrank;

  /* init */

  MPI_Init (&argc, &argv);

  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

  if (myrank == 0)
  {
    /* set nt */

    if (argc > 1) nt = atoi (argv[1]);
    else nt = 10000;

    /* buffers */

    size = 4*nt;

    for (i = 0; i < 3; i ++)
    {
      ERRMEM (t[0][i] = malloc (sizeof(REAL[size])));
      ERRMEM (t[1][i] = malloc (sizeof(REAL[size])));
      ERRMEM (t[2][i] = malloc (sizeof(REAL[size])));
      ERRMEM (v[i] = malloc (sizeof(REAL[size])));
    }
    ERRMEM (rank = malloc (sizeof(int[size])));
    ERRMEM (tid = malloc (sizeof(unsigned int[size])));

    /* generate triangles and velocities */

    generate_triangles_and_velocities (lo, hi, nt, t, v, tid);
  }
  else
  {
    /* set nt */

    nt = 0;

    /* buffers */

    if (argc > 1) size = atoi (argv[1])*4;
    else size = 10000*4;

    for (i = 0; i < 3; i ++)
    {
      ERRMEM (t[0][i] = malloc (sizeof(REAL[size])));
      ERRMEM (t[1][i] = malloc (sizeof(REAL[size])));
      ERRMEM (t[2][i] = malloc (sizeof(REAL[size])));
      ERRMEM (v[i] = malloc (sizeof(REAL[size])));
    }
    ERRMEM (rank = malloc (sizeof(int[size])));
  }

  /* create load balancer */

  struct loba *lb = loba_create (ZOLTAN_RCB);

  /* perform time stepping */

  REAL step = 1E-3, time;

  for (time = 0.0; time < 1.0; time += step)
  {
    loba_balance (lb, nt, t[0], tid, 1.1, rank);

    migrate_triangles (size, nt, t, v, rank);

    integrate_triangles (step, lo, hi, nt, t, v);
  }

  /* finalise */

  loba_destroy (lb);

  for (i = 0; i < 3; i ++)
  {
    free (t[0][i]);
    free (t[1][i]);
    free (t[2][i]);
    free (v[i]);
  }
  free (rank);

  MPI_Finalize ();

  return 0;
}
