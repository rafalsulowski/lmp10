#include "points.h"
#include "splines.h"
#include "makespl.h"

#define DEBUG 1

//#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

char *usage =
    "Usage: %s -s spline-file [-p points-file] [ -g gnuplot-file [-f from_x -t to_x -n n_points ] ]\n"
    "            if points-file is given then\n"
    "               reads discrete 2D points from points-file\n"
    "               writes spline approximation to spline-file\n"
    "               - number of points should be >= 4\n"
    "            else (points-file not given)\n"
    "               reads spline from spline-file\n"
    "            endfi\n"
    "            if gnuplot-file is given then\n"
    "               makes table of n_points within <from_x,to_x> range\n"
    "               - from_x defaults to x-coordinate of the first point in points-file,\n"
    "               - to_x defaults to x-coordinate of the last point\n"
    "               - n_points defaults to 100\n"
    "               - n_points must be > 1\n"
    "            endif\n";

int main(int argc, char **argv)
{
  int opt;
  char *inp = (char *)malloc(sizeof(char) * 50);
  char *out = (char *)malloc(sizeof(char) * 50);
  char *gpt = (char *)malloc(sizeof(char) * 50);
  double fromX = 0;
  double toX = 0;
  int n = 100;
  char *progname = argv[0];

  points_t pts;
  spline_t spl;

  pts.n = 0;
  spl.n = 0;

  /*//process options, save user choices
  while ((opt = getopt(argc, argv, "p:s:g:f:t:n:")) != -1)
  {
    switch (opt)
    {
    case 'p':
      inp = optarg;
      break;
    case 's':
      out = optarg;
      break;
    case 'g':
      gpt = optarg;
      break;
    case 'f':
      fromX = atof(optarg);
      break;
    case 't':
      toX = atof(optarg);
      break;
    case 'n':
      n = atoi(optarg);
      break;
    default:
      fprintf(stderr, usage, progname);
      exit(EXIT_FAILURE);
    }
  }
  if (optind < argc)
  {
    fprintf(stderr, "\nBad parameters!\n");
    for (; optind < argc; optind++)
      fprintf(stderr, "\t\"%s\"\n", argv[optind]);
    fprintf(stderr, "\n");
    fprintf(stderr, usage, progname);
    exit(EXIT_FAILURE);
  }*/

  printf("\n%s\n\n", usage);
  strcpy(inp, argv[4]);
  strcpy(out, argv[2]);
  strcpy(gpt, argv[6]);
  fromX = atof(argv[8]);
  toX = atof(argv[10]);
  n = atoi(argv[12]);

  /* if points-file was given, then read points, generate spline, save it to file */
  if (inp != NULL)
  {
    FILE *OUT = NULL; /* we shall open it later, when we shall get points */

    FILE *in = fopen(inp, "r");
    if (in == NULL)
    {
      fprintf(stderr, "%s: can not read points file: %s\n\n", argv[0], inp);
      exit(EXIT_FAILURE);
    }

    if (read_pts_failed(in, &pts))
    {
      fprintf(stderr, "%s: bad contents of points file: %s\n\n", argv[0],
              inp);
      exit(EXIT_FAILURE);
    }
    else
      fclose(in);

    OUT = fopen(out, "w");
    if (OUT == NULL)
    {
      fprintf(stderr, "%s: can not write spline file: %s\n\n", argv[0], OUT);
      exit(EXIT_FAILURE);
    }

    //tworzenie pliku spl
    //wlasny aproksymator!
    make_spl_4(&pts, &spl);

    if (spl.n > 0)
      write_spl(&spl, OUT);

    fclose(OUT);
  }
  else if (out != NULL)
  { /* if point-file was NOT given, try to read splines from a file */
    FILE *splf = fopen(out, "r");
    if (splf == NULL)
    {
      fprintf(stderr, "%s: can not read spline file: %s\n\n", argv[0], inp);
      exit(EXIT_FAILURE);
    }
    if (read_spl(splf, &spl))
    {
      fprintf(stderr, "%s: bad contents of spline file: %s\n\n", argv[0],
              inp);
      exit(EXIT_FAILURE);
    }
  }
  else
  { /* ponts were not given nor spline was given -> it is an error */
    fprintf(stderr, usage, argv[0]);
    exit(EXIT_FAILURE);
  }

  if (spl.n < 1)
  { /* check if there is a valid spline */
    fprintf(stderr, "%s: bad spline: n=%d\n\n", argv[0], spl.n);
    exit(EXIT_FAILURE);
  }

  /* check if plot was requested and generate it if yes */
  if (gpt != NULL && n > 1)
  {
    FILE *gpf = fopen(gpt, "w");
    int i;
    double dx;
    if (fromX == 0 && toX == 0)
    { /* calculate plot range if it was not specified */
      if (pts.n > 1)
      {
        fromX = pts.x[0];
        toX = pts.x[pts.n - 1];
      }
      else if (spl.n > 1)
      {
        fromX = spl.x[0];
        toX = spl.x[spl.n - 1];
      }
      else
      {
        fromX = 0;
        toX = 1;
      }
    }
    dx = (toX - fromX) / (n - 1); //podzialka, przeskok co jaka wartosc

    if (gpf == NULL)
    {
      fprintf(stderr, "%s: can not write gnuplot file: %s\n\n", argv[0],
              gpt);
      exit(EXIT_FAILURE);
    }

    for (i = 0; i < n; i++)
      fprintf(gpf, "%g %g\n", fromX + i * dx,
              value_spl(&spl, fromX + i * dx));

    fclose(gpf);
  }

  return 0;
}
