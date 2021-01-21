#include "piv_ge_solver.h"
#include <stdlib.h>

//rozwiazuje uklad metoda eliminacji gaussa
int piv_ge_solver(matrix_t *eqs)
{
  if (eqs != NULL)
  {
    pivot_ge_in_situ_matrix(eqs);
    if (bs_matrix(eqs) == 0) //wsteczne podstawienie
    {
      return 0; //succes
    }
    else
    {
      return 1; //failure
    }
  }
  else
    return 1; //failure
}
