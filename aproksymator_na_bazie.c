#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE
*/

/*
 * Funkcje bazowe: n - stopein wielomianu (baza) a,b - granice przedzialu aproksymacji 
 * i - numer funkcji x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
 */
double
fi(double a, double b, int n, int i, double x)
{
	double h = (b - a) / (n - 1); //co ile rosnie x
	double h3 = h * h * h;		  //h^3
	int hi[5] = {i - 2, i - 1, i, i + 1, i + 2};
	double hx[5];
	int j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j]; //wartosc wspolrzednej x

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 1 / h3 * (x - hx[0]) * (x - hx[0]) * (x - hx[0]);
	else if (x > hx[1] && x <= hx[2])
		return 1 / h3 * (h3 + 3 * h * h * (x - hx[1]) + 3 * h * (x - hx[1]) * (x - hx[1]) - 3 * (x - hx[1]) * (x - hx[1]) * (x - hx[1]));
	else if (x > hx[2] && x <= hx[3])
		return 1 / h3 * (h3 + 3 * h * h * (hx[3] - x) + 3 * h * (hx[3] - x) * (hx[3] - x) - 3 * (hx[3] - x) * (hx[3] - x) * (hx[3] - x));
	else /* if (x > hx[3]) && (x <= hx[4]) */
		return 1 / h3 * (hx[4] - x) * (hx[4] - x) * (hx[4] - x);
}

/* Pierwsza pochodna fi */
double
dfi(double a, double b, int n, int i, double x)
{
	double h = (b - a) / (n - 1);
	double h3 = h * h * h;
	int hi[5] = {i - 2, i - 1, i, i + 1, i + 2};
	double hx[5];
	int j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 3 / h3 * (x - hx[0]) * (x - hx[0]);
	else if (x > hx[1] && x <= hx[2])
		return 1 / h3 * (3 * h * h + 6 * h * (x - hx[1]) - 9 * (x - hx[1]) * (x - hx[1]));
	else if (x > hx[2] && x <= hx[3])
		return 1 / h3 * (-3 * h * h - 6 * h * (hx[3] - x) + 9 * (hx[3] - x) * (hx[3] - x));
	else /* if (x > hx[3]) && (x <= hx[4]) */
		return -3 / h3 * (hx[4] - x) * (hx[4] - x);
}

/* Druga pochodna fi */
double
d2fi(double a, double b, int n, int i, double x)
{
	double h = (b - a) / (n - 1);
	double h3 = h * h * h;
	int hi[5] = {i - 2, i - 1, i, i + 1, i + 2};
	double hx[5];
	int j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 6 / h3 * (x - hx[0]);
	else if (x > hx[1] && x <= hx[2])
		return 1 / h3 * (6 * h - 18 * (x - hx[1]));
	else if (x > hx[2] && x <= hx[3])
		return 1 / h3 * (6 * h - 18 * (hx[3] - x));
	else /* if (x > hx[3]) && (x <= hx[4]) */
		return 6 / h3 * (hx[4] - x);
}

/* Trzecia pochodna fi */
double
d3fi(double a, double b, int n, int i, double x)
{
	double h = (b - a) / (n - 1);
	double h3 = h * h * h;
	int hi[5] = {i - 2, i - 1, i, i + 1, i + 2};
	double hx[5];
	int j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 6 / h3;
	else if (x > hx[1] && x <= hx[2])
		return -18 / h3;
	else if (x > hx[2] && x <= hx[3])
		return 18 / h3;
	else /* if (x > hx[3]) && (x <= hx[4]) */
		return -6 / h3;
}

/* Pomocnicza f. do rysowania bazy */
double
xfi(double a, double b, int n, int i, FILE *out)
{
	double h = (b - a) / (n - 1);
	double h3 = h * h * h;
	int hi[5] = {i - 2, i - 1, i, i + 1, i + 2};
	double hx[5];
	int j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	fprintf(out, "# nb=%d, i=%d: hi=[", n, i);
	for (j = 0; j < 5; j++)
		fprintf(out, " %d", hi[j]);
	fprintf(out, "] hx=[");
	for (j = 0; j < 5; j++)
		fprintf(out, " %g", hx[j]);
	fprintf(out, "]\n");
}

void make_spl(points_t *pts, spline_t *spl)
{ //pts->n - 1 zwraca liczbe punktow wczytanych z pliku spl

	matrix_t *eqs = NULL;
	double *x = pts->x;		  //wskazuje na tablice wspolrzednych x o rozmiarze pts->n
	double *y = pts->y;		  //wskazuje na tablice wspolrzednych y o rozmiarze pts->n
	double a = x[0];		  //poczatek przedzialu
	double b = x[pts->n - 1]; //koniec przedzialu
	int i, j, k;			  //zmienne pomocnicze

	//ustalamy rozmiar macierzy czyli stopien aproksymacji
	int nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
	char *nbEnv = getenv("APPROX_BASE_SIZE"); //pobiera rozmiar bazy

	if (nbEnv != NULL && atoi(nbEnv) > 0)
		nb = atoi(nbEnv); //ustawia ten rozmiar

	//tworzymy macierz o rozm 5x6 dla aproksymacji o bazie wiel. 4
	eqs = make_matrix(nb, nb + 1);

#ifdef DEBUG
#define TESTBASE 500
	{
		FILE *tst = fopen("debug_base_plot.txt", "w");
		double dx = (b - a) / (TESTBASE - 1);
		for (j = 0; j < nb; j++)
			xfi(a, b, nb, j, tst);
		for (i = 0; i < TESTBASE; i++)
		{
			fprintf(tst, "%g", a + i * dx);
			for (j = 0; j < nb; j++)
			{
				fprintf(tst, " %g", fi(a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", dfi(a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", d2fi(a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", d3fi(a, b, nb, j, a + i * dx));
			}
			fprintf(tst, "\n");
		}
		fclose(tst);
	}
#endif

	//uzupelnia macierz wartosciami (sumy szeregow)
	for (j = 0; j < nb; j++)
	{
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, fi(a, b, nb, i, x[k]) * fi(a, b, nb, j, x[k]));
		//eqs - macierz, i j indeksy, fi() * fi() wartosc

		//n - ilosc wczytanych punktow
		for (k = 0; k < pts->n; k++) //dodaje prawe strony macierzy, czyli macierz B
			add_to_entry_matrix(eqs, j, nb, y[k] * fi(a, b, nb, j, x[k]));
		//eqs - macierz, j indeks wiersza, nb - indeks kolumny (kolumna z prawymi stronami)
		//y[k] - wartosc wspolrzednej y danego punktu pomnozona razy wartosc funkcji
	}

#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	//rozwiazanie ukladu rownan w postaci macierzowej metoda eliminacji gaussa i wstecznego podstawienia
	if (piv_ge_solver(eqs))
	{
		spl->n = 0;
		return;
	}

	//terazm mamy juz odpowiedzi (wartosci wspolczynikow a0, a1, ...) w ostatniej kolumnie macierzy esq

#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	if (alloc_spl(spl, nb) == 0) //alokuje mniejsce
	{
		for (i = 0; i < spl->n; i++)
		{
			double xx = spl->x[i] = a + i * (b - a) / (spl->n - 1);
			xx += 10.0 * DBL_EPSILON; // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++)
			{
				double ck = get_entry_matrix(eqs, k, nb); //pobiera element o indeksie k,nb czyli element z wyliczoną wartoscią współczynika a0, a1 itd...
				spl->f[i] += ck * fi(a, b, nb, k, xx);
				spl->f1[i] += ck * dfi(a, b, nb, k, xx);
				spl->f2[i] += ck * d2fi(a, b, nb, k, xx);
				spl->f3[i] += ck * d3fi(a, b, nb, k, xx);
			}
		}
	}

#ifdef DEBUG
	{
		FILE *tst = fopen("debug_spline_plot.txt", "w");
		double dx = (b - a) / (TESTBASE - 1);
		for (i = 0; i < TESTBASE; i++)
		{
			double yi = 0;
			double dyi = 0;
			double d2yi = 0;
			double d3yi = 0;
			double xi = a + i * dx;
			for (k = 0; k < nb; k++)
			{
				yi += get_entry_matrix(eqs, k, nb) * fi(a, b, nb, k, xi);
				dyi += get_entry_matrix(eqs, k, nb) * dfi(a, b, nb, k, xi);
				d2yi += get_entry_matrix(eqs, k, nb) * d2fi(a, b, nb, k, xi);
				d3yi += get_entry_matrix(eqs, k, nb) * d3fi(a, b, nb, k, xi);
			}
			fprintf(tst, "%g %g %g %g %g\n", xi, yi, dyi, d2yi, d3yi);
		}
		fclose(tst);
	}
#endif
}

double potega(double x, int n)
{ // Potęga x^n
	double sum = 1;
	for (int i = 0; i < n; i++)
	{
		sum *= x;
	}
	return sum;
}

void make_spl_4(points_t *pts, spline_t *spl)
{ //pts->n - 1 zwraca liczbe punktow wczytanych z pliku spl

	matrix_t *eqs = NULL;
	double *x = pts->x;		  //wskazuje na tablice wspolrzednych x o rozmiarze pts->n
	double *y = pts->y;		  //wskazuje na tablice wspolrzednych y o rozmiarze pts->n
	double a = x[0];		  //poczatek przedzialu
	double b = x[pts->n - 1]; //koniec przedzialu
	int i, j, k;			  //zmienne pomocnicze
	double s = 0;

	//ustalamy rozmiar macierzy czyli stopien aproksymacji
	int nb = 4;
	char *nbEnv = getenv("APPROX_BASE_SIZE"); //pobiera rozmiar bazy

	if (nbEnv != NULL && atoi(nbEnv) > 0)
		nb = atoi(nbEnv); //ustawia ten rozmiar

	//tworzymy macierz o rozm 5x6 dla aproksymacji o bazie wiel. 4
	nb = 4; //wiel stopnia 4
	eqs = make_matrix(nb, nb + 1);

	double xsum[10]; // tablica na potęgi xi^k
	double A[5][5];	 // Macierz A
	double B[5];	 // Macierz b
	double Xodp[5];	 //macierz na rozwiazania

	// uzupełnienie xsum
	for (int i = 0; i < 10; i++)
	{
		double s = 0;
		for (int j = 0; j < 30; j++)
		{
			s += potega(*(x + j), i);
		}
		xsum[i] = s;
	}

	// Uzupełnienie macierzy A i b
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			A[i][j] = xsum[i + j];
		}
		double s = 0;
		for (int j = 0; j < 30; j++)
		{
			s += potega(*(x + j), i) * *(y + j);
		}
		B[i] = s;
	}

	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			printf("%lf , ", A[i][j]);
		}
		printf("\n");
	}

	printf("\n");
	for (int i = 0; i < 5; i++)
		printf("%lf\n", B[i]);

	//uzupelnia macierz wartosciami (sumy szeregow)
	for (j = 0; j < nb; j++)
	{
		for (i = 0; i < nb; i++)
			add_to_entry_matrix(eqs, j, i, A[i][j]);
		//eqs - macierz, i -> (id kolumny) j -> (id wiersza)

		//n - ilosc wczytanych punktow
		for (k = 0; k < pts->n; k++) //dodaje prawe strony macierzy, czyli macierz B
			add_to_entry_matrix(eqs, j, nb, B[i]);
		//eqs - macierz, j indeks wiersza, nb - indeks kolumny (kolumna z prawymi stronami)
		//y[k] - wartosc wspolrzednej y danego punktu pomnozona razy wartosc funkcji
	}

#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	//rozwiazanie ukladu rownan w postaci macierzowej metoda eliminacji gaussa i wstecznego podstawienia
	if (piv_ge_solver(eqs))
	{
		spl->n = 0;
		return;
	}

	//terazm mamy juz odpowiedzi (wartosci wspolczynikow a0, a1, ...) w ostatniej kolumnie macierzy esq

#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	if (alloc_spl(spl, nb) == 0) //alokuje mniejsce
	{
		for (i = 0; i < spl->n; i++)
		{
			double xx = spl->x[i] = a + i * (b - a) / (spl->n - 1);
			xx += 10.0 * DBL_EPSILON; // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++)
			{
				double ck = get_entry_matrix(eqs, k, nb); //pobiera element o indeksie k,nb czyli element z wyliczoną wartoscią współczynika a0, a1 itd...
				spl->f[i] += ck * fi(a, b, nb, k, xx);
				spl->f1[i] += ck * dfi(a, b, nb, k, xx);
				spl->f2[i] += ck * d2fi(a, b, nb, k, xx);
				spl->f3[i] += ck * d3fi(a, b, nb, k, xx);
			}
		}
	}

#ifdef DEBUG
	{
		FILE *tst = fopen("debug_spline_plot.txt", "w");
		double dx = (b - a) / (TESTBASE - 1);
		for (i = 0; i < TESTBASE; i++)
		{
			double yi = 0;
			double dyi = 0;
			double d2yi = 0;
			double d3yi = 0;
			double xi = a + i * dx;
			for (k = 0; k < nb; k++)
			{
				yi += get_entry_matrix(eqs, k, nb) * fi(a, b, nb, k, xi);
				dyi += get_entry_matrix(eqs, k, nb) * dfi(a, b, nb, k, xi);
				d2yi += get_entry_matrix(eqs, k, nb) * d2fi(a, b, nb, k, xi);
				d3yi += get_entry_matrix(eqs, k, nb) * d3fi(a, b, nb, k, xi);
			}
			fprintf(tst, "%g %g %g %g %g\n", xi, yi, dyi, d2yi, d3yi);
		}
		fclose(tst);
	}
#endif
}
