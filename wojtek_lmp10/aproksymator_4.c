#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#define DEBUG

#define TESTBASE 500

double potega(double x, int n)
{ // Potęga x^n
	double sum = 1;
	for (int i = 0; i < n; i++)
	{
		sum *= x;
	}
	return sum;
}

void make_xsum(double t[], double *x, int n)
{
	double s;
	for (int i = 0; i < 10; i++)
	{
		s = 0;
		for (int j = 0; j < n; j++)
		{
			s += potega(*(x + j), i);
		}
		t[i] = s;
	}
}

void make_AB(double a[][5], double b[], double *x, double *y, double xsum[])
{
	double s = 0;
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
			a[i][j] = xsum[i + j];
		s = 0;
		for (int j = 0; j < 30; j++)
		{
			s += potega(*(x + j), i) * *(y + j);
		}
		b[i] = s;
	}
}

void make_spl_4(points_t *pts, spline_t *spl, double *Tab_A)
{
	//pts->n - 1 zwraca liczbe punktow wczytanych z pliku spl

	matrix_t *eqs = NULL;
	double *x = pts->x; //wskazuje na tablice wspolrzednych x o rozmiarze pts->n
	double *y = pts->y; //wskazuje na tablice wspolrzednych y o rozmiarze pts->n
	int i, j, k;		//zmienne pomocnicze
	int num_p = pts->n;
	double s = 0;

	//ustalamy rozmiar macierzy czyli stopien aproksymacji
	int nb = 5;								  // wielomian stopnia 4
	char *nbEnv = getenv("APPROX_BASE_SIZE"); //pobiera rozmiar bazy

	if (nbEnv != NULL && atoi(nbEnv) > 0)
		nb = atoi(nbEnv); //ustawia ten rozmiar

	//tworzymy macierz o rozm 5x6 dla aproksymacji o bazie wiel. 4
	eqs = make_matrix(nb, nb + 1);

	double xsum[10]; // tablica na potęgi xi^k
	double A[5][5];	 // Macierz A
	double B[5];	 // Macierz b

	make_xsum(xsum, x, num_p); // uzupełnienie xsum
	make_AB(A, B, x, y, xsum); // Uzupełnienie macierzy A i Bi

	/* Wypisanie Macierzy A i B */
	/*
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
	*/
	//uzupelnia macierz wartosciami (sumy szeregow)
	for (j = 0; j < nb; j++)
	{
		for (i = 0; i < nb; i++)
			add_to_entry_matrix(eqs, j, i, A[i][j]);
		//eqs - macierz, i -> (id kolumny) j -> (id wiersza)

		//n - ilosc wczytanych punktow
		add_to_entry_matrix(eqs, j, nb, B[j]);
		//eqs - macierz, j indeks wiersza, nb - indeks kolumny (kolumna z prawymi stronami)
		//y[k] - wartosc wspolrzednej y danego punktu pomnozona razy wartosc funkcji
	}
	write_matrix(eqs, stdout);

#ifdef DEBUG
	//write_matrix(eqs, stdout);
#endif

	//rozwiazanie ukladu rownan w postaci macierzowej metoda eliminacji gaussa i wstecznego podstawienia
	if (piv_ge_solver(eqs))
	{
		spl->n = 0;
		return;
	}

	//przypisanie wartosci wspolczynikow do zewnetrznej tablicy
	for (int i = 0; i < 5; i++)
	{
		*(Tab_A + i) = get_entry_matrix(eqs, i, nb);
	}

	//terazm mamy juz odpowiedzi (wartosci wspolczynikow a0, a1, ...) w ostatniej kolumnie macierzy esq
	write_matrix(eqs, stdout);
}
