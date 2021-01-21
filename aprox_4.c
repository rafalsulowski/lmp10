#include "makespl.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

double potega(double x, int n)
{ // Potęga x^n
    double sum = 1;
    for (int i = 0; i < n; i++)
    {
        sum *= x;
    }
    return sum;
}

double ValuePolymainFromX(double *A, double x)
{
    double sol = 0;
    double ax = *(A + 1) * x;              //x^1 * a1  - skladnik wielomianu
    double ax2 = *(A + 2) * x * x;         //x^2 * a2  - skladnik wielomianu
    double ax3 = *(A + 3) * x * x * x;     //x^3 * a3  - skladnik wielomianu
    double ax4 = *(A + 4) * x * x * x * x; //x^4 * a4  - skladnik wielomianu

    sol = *(A + 0) + ax + ax2 + ax3 + ax4; //wylicza wartosc wielomianu dla danego argumentu sumujac skladniki

    return sol;
}

void make_spl_4(points_t *pts, spline_t *spl)
{ //pts->n - 1 zwraca liczbe punktow wczytanych z pliku spl

    matrix_t *eqs = NULL;
    double *x = pts->x;       //wskazuje na tablice wspolrzednych x o rozmiarze pts->n
    double *y = pts->y;       //wskazuje na tablice wspolrzednych y o rozmiarze pts->n
    double a = x[0];          //poczatek przedzialu
    double b = x[pts->n - 1]; //koniec przedzialu
    int j, k;                 //zmienne pomocnicze
    double s = 0;

    //ustalamy rozmiar macierzy czyli stopien aproksymacji
    int nb = 5;
    char *nbEnv = getenv("APPROX_BASE_SIZE"); //pobiera rozmiar bazy

    if (nbEnv != NULL && atoi(nbEnv) > 0)
        nb = atoi(nbEnv); //ustawia ten rozmiar

    //tworzymy macierz o rozm 5x6 dla aproksymacji o bazie wiel. 4
    nb = 5; //wiel stopnia 4 (razem z wyrazem wolnym )
    eqs = make_matrix(nb, nb + 1);

    double xsum[10]; // tablica na potęgi xi^k
    double A[5][5];  // Macierz A
    double B[5];     // Macierz b

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
        //for (k = 0; k < pts->n; k++) //dodaje prawe strony macierzy, czyli macierz B
        add_to_entry_matrix(eqs, j, nb, B[j]);
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

    ///////////////////////////////////////////////

    int NumberOfPoints = pts->n; //liczba wczytanych punktow
    double argument = a;         //W(argument), bedzie sie on zwiekszal z kazda iteracja o tyle: (b - a) / NumberOfPoints

    //dla latwiejszego przekazywania wielomianu do funkcji tworze tablice wspolczynikow:
    double *A_wsp = (double *)malloc(sizeof(double) * 5);
    if (A_wsp)
        for (int i = 0; i < 5; i++)
            *(A_wsp + i) = get_entry_matrix(eqs, i, nb); //pobiera wartosc wspolczynnika ai
    else                                                 //obsluga bledu
    {
        printf("Nie udalo sie zaalokowac miejsca na tablice wspolczynnikow!\n");
        exit(1);
    }

    //allokuje miejsce dla punktow wynikowych
    //spl->x = (double *)malloc(NumberofPoints * sizeof(double)); //wlsciwie to chyba nie trzeba alokowac miejsca na argumenty
    spl->y = (double *)malloc((NumberofPoints + 1) * sizeof(double));
    spl->y + NumberofPoints + 1 = NULL;

    if (spl->y)
    {
        for (int i = 0; i < NumberOfPoints; i++)
        { //przypisuje wartoci y = W(argument)
            spl->y + i = ValuePolymainFromX(A_wsp, argument);
            argument += (b - a) / NumberOfPoints; //zwiekszam argument
        }
    }
    else //obsluga bledu
    {
        printf("Nie udalo sie zaalokowac miejsca na wynikowe punkty!\n");
        exit(1);
    }

    //////////////////////////////////////////////

    /*if (alloc_spl(spl, nb) == 0) //alokuje mniejsce
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
    }*/

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
