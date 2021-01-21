CC = GCC

aprox: main.o splines.o points.o aprox_4.o gaus/libge.a
	gcc -o aprox -Wall main.o splines.o points.o aprox_4.o -L gaus -l ge

intrp: main.o splines.o points.o interpolator.o gaus/libge.a
	gcc -o intrp  main.o splines.o points.o interpolator.o -L gaus -l ge

prosta: main.o splines.o points.o prosta.o
	gcc -o prosta  main.o splines.o points.o prosta.o	

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h
	gcc -I gaus -c aproksymator_na_bazie.c

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	gcc -I gaus -c interpolator.c

.PHONY: clean 

clean:
	-rm *.o aprox intrp prosta
