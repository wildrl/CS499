TARGETS=solve_dpend dpend_mpfr

all: $(TARGETS)

orig: solve_dpend dpend_mpfr
	./solve_dpend 0.0 10.0 90.0 0.00 -10.0 0.0 1000 > outfile-original.txt

mpfr: dpend_mpfr
	./dpend_mpfr 0.0 10.0 90.0 0.00 -10.0 0.0 1000 8
	./dpend_mpfr 0.0 10.0 90.0 0.00 -10.0 0.0 1000 16
	./dpend_mpfr 0.0 10.0 90.0 0.00 -10.0 0.0 1000 32
	./dpend_mpfr 0.0 10.0 90.0 0.00 -10.0 0.0 1000 64
	./dpend_mpfr 0.0 10.0 90.0 0.00 -10.0 0.0 1000 128

solve_dpend: solve_dpend.c
	$(CC) -g -O3 -Wall -std=c99 -o $@ $< -lm

dpend_mpfr: dpend_mpfr.c
	$(CC) -g -O3 -Wall -std=c99 -o $@ $< -lgmp -lmpfr

clean:
	rm -f $(TARGETS)

