TARGETS=solve_dpend dpend_mpfr dpend_out dpend_math

all: $(TARGETS)

orig: solve_dpend dpend_mpfr
	./solve_dpend 0.0 10.0 180.0 0.00 -10 0.0 2000 > outfile-original.txt

mpfr: dpend_mpfr
	./dpend_mpfr 0.0 10.0 180 0 -10 40 2000

solve_dpend: solve_dpend.c
	$(CC) -g -O3 -Wall -std=c99 -o $@ $< -lm

dpend_mpfr: dpend_out.c dpend_math.c dpend_mpfr.c
	$(CC) -g -O3 -Wall -std=c99 -D_BSD_SOURCE -o $@ $^ -lgmp -lmpfr

clean:
	rm -f $(TARGETS)
