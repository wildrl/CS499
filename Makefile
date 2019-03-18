TARGETS=solve_dpend dpend_mpfr dpend_out dpend_math

all: $(TARGETS)

orig: solve_dpend dpend_mpfr
	./solve_dpend 0.0 10.0 180.0 0.00 -10 0.0 2000 > outfile-original.txt

mpfr: dpend_mpfr
<<<<<<< HEAD
	./dpend_mpfr 0 100 90 0 90 0
=======
	./dpend_mpfr 0 50 90 0 90 0
>>>>>>> 422e65c398c33585b31111b2e18410e890a96293

solve_dpend: solve_dpend.c
	$(CC) -g -O3 -Wall -std=c99 -o $@ $< -lm

dpend_mpfr: dpend_out.c dpend_math.c dpend_mpfr.c
	$(CC) -g -O3 -Wall -std=c99 -D_BSD_SOURCE -o $@ $^ -lgmp -lmpfr

clean:
	rm -f $(TARGETS)
