TARGETS=solve_dpend dpend_mpfr dpend_out dpend_out_hdr

all: $(TARGETS)

orig: solve_dpend dpend_mpfr
	./solve_dpend 0.0 10.0 90.0 0.00 -10.0 0.0 1000 > outfile-original.txt

mpfr: dpend_mpfr
	./dpend_mpfr 0.0 20.0 180.45 0 16.25 0 8000

solve_dpend: solve_dpend.c
	$(CC) -g -O3 -Wall -std=c99 -o $@ $< -lm

dpend_mpfr: dpend_out.c dpend_mpfr.c
	$(CC) -g -O3 -Wall -std=c99 -D_BSD_SOURCE -o $@ $^ -lgmp -lmpfr

clean:
	rm -f $(TARGETS)

