TARGETS=solve_dpend dpend_mpfr

all: $(TARGETS)

demo: solve_dpend dpend_mpfr
	./solve_dpend 0.0 10.0 90.0 0.00 -10.0 0.0 1000 > outfile-original.txt
	./dpend_mpfr  0.0 10.0 90.0 0.00 -10.0 0.0 1000 64 > outfile-mpfr.txt

alldata: dpend_mpfr
	./dpend_mpfr 0.0 10.0 90.0 0.00 -10.0 0.0 1000 8 > mpfr8.txt
	./dpend_mpfr 0.0 10.0 90.0 0.00 -10.0 0.0 1000 16 > mpfr16.txt
	./dpend_mpfr 0.0 10.0 90.0 0.00 -10.0 0.0 1000 32 > mpfr32.txt
	./dpend_mpfr 0.0 10.0 90.0 0.00 -10.0 0.0 1000 64 > mpfr64.txt

solve_dpend: solve_dpend.c
	$(CC) -g -O3 -Wall -std=c99 -o $@ $< -lm

dpend_mpfr: dpend_mpfr.c
	$(CC) -g -O3 -Wall -std=c99 -o $@ $< -lgmp -lmpfr

clean:
	rm -f $(TARGETS)

