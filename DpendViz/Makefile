TARGETS=dpend_viz

all: $(TARGETS)

demo: dpend_viz
        ./dpend_viz 7
       
dpend_viz: dpend_viz.py
        $(python) $@

clean:
	rm -f $(TARGETS)
