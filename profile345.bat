python -m cProfile -o profile345.pstats oscillators345.py
gprof2dot -f pstats profile345.pstats | dot -Tpdf -o graphics/profile345.pdf