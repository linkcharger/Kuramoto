python -m cProfile -o profile12.pstats oscillators12.py
gprof2dot -f pstats profile12.pstats | dot -Tpdf -o graphics/profile12.pdf