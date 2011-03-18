BINPATH = ../../../bin
PRISM_COORDS = 0 0 0
ORDER = 2/2/2

comparison.png: tessgrav.txt prismgrav.txt compare.py
	./compare.py

prismgrav.txt: prismgrav_local.txt local2global.py
	./local2global.py ${PRISM_COORDS} > $@

tessgrav.txt: tess.txt points.py
	./points.py | ${BINPATH}/tessgz $< -o${ORDER} | ${BINPATH}/tessgxx $< -o${ORDER} | ${BINPATH}/tessgxy $< -o${ORDER} | ${BINPATH}/tessgxz $< -o${ORDER} | ${BINPATH}/tessgyy $< -o${ORDER} | ${BINPATH}/tessgyz $< -o${ORDER} | ${BINPATH}/tessgzz $< -o${ORDER} > $@

prismgrav_local.txt: prism.txt points.py sph2prism.py
	./points.py | ./sph2prism.py ${PRISM_COORDS} | ${BINPATH}/prismgx $< | ${BINPATH}/prismgy $< | ${BINPATH}/prismgz $< | ${BINPATH}/prismgxx $< | ${BINPATH}/prismgxy $< | ${BINPATH}/prismgxz $< | ${BINPATH}/prismgyy $< | ${BINPATH}/prismgyz $< | ${BINPATH}/prismgzz $< > $@

prism.txt: tess.txt
	${BINPATH}/tess2prism < $< > $@
