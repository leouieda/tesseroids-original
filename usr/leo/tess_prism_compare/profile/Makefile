BINPATH = ../../../../bin
ORDER = -o2/2/2


comparison.png: compare.py prismgrav.txt tessgrav.txt tess.txt
	./$< prismgrav.txt tessgrav.txt tess.txt
	
prismgrav.txt: ../local2global.py prismgrav_local.txt tessgrav.txt
	./$< prismgrav_local.txt tessgrav.txt > $@

tessgrav.txt: tess.txt points.py
	./points.py $< | ${BINPATH}/tesspot $< ${ORDER} | ${BINPATH}/tessgz $< ${ORDER} | ${BINPATH}/tessgxx $< ${ORDER} | ${BINPATH}/tessgxy $< ${ORDER} | ${BINPATH}/tessgxz $< ${ORDER} | ${BINPATH}/tessgyy $< ${ORDER} | ${BINPATH}/tessgyz $< ${ORDER} | ${BINPATH}/tessgzz $< ${ORDER} > $@

prismgrav_local.txt: prism.txt points.py ../sph2prism.py
	./points.py tess.txt | ./../sph2prism.py | ${BINPATH}/prismpot $< | ${BINPATH}/prismgx $< | ${BINPATH}/prismgy $< | ${BINPATH}/prismgz $< | ${BINPATH}/prismgxx $< | ${BINPATH}/prismgxy $< | ${BINPATH}/prismgxz $< | ${BINPATH}/prismgyy $< | ${BINPATH}/prismgyz $< | ${BINPATH}/prismgzz $< > $@

prism.txt: tess.txt
	${BINPATH}/tess2prism < $< > $@

clean:
	rm prism.txt tessgrav.txt prismgrav_local.txt prismgrav.txt
