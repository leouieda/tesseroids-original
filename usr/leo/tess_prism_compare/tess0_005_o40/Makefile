BINPATH = ../../../../bin
ORDER = -o2/2/2 -a
HIGHORDER = -o40/40/40 -a


comparison.png: ../compare.py differences.txt tess.txt
	./$< differences.txt tess.txt

differences.txt: ../calc_diff.py tessgrav.txt tessgrav_highorder.txt 
	./$< tessgrav.txt tessgrav_highorder.txt > $@

tessgrav.txt: tess.txt ../points.py
	./../points.py $< | ${BINPATH}/tesspot $< ${ORDER} | ${BINPATH}/tessgz $< ${ORDER} | ${BINPATH}/tessgxx $< ${ORDER} | ${BINPATH}/tessgxy $< ${ORDER} | ${BINPATH}/tessgxz $< ${ORDER} | ${BINPATH}/tessgyy $< ${ORDER} | ${BINPATH}/tessgyz $< ${ORDER} | ${BINPATH}/tessgzz $< ${ORDER} > $@

tessgrav_highorder.txt: tess.txt ../points.py
	./../points.py $< | ${BINPATH}/tesspot $< ${HIGHORDER} | ${BINPATH}/tessgz $< ${HIGHORDER} | ${BINPATH}/tessgxx $< ${HIGHORDER} | ${BINPATH}/tessgxy $< ${HIGHORDER} | ${BINPATH}/tessgxz $< ${HIGHORDER} | ${BINPATH}/tessgyy $< ${HIGHORDER} | ${BINPATH}/tessgyz $< ${HIGHORDER} | ${BINPATH}/tessgzz $< ${HIGHORDER} > $@
