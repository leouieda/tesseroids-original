BINPATH = ../../bin

all: parana-ggt-10sec.png parana-dem-10sec-density.png parana-dem-10sec.png


parana-dem-10sec.png: plot_dem.py parana-dem-10sec.xyz
	@echo "\nPLOTTING DEM\n"
	python $<

parana-dem-10sec-density.txt: dem_density.py parana-dem-10sec.xyz
	@echo "\nASSIGNING DENSITY VALUES TO DEM\n"
	python $< > $@

parana-dem-10sec-density.png: plot_density.py parana-dem-10sec-density.txt
	@echo "\nPLOTTING DENSITY\n"
	python $<

dem-tess-model.txt: parana-dem-10sec-density.txt
	@echo "\nCONVERTING DEM TO TESSEROIDS\n"
	@${BINPATH}/tessmodgen -v -s0.166666667/0.166666667 -z0 < $< > $@

parana-ggt-10sec.txt: dem-tess-model.txt
	@echo "\nCALCULATING GGT OF DEM MODEL AT GOCE HEIGHT\n"
	@(time ${BINPATH}/tessgrd -r-60/-45/-30/-15 -b100/100 -z250e03 | ${BINPATH}/tessgxx $< -a | \
	${BINPATH}/tessgxy $< -a | ${BINPATH}/tessgxz $< -a | ${BINPATH}/tessgyy $< -a | \
	${BINPATH}/tessgyz $< -a | ${BINPATH}/tessgzz $< -a -v > $@) 2> ggt_calc_time.txt
	
parana-ggt-10sec.png: plot_ggt.py parana-ggt-10sec.txt
	@echo "\nPLOTTING GGT\n"
	python $<
	
