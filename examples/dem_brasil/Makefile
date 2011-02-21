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

parana-ggt-10sec.txt: dem-tess-model.txt calc_ggt.sh 
	@echo "\nCALCULATING GGT OF DEM MODEL AT GOCE HEIGHT\n"
	@/usr/bin/time -v -o ggt_exec_time.txt ./calc_ggt.sh $< $@
	
parana-ggt-10sec.png: plot_ggt.py parana-ggt-10sec.txt
	@echo "\nPLOTTING GGT\n"
	python $<
	
