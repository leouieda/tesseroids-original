BINPATH = ../../../bin

all: time_normal.txt

model.txt: 
	${BINPATH}/tessgrd -r0/359/-80/80 -b360/161 -z-10e03 -v -lmodgrid.log | ${BINPATH}/tessmodgen -s1/1 -d-1e03 -z0 -lmodgen.log > $@

grid.txt:
	${BINPATH}/tessgrd -r0/359/-90/90 -b360/181 -z50e03 -v -lgrid.log > $@

time_normal.txt: time_normal.sh model.txt grid.txt
	./$<

time_adapt.txt: time_adapt.sh model.txt grid.txt
	./$<

