BINPATH = ../../../bin

all: time_normal.txt time_adapt.txt

model.txt: 
	@echo "\nGENERATING MODEL\n"
	${BINPATH}/tessgrd -r-60/60/-60/60 -b121/121 -z-10e03 -v -lmodgrid.log | ${BINPATH}/tessmodgen -s1/1 -d-1e03 -z0 -lmodgen.log > $@

grid.txt:
	@echo "\nGENERATING GRID\n"
	${BINPATH}/tessgrd -r-80/80/-80/80 -b160/160 -z50e03 -v -lgrid.log > $@

time_normal.txt: time_normal.sh model.txt grid.txt
	@echo "\nTIMING NORMAL EXECUTION\n"
	./$< -o5/5/2

time_adapt.txt: time_adapt.sh model.txt grid.txt
	@echo "\nTIMING ADAPTATIVE EXECUTION\n"
	./$<

clean:
	rm *.log model.txt grid.txt normal_data.txt adapt_data.txt
