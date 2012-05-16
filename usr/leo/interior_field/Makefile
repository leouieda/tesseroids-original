BIN=../../../bin

all: Figure1.png Figure2.png

model.txt:
	echo "-0.5 0.5 -0.5 0.5 0 -5000 1000" > $@
	
output.txt: points.py model.txt
	python $< | \
	${BIN}/tessgxx model.txt | \
	${BIN}/tessgyy model.txt | \
	${BIN}/tessgzz model.txt -v > $@

Figure1.png: plot.py output.txt
	python $< output.txt
	mv output.png $@
	
output-a.txt: points.py model.txt
	python $< | \
	${BIN}/tessgxx model.txt -a | \
	${BIN}/tessgyy model.txt -a | \
	${BIN}/tessgzz model.txt -a -v > $@

Figure2.png: plot.py output-a.txt
	python $< output-a.txt
	mv output.png $@

clean:
	rm -f model.txt output*.txt Figure*.png
