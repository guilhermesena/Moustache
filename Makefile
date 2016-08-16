run: 
	@if [ -f bin/moustache ] && [ -f ${f} ]; then mkdir -p out; ./bin/moustache <${f} >out/`basename $f`; fi
all: 
	@mkdir -p bin
	@g++ -c src/*.cpp -std=c++11
	@g++ -o bin/moustache moustache.o cellset.o cell.o
	@rm *.o




