FLAGS = -c

run: 
	@if [ -f bin/moustache ] && [ -f ${f} ]; then mkdir -p out; ./bin/moustache <${f} >out/`basename $f`; fi

debug: FLAGS += -DDEBUG
debug: all

all: 
	@mkdir -p bin
	@g++ ${FLAGS} src/*.cpp -std=c++11
	@g++ -o bin/moustache moustache.o cellset.o cell.o
	@rm *.o




