all : moustache

run : 
	if [ -f bin/moustache ] && [ -f ${f} ]; then mkdir -p out; ./bin/moustache <${f} >out/`basename $f`; fi

moustache : 
	mkdir -p bin
	g++ -o bin/moustache src/moustache.cpp




