all : moustache

run : 
	if [ -f bin/moustache ] && [ -f in/${FILE} ]; then ./bin/moustache <in/${FILE}; fi

moustache : 
	mkdir -p bin
	g++ -o bin/moustache src/moustache.cpp




