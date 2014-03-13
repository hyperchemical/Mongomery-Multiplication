#Compiler Settings
PATH := /usr/um/gcc-4.7.0/bin:$(PATH)
LD_LIBRARY_PATH := /usr/um/gcc-4.7.0./lib64
LD_RUN_PATH := /usr/um/gcc-4.7.0/lib64

#Shortcuts
SFLAG = -std=c++11 -O3
TFLAG = 
DFLAG = -Wall -Wextra -pedantic -Wvla -g -ggdb  -pg
FILE = modex.cpp uberzahl.cpp
OUT= -o modex

#Submit
normal:
	g++ $(TFLAG) $(SFLAG) $(FILE) $(OUT)

clean:
	rm -rf *o modex
