all: tags ut

CPPFLAGS = -Wall -lm -DNDEBUG

testMatrix: testMatrix.cpp matrix.h
	apt-get install libunittest++-dev
	g++ $< $(CPPFLAGS) -lunittest++ -g -o $@

ut: testMatrix
	./testMatrix

tags: matrix.h image2D.h
	ctags $^

html: matrix.h image2D.h Doxyfile
	doxygen

clean:
	-rm -rf html *.o testMatrix

.PHONY: all clean ut
