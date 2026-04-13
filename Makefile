all: 
	cd build && make

clean:
	cd build && make clean

cleanobj:
	cd build && make cleanobj

dist: 
	make clean
	mkdir basic_nm_codes
	cp -r LICENSE Makefile README build include lib src basic_nm_codes
	tar -czvf basic_nm_codes.tar.gz basic_nm_codes
	rm -rf basic_nm_codes

.PHONY: distclean
distclean:
	rm -f basic_nm_codes.tar.gz
