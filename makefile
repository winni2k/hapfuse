all:
	make -C src

debug:
	make clean -C src
	make debug -C src

test:
	make test -C src

clean:
	make clean -C src

oxford:
	make oxford -C src
