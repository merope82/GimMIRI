.PHONY: all clean tgz

all:
	@make --no-print-directory -C src

clean:
	@make --no-print-directory -C src clean

nuke:
	@make --no-print-directory -C src nuke

tgz:
	tar cf GimMIRI.tar ../GimMIRI
	gzip GimMIRI.tar
                                              
                                       