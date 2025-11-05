grit: grit.o perm.o perm_etc.o chain.o path.o gamma_dev.o path_tree_etc.o \
	 		histogram.o get_CPU_time.o rngs.o signs.o gr_misc.o path_etc.o run_chains.o
	gcc   -Wall -g -O3  grit.o perm.o perm_etc.o chain.o run_chains.o path.o \
			gamma_dev.o path_tree_etc.o histogram.o get_CPU_time.o rngs.o signs.o gr_misc.o path_etc.o -o grit  -lm -lmcheck

grit.o: grit.c grit.h perm.h perm_etc.h structs.h chain.h rngs.h run_chains.h path_etc.h path_tree_etc.h histogram.h params.h
	gcc -c   -Wall -g -O3 grit.c

run_chains.o: run_chains.c grit.h perm.h perm_etc.h path.h exts.h structs.h chain.h \
			run_chains.h gamma_dev.h rngs.h signs.h gr_misc.h path_etc.h path_tree_etc.h histogram.h params.h
	gcc -c   -Wall -g -O3 run_chains.c

chain.o: chain.c grit.h perm.h perm_etc.h grit.h path.h path_tree_etc.h histogram.h exts.h structs.h \
			chain.h gamma_dev.h rngs.h signs.h gr_misc.h path_etc.h params.h
	gcc -c   -Wall -g -O3 chain.c 

perm.o: perm.c perm.h perm_etc.h exts.h structs.h grit.h rngs.h params.h gr_misc.h
	gcc -c   -Wall -g -O3 perm.c  

perm_etc.o: perm_etc.c perm_etc.h exts.h structs.h grit.h params.h
	gcc -c   -Wall -g -O3 perm_etc.c  

path.o: path.c perm.h perm_etc.h exts.h structs.h path.h grit.h rngs.h signs.h gr_misc.h path_etc.h params.h
	gcc -c  -Wall -g -O3 path.c 

path_etc.o: path_etc.c perm.h perm_etc.h exts.h structs.h path.h grit.h path_etc.h gr_misc.h chain.h params.h
	gcc -c  -Wall -g -O3 path_etc.c 

signs.o: signs.c grit.h structs.h exts.h signs.h rngs.h gr_misc.h perm.h params.h
	gcc -c  -Wall -g -O3 signs.c

gamma_dev.o: gamma_dev.c gamma_dev.h grit.h exts.h rngs.h
	gcc -c  -Wall -g -O3 gamma_dev.c

path_tree_etc.o: path_tree_etc.c structs.h grit.h exts.h path.h path_tree_etc.h histogram.h gr_misc.h
	gcc -c  -Wall -g -O3 path_tree_etc.c

histogram.o: histogram.c structs.h grit.h exts.h path.h histogram.h gr_misc.h
	gcc -c  -Wall -g -O3 histogram.c

get_CPU_time.o: get_CPU_time.c
	gcc -c  -Wall -g -O3 get_CPU_time.c 

rngs.o: rngs.c rngs.h
	gcc -c  -Wall -g -O3 rngs.c

gr_misc.o: gr_misc.c gr_misc.h
	gcc -c  -Wall -g -O3 gr_misc.c
