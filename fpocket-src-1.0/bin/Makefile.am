# Compiler and linkers:
CFLAGS = -O
LDFLAGS = -fno-underscoring -lm
 
# Targets: bin means that programs will be installed in /usr/local/bin
# I guess the use of --prefix PATH when launching ./configure changes
# the prefix path from /usr/local to PATH
bin_PROGRAMS = fpocket tpocket dpocket
 
# Dependances
fpocket_SOURCES =   ../src/fpmain.c ../src/qhull/qvoronoi.c ../src/qhull/user.c ../src/qhull/global.c \
					../src/qhull/io.c ../src/qhull/geom2.c ../src/qhull/poly2.c \
					../src/qhull/merge.c ../src/qhull/qhull.c ../src/qhull/geom.c  \
					../src/qhull/poly.c ../src/qhull/qset.c ../src/qhull/mem.c \
					../src/qhull/stat.c \
					 ../src/psorting.c ../src/pscoring.c  \
                    ../src/utils.c ../src/pertable.c ../src/memhandler.c \
                    ../src/voronoi.c ../src/sort.c ../src/calc.c \
                    ../src/writepdb.c ../src/rpdb.c ../src/tparams.c \
                    ../src/fparams.c ../src/pocket.c ../src/refine.c \
                    ../src/descriptors.c ../src/cluster.c ../src/aa.c \
                    ../src/fpocket.c ../src/write_visu.c  ../src/fpout.c \
                    ../src/atom.c ../src/writepocket.c ../src/voronoi_lst.c 
 
tpocket_SOURCES =   ../src/qhull/qvoronoi.c ../src/qhull/user.c ../src/qhull/global.c \
					../src/qhull/io.c ../src/qhull/geom2.c ../src/qhull/poly2.c \
					../src/qhull/merge.c ../src/qhull/qhull.c ../src/qhull/geom.c  \
					../src/qhull/poly.c ../src/qhull/qset.c ../src/qhull/mem.c \
					../src/qhull/stat.c \
					../src/tpmain.c ../src/psorting.c ../src/pscoring.c ../src/utils.c \
                    ../src/pertable.c ../src/memhandler.c ../src/voronoi.c ../src/sort.c \
                    ../src/calc.c ../src/writepdb.c ../src/rpdb.c ../src/tparams.c \
                    ../src/fparams.c ../src/pocket.c ../src/refine.c ../src/tpocket.c  \
                    ../src/descriptors.c ../src/cluster.c ../src/aa.c ../src/fpocket.c \
                    ../src/write_visu.c ../src/fpout.c ../src/atom.c ../src/writepocket.c \
                    ../src/voronoi_lst.c ../src/neighbor.c
 
dpocket_SOURCES =  	../src/qhull/qvoronoi.c ../src/qhull/user.c ../src/qhull/global.c \
					../src/qhull/io.c ../src/qhull/geom2.c ../src/qhull/poly2.c \
					../src/qhull/merge.c ../src/qhull/qhull.c ../src/qhull/geom.c  \
					../src/qhull/poly.c ../src/qhull/qset.c ../src/qhull/mem.c \
					../src/qhull/stat.c \
					../src/dpmain.c ../src/psorting.c ../src/pscoring.c ../src/dpocket.c \
                    ../src/dparams.c  ../src/voronoi.c ../src/sort.c  ../src/rpdb.c \
                    ../src/descriptors.c ../src/neighbor.c ../src/atom.c ../src/aa.c \
                    ../src/pertable.c ../src/calc.c ../src/utils.c ../src/writepdb.c \
                    ../src/memhandler.c ../src/pocket.c ../src/refine.c ../src/cluster.c \
                    ../src/fparams.c ../src/voronoi_lst.c ../src/fpocket.c
					


