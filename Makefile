INCLUDE_PATH= -I -Wno-deprecated-declarations -fopenmp  -I.  -I/usr/local/openfpm_numerics/include -I/usr/local/openfpm_pdata/include/config -I/usr/local/openfpm_pdata/include -I/usr/local/openfpm_data/include -I/usr/local/openfpm_vcluster/include -I/usr/local/openfpm_io/include -I/usr/local/openfpm_devices/include -I/root/VCDEVEL/include  -I/root/METIS/include -I/root/PARMETIS/include -I/root/BOOST/include -I/root/HDF5/include -I/root/LIBHILBERT/include   -I/root/PETSC/include -I/root/OPENBLAS/include -I/root/SUITESPARSE/include -I/root/EIGEN -I/root/BLITZ/include -I/root/ALGOIM/include
LIBS_PATH= -L/usr/local/openfpm_devices/lib -L/usr/local/openfpm_pdata/lib  -L/usr/local/openfpm_vcluster/lib -L/root/VCDEVEL/lib  -L/root/METIS/lib -L/root/PARMETIS/lib  -L/root/BOOST/lib -L/root/HDF5/lib -L/root/LIBHILBERT/lib   -L/root/PETSC/lib -L/root/OPENBLAS/lib -L/root/SUITESPARSE/lib 
LIBS=-fopenmp  /usr/lib/gcc/x86_64-linux-gnu/9/libgomp.so -lvcluster -lofpm_pdata -lofpmmemory -lparmetis -lmetis -lboost_iostreams -lboost_program_options -lhdf5 -llibhilbert -lVc     -lrt -lpetsc -lopenblas -lumfpack -lamd -lbtf -lcamd -lccolamd -lcholmod -lcolamd -lcxsparse -lklu -ldl -lrbio -lspqr -lsuitesparseconfig -ldl -lboost_filesystem -lboost_system -L/root/BOOST/lib -lboost_context
LIBS_NVCC=-Xcompiler=-fopenmp  -lvcluster -lofpm_pdata -lofpmmemory -lparmetis -lmetis -lboost_iostreams -lboost_program_options -lhdf5 -llibhilbert -lVc     -lrt -lpetsc -lopenblas -lumfpack -lamd -lbtf -lcamd -lccolamd -lcholmod -lcolamd -lcxsparse -lklu -ldl -lrbio -lspqr -lsuitesparseconfig -ldl -lboost_filesystem -lboost_system -L/root/BOOST/lib -lboost_context
INCLUDE_PATH_NVCC=-Xcompiler=-Wno-deprecated-declarations -Xcompiler=-fopenmp   -I. -I/usr/local/openfpm_numerics/include -I/usr/local/openfpm_pdata/include/config -I/usr/local/openfpm_pdata/include -I/usr/local/openfpm_data/include -I/usr/local/openfpm_vcluster/include -I/usr/local/openfpm_io/include -I/usr/local/openfpm_devices/include -I/root/METIS/include -I/root/PARMETIS/include -I/root/BOOST/include -I/root/HDF5/include -I/root/LIBHILBERT/include   -I/root/PETSC/include -I/root/OPENBLAS/include -I/root/SUITESPARSE/include -I/root/EIGEN -I/root/BLITZ/include -I/root/ALGOIM/include

CC=mpic++

LDIR =
OPT=

OBJ = main.o

sph_dlb:
sph_dlb_test: OPT += -DTEST_RUN
sph_dlb_test: sph_dlb


%.o: %.cpp
	$(CC) -O3 -g  $(OPT) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

sph_dlb: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: sph_dlb

run: sph_dlb_test
	mpirun --allow-run-as-root --oversubscribe -np 4 ./sph_dlb

.PHONY: clean all run

clean:
	rm -f *logfile.txt
	rm -f *.o *~ core sph_dlb log
	cd  output/
	rm -rf output/*.vtp output/*.pvtp

