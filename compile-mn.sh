module purge
module load impi intel mkl gsl fftw

mpiicc -O3 -funroll-loops -march=native -mfpmath=sse -O3 -I../../shared/include/ $FFTW_INCL src/compute-nonEqOp-KPM+FFT.cpp  -lgsl -lgslcblas $FFTW_LDFLAGS -lm -o compute-nonEqOp-KPM+FFT
mpiicc -O3 -funroll-loops -march=native -mfpmath=sse -O3 -I../../shared/include/ $FFTW_INCL src/compute-DOS+FFT.cpp -lgsl -lgslcblas $FFTW_LDFLAGS -lm -o compute-DOS+FFT

