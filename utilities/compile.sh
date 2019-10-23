nvcc -DNEIGEN -c src/full_spectral_sum.cu -Iinclude -I../include -o obj/full_spectral_sum.o
nvcc -DNEIGEN -c -x cu src/spectral_conductivity.cpp -Iinclude -I../include -o obj/spectral_conductivity.o
nvcc -DNEIGEN obj/full_spectral_sum.o  obj/spectral_conductivity.o -o main
