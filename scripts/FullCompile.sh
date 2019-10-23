

#El siguiente es un script diseñado para compilar de manera eficiente mi proyecto de tesis doctoral
# autor: Jose Hugo Garcia Aguila phD student at Universidade Federal do Rio de Janeiro
#**************Se definen los paths principales de mi proyecto**************#
COMPILEROOT=$(pwd)
INCLUDEPATH=$COMPILEROOT/include
SOURCEPATH=$COMPILEROOT/src
LIBPATH=$COMPILEROOT/lib
BINPATH=$COMPILEROOT/bin
#**************Se definen las opciones los compiladores*********************#
MAINFLAG="-I"$COMPILEROOT"/include"  #-lgraphene_sparse_cuda -lkpm_sparse_cuda"  
#************Se compilan las librerias asociadas al proyecyo****************#
#************Se compilan los programas principales**************************
PCID=`uname -n`
NAME=$1
nvcc -o $BINPATH/$NAME$PCID $SOURCEPATH/$NAME.cu $MAINFLAG -lcudart -lcublas -lcusparse -arch=sm_20    #-ccbin=icc
