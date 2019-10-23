

#El siguiente es un script diseñado para compilar de manera eficiente mi proyecto de tesis doctoral
# autor: Jose Hugo Garcia Aguila phD student at Universidade Federal do Rio de Janeiro
#**************Se definen los paths principales de mi proyecto**************#
COMPILEROOT=~/Dropbox/GraphenePerioBC
INCLUDEPATH=$COMPILEROOT/include
SOURCEPATH=$COMPILEROOT/src
LIBPATH=$COMPILEROOT/lib
BINPATH=$COMPILEROOT/bin
#**************Se definen las opciones los compiladores*********************#
MAINFLAG="-I"$COMPILEROOT"/include"  #-lgraphene_sparse_cuda -lkpm_sparse_cuda"  
#************Se compilan las librerias asociadas al proyecyo****************#
#************Se compilan los programas principales**************************
PCID=`uname -n`
#NAME=GrapheneConductivityVacancies
#NAME=PureGrapheneConductivity
#NAME=PureGrapheneConductivityStreda
#NAME=HaldaneHoneyCombStreda
#NAME=HaldaneHoneyComb
#NAME=HaldaneHoneyComb
NAME=HaldaneVacancy
#NAME=ExtremaTestMagneticField
#NAME=FilterOperator
nvcc -o $BINPATH/$NAME$PCID $SOURCEPATH/$NAME.cu $MAINFLAG -lcudart -lcublas -lcusparse -arch=sm_20  
#icc -o $BINPATH/$NAME$PCID $SOURCEPATH/$NAME.cu $MAINFLAG -lcudart -lcublas -lcusparse -arch=sm_20  -I${MKLROOT}/include  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm  --compiler-bindir=icc
