INCLUDE="-Iinclude -I../../../include/"
LDLIB="-L../../../lib/"
LIB="../../../lib/parser_util.o ../../../lib/mat_op.o ../../../lib/add_k_dep.o ../../../lib/write_operator.o" 
#g++ -O3 -Iinclude -I../../kpm-non-equilibrium/shared/include/ src/graphene+disorder.cpp -o graphene+disorder
#g++ -O3 -Iinclude -I../../kpm-non-equilibrium/shared/include/ src/graphene+soc+disorder.cpp -o  graphene+soc+disorder
#g++ -O3 -Iinclude -I../../kpm-non-equilibrium/shared/include/ src/graphene+disorder-projGen.cpp -o graphene+disorder-projGen 
#g++ -O3 -Iinclude -I../../kpm-non-equilibrium/shared/include/ src/graphene+soc+kspace.cpp -o graphene+soc+kspace

icpc -O3 $INCLUDE  src/gr+tmd.cpp  -o gr+tmd  $LIB

