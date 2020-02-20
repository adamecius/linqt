BUILD_DIR="../../../build"
printf "This example assumes the package is built on ${BUILD_DIR} directory \n\n"


SYSTEM_LABEL="graphene_example"

printf "(1) Run the python script ${SYSTEM_LABEL}.ipynb for creating the wannier files \n\n"
ipython ${SYSTEM_LABEL}.ipynb


DIMX=100
DIMY=100
DIMZ=1
printf "(2) Use wannier2sparse for expanding wannier system with label $WANNIER_LABEL to sparse matrices in CSR\n"
printf "into a supercell of the system  with dimensions $DIMX $DIMY $DIMZ in sparse CSR format\n\n"
$BUILD_DIR/wannier2sparse ${SYSTEM_LABEL} $DIMX $DIMY $DIMZ

printf "Moving these CSR matrices to the directory operators \n"
mkdir -p operators
mv *.CSR operators/

NUM_MOMS=100
SCAL_FACTOR=0.1;
SHIFT_FACTOR=0.0
printf "(3) Use inline_compute-kpm-nonEqOp for computing $NUM_MOMS chebyshevs moments \n\n"
$BUILD_DIR/inline_compute-kpm-nonEqOp $SYSTEM_LABEL VX   VX $NUM_MOMS $SCAL_FACTOR $SHIFT_FACTOR
$BUILD_DIR/inline_compute-kpm-nonEqOp $SYSTEM_LABEL VYSZ VX $NUM_MOMS $SCAL_FACTOR $SHIFT_FACTOR


#printf "(4) Use inline_kuboBastinFromChebmom for computing XX AND XYSZ conductivities \n\n"
$BUILD_DIR/inline_kuboBastinFromChebmom NonEqOpVXgraphene_exampleKPM_M100x100RV1.chebmom2D  $NUM_MOMS
$BUILD_DIR/inline_kuboBastinFromChebmom NonEqOpVYSZgraphene_exampleKPM_M100x100RV1.chebmom2D  $NUM_MOMS


printf "FINISHED EXAMPLE \n"

