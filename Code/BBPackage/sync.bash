#!/bin/bash
# Copy the necessary files to the BayesLogit Package directory.

rsyncit="rsync -Crvzut --exclude-from=$HOME/.rsync-exclude $@"

BASE=$HOME/RPackage/BayesBridge/Code
MYLIB=$HOME/Code

BBDIR=$HOME/RPackage/BayesBridge/Code/BBPackage/BayesBridge

# CPP files.
$rsyncit $BASE/C/BridgeRegression.h    $BBDIR/src/
$rsyncit $BASE/C/BridgeWrapper.h       $BBDIR/src/
$rsyncit $BASE/C/BridgeWrapper.cpp     $BBDIR/src/
$rsyncit $BASE/C/retstable.h           $BBDIR/src/
$rsyncit $BASE/C/retstable.c           $BBDIR/src/
# $rsyncit $BASE/C/HmcSampler.cpp        $BBDIR/src/
# $rsyncit $BASE/C/HmcSampler.h          $BBDIR/src/

$rsyncit $BASE/C/magnet/ $BBDIR/inst/include/magnet/

$rsyncit $MYLIB/Matrix/Matrix.h            $BBDIR/src/
$rsyncit $MYLIB/Matrix/MatrixFrame.h       $BBDIR/src/
$rsyncit $MYLIB/Matrix/Matrix.cpp          $BBDIR/src/
$rsyncit $MYLIB/Matrix/MatrixFrame.cpp     $BBDIR/src/
$rsyncit $MYLIB/RNG/RNG.hpp                $BBDIR/src/RNG.h
$rsyncit $MYLIB/RNG/RRNG.hpp               $BBDIR/src/RRNG.h
$rsyncit $MYLIB/RNG/RNG.cpp                $BBDIR/src/
$rsyncit $MYLIB/RNG/RRNG.cpp               $BBDIR/src/
# $rsyncit $MYLIB/RandomVariates/Normal.hpp  $BBDIR/src/Normal.h



# R files.
$rsyncit $BASE/C/BridgeWrapper.R  $BBDIR/R/
$rsyncit $BASE/R/BridgeTMix.R     $BBDIR/R/
$rsyncit $BASE/R/BridgeNMix.R     $BBDIR/R/
$rsyncit $BASE/R/bridge-trace.R   $BBDIR/R/

# Data files.
$rsyncit $BASE/C/diabetes.RData  $BBDIR/data/

sed -i~ s/\"Bridge\"/\"BayesBridge\"/ $BBDIR/R/BridgeWrapper.R

# sed -i~ s/retstable\.c/retstable\.h/ $BBDIR/src/BridgeRegression.h
# sed -i~ s/quartic\.h/quartic\.hpp/ $BBDIR/src/HmcSampler.cpp
sed -i~ s/\.hpp/\.h/ $BBDIR/src/*.h
sed -i~ s/\.hpp/\.h/ $BBDIR/src/*.cpp

# Change to Rprintf
sed -i~ s/fprintf\(stderr,/printf\(/g $BBDIR/src/*.[ch]
sed -i~ s/fprintf\(stderr,/printf\(/g $BBDIR/src/*.cpp

# There must be a better way.
sed -i~ -e 's/Rprintf/printf/g' $BBDIR/src/*.[ch]
sed -i~ -e 's/Rprintf/printf/g' $BBDIR/src/*.cpp

sed -i~ -e 's/printf/Rprintf/g' $BBDIR/src/*.[ch]
sed -i~ -e 's/printf/Rprintf/g' $BBDIR/src/*.cpp

sed -i~ -e 's/source/\#source/' $BBDIR/R/BridgeNMix.R