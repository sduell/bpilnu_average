#!/bin/sh

code=plotting_code/LikelihoodScanPlots.C
exec=plotting_code/plot_likelihood_scans.exe
rm -f $exec

CXX=$($ROOTSYS/bin/root-config --cxx)
flags="$($ROOTSYS/bin/root-config --cflags --glibs) -lTreePlayer -Iatlasstyle -lRooFitCore -lRooFit -lRooStats"

echo "Compilation flags"
echo $flags
echo "------------"
success=no
$CXX $flags $code -o $exec && {
    echo ; echo "  Successful compilation.";
    echo "  Produced:"; echo "    $exec" ; echo
    success=yes
} || {
    echo ; echo "Compilation failed :("
    echo "$PWD"
}
echo "------------"

vars="Inclusive pTH yAbsH Njets pTj1"
[[ $success = yes ]] && {
    echo ; echo "  Running:"
    ./$exec $vars

    echo ; echo "------------"
}
