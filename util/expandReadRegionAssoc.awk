#!/usr/bin/awk -f

BEGIN {
    if(ARGC < 3) {
        print "Please provide an edge-candidates and an expanded fragment-edge-associations" > "/dev/stderr"
        exit ++ExitCode;
    }
    FS="\t"
    OFS="\t"
}

(ARGIND == 1) {
    HReg[$1]=$2;
    VReg[$1]=$3;
    next
}

{
    print $1,HReg[$2],192
    print $1,VReg[$2],192
}

END {
    if(ExitCode) {exit ExitCode; }
}
