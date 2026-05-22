#!/bin/bash

RED="\033[0;31m"
GREEN="\033[0;32m"
NC="\033[0m"

TestDirBase="testfiles"

IsolateKmerLen=18;


function main {
    if [ "$#" -lt 2 ]; then
        >&2 echo -e "Usage: $(basename "$0") ReadInfo.tsv workDir [ test1=all ... ]\n" \
                    "\tReadInfo is tab sep with headers: ReadID, Parity, Length, BPPos, HostLen, and VirusLen\b" \
                    "\tworkDir is a ViroJoin output directory, some intermediate files may be created at workdir/$TestDirBase\n" \
                    "\tall tests are:\n" \
                    "\t\tisolation - All reads with at least $IsolateKmerLen viral bp are retained\n";
        exit 1;
    fi
    readInfoFile=$1; shift;
    workingDir=$1; shift;
    #Determine the requested tests
    bAll=0;
    [ "$#" -lt 1 ] && bAll=1
    declare -A testSet;
    for testNm in "$@"; do
        [ "$testNm" == "all" ] && bAll=1;
        testSet["$testNm"]=1;
    done
    #Perform requested Tests 
    bFail=0;
    if [[ $bAll == 1 || -n "${testSet["isolation"]}" ]]; then
        RunTest "isolation" "$readInfoFile" "$workingDir" || ((bFail++));
    fi
    #Report the overall result
    finalRes="${GREEN} All Requested Tests Passed${NC}"
    [ "$bFail" -gt 0 ] && finalRes="${RED} Some Requested Tests Failed${NC}";
    >&2 echo -e "$finalRes"
}

#Requires testDir to exist
function NoteRetainedReads {
    workDir=$1; shift
    outFile=$1; shift
    fwd="$workDir/bam_0/retained-pairs_1.fq";
    rev="$workDir/bam_0/retained-pairs_2.fq";
    awk -v rev="$rev" '
        (FNR % 4 != 1){getline < rev; next}
        {
            id1=substr($1,2);
            getline < rev    
            id2=substr($1,2);
            print id1"\t"id2;
        } 
    ' "$fwd" >| "$outFile"
    echo "$outFile"
}

function RunTest {
    key=$1; shift;
    infoFile=$1; shift;
    workDir=$1; shift
    testDir="$workDir/$TestDirBase"
    mkdir -p "$testDir"
    retVal=1;
    resStr="${RED}FAIL";
    if msg=" : $("test_${key}" "$infoFile" "$workDir" "$testDir")"; then
        resStr="${GREEN}PASS"   
        retVal=0;
    fi
    >&2 echo -e "[$resStr${NC}] ${key}${msg}"
    return "$retVal"
}


function test_isolation {
    infoFile=$1; shift
    workDir=$1; shift
    testDir=$1; shift
    retReadsFile="$testDir/retained.pairs.tab"
    log="$testDir/isolate.log"
    NoteRetainedReads "$workDir" "$retReadsFile"
    awk -v minVL="$IsolateKmerLen" -v lf="$log" '
        function failure(msg) {
            print msg ", see", lf
            print msg > lf
            print FNR ": " $0 > lf
            exit 1;
        }
        (ARGIND == 1){
            if($1 != $2) { failure("Pairing Failure"); } 
            if(++InSet[$1] > 1){ failure("Duplicate Read ID"); }
            next
        }
        # process read info
        (FNR == 1){ #header line
            for(i=1;i<=NF;i++){ colIdx[$i]=i; }
            next
        }
        #Skip low virus reads
        ($colIdx["VirusLen"] < minVL){ next; }
        (!InSet[$colIdx["ReadID"]]) { failure("Missing Read"); }
    ' "$retReadsFile" "$infoFile" 
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then
    main "$@"
fi

#samtools view bam_0/all_alignments.ns.bam | awk '($1 ~ /_[LR]_[12]/){n=split($1,a,"_"); str=a[1]; for(i=2;i<=n-2;i++){str=str"_"a[i]}; $1=str} {id=$1} (!Seen[id]++){print id;}' >| mappedIDs.list
