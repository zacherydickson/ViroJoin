#!/bin/bash

RED="\033[0;31m"
GREEN="\033[0;32m"
YELLOW="\033[0;33m"
NC="\033[0m"

TestDirBase="testfiles"

IsolateKmerLen=18;
MinClipLen=20;


function main {
    if [ "$#" -lt 2 ]; then
        >&2 echo -e "Usage: $(basename "$0") ReadInfo.tsv workDir [ test1=all ... ]\n" \
                    "\tReadInfo is tab sep with headers: ReadID, Parity, Length, BPPos, HostLen," \
                    "\t VirusLen, MaxViralNucProp, Config\b" \
                    "\tworkDir is a ViroJoin output directory, some intermediate files may be created at workdir/$TestDirBase\n" \
                    "\tall tests are:\n" \
                    "\t\tisolation - All reads with at least $IsolateKmerLen viral bp are retained\n" \
                    "\t\tmapping - All relevant isolated reads map\n" \
                    "\t\tenumerate_chimeras - All reads at least $MinClipLen viral bp have a junction\n";
        exit 1;
    fi
    readInfoFile=$1; shift;
    workingDir=$1; shift;
    testDir="$workingDir/$TestDirBase"
    [ -d "$testDir" ] && 
        >&2 echo -e "[${YELLOW}WARNING${NC}] $testDir already exists: Test results may not be up to date";
    mkdir -p "$testDir"
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
    nTests=0;
    if [[ $bAll == 1 || -n "${testSet["isolation"]}" ]]; then
        ((nTests++));
        RunTest "isolation" "$readInfoFile" "$workingDir" || ((bFail++));
    fi
    if [[ $bAll == 1 || -n "${testSet["mapping"]}" ]]; then
        ((nTests++));
        RunTest "mapping" "$readInfoFile" "$workingDir" || ((bFail++));
    fi
    if [[ $bAll == 1 || -n "${testSet["enumerate_chimeras"]}" ]]; then
        ((nTests++));
        RunTest "enumerate_chimeras" "$readInfoFile" "$workingDir" || ((bFail++));
    fi
    if [ "$nTests" == 0 ]; then
        >&2 echo -e "[${YELLOW}WARNING${NC}] No implemented tests requested";
    fi
    #Report the overall result
    finalRes="${GREEN} All Requested Tests Passed${NC} ($nTests/$nTests)"
    [ "$bFail" -gt 0 ] && finalRes="${RED} Some Requested Tests Failed${NC} ($bFail/$nTests)";
    >&2 echo -e "$finalRes"
}


function NoteMappedReads {
    workDir=$1; shift
    outFile=$1; shift
    bam="$workDir/bam_0/all_alignments.ns.bam"
    [ -s "$bam" ] ||
        { echo "all_alignments.ns is missing or empty"; return 1; }
    samtools view "$bam" | 
        awk '
            ($1 ~ /_[LR]_[12]/){
                n=split($1,a,"_");
                str=a[1];
                for(i=2;i<=n-2;i++){
                    str=str"_"a[i]
                };
                $1=str
            }
            {id=$1}
            (!Seen[id]++){print id; }
        ' >| "$outFile"
}

function NoteRetainedReads {
    workDir=$1; shift
    outFile=$1; shift
    fwd="$workDir/bam_0/retained-pairs_1.fq";
    rev="$workDir/bam_0/retained-pairs_2.fq";
    [ -s "$fwd" ] && [ -s "$rev" ] || 
        { echo "retained-pairs_[12].fq missing or empty"; return 1; }
    awk -v rev="$rev" '
        (FNR % 4 != 1){getline < rev; next}
        {
            id1=substr($1,2);
            getline < rev    
            id2=substr($1,2);
            print id1"\t"id2;
        } 
    ' "$fwd" >| "$outFile"
}

function RunTest {
    key=$1; shift;
    infoFile=$1; shift;
    workDir=$1; shift
    testDir="$workDir/$TestDirBase"
    mkdir -p "$testDir"
    retVal=1;
    resStr="${RED}FAIL";
    if msg="$("test_${key}" "$infoFile" "$workDir" "$testDir")"; then
        resStr="${GREEN}PASS"   
        retVal=0;
    else
        msg=" : $msg"
    fi
    >&2 echo -e "[$resStr${NC}] ${key}${msg}"
    return "$retVal"
}

function test_enumerate_chimeras {
    infoFile=$1; shift
    workDir=$1; shift
    testDir=$1; shift
    mapReadsFile="$testDir/mapped.list"
    resFile="$workDir/junction-candidates.bedpe"
    log="$testDir/enumerate_chimeras.log"
    [ -s "$mapReadsFile" ] ||
        NoteMappedReads "$workDir" "$mapReadsFile" ||
        return 1;
    [ -s "$resFile" ] ||
        { echo "$resFile is missing or empty"; return 1; }
    awk -v minClip="$MinClipLen" -v lf="$log" '
        function failure(msg) {
            print msg ", see", lf
            print msg > lf
            print FNR ": " $0 > lf
            exit 1;
        }
        function warning(msg) {
            if(!bWarned){
                print " Warning - " msg ", see", lf
                bWarned=1;
            }
            print msg > lf
            print FNR ": " $0 > lf
        }
        (ARGIND == 1) {
            config = "H" $9 "V" $10;
            ObsConfig[$7,config] = 1
            next;
        }
        (ARGIND == 2){ InMapSet[$1]=1; next; }
        # process read info
        (FNR == 1){ #header line
            for(i=1;i<=NF;i++){ colIdx[$i]=i; }
            next
        }
        #Skip reads which did not map
        { id=$colIdx["ReadID"]; }
        (!InMapSet[id]) { next; }
        { config=$colIdx["Config"]; }
        #Combine Read Pair Info
        {
            hL = $colIdx["HostLen"];
            vL = $colIdx["VirusLen"];
            getline;
            hL = ($colIdx["HostLen"] > hL) ? $colIdx["HostLen"] : hL;
            vL = ($colIdx["VirusLen"] > vL) ? $colIdx["VirusLen"] : vL;
        }
        # Skip reads not expected to support a chimera
        ((vL < minClip || hL < minClip)) { 
            # warn if it was observed anyway
            if(ObsConfig[id,config]) { warning("Extra Chimera"); }
            next;   
        }
        (!ObsConfig[id,config]) { failure("Missing Chimera"); }

    ' "$resFile" "$mapReadsFile" "$infoFile"
}

#Ensure that the mapped reads are properly paired,
#   there are no duplicate read ids, and
#   that the reads with sufficient viral sequence are retained
function test_isolation {
    infoFile=$1; shift
    workDir=$1; shift
    testDir=$1; shift
    retReadsFile="$testDir/retained.pairs.tab"
    log="$testDir/isolate.log"
    [ -s "$retReadsFile" ] ||
        NoteRetainedReads "$workDir" "$retReadsFile" ||
        return 1;
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

function test_mapping {
    infoFile=$1; shift
    workDir=$1; shift
    testDir=$1; shift
    retReadsFile="$testDir/retained.pairs.tab"
    mapReadsFile="$testDir/mapped.list"
    log="$testDir/mapping.log"
    [ -s "$retReadsFile" ] ||
        NoteRetainedReads "$workDir" "$retReadsFile" ||
        return 1;
    [ -s "$mapReadsFile" ] ||
        NoteMappedReads "$workDir" "$mapReadsFile" ||
        return 1;
    awk -v lf="$log" '
        function failure(msg) {
            print msg ", see", lf
            print msg > lf
            print FNR ": " $0 > lf
            exit 1;
        }
        (ARGIND == 1){ InMapSet[$1]=1; next; }
        (ARGIND == 2){ InIsoSet[$1]=1; next; }
        # process read info
        (FNR == 1){ #header line
            for(i=1;i<=NF;i++){ colIdx[$i]=i; }
            next
        }
        { id = $colIdx["ReadID"]; }
        (InIsoSet[id] && !InMapSet[id]) { failure("Missing Read"); }
    ' "$mapReadsFile" "$retReadsFile" "$infoFile"
    return 0;
}

if [ "${BASH_SOURCE[0]}" == "${0}" ]; then
    main "$@"
fi

