#!/usr/bin/awk -f

function output(    fragGrpStr, fragIdxList, nFrag, grpIdx, a, line){
    nFrag=asorti(InGroup,fragIdxList);
    #Build a key for the particular set of fragment Indexes
    fragGrpStr=fragIdxList[1]; 
    for(i=2;i<=nFrag;i++){
        fragGrpStr=fragGrpStr","fragIdxList[i]
    }
    #If the key has not been seen, assign a fragment group ID
    if(!GroupIdx[fragGrpStr]) {
        GroupIdx[fragGrpStr]=NextGroupID++
        str=GroupIdx[fragGrpStr] "\t" nFrag
        for(i=1;i<=nFrag;i++){
            str = str "\t" fragIdxList[i]
        }
        print str > fragmentGroupFile
    }
    grpIdx=GroupIdx[fragGrpStr]
    split(lastline,a,"\t")
    a[7] = grpIdx
    line=a[1];
    for(i=2;i<=12;i++) {
        line = line "\t" a[i]
    }
    print line
}

BEGIN {
#    if(ARGC < 2) {
#        print "Usage: compressCandidates.awk sorted-junction-candidates.bedpe > compressed\n\tAlso creates fragment-index.tsv and fragment-group.tsv" > "/dev/stderr";
#        exit (ExitCode = 1);
#    }
    OFS="\t"
    if(!outDir) {outDir = "."}
    NextFragmentID=0;
    NextGroupID=0;
    fragmentIndexFile=outDir"/fragment-index.tsv"
    print "Name\tIndex" > fragmentIndexFile
    fragmentGroupFile=outDir"/fragment-group.tsv"
    print "GroupIndex\tNFragments\tFragmentIndex1\t..." >  fragmentGroupFile
}

{ # Obtain the fragment index for the bedpe entry, and then remove the name from the entry
    fragName = $7;
    $7="dummy"
    #If this fragName does not have an id, assign one and record this assignment
    if(!FragIdx[fragName]){
        FragIdx[fragName]=NextFragmentID++;
        print fragName"\t"FragIdx[fragName] > fragmentIndexFile
    }
    fragID=FragIdx[fragName]; 
}

($0 != lastline){ #Check if we have moved to a new block of entries
    if(lastline){output()};
    delete InGroup;
    lastline=$0
}

{ #Record this fragment as being associated with the current entry
    InGroup[fragID]=1;
}

