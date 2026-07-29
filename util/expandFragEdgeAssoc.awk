#!/usr/bin/awk -f

BEGIN {
    if(ARGC < 2){
        print "Provide a fragment-edge-association.tab file" > "/dev/stderr"
        exit ++ExitCode;
    }
    FS ="\t"
    OFS="\t"
    print "Loading" > "/dev/stderr";
}

{ #Load the associations
    EdgeForFrag[$1,++NEdgeByFrag[$1]]=$2;
    FragInEdge[$2,++NFragInEdge[$2]]=$1;
    AdjMat[$2,$2]=1; #Make sure all edges are adjacent to themselves
    printf("%d\r",FNR) > "/dev/stderr";
}


END {
    if(ExitCode){exit ExitCode;}
    #Build an adjacency matrix for edges (edges are nodes in a graph, connected if they share a read)
    n=0;
    print "Building AdjMat" > "/dev/stderr"
    for(frag in NEdgeByFrag){
        for(i=1;i<NEdgeByFrag[frag];i++){
            edge_i=EdgeForFrag[frag,i];
            for(j=1;j<=NEdgeByFrag[frag];j++){
                edge_j =EdgeForFrag[frag,j]
                AdjMat[edge_i,edge_j]=1
                AdjMat[edge_j,edge_i]=1
            }
        }
        printf("%d\r", ++n) > "/dev/stderr";
    }
    print "Building Edges" > "/dev/stderr"
    #For each edge, associate all unique fragments from edges adjacent to this edge (should include itself)
    n=0;
    for(edge_i in NFragInEdge){
        #Get all unique fragments from adjacent edges
        delete fragSet;
        for(edge_j in NFragInEdge){
            if(!AdjMat[edge_i,edge_j]) { continue; }
            for(i=1;i<=NFragInEdge[edge_j];i++){
                frag = FragInEdge[edge_j,i];
                fragSet[frag] = 1;
            }
        }
        #Output the result
        for(frag in fragSet){
            print frag, edge_i, 0
        }
        printf("%d\r", ++n) > "/dev/stderr";
    }
    print "Done" > "/dev/stderr"
}

