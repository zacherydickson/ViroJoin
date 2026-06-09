#include <algorithm>
#include <iostream>
#include <functional>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "BreakpointGraph.h"
#include "ChimericFragment.h"
#include <cptl_stl.h>
#include "RegionGraph.h"
#include "igraph/igraph.h"
#include "config.h"
#include "utils.h"

//==== TYPE DECLARATIONS


typedef std::unordered_map<std::string,CBPGraph> BPGraphMap_t;
//Type to contain a label for a host-virus contig (chr and strand) pair
typedef std::pair<std::string,std::string> StrandLabelPair_t;

struct StrandLabelPair_HashFunctor {
    std::size_t operator() (const StrandLabelPair_t & pair) const {
        return std::hash<std::string>{}(pair.first + pair.second);
    }
};

struct StrandLabelPair_EqualFunctor {
    bool operator() (const StrandLabelPair_t & a, const StrandLabelPair_t & b) const {
        return (a.first + a.second) == (b.first + b.second);
    }
};

//Type to contain a vector of Chimeric Fragments
typedef std::vector<ChimericFragment_t> ChimericFragmentVec_t;
//Type to associate a vector of chimeric fragments with a host-virus contig pair
typedef std::unordered_map< StrandLabelPair_t, ChimericFragmentVec_t,
                            StrandLabelPair_HashFunctor,
                            StrandLabelPair_EqualFunctor> ChimericFragmentVecMap_t;

//Type to contain a vector of CBPGraphs (post decomposition)
typedef std::vector<CBPGraph> BPGraphVec_t;
//Type to contain a host-virus pair of Graph Vectors
typedef std::pair<BPGraphVec_t,BPGraphVec_t> BPGraphVecPair_t;
//Type to map host-virus contig labels to corresponding graph vector pairs
typedef std::unordered_map<StrandLabelPair_t,BPGraphVecPair_t,
                            StrandLabelPair_HashFunctor,
                            StrandLabelPair_EqualFunctor> BPGraphVecPairMap_t;

//typedef std::unordered_map<StrandLabelPair_t,CRegionGraph,
//                            StrandLabelPair_HashFunctor,
//                            StrandLabelPair_EqualFunctor> PairedRegGraphMap_t;
typedef std::vector<CRegionGraph> RegGraphVec_t;

//==== GLOBAL VARIABLE DECLARATIONS

static size_t MinimumReads = 4;
static int SplitBonus = 1;
int MaxInsertSize;
int ReadLength;
int UpstreamSize = 5;
double SplitFactor = 2.0;
std::unordered_set<std::string> VirusNameSet;

//==== FUNCTION DECLARATIONS

bool ClusterBPGraph(CBPGraph & graph);
bool ClusterBPGraphs(BPGraphMap_t & graphMap);
bool ClusterBPGraphVecPairMap(BPGraphVecPairMap_t & graphVecPairMap);
BPGraphVecPair_t ConstructBPGraphVecPair(const ChimericFragmentVec_t & fragVec);
BPGraphVecPairMap_t ConstructBPGraphVecPairMap(
        const ChimericFragmentVecMap_t & fvMap );
//void ConnectBPGraphs(BPGraphMap_t & graphMap);
CRegionGraph ConstructRegionGraph(const BPGraphVecPair_t & graphMap);
RegGraphVec_t ConstructRegionGraphVec(
        const BPGraphVecPairMap_t & graphVecPairMap);
CRegionGraph ConstructAndFilterRegionGraph(
        const BPGraphVecPair_t & graphVecPair);
//PairedRegGraphMap_t ConstructPairedRegionGraphMap(
//        const BPGraphVecPairMap_t & graphVecPairMap);
BPGraphMap_t DecomposeBPGraph(  const std::string & baseLabel,
                                CBPGraph & graph);
BPGraphVec_t DecomposeBPGraph(CBPGraph & graph);
bool DecomposeBPGraphs(BPGraphMap_t & graphMap);
bool DecomposeBPGraphs(BPGraphMap_t & graphMap, ctpl::thread_pool & threadPool);
bool DecomposeBPGraphVec(BPGraphVec_t & graphVec);
bool FilterBPGraph(CBPGraph & graph);
bool FilterBPGraphs(BPGraphMap_t & graphMap);
bool FilterBPGraphVec(BPGraphVec_t & graphVec);
bool FilterBPGraphVecPairMap(BPGraphVecPairMap_t & graphVecPairMap);
bool FilterAndClusterBPGraph(CBPGraph & graph);
RegGraphVec_t FragmentMapToRegionGraphVec(const ChimericFragmentVecMap_t & fvMap);
RegGraphVec_t FragmentMapToRegionGraphVec(size_t nThread, 
                                    const ChimericFragmentVecMap_t & fvMap);
CRegionGraph FragmentsToRegionGraph(const ChimericFragmentVec_t & fragVec);
bool OperateOnBPGraphVecPairMap(BPGraphVecPairMap_t & graphVecPairMap,
                                const std::string & operationName,
                                std::function<bool(BPGraphVec_t &)> operation);
bool ProcessBPGraphs(BPGraphMap_t & graphMap);
bool ProcessBPGraphs(size_t nThread, BPGraphMap_t & graphMap);
BPGraphVecPairMap_t ProcessFragments(const ChimericFragmentVecMap_t & fvMap);
BPGraphVecPairMap_t ProcessFragments(size_t nThread, const ChimericFragmentVecMap_t & fvMap);
#ifndef NDEBUG
void OutputDebugBPGraph(const BPGraphMap_t & graphMap,
                        std::string BPAdjFileName,
                        std::string BPVertFileName);
void OutputDebugRegGraph(   const CRegionGraph & regGraph,
                            std::string RegAdjFileName,
                            std::string RegVertFileName);
#endif //NDEBUG
void OutputResults( const CRegionGraph & regGraph,
                    const std::string & regFileName,
                    const std::string & edgeFileName,
                    const std::string & assocFileName);

BPGraphMap_t LoadBPGraphs(const std::string & fname,bool bOnline = true);
ChimericFragmentVecMap_t LoadFragments(const std::string & fname);
std::string to_bed(CRegionGraph::VertexProps);

//==== MAIN

//Parses candidate junctions into a graph of breakpoint positions in each contig
// Breakpoints defined by fragments with the same proximal and distal positions
//  are counted as one breakpoint position (all such fragments for a fragment group)
// Those are then clustered into regions by identifying maximal cliques
// Cliques are filtered for minimum size
// The regions are then used to construct a bipartite graph with edges between host and
//  viral regions which have fragments associated with both
// The edges are filtered for minimum support
// Regions with edges are output to a bed file, with the regionID as the name field
// Edges are output to a tab delim file:
//  fragmentName, regionIDHost, regionIDVirus, edgeID, grpIDList
//  each regIDHost- regionIDVirus combo is associted with one edgeID
//  a given fragment may occur in multiple fragment groups associated with an edge
//  the amount of support for an edge is then the number of unique grpIDs associated with that edgeID
//Inputs - the virus ref name (to differentiate host and viral references)
//       - the workspace (to find the stats file)
//       - the working dir (for output and to find the config and junction candidates)
//Output - A bed file defining the regions associated with any edges
//       - A tab delim file specifying the edges:
//          edgeID, hostRegID, virusRegID, edgeMetadata(NFragGrp)
//       - A tab delim file associating fragments with edges
//          fragmentName, edgeID , fragmentGrpInEdgeID
int main(int argc, char* argv[]) {
    //#Parse Inputs
    std::string virus_ref_fname = argv[1];
    std::string workdir = argv[2];
    std::string workspace = argv[3];

    //## Files to be used from the workdir
    std::string stats_file_name = workspace + "/stats.txt";
    std::string candidate_file_name = workdir + "/junction-candidates.bedpe";
    std::string config_file_name = workdir + "/config.txt";
    //## Output Files
    std::string regFileName = workdir + "/region-candidates.bed";
    std::string edgeFileName = workdir + "/edge-candidates.tab";
    std::string assocFileName = workdir + "/fragment-edge-associations.tab";
#ifndef NDEBUG
    //##Debug Files
    std::string BPAdjFileName = workdir + "/BPadj.tab";
    std::string BPVertFileName = workdir + "/BPvertex.tab";
    std::string RegAdjFileName = workdir + "/Regadj.tab";
    std::string RegVertFileName = workdir + "/Regvertex.tab";
#endif //NDEBUG

    LoadVirusNames(virus_ref_fname,VirusNameSet);

    MaxInsertSize = parse_stats(stats_file_name).max_is;
    ReadLength = parse_config(config_file_name).read_len;
    size_t nThread = parse_config(config_file_name).threads;

    igraph_setup();

    ChimericFragmentVecMap_t fragVecMap = LoadFragments(candidate_file_name);

    RegGraphVec_t regGraphVec = (nThread == 1) ? 
                                FragmentMapToRegionGraphVec(fragVecMap) :
                                FragmentMapToRegionGraphVec(nThread, fragVecMap);
    fprintf(stderr,"Merging Paired Region Graphs ...\n");
    CRegionGraph regGraph = CRegionGraph::merge_graphs(regGraphVec);
    regGraph.ensureConstructed(); //Explicit Call 
    regGraph.mergeUninformitiveOverlap();
    regGraph.filterEdges(MinimumReads,SplitBonus);
    fprintf(stderr,
            "After merging and filtering, The region graph contains %d edges between %d regions  ...\n",
            regGraph.ecount(),regGraph.vcount());
    
    OutputResults(regGraph,regFileName,edgeFileName,assocFileName);
    //TODO Fix the giant region bug
           
    fprintf(stderr,"Done - enumerate_edges\n");
}

//==== FUNCTION DEFINITIONS

//Identifies all maximal cliques within a graph
//Inputs    - a graph to filter
//Output    - true if the graph still has sufficient support, false otherwise
bool ClusterBPGraph(CBPGraph & graph) {
    return graph.maximalCliques(MinimumReads, SplitBonus);
}

//Filters each graph in a graph map
//Inputs    - a graph map containing graphs to filter
//Output    - None, modifies the input
bool ClusterBPGraphs(BPGraphMap_t & graphMap) {
    fprintf(stderr,"Clustering fragments within graphs ...\n");
    //int counter = 0;
    for(auto it = graphMap.begin(); it != graphMap.end(); ){
        fprintf(stderr,"Clustering graph with %d nodes and %d edges ...\n",it->second.vcount(),it->second.ecount());
        if(ClusterBPGraph(it->second)){ 
            it++;
        } else {
            it = graphMap.erase(it);
        }
    }
    fprintf(stderr,"After clustering %lu graphs remain\n",graphMap.size());
    return (graphMap.size() > 0);
}

bool ClusterBPGraphVec(BPGraphVec_t & graphVec) {
    auto removeIt = std::remove_if( graphVec.begin(),graphVec.end(),
                                    [](CBPGraph & graph) {
                                        return !ClusterBPGraph(graph);
                                    } );
    graphVec.erase(removeIt,graphVec.end());
    return !graphVec.empty();
}

//Handle construction of a pair of BP Graph Vectors (with one element each)
// from a vector of chimeric fragments
//The first Vector refers to host regions and
//  the second vector refers to viral regions
//Construction is only performed if there are sufficient fragments
//Input - A vector of chimeric fragments
//Output - A pair of CBPGraph vectors, if there were enough fragments
//          each vector will have one element
//        Otherwise both vectors will be empty
BPGraphVecPair_t ConstructBPGraphVecPair(const ChimericFragmentVec_t & fragVec) {
    BPGraphVecPair_t gvPair;
    //If there are not enough fragments for this contig pair, then there
    //will not be a valid breakpoint
    if(fragVec.empty() ||fragVec.size() < MinimumReads - SplitBonus) {
        return gvPair;
    }
    const ChimericFragment_t & front = fragVec.front();
    for( ChimericFragment_t::IV_IDX ivIdx :
            {ChimericFragment_t::IV1, ChimericFragment_t::IV2} ) 
    {
        //Construct the empty graph for this interval
        BPGraphVec_t * vec_ptr = (ivIdx == ChimericFragment_t::IV1) ?
                                &gvPair.first : &gvPair.second;
        vec_ptr->emplace_back(front.getChr(ivIdx), front.opens_left(ivIdx),
                                UpstreamSize, ReadLength,
                                MaxInsertSize, SplitFactor);
        CBPGraph & graph = vec_ptr->back();
        //Add each fragment to the graph
        for(const ChimericFragment_t & frag : fragVec){
           graph.addOrUpdateVertex( frag.proximal_pos(ivIdx),
                                    frag.distal_pos(ivIdx),
                                    frag.is_split(ivIdx),
                                    frag.getName(),false);
        }
        //Explicity request edges be constructed now
        graph.ensureConstructed();
    }
    return gvPair;
}

//Handle construction of a map associating a host-virus contig pair to a host-virus
// graph vector pair; each vector will have just one element
BPGraphVecPairMap_t ConstructBPGraphVecPairMap(
        const ChimericFragmentVecMap_t & fvMap )
{
    fprintf(stderr, "Constructing Breakpoint Graph Pairs ...\n");
    BPGraphVecPairMap_t gvPairMap;
    for (const auto & label_fv_pair : fvMap) {
        BPGraphVecPair_t gvPair = ConstructBPGraphVecPair(label_fv_pair.second);
        //Skip low support pairs
        if(gvPair.first.empty()) { continue; }
        gvPairMap.insert({label_fv_pair.first,std::move(gvPair)});
    }
    fprintf(stderr, "Constructed %lu Breakpoint Graph Pairs\n",gvPairMap.size());
    return gvPairMap;
}

//void ConnectBPGraphs(BPGraphMap_t & graphMap) {
//    fprintf(stderr,"Constructing edges in the graphs ...\n");
//    for(auto & pair : graphMap){
//        fprintf(stderr,"Constructing edges for the graph on %s ... %-10s\r",pair.first.c_str(),"");
//        pair.second.constructEdges();
//    }
//    fprintf(stderr,"\n");
//}

//Takes the cliques generated in the graph map and builds regions from them
//Which are then placed into a bipartite graph of host and viral regions
CRegionGraph ConstructRegionGraph(const BPGraphVecPair_t & graphVecPair) {
    //fprintf(stderr, "Constructing Region Graph ... \n");
    CRegionGraph regGraph;
    for(const BPGraphVec_t * vec_ptr :
            {&graphVecPair.first,&graphVecPair.second} )
    {
        std::string chr = vec_ptr->front().get_chromosome();
        bool opensLeft = vec_ptr->front().opens_left();
        bool isHost = (vec_ptr == &graphVecPair.second);
        for(const CBPGraph & graph : *vec_ptr){
            std::map<size_t,CRegionGraph::VertexProps> regionPropMap;
            //Construct the Regions from the cliques in the graph
            for(igraph_int_t id = 0; id < graph.vcount(); id++){
                CBPGraph::VertexProps fragGrpProp = graph.get_vertex_properties(id);
                //Construct a fragment Group String compatible with CRegionGraphs
                std::string fragGrpStr = strjoin(   fragGrpProp.assocFragments.begin(),
                                                    fragGrpProp.assocFragments.end(),
                                                    CRegionGraph::DupDelim);
                size_t left = std::min(fragGrpProp.proximalPos,fragGrpProp.distalPos);
                size_t right = std::max(fragGrpProp.proximalPos,fragGrpProp.distalPos);
                //Iterate over cliques for this id
                for( int cID : fragGrpProp.cliques) {
                    if(!regionPropMap.count(cID)){
                        CRegionGraph::VertexProps regProp = {   cID, chr, opensLeft,
                                                                false, isHost, 
                                                                size_t(~0),
                                                                0, {}};
                        regionPropMap.emplace(cID,regProp);
                    }
                    CRegionGraph::VertexProps & regProp = regionPropMap.at(cID);
                    regProp.fromSplit |= fragGrpProp.isSplit;
                    regProp.assocFragGroups.push_back(fragGrpStr);
                    if(left < regProp.left) { regProp.left = left; }
                    if(right > regProp.right) { regProp.right = right; }
                }
            }
            //Add the regions to the region graph
            for(const auto & pair : regionPropMap) {
                //std::cerr << "Vertex " << regGraph.vcount() << "\n";
                //std::cerr << strjoin(pair.second.assocFragGroups.begin(),pair.second.assocFragGroups.end(),'\t') << "\n";
                regGraph.addOrUpdateVertex( pair.second);
            }
        }
    }
    ////Explicitly request construction of edges
    regGraph.ensureConstructed();
    //fprintf(stderr, "Region graph with %d regions and %d edges created\n",regGraph.vcount(),regGraph.ecount());
    return regGraph;
}


CRegionGraph ConstructAndFilterRegionGraph(
        const BPGraphVecPair_t & graphVecPair)
{
    CRegionGraph regGraph = ConstructRegionGraph(graphVecPair);
    //fprintf(stderr,"Merging uninformatively different overlapping regions ...\n");
    regGraph.mergeUninformitiveOverlap();
    //fprintf(stderr,"After merging, %d regions and %d edges remain ...\n", regGraph.vcount(), regGraph.ecount());
    //fprintf(stderr,"Filtering low support edges ...\n");
    regGraph.filterEdges(MinimumReads,SplitBonus);
    //fprintf(stderr,"After filtering, %d edges remain\n",regGraph.ecount());
    return regGraph;
}


RegGraphVec_t ConstructRegionGraphVec(
        const BPGraphVecPairMap_t & graphVecPairMap)
{
    fprintf(stderr,"Constructing Paired Region Graphs ...\n");
    RegGraphVec_t prgMap;
    size_t regCounter = 0;
    size_t edgeCounter = 0;
    for(const auto & pair : graphVecPairMap) {
        CRegionGraph regGraph  = ConstructRegionGraph(pair.second);//ConstructAndFilterRegionGraph(pair.second);
        //Skip graphs with no edges
        if(regGraph.ecount() == 0) { continue; }
        regCounter += regGraph.vcount();
        edgeCounter += regGraph.ecount();
        prgMap.push_back(std::move(regGraph));
    }
    fprintf(stderr,
            "Constructed %lu edges between %lu regions across %lu paired region graphs\n",
            edgeCounter,regCounter,prgMap.size());
    return prgMap;
}

BPGraphMap_t DecomposeBPGraph(  const std::string & baseLabel,
                                CBPGraph & graph)
{
    BPGraphMap_t res;
    //std::cerr << "\tStart Decompose for " << baseLabel << "\n";
    std::vector<CBPGraph> resVec = graph.decompose(int(MinimumReads - SplitBonus));
    for(size_t i = 0; i < resVec.size(); i++){
        std::string label = baseLabel + "_" + std::to_string(i);
        res.emplace(label, std::move(resVec[i]));
    }
    //std::cerr << "\tEnd Decompose for " << baseLabel << "\n";
    return res;
}

BPGraphVec_t DecomposeBPGraph(CBPGraph & graph)
{
    return graph.decompose(int(MinimumReads - SplitBonus));
}

//Separates each graph in the graph map into separate graphs as connected
//components, only components with enough vertexes are retained
//Output - true if there are still graphs remaining, false otherwise
bool DecomposeBPGraphs(BPGraphMap_t & graphMap) {
    fprintf(stderr,"Decomposing graphs into connected components...\n");
    BPGraphMap_t tmp;
    while(graphMap.size()) {
        auto it = graphMap.begin();
        tmp.merge(DecomposeBPGraph(it->first,it->second));
        graphMap.erase(it);
    }
    std::swap(tmp,graphMap);
    fprintf(stderr,"After decomposition, there are %lu graphs...\n",graphMap.size());
    return (graphMap.size() > 0);
}

bool DecomposeBPGraphs(BPGraphMap_t & graphMap, ctpl::thread_pool & threadPool) {
    fprintf(stderr,"Decomposing Breakpoint Graphs ...\n");
    BPGraphMap_t tmp;
    std::vector<std::future<BPGraphMap_t>> decompFutureVec;
    for(auto it = graphMap.begin(); it != graphMap.end(); it++) {
        std::future<BPGraphMap_t> future = threadPool.push(
                [it](int id) { return DecomposeBPGraph(it->first,it->second); } ); 
        decompFutureVec.push_back(std::move(future));
    }
    //Retain the split graphs
    for(auto & future : decompFutureVec){
        tmp.merge(future.get());
    }
    fprintf(stderr,"After decomposition, there are %lu graphs ...\n",tmp.size());
    if(!tmp.size()){ return false; }
    std::swap(tmp,graphMap);
    return true;
}

bool DecomposeBPGraphVec(BPGraphVec_t & graphVec) {
   BPGraphVec_t all_res;
   for(CBPGraph & graph : graphVec){
       //Break the graph into connected components of sufficient size
       BPGraphVec_t local_res = DecomposeBPGraph(graph);
       all_res.insert( all_res.end(),
                       std::make_move_iterator(local_res.begin()),
                       std::make_move_iterator(local_res.end()) );
       //locar_res goes out of scope and doesn't need to be erased
   }
   graphVec = std::move(all_res);
   return !graphVec.empty();
}

//bool DecomposeBPGraphVecPairMap(BPGraphVecPairMap_t & graphVecPairMap) {
//    fprintf(stderr,"Decomposing Breakpoint Graphs ...\n");
//    size_t hostCounter = 0;
//    size_t virusCounter = 0;
//    //for( auto & pair : graphVecPairMap){
//    for( auto it = graphVecPairMap.begin(); it != graphVecPairMap.end(); ) {
//        //const StrandLabelPair_t & label = pair.first;
//        BPGraphVecPair_t & graphVecPair = it->second;
//        for(BPGraphVec_t * vec_ptr :
//                {&graphVecPair.first, &graphVecPair.second})
//        {
//            BPGraphVec_t all_res;
//            for(CBPGraph & graph : *vec_ptr){
//                //Break the graph into connected components of sufficient size
//                BPGraphVec_t local_res = DecomposeBPGraph(graph);
//                all_res.insert( all_res.end(),
//                                std::make_move_iterator(local_res.begin()),
//                                std::make_move_iterator(local_res.end()) );
//                //locar_res goes out of scope and doesn't need to be erased
//            }
//            *vec_ptr = std::move(all_res);
//        }
//        //Filter pairs where either the host or viral graph has no connected components of sufficient size
//        if(graphVecPair.first.empty() || graphVecPair.second.empty()){
//            it = graphVecPairMap.erase(it);
//        } else {
//            it++;
//            hostCounter += graphVecPair.first.size();
//            virusCounter += graphVecPair.second.size();
//        }
//    }
//    fprintf(stderr,
//            "After decomposition, there are %lu host graphs and %lu viral graphs across %lu graph pairs ...\n",
//            hostCounter, virusCounter,graphVecPairMap.size() );
//    return !graphVecPairMap.empty(); 
//}

//Given a graph, filters undersupported vertexes, and if sufficient support remains
// then it will cluster vertexes
//Inputs - an arbitrary id for the function call
//       a CBPGraph object on which to operate
//Output - true if there are cliques with sufficent support, false otherwise
bool FilterAndClusterBPGraph(CBPGraph & graph) {
    if(!FilterBPGraph(graph)){ return false; }
    //return ClusterBPGraph(graph);
    return true;
}



//Ensures that every vertex within a graph has sufficient edges to contribute to 
//  a valid region, then ensures each graph has sufficient nodes to contibute to
//  a valid region
//Inputs    - a graph to filter
//Output    - true if the graph still has sufficient support, false otherwise
bool FilterBPGraph(CBPGraph & graph) {
    graph.filterVertices(MinimumReads, SplitBonus);
    if(size_t(graph.vcount() + SplitBonus) < MinimumReads){
        return false;
    }
    return true;
}

//Filters each graph in a graph map
//Inputs    - a graph map containing graphs to filter
//Output    - None, modifies the input
bool FilterBPGraphs(BPGraphMap_t & graphMap) {
    fprintf(stderr,"Filtering Initial BP Graphs ...\n");
    for(auto it = graphMap.begin(); it != graphMap.end(); ){
        if(FilterBPGraph(it->second)){ 
            it++;
        } else {
            it = graphMap.erase(it);
        }
    }
    fprintf(stderr,"Filtered Down to %lu graphs\n",graphMap.size());
    return (graphMap.size() > 0);
}

//Removes any graphs with insufficient verticies
//Output - true if there is at least one remaining graph, false otherwise
bool FilterBPGraphVec(BPGraphVec_t & graphVec) {
    auto removeIt = std::remove_if( graphVec.begin(),graphVec.end(),
                                    [](CBPGraph & graph) {
                                        return !FilterBPGraph(graph);
                                    } );
    graphVec.erase(removeIt,graphVec.end());
    return !graphVec.empty();
}


BPGraphMap_t LoadBPGraphs(const std::string & fname,bool bOnline) {
    //std::cerr << "Loading graphs ..." << "\n";
    fprintf(stderr,"Loading Graphs ...\n");
    std::ifstream in(fname);
    BPGraphMap_t graphByContig;
    std::string bedpeStr;
    size_t counter = 0;
    while(getline(in,bedpeStr)){
        if(++counter % 10000 == 1){
            fprintf(stderr,"At least %lu fragments loaded\r",counter);
        }
        ChimericFragment_t frag = ChimericFragment_t::from_bedpe(bedpeStr);
        //Each fragment implies two regions: The Host and viral
        for( ChimericFragment_t::IV_IDX ivIdx :
                {ChimericFragment_t::IV1, ChimericFragment_t::IV2} ) 
        {
            std::string contig =    frag.getChr(ivIdx) +
                                    ((frag.opens_left(ivIdx)) ? "L" : "R");
            if(!graphByContig.count(contig)){
                CBPGraph graph( frag.getChr(ivIdx), frag.opens_left(ivIdx),
                                UpstreamSize, ReadLength,
                                MaxInsertSize, SplitFactor);
               graphByContig.insert({contig, std::move(graph)}); 
            }
            CBPGraph & graph = graphByContig.at(contig);
            graph.addOrUpdateVertex( frag.proximal_pos(ivIdx),
                                     frag.distal_pos(ivIdx),
                                     frag.is_split(ivIdx),
                                     frag.getName(),bOnline);
        }
    }
    
    fprintf(stderr,"\nLoaded %lu graphs\n",graphByContig.size());
    return graphByContig;
}



ChimericFragmentVecMap_t LoadFragments(const std::string & fname) {
    fprintf(stderr,"Loading Fragments ...\n");
    ChimericFragmentVecMap_t fvMap;
    std::ifstream in(fname);
    std::string bedpeStr;
    size_t counter = 0;
    while(getline(in,bedpeStr)){
        if(++counter % 10000 == 1){
            fprintf(stderr,"At least %lu fragments loaded\r",counter);
        }
        //Construct Fragment
        ChimericFragment_t frag = ChimericFragment_t::from_bedpe(bedpeStr);
        //Construct Label
        StrandLabelPair_t label;
        for( ChimericFragment_t::IV_IDX ivIdx :
                {ChimericFragment_t::IV1, ChimericFragment_t::IV2} ) 
        {
            std::string * contig_ptr = (ivIdx == ChimericFragment_t::IV1) ?
                                    &(label.first) : &(label.second);
            *contig_ptr =   frag.getChr(ivIdx) +
                            ((frag.opens_left(ivIdx)) ? "L" : "R");
        }
        //Add the fragment to the appropriate vector
        if(!fvMap.count(label)){
            fvMap.emplace(label,ChimericFragmentVec_t());
        }
        fvMap.at(label).push_back(frag);
    }
    fprintf(stderr,"Loaded %lu fragments across %lu host-virus contig pairs\n",
            counter,fvMap.size());
    return fvMap;
}

//Single threaded implementation which processes from raw chimeic fragments 
// through to a final merged region graph
RegGraphVec_t FragmentMapToRegionGraphVec(const ChimericFragmentVecMap_t & fvMap) {
     BPGraphVecPairMap_t graphVecPairMap =  ProcessFragments(fvMap);
    //REGION GRAPH TO ID EDGES
    return ConstructRegionGraphVec(graphVecPairMap);
}


//Multithreaded threaded implementation which processes from raw chimeic fragments 
// through to a final merged region graph
RegGraphVec_t FragmentMapToRegionGraphVec(size_t nThread, 
                                    const ChimericFragmentVecMap_t & fvMap)
{
    fprintf(stderr,"Constructing Region Graphs from Fragments ...\n");
    ctpl::thread_pool threadPool(nThread);
    RegGraphVec_t regGraphVec;
    std::vector<std::future<CRegionGraph>> futureVec;
    for(const auto & pair : fvMap){
        std::future<CRegionGraph> future = threadPool.push(
                [&threadPool,&pair](int id) {
                    return FragmentsToRegionGraph(pair.second);
                } );
        futureVec.push_back(std::move(future));
    }
    size_t ecount =0;
    size_t vcount =0;
    size_t counter = 0;
    for(auto & future : futureVec){
        CRegionGraph reg = future.get();
        fprintf(stderr,"Completed at least %ld of %ld graphs%-10s\r",++counter,futureVec.size(),"");
        if(reg.ecount()){
            ecount += reg.ecount();
            vcount += reg.vcount();
            regGraphVec.push_back(std::move(reg));
        }
    }
    fprintf(stderr,
            "Constructed %ld disjoint graphs with %ld edges between %ld regions\n",
            regGraphVec.size(), ecount, vcount );
    return regGraphVec;
}


CRegionGraph FragmentsToRegionGraph(const ChimericFragmentVec_t & fragVec) {
    BPGraphVecPair_t gvPair = ConstructBPGraphVecPair(fragVec);
    if(gvPair.first.empty()) { return CRegionGraph(); }
    if( !DecomposeBPGraphVec(gvPair.first) ||
        !DecomposeBPGraphVec(gvPair.second) )
    {
        return CRegionGraph();
    }
    if( !FilterBPGraphVec(gvPair.first) ||
        !FilterBPGraphVec(gvPair.second) )
    {
        return CRegionGraph();
    }
    if( !ClusterBPGraphVec(gvPair.first) ||
        !ClusterBPGraphVec(gvPair.second) )
    {
        return CRegionGraph();
    }
    return ConstructRegionGraph(gvPair);//ConstructAndFilterRegionGraph(gvPair);
}

bool OperateOnBPGraphVecPairMap(BPGraphVecPairMap_t & graphVecPairMap,
                                const std::string & operationName,
                                std::function<bool(BPGraphVec_t &)> operation)
{
    fprintf(stderr,"Performing %s on Breakpoint Graph Pairs ...\n",
            operationName.c_str());
    size_t hostCounter = 0;
    size_t virusCounter = 0;
    size_t hostFragmentCounter = 0;
    size_t virusFragmentCounter = 0;
    size_t maxHFrag = 0;
    size_t maxVFrag = 0;
    //for( auto & pair : graphVecPairMap){
    for( auto it = graphVecPairMap.begin(); it != graphVecPairMap.end(); ) {
        //const StrandLabelPair_t & label = pair.first;
        BPGraphVecPair_t & graphVecPair = it->second;
        //Filter pairs where either the operation fails for either the host or viral graphs
        if(!operation(graphVecPair.first) || !operation(graphVecPair.second)){
            it = graphVecPairMap.erase(it);
        } else {
            it++;
            hostCounter += graphVecPair.first.size();
            virusCounter += graphVecPair.second.size();
            for(auto & g : graphVecPair.first){
                hostFragmentCounter += g.vcount();
                if(size_t(g.vcount()) > maxHFrag){ maxHFrag = g.vcount(); }
            }
            for(auto & g : graphVecPair.second){
                virusFragmentCounter += g.vcount();
                if(size_t(g.vcount()) > maxVFrag){ maxVFrag = g.vcount(); }
            }
        }
    }
    fprintf(stderr,
            "After %s, there are (%lu , %lu, %lu) host and (%lu , %lu, %lu) viral (f,g,max(f/g)) across %lu graph pairs ...\n",
            operationName.c_str(),hostFragmentCounter,hostCounter, maxHFrag, virusFragmentCounter, virusCounter, maxVFrag, graphVecPairMap.size() );
    return !graphVecPairMap.empty();
}



#ifndef NDEBUG
void OutputDebugBPGraph(const BPGraphMap_t & graphMap,
                        std::string BPAdjFileName,
                        std::string BPVertFileName)
{ 
    FILE* adjFile_ptr = fopen(BPAdjFileName.c_str(), "w");
    std::ofstream vertFile(BPVertFileName);
    for(auto & pair : graphMap){
        fprintf(adjFile_ptr,"===%s\n",pair.first.c_str());
        pair.second.write_edgelist(adjFile_ptr);
        vertFile << "===" << pair.first << "\n";
        for(int id = 0; id < pair.second.vcount();id++){
            vertFile << pair.second.get_vertex_properties(id).to_string() <<
                        "\n";
        }
    }
    fclose(adjFile_ptr);
}

void OutputDebugRegGraph(   const CRegionGraph & regGraph,
                            std::string RegAdjFileName,
                            std::string RegVertFileName)
{
    FILE* adjFile_ptr = fopen(RegAdjFileName.c_str(), "w");
    regGraph.write_edgelist(adjFile_ptr);
    fclose(adjFile_ptr);
    std::ofstream vertFile(RegVertFileName);
    for(int id = 0; id < regGraph.vcount();id++){
        vertFile << regGraph.get_vertex_properties(id).to_string() << "\n";
    }
}

#endif //NDEBUG


void OutputResults( const CRegionGraph & regGraph,
                    const std::string & regFileName,
                    const std::string & edgeFileName,
                    const std::string & assocFileName)
{
    fprintf(stderr,"Outputting Results ... \n");
        // Open output stream 
    std::ofstream regionBedFile(regFileName);
    std::ofstream edgeTabFile(edgeFileName);
    std::ofstream assocTabFile(assocFileName);
    // Track unique regions
    size_t nAssoc = 0;
    std::set<std::string> seenRegions;
    std::set<std::string> seenFragments;
    for(int eid = 0; eid < regGraph.ecount();eid++){
        CRegionGraph::EdgeProps prop = regGraph.get_edge_properties(eid);
        //Output the fragment-edge associations
        for(size_t fragGrpIdx = 0; fragGrpIdx < prop.assocFragGroups.size(); fragGrpIdx++) {
            const std::string & fragGrp = prop.assocFragGroups[fragGrpIdx];
            std::vector<std::string> fragNameVec = strsplit(fragGrp,CRegionGraph::DupDelim);
            for(const std::string & fragName : fragNameVec){
                assocTabFile << fragName << "\t" << eid << "\t" << fragGrpIdx << "\n";
                seenFragments.insert(fragName);
                nAssoc++;
            }
        }
        //Get Endpoints of the edge, and the region defining information
        std::pair<igraph_int_t,igraph_int_t> endpoints =
            regGraph.get_edge_endpoints(eid);
        //Output the edge information
        edgeTabFile << eid << "\t" << endpoints.first << "\t" << endpoints.second << "\t" << prop.weight << "\n";
        //Output each region
        for(auto & prop : {  regGraph.get_vertex_properties(endpoints.first),
                            regGraph.get_vertex_properties(endpoints.second)} )
        {
            std::string regStr = to_bed(prop);
            auto res = seenRegions.insert(regStr);
            if(res.second){
                regionBedFile << regStr << "\n";
            }
        }
    }
    fprintf(stderr,"Wrote %lu Unique Regions to %s\n",seenRegions.size(),regFileName.c_str());
    fprintf(stderr,"Wrote %d Edges to %s\n",regGraph.ecount(),edgeFileName.c_str());
    fprintf(stderr,"Wrote %lu associations to %lu unique fragments to %s\n",nAssoc,seenFragments.size(),edgeFileName.c_str());
}


//Single theaded processing of BPGraphs
bool ProcessBPGraphs(BPGraphMap_t & graphMap) {
    if(!DecomposeBPGraphs(graphMap)) { return false ; }
    if(!FilterBPGraphs(graphMap)) { return false;}
    //return ClusterBPGraphs(graphMap);
    return true;
}

//Multithreaded processing of BP graphs to identify regions
bool ProcessBPGraphs(size_t nThread, BPGraphMap_t & graphMap)
{
    ctpl::thread_pool threadPool(nThread);
    //In parallel Decompose each graph, and stop if there are no graphs left
    if(!DecomposeBPGraphs(graphMap,threadPool) ) { return false;}
    fprintf(stderr,"Filtering and Clustering Breakpoint Graphs ...\n");
    //Launch processes for each graph
    std::map<std::string,std::future<bool>> futureMap;
    for(auto it = graphMap.begin(); it != graphMap.end(); it++){
        std::future<bool> future = threadPool.push(
                [it](int id) {return FilterAndClusterBPGraph(it->second);} );
        futureMap.insert({it->first,std::move(future)});
    }
    //Note which graphs have insufficient support
    std::vector<std::string> toFilter;
    size_t counter = 0;
    for(auto & pair : futureMap){
        fprintf(stderr,"Waiting for %s (%lu/%lu) ... %-10s\r",
                pair.first.c_str(),counter++,futureMap.size(),"");
        if(!pair.second.get()){
            toFilter.push_back(pair.first);
        }
    }
    //Remove the low support graphs AFTER all processing is done
    for(const std::string & contig : toFilter){
        graphMap.erase(contig);
    }
    fprintf(stderr,"After clustering, %lu Breakpoint Graphs remain\n", graphMap.size());
    return (graphMap.size() > 0);
}


//Single threaded version processing fragments into BPGraphs
BPGraphVecPairMap_t ProcessFragments(const ChimericFragmentVecMap_t & fvMap) {
    BPGraphVecPairMap_t graphVecPairMap = ConstructBPGraphVecPairMap(fvMap);
    if(!OperateOnBPGraphVecPairMap( graphVecPairMap,"Decomposition",
                                    DecomposeBPGraphVec) )
    {
        return graphVecPairMap;
    }
    if(!OperateOnBPGraphVecPairMap( graphVecPairMap,"Filtering",
                                    FilterBPGraphVec)) {
        return graphVecPairMap;
    }
    if(!OperateOnBPGraphVecPairMap( graphVecPairMap,"Clustering",
                                    ClusterBPGraphVec)) {
        return graphVecPairMap;
    }
    return graphVecPairMap;
}

//Takes a set of vertex properties defining a region, and constructs a bed formated string
std::string to_bed(CRegionGraph::VertexProps props) {
    std::string bed = props.chromosome;
    bed += '\t' + std::to_string(props.left);
    bed += '\t' + std::to_string(props.right);
    bed += '\t' + std::to_string(props.id);
    uint16_t flag = (ChimericFragment_t::HAS_INTERVAL);
    flag |= (props.opensLeft) ? ChimericFragment_t::OPENS_LEFT : 0;
    flag |= (props.fromSplit) ? ChimericFragment_t::IS_SPLIT : 0;
//    flag |= (props.fromSplit) ? ChimericFragment_t::IS_SPLIT : 0;
    bed += '\t' + std::to_string(int(flag));
    bed += '\t' + std::string((props.opensLeft == props.isHost) ? "-" : "+");
    return bed;
}

