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
#include "htslib/sam.h"

//==== TYPE DECLARATIONS

struct ReadRegionAssoc_t {
    std::string readName;
    uint16_t flag;
    igraph_int_t regionId;
    std::string to_string() const {
        return  readName + "\t" + std::to_string(regionId) + "\t" +
                std::to_string(int(flag));
    }

};

typedef std::vector<ReadRegionAssoc_t> ReadRegionAssocVec_t;

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

//Type for looking up all fragments by name
typedef std::unordered_map<std::string,std::vector<const ChimericFragment_t*>>
        ChimericFragmentIndex_t;

//Type to contain a vector of CBPGraphs (post decomposition)
typedef std::vector<CBPGraph> BPGraphVec_t;
//Type to contain a host-virus pair of Graph Vectors
typedef std::pair<BPGraphVec_t,BPGraphVec_t> BPGraphVecPair_t;
//Type to map host-virus contig labels to corresponding graph vector pairs
typedef std::unordered_map<StrandLabelPair_t,BPGraphVecPair_t,
                            StrandLabelPair_HashFunctor,
                            StrandLabelPair_EqualFunctor> BPGraphVecPairMap_t;

typedef std::vector<CRegionGraph> RegGraphVec_t;

//==== GLOBAL VARIABLE DECLARATIONS

int             MaxInsertSize;
const size_t    MinimumReads = 4;
int             ReadLength;
const int       SplitBonus = 1;
const double    SplitFactor = 2.0;
int             UpstreamSize = 5;
bool	        BCliqueClustering = false;

//==== FUNCTION DECLARATIONS

//Retval                    Function Name
bool                        ChimericFragmentOverlapsRegion(
                                const ChimericFragment_t & frag,
                                ChimericFragment_t::IV_IDX ivIdx,
                                const CRegionGraph::VertexProps & regProp);
bool                        ClusterBPGraph(CBPGraph & graph);
BPGraphVecPair_t            ConstructBPGraphVecPair(
                                const ChimericFragmentVec_t & fragVec);
BPGraphVecPairMap_t         ConstructBPGraphVecPairMap(
                                const ChimericFragmentVecMap_t & fvMap );
ReadRegionAssocVec_t        ConstructReadRegionAssociations(
                                const CRegionGraph & regGraph,
                                const ChimericFragmentIndex_t & cfIndex);
CRegionGraph                ConstructRegionGraph(
                                const BPGraphVecPair_t & graphMap);
RegGraphVec_t               ConstructRegionGraphVec(
                                const BPGraphVecPairMap_t & graphVecPairMap);
BPGraphVec_t                DecomposeBPGraph(CBPGraph & graph);
bool                        DecomposeBPGraphVec(BPGraphVec_t & graphVec);
bool                        FilterBPGraph(CBPGraph & graph);
bool                        FilterBPGraphVec(BPGraphVec_t & graphVec);
RegGraphVec_t               FragmentMapToRegionGraphVec(
                                const ChimericFragmentVecMap_t & fvMap);
RegGraphVec_t               FragmentMapToRegionGraphVec(size_t nThread,
                                const ChimericFragmentVecMap_t & fvMap);
CRegionGraph                FragmentsToRegionGraph(
                                const ChimericFragmentVec_t & fragVec);
ChimericFragmentIndex_t     IndexChimericFragments(
                                const ChimericFragmentVecMap_t & fvMap);
ChimericFragmentVecMap_t    LoadFragments(const std::string & fname);
bool                        OperateOnBPGraphVecPairMap(
                                BPGraphVecPairMap_t & graphVecPairMap,
                                const std::string & operationName,
                                std::function<bool(BPGraphVec_t &)> operation);
BPGraphVecPairMap_t         ProcessFragments(const ChimericFragmentVecMap_t & fvMap);
void                        OutputResults(
                                const CRegionGraph & regGraph,
                                const ReadRegionAssocVec_t & rrAssocVec,
                                const std::string & regFileName,
                                const std::string & edgeFileName,
                                const std::string & assocFileName,
                                const std::string & regAssocFileName);
std::string                 to_bed(CRegionGraph::VertexProps);

//==== MAIN

//Parses candidate junctions into a pair of graph of breakpoint positions for each
// host-virla strand pair
// Breakpoints defined by fragments with the same proximal and distal positions
//  are counted as one breakpoint position (all such fragments form a fragment group)
// Those are then clustered into regions by identifying maximal cliques
// Cliques are filtered for minimum size
// The regions are then used to construct a bipartite graph with edges between host and
//  viral regions which have fragments associated with both
// The edges are filtered for minimum support
// Regions with edges are output to a bed file, with the regionID as the name field
//  The score field is a decimal flag:
//      0x1	1	Set if the region opens left
//      0x2	2	Set if the region was supported by split reads
//      0x4	4	Always Set
//      0x80	128	Set if the region is a host region
// Edges are output to a tab delim file:
//  edgeID, regionIDHost, regionIDVirus,Support
// Fragment-edge associations are output to a tab delim file:
//  fragment name, edgeID, grpID 
//  a given fragment may occur in multiple fragment groups associated with an edge
//  the amount of support for an edge is then the number of unique grpIDs associated with that edgeID
// Read-region associations are output to a tab delim file:
//  readName, regionID, flag
//  Flag is a decimal flag: (Like BAM flag)
//  0x40    64  READ1
//  0x80    128 READ2
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
    std::string workdir = argv[1];
    std::string workspace = argv[2];

    //## Files to be used from the workdir
    std::string stats_file_name = workspace + "/stats.txt";
    std::string candidate_file_name = workdir + "/junction-candidates.bedpe";
    std::string config_file_name = workdir + "/config.txt";
    //## Output Files
    std::string regFileName = workdir + "/region-candidates.bed";
    std::string edgeFileName = workdir + "/edge-candidates.tab";
    std::string assocFileName = workdir + "/fragment-edge-associations.tab";
    std::string regAssocFileName = workdir + "/read-region-associations.tab";

    //Load global variables
    MaxInsertSize = parse_stats(stats_file_name).max_is;
    auto config = parse_config(config_file_name);
    ReadLength = config.read_len;
    size_t nThread = config.threads;
    BCliqueClustering = config.clique;

    //Initialize igraph
    igraph_setup();

    //Load Fragment Data
    ChimericFragmentVecMap_t fragVecMap = LoadFragments(candidate_file_name);

    //Construct Region graphs for each host-virus strand combo
    RegGraphVec_t regGraphVec = (nThread == 1) ? 
                                FragmentMapToRegionGraphVec(fragVecMap) :
                                FragmentMapToRegionGraphVec(nThread, fragVecMap);
    //Construct one joint genome wide region graph
    fprintf(stderr,"Merging Paired Region Graphs ...\n");
    CRegionGraph regGraph = CRegionGraph::merge_graphs(regGraphVec);
    regGraph.ensureConstructed(); //Explicit Call 
    regGraph.mergeUninformitiveOverlap();
    regGraph.filterEdges(MinimumReads,SplitBonus);
    fprintf(stderr,
            "After merging and filtering, The region graph contains %d edges between %d regions  ...\n",
            regGraph.ecount(),regGraph.vcount());

    //TODO: It can occur where two indistinguishable regions with seprate ids are created
    //  this can cause some downstream issues (mostly handled), with the biggest
    //  extant problem being bloat from duplicate entries. 
    //  The todo list item is to either:
    //      ID why this can happen and eliminate it
    //      Identify that it has happened and collapse the duplicate regions together
    //      updating associated information in edges as well

    OutputResults(  regGraph, ConstructReadRegionAssociations(regGraph,
                                IndexChimericFragments(fragVecMap) ),
                    regFileName,edgeFileName,assocFileName,regAssocFileName);
           
    fprintf(stderr,"Done - enumerate_edges\n");
}

//==== FUNCTION DEFINITIONS

//Given a chimeric fragment, an interval to look at, and a region
// tests if the corresponding interval and the region overlap
//Output - true if overlapping, false otherwise
bool ChimericFragmentOverlapsRegion(const ChimericFragment_t & frag,
                                    ChimericFragment_t::IV_IDX ivIdx,
                                    const CRegionGraph::VertexProps & regProp)
{
    if(frag.getChr(ivIdx) != regProp.chromosome) { return false; }
    if(frag.opens_left(ivIdx) != regProp.opensLeft) { return false; }
    size_t minRight = std::min(frag.getEnd(ivIdx),regProp.right);
    size_t maxLeft = std::max(frag.getOffset(ivIdx),regProp.left);
    return (minRight >= maxLeft);
}

//Identifies all maximal cliques within a graph
//Inputs    - a graph to filter
//Output    - true if the graph still has sufficient support, false otherwise
bool ClusterBPGraph(CBPGraph & graph) {
    if(BCliqueClustering){
        return graph.maximalCliques(MinimumReads, SplitBonus);
    } else { //Otherwise go with connected components
        return true;
    }
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

//Given a Region graph with edges which associate chimeric fragment names between regions
//  and an index of all chimeric fragments with the same name, uses the FROM_R1 and FROM_R2
//  bits in the Chimeric fragment index to assocociate particular reads with with 
//  particular regions
//Inputs - a regionGraph with edges with host and viral endpoints
//       - a mapping from fragment names to chimeric fragment object pointers
//Output - a vector or ReadRegionAssoc Options 
ReadRegionAssocVec_t ConstructReadRegionAssociations(
                        const CRegionGraph & regGraph,
                        const ChimericFragmentIndex_t & cfIndex)
{
    
    std::vector<ReadRegionAssoc_t> rrAssocVec;
    std::unordered_set<std::string> seenRR;
    for(int eid = 0; eid < regGraph.ecount();eid++){
        CRegionGraph::EdgeProps prop = regGraph.get_edge_properties(eid);
        std::unordered_set<std::string> seenFrags;
        for(const std::string & fragGroup : prop.assocFragGroups){
            for(const std::string & fName : 
                    strsplit(fragGroup,CRegionGraph::DupDelim))
            {
                //Skip fragments already processed for this edge
                if(!seenFrags.insert(fName).second) { continue; }
                if(!cfIndex.count(fName)){
                    throw std::runtime_error(
                            "Encountered an fragment Name " + fName +
                            " in the region graph not in the junction candidates");
                }
                //Cache for region properties
                std::unordered_map<igraph_int_t,CRegionGraph::VertexProps>
                    vPropMap;
                for(igraph_int_t rid : {    prop.endpoints.first,
                                            prop.endpoints.second } )
                {
                    if(!vPropMap.count(rid)){
                        vPropMap[rid] = regGraph.get_vertex_properties(rid);
                    }
                    ChimericFragment_t::IV_IDX ivIdx =
                        (rid == prop.endpoints.first) ?
                            ChimericFragment_t::IV1 : ChimericFragment_t::IV2;
                    //Skip read region associations already observed
                    if(!seenRR.insert(fName + std::to_string(rid)).second) {
                        continue;
                    }
                    ReadRegionAssoc_t rrAssoc = {fName,0,rid};
                    //Iterate over fragments with this name and find which (if any) overlap this region
                    for(const ChimericFragment_t * fragPtr : cfIndex.at(fName)){
                        if(ChimericFragmentOverlapsRegion(  *fragPtr,ivIdx,
                                                            vPropMap[rid]) )
                        {
                            if(fragPtr->fromR1(ivIdx)) {
                                rrAssoc.flag |= BAM_FREAD1;
                            }
                            if(fragPtr->fromR2(ivIdx)) {
                                rrAssoc.flag |= BAM_FREAD2;
                            }
                        }
                    }
                    if(rrAssoc.flag == 0){
                        throw std::runtime_error(
                                "Encountered a region (" +
                                vPropMap[rid].to_string() +
                                ") associated non-overlapping junction candidates (" +
                                fName + ")");
                    }
                    rrAssocVec.push_back(rrAssoc);
                }
            }
        }
    }
    return rrAssocVec;
}

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
        bool isHost = (vec_ptr == &graphVecPair.first);
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
                regGraph.addOrUpdateVertex( pair.second);
            }
        }
    }
    ////Explicitly request construction of edges
    regGraph.ensureConstructed();
    //fprintf(stderr, "Region graph with %d regions and %d edges created\n",regGraph.vcount(),regGraph.ecount());
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
        CRegionGraph regGraph  = ConstructRegionGraph(pair.second);
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

BPGraphVec_t DecomposeBPGraph(CBPGraph & graph)
{
    return graph.decompose(int(MinimumReads - SplitBonus));
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


ChimericFragmentIndex_t IndexChimericFragments(
                                const ChimericFragmentVecMap_t & fvMap)
{
    ChimericFragmentIndex_t fragVecByNameMap;
    for(const auto & pair : fvMap){
        //const StrandLabelPair_t & label = pair.first;
        const ChimericFragmentVec_t & fv = pair.second;
        for(const ChimericFragment_t & frag : fv){
            const std::string & fName = frag.getName();
            fragVecByNameMap[fName].push_back(&frag);
        }
    }
    return fragVecByNameMap;
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
        //Add the fragment to the appropriate ve/ctor
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
    //Construct breakpoint graphs
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
    return ConstructRegionGraph(gvPair);
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

void OutputResults( const CRegionGraph & regGraph,
                    const ReadRegionAssocVec_t & rrAssocVec,
                    const std::string & regFileName,
                    const std::string & edgeFileName,
                    const std::string & assocFileName,
                    const std::string & regAssocFileName)
{
    fprintf(stderr,"Outputting Results ... \n");
        // Open output stream 
    std::ofstream regionBedFile(regFileName);
    std::ofstream edgeTabFile(edgeFileName);
    std::ofstream assocTabFile(assocFileName);
    std::ofstream rrAssocTabFile(regAssocFileName);
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
    //Output the read region associations
    std::unordered_set<std::string> seenFrag;
    for(const ReadRegionAssoc_t & rrAssoc : rrAssocVec){
        rrAssocTabFile << rrAssoc.to_string() << "\n";
        if(rrAssoc.flag & BAM_FREAD1){
            seenFrag.insert(rrAssoc.readName + "R1");
        }
        if(rrAssoc.flag & BAM_FREAD1){
            seenFrag.insert(rrAssoc.readName + "R2");
        }
    }
    fprintf(stderr,"Wrote %lu Unique Regions to %s\n",seenRegions.size(),
            regFileName.c_str() );
    fprintf(stderr,"Wrote %d Edges to %s\n",regGraph.ecount(),
            edgeFileName.c_str() );
    fprintf(stderr,
            "Wrote %lu edge associations to %lu unique fragments to %s\n",
            nAssoc, seenFragments.size(), assocFileName.c_str() );
    fprintf(stderr,"Wrote %lu region associations to %lu unique reads to %s\n",
            rrAssocVec.size(),seenFrag.size(),regAssocFileName.c_str());
    
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
    if(props.isHost){ flag |= (1 << ChimericFragment_t::FLAG_BITS); };
//    flag |= (props.fromSplit) ? ChimericFragment_t::IS_SPLIT : 0;
    bed += '\t' + std::to_string(int(flag));
    bed += '\t' + std::string((props.opensLeft == props.isHost) ? "-" : "+");
    return bed;
}

