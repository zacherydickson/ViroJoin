#include <algorithm>
#include <iostream>
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

struct jRegLabel_t {
    std::string chr;
    size_t distPos;
    size_t proxPos;
    bool opensLeft;
    bool isSplit;
    int compare(const jRegLabel_t & other, size_t dist = 0) const {
	if(this->chr != other.chr) {
	    return (this->chr < other.chr) ? -1 : 1;
	}
	if(this->opensLeft != other.opensLeft) {
	    return (this->opensLeft) ? 1 : -1;
	}
	if(this->proxPos > other.proxPos && this->proxPos - other.proxPos > dist) return 1;
	if(this->proxPos < other.proxPos && other.proxPos - this->proxPos > dist) return -1;
	if(this->distPos > other.distPos && this->distPos - other.distPos > dist) return 1;
	if(this->distPos < other.distPos && other.distPos - this->distPos > dist) return -1;
        if(this->isSplit != other.isSplit){
            return (this->isSplit) ? 1 : -1;
        }
	return 0;
    }
};

struct jRegLabel_EqFunctor {
    bool operator()(const jRegLabel_t & a, const jRegLabel_t & b) const {
	return (a.compare(b) == 0);
    }
};

struct jRegLabel_HashFunctor {
    size_t operator()(const jRegLabel_t & a) const {
	return std::hash<std::string>{}(
                a.chr + std::to_string(a.opensLeft << 1 | a.isSplit) +
                std::to_string(a.proxPos) + std::to_string(a.distPos)
        );
    }
};

struct junctionRegion_t {
    junctionRegion_t() : left(0), right(0), nSplit(0) {}
    junctionRegion_t(size_t pos)
	: left(pos), right(pos), nSplit(0) {}
    junctionRegion_t(const junctionRegion_t & other)
	:   left(other.left), right(other.right), nSplit(other.nSplit),
	    QNameSet(other.QNameSet) {}
    size_t left;
    size_t right;
    size_t nSplit;
    std::unordered_set<std::string> QNameSet;
};


typedef std::vector<jRegLabel_t> jRegLabelVector_t;
typedef std::unordered_map< jRegLabel_t,size_t,jRegLabel_HashFunctor,
			    jRegLabel_EqFunctor>
	    jRegLabelCount_t;
typedef std::unordered_map< jRegLabel_t,bool,jRegLabel_HashFunctor,
			    jRegLabel_EqFunctor>
	    jRegSplitStatus_t;
typedef std::unordered_map< jRegLabel_t,junctionRegion_t,jRegLabel_HashFunctor,
			    jRegLabel_EqFunctor>
	    jRegMap_t;
typedef std::unordered_set< jRegLabel_t, jRegLabel_HashFunctor,
			    jRegLabel_EqFunctor>
	    jRegLabelSet_t;
typedef std::unordered_map< jRegLabel_t, jRegLabelSet_t,
			    jRegLabel_HashFunctor,jRegLabel_EqFunctor> 
	    BestJRegSetMap_t;
typedef std::unordered_map< jRegLabel_t, jRegLabelSet_t,
			    jRegLabel_HashFunctor,jRegLabel_EqFunctor> 
	    MututalJRegSetMap_t;

typedef std::unordered_map<std::string,CBPGraph> BPGraphMap_t;

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
//void ConnectBPGraphs(BPGraphMap_t & graphMap);
CRegionGraph ConstructRegionGraph(const BPGraphMap_t & graphMap);
CRegionGraph ConstructAndFilterRegionGraph(const BPGraphMap_t & graphMap);
BPGraphMap_t DecomposeBPGraph(  const std::string & baseLabel,
                                CBPGraph & graph);
bool DecomposeBPGraphs(BPGraphMap_t & graphMap);
bool DecomposeBPGraphs(BPGraphMap_t & graphMap, ctpl::thread_pool & threadPool);
bool FilterBPGraph(CBPGraph & graph);
bool FilterBPGraphs(BPGraphMap_t & graphMap);
bool FilterAndClusterBPGraph(CBPGraph & graph);
bool ProcessBPGraphs(BPGraphMap_t & graphMap);
bool ProcessBPGraphs(size_t nThread, BPGraphMap_t & graphMap);
void FilterRegions(jRegMap_t & regionMap);
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

    jRegLabelVector_t labelVec;
    jRegLabelCount_t labelCount;

    igraph_setup();
    //BREAKPOINT GRAPH TO ID REGIONS
    bool bOnline = false;
    BPGraphMap_t graphMap = LoadBPGraphs(candidate_file_name,bOnline);

    //std::cerr
    //for(auto it = graphMap.begin(); it != graphMap.end();){
    //    if(it->first == "NC_007605.1R"){
    //        it++;
    //    } else {
    //        it = graphMap.erase(it);
    //    }
    //}

    if(nThread == 1){ //Single Threaded version - maybe avoid some overhead
        ProcessBPGraphs(graphMap);
    } else { //MultiThreaded Version
        ProcessBPGraphs(nThread,graphMap);
    }
#ifndef NDEBUG
    OutputDebugBPGraph(graphMap,BPAdjFileName,BPVertFileName);
#endif //NDEBUG
    //REGION GRAPH TO ID EDGES
    CRegionGraph regGraph  = ConstructAndFilterRegionGraph(graphMap);
#ifndef NDEBUG
    OutputDebugRegGraph(regGraph,RegAdjFileName,RegVertFileName);
#endif //NDEBUG
    OutputResults(regGraph,regFileName,edgeFileName,assocFileName);
           
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
CRegionGraph ConstructRegionGraph(const BPGraphMap_t & graphMap) {
    fprintf(stderr, "Constructing Region Graph ... \n");
    CRegionGraph regGraph;
    for( const auto & pair : graphMap ){
        const CBPGraph & graph = pair.second;
        std::string chr = graph.get_chromosome();
        bool opensLeft = graph.opens_left();
        bool isHost = !VirusNameSet.count(chr);
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
    fprintf(stderr, "Region graph with %d regions and %d edges created\n",regGraph.vcount(),regGraph.ecount());
    return regGraph;
}


CRegionGraph ConstructAndFilterRegionGraph(const BPGraphMap_t & graphMap) {
    CRegionGraph regGraph = ConstructRegionGraph(graphMap);
    fprintf(stderr,"Merging uninformatively different overlapping regions ...\n");
    regGraph.mergeUninformitiveOverlap();
    fprintf(stderr,"After merging, %d regions and %d edges remain ...\n", regGraph.vcount(), regGraph.ecount());
    fprintf(stderr,"Filtering low support edges ...\n");
    regGraph.filterEdges(MinimumReads,SplitBonus);
    fprintf(stderr,"After filtering, %d edges remain\n",regGraph.ecount());
    return regGraph;
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

//Outputs regions which have enough reads assigned,
//  enough is defined as the minimum reads
//  minus the split bonus if any split reads are present
//After this check, removes reads which now map to only host or only virus
//This pair of filters is repeated until no filtering occurs
//Then the final set of junctions is written
//Inputs - a mapping of region labels to regions
//	 - also uses global min reads and split bonus, and viral names
//Output - None, modifes the regionMap
void FilterRegions(jRegMap_t & regionMap){
    //FUTURE:
    //Current implementation does a lot of unecessary filtering
    //One could map each qname to a region and vice versa so that the
    //on each iteration only regions which may have changed are checked
    //Cost - extra complexity and memory
    bool bFilter;
    size_t filterRound= 0;
    do {
	fprintf(stderr,"Filtering Regions - Round %lu ...\n",filterRound++);
	bFilter=false;
	//First Pass Remove regions with too few reads
	//Containers tracking the number of host/viral regions a particular
    	//qname is present in After the first filter pass
    	std::unordered_map<std::string,size_t> hostCount;
    	std::unordered_map<std::string,size_t> viralCount;
    	for( auto it = regionMap.begin(); it != regionMap.end();){
    	    const jRegLabel_t & label = it->first;
    	    const junctionRegion_t & reg = it->second;
    	    size_t effectiveReads = reg.QNameSet.size();
    	    if(reg.nSplit) effectiveReads += SplitBonus;
    	    if(effectiveReads < MinimumReads){ //Fails filter remove
		bFilter=true;
		it = regionMap.erase(it);
    	    } else { // Keep the region
		std::unordered_map<std::string,size_t> * pCountObj =
		    (VirusNameSet.count(label.chr)) ? &viralCount : &hostCount;
		//Increment the host/virus reg count for the qname
		for(const std::string & qname : reg.QNameSet){
		    if(!pCountObj->count(qname)){
			(*pCountObj)[qname] = 0;
		    }
		    (*pCountObj)[qname]++;
		}
    	        it++;
    	    }
    	}
	//If no regions were removed, second pass is unecessary
	if(!bFilter && filterRound > 1) continue;
	fprintf(stderr,"\tAfter 1st Pass: %lu Regions remain\n",regionMap.size());
	//Second Pass - Remove qnames which now map to only host or 
	//If no qnames get removed the next iteration won't remove any regions
	bFilter=false; 
	for( auto & pair : regionMap){
	    for (   auto it = pair.second.QNameSet.begin();
		    it != pair.second.QNameSet.end(); )
	    {
		//If the qname is associated with both a host and virus
		//region, it may stay
		if(hostCount[*it] && viralCount[*it]){
		    it++;
		} else { // erase the non-junction qname
		    bFilter = true;
		    //Update the split read count for the region if
		    //necessary
		    bool bSplit = ((*it)[it->length()-2] == '_');
		    if(bSplit){
			//The weird syntax is to prevent underflow in a
			// case which shouldn't happen
			pair.second.nSplit += (pair.second.nSplit) ? -1 : 0;
		    }
		    it = pair.second.QNameSet.erase(it);
		}
	    }
	}
    } while(bFilter);

    fprintf(stderr,"Filtered Down to %lu Regions\n",regionMap.size());
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

