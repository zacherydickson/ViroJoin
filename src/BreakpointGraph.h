#ifndef BREAKPOINT_GRAPH_H
#define BREAKPOINT_GRAPH_H

//This code was originally generated with Google Gemini 3.5 Flash (Extended)
// on May 27th 2026, with the following prompt:
// =============================================================================
/*
I am using the igraph library in my c++ project.
I would like to write a class `CBPGraph` using igraph functionality to handle the construction and destruction of the graph. The graph should ultimately be compatible with the libleidenalg library.

I would like the graph to have a string attribute "Chromosome", and a boolean attribute "OpensLeft"
Each vertex should have the following attributes: int "ProximalPos", int "DistalPos", bool "IsSplit"; string "assocFragments".
Each edge should have the following attributes: Numeric weight.

I will also want some extra functionality, to check if the combo of distal and proximal positions of a new vertex is unique amongst verticies. If so a new vertex can be added, otherwise the assocFragments information should be added to the existing vertex. The combined IsSplit attribute is true if either vertex is true.

When a vertex is added, each pre-existing vertex should be checked to see if an edge should be created.

Whether an edge is created between verticies and the weight depends on the graph attribute OpensLeft, The vertex attributes IsSplit, and the proximal positions.
Consider an asymmetrical window around the proximal position. An edge is formed if the windows for both verticies overlap. The OpensLeft attribute determines which side of the window is up- vs downstream. The upstream side of the window is smaller.  Upstream is numerically smaller if false, or numerically higher direction if true.  IsSplit determines the size of the window. If true then the window extends downstream a smaller ReadLen distance,otherwise the downstream extent is a larger MaxInsertSize distance. The upstream distance is always a fixed value UpstreamDist.

For each vertex we can determine the proportion of that vertex's window which is overlapped by the other window. The weight for the edge is the sum sum of these two proportions multiplied by some value SplitFactor for each of the vertexes with the is split attribute. As an example consider an Unsplit vertex upstream of a split vertex. 11% of the former's window is covered by the latter, and 90% of the latter's window is covered by the former. The Weight of the edge between these two would be (0.11 + 0.9)  * 2 = 2.02.

Could you write c++ code that matches this specification?
*/
//==============================================================================
//Since modified for style and readability


#include <algorithm>
#include "igraph/igraph.h"
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <set>
#include <string>
#include "str_utils.h"
#include <unordered_set>
#include <vector>


//Note: Current implementation stores all graph and vertex attributes twice
//  in the object, and in the underyling graph
class CBPGraph {
public:
    enum GRAPH_STATES {
        OWNS_GRAPH = 0x1,
        VALID_LOOKUP = 0x2,
        VALID_EDGES = 0x4,
        VALID_CLIQUES = 0x8,
    };
// Structure to cache vertex properties internally for quick lookup and manipulation
struct VertexProps {
    igraph_integer_t id;
    int proximalPos;
    int distalPos;
    bool isSplit;
    std::vector<int> cliques;
    std::vector<std::string> assocFragments;
    std::string to_string() {
        return  std::to_string(id) + ") " +
                std::to_string(distalPos) + "-" + std::to_string(proximalPos) + 
                std::string((isSplit) ? "|" : ">") + "\t" +
                to_strjoin(cliques.begin(),cliques.end(),',') + "\t" +
                strjoin(assocFragments.begin(),assocFragments.end(),',');
    }
};
    //Members
public:
    static const char DupDelim = 29;
protected:
    igraph_t graph;
    uint8_t flag;
    // Internal caches for O(log N) vertex uniqueness checks and property tracking
    std::map<std::pair<int, int>, igraph_integer_t> vertex_lookup;
    //std::vector<VertexProps> vertices;
    // Window and weight configuration parameters
    int upstreamDist;
    int readLen;
    int maxInsertSize;
    double splitFactor;
    
//Con-/Destruction
public:
    CBPGraph() = delete;
    CBPGraph(const std::string& chrom, bool opensLeftVal, 
             int upsDist, int rLen, int maxInsert, double sFactor);
    CBPGraph(   igraph_t && graph, const std::string & chrom, bool opensLeftVal,
                int upsDist, int rLen, int maxInsert, int sFactor,
                bool hasEdges, bool hasCliques);
    CBPGraph(igraph_t && graph, const CBPGraph & parent);
    ~CBPGraph() { if(flag & OWNS_GRAPH) {igraph_destroy(&graph); } }
    // Delete copy semantics to prevent double-freeing the underlying igraph_t resource
    CBPGraph(const CBPGraph&) = delete;
    CBPGraph& operator=(const CBPGraph&) = delete;
    //Move semantics
    CBPGraph(CBPGraph&& other);
    CBPGraph& operator=(CBPGraph&& other);
//Accessors
public:
    std::string get_chromosome() const { return GAS(&graph,"Chromosome"); }
    igraph_t* get_igraph() { assertOwnership(); return &graph; }
    bool opens_left() const {return GAB(&graph,"OpensLeft");}
    VertexProps get_vertex_properties(int id) const;
    int vcount() const { return igraph_vcount(&graph); }
    int ecount() const { return igraph_ecount(&graph); }
//Methods:
public:
    void addOrUpdateVertex( int proximalPos, int distalPos, bool isSplit,
                            const std::string & assocFragments,
                            bool bOnline = true);
    std::vector<CBPGraph> decompose(int minVertex);
    void ensureConstructed();
    void filterVertices( double minDegree, double splitBonus);
    bool maximalCliques(  double minVertex, double splitBonus);
    void write_edgelist(FILE * outstream) const {
        igraph_write_graph_edgelist(&graph,outstream);
    }
private:
    void assertOwnership() const;
    void BronKerbosh2 ( std::set<igraph_int_t> R,
                        std::set<igraph_int_t> P,
                        std::set<igraph_int_t> X,
                        std::vector<std::set<igraph_int_t>> & res) const;
    bool checkAndCreateEdge(igraph_integer_t v1_id, igraph_integer_t v2_id);
    bool checkAndCreateEdge(igraph_integer_t v1_id, igraph_integer_t v2_id,
                            double s1, double e1, double s2, double e2);
    void constructEdges();
    void ensureValidLookup();
    static bool fragsets_are_comparable(    std::vector<std::string> fragVec1,
                                            std::vector<std::string> fragVec2);
    void getWindow(int proxPos, bool isSplit, double& start, double& end) const;
    void init_attribute_table();
    void removeSharedFragEdges(std::string frag, igraph_int_t vid);
    igraph_int_t selectPivot(   const std::set<igraph_int_t> & P,
                                std::set<igraph_int_t> &symDiff) const;
    bool vertexesHaveIndependentSupport(    igraph_int_t vid1,
                                            igraph_int_t vid2 ) const 
    {
        return  !CBPGraph::fragsets_are_comparable(
                    this->get_vertex_properties(vid1).assocFragments,
                    this->get_vertex_properties(vid1).assocFragments
                );
    }
    void weightEdge(    igraph_int_t eid, igraph_int_t v1_id, igraph_int_t v2_id,
                        double s1, double e1, double s2,
                        double e2, double eMin, double sMax);
};

//DEFINITIONS

//Constructor
CBPGraph::CBPGraph(const std::string& chrom, bool opensLeftVal, 
             int upsDist, int rLen, int maxInsert, double sFactor)
        : flag(OWNS_GRAPH | VALID_LOOKUP | VALID_EDGES),
          upstreamDist(upsDist), readLen(rLen), maxInsertSize(maxInsert),
          splitFactor(sFactor)
{
    init_attribute_table();
    // Initialize an undirected graph
    if (igraph_empty(&graph, 0, IGRAPH_UNDIRECTED) != IGRAPH_SUCCESS) {
        throw std::runtime_error("Failed to initialize igraph object.");
    }
    // Set graph-level attributes
    SETGAS(&graph, "Chromosome", chrom.c_str());
    SETGAB(&graph, "OpensLeft", opensLeftVal);
}


//Pre-constructed graph constructor
CBPGraph::CBPGraph( igraph_t && graph, const std::string & chrom,
                    bool opensLeftVal, int upsDist, int rLen, int maxInsert,
                    int sFactor, bool hasEdges, bool hasCliques) :
    graph(std::move(graph)), flag(OWNS_GRAPH), upstreamDist(upsDist),
    readLen(rLen), maxInsertSize(maxInsert),splitFactor(sFactor)
{
    if(hasEdges) { flag |= VALID_EDGES; }
    if(hasCliques) { flag |= VALID_CLIQUES; }
}

//Child graph constructor
CBPGraph::CBPGraph(igraph_t && graph, const CBPGraph & parent) :
    CBPGraph(   std::move(graph), parent.get_chromosome(), parent.opens_left(),
                parent.upstreamDist, parent.readLen, parent.maxInsertSize,
                parent.splitFactor, parent.flag & VALID_EDGES,
                parent.flag & VALID_CLIQUES)
{
}

//Move Constructor
CBPGraph::CBPGraph(CBPGraph&& other) :
    graph(std::move(other.graph)),
    flag(other.flag),
    vertex_lookup(std::move(other.vertex_lookup)),
    upstreamDist(other.upstreamDist),
    readLen(other.readLen),
    maxInsertSize(other.maxInsertSize),
    splitFactor(other.splitFactor)
{
    other.flag &= ~OWNS_GRAPH;
}

//Move Assignment Operator
CBPGraph& CBPGraph::operator=(CBPGraph&& other) {
    if(this == &other) { return *this; }
    if(flag & OWNS_GRAPH) { igraph_destroy(&graph); }
    graph = std::move(other.graph);
    flag = other.flag;
    vertex_lookup = std::move(other.vertex_lookup);
    upstreamDist = other.upstreamDist;
    readLen = other.readLen;
    maxInsertSize = other.maxInsertSize;
    splitFactor = other.splitFactor;
    other.flag &= ~OWNS_GRAPH;
    return *this;
}


//Checks that this object owns its underlying graph and can make changes
void CBPGraph::assertOwnership() const { 
    if(!(flag & OWNS_GRAPH)) {
        throw std::logic_error("Attempt to call non-const function from moved graph");
    }
}


/**
     * Adds a vertex if the (proximalPos, distalPos) pair is unique.
     * If it already exists, merges attributes with the existing vertex.
     */
void CBPGraph::addOrUpdateVertex(   int proximalPos, int distalPos, bool isSplit,
                                    const std::string& assocFragment, bool bOnline)
{
    //Skip empty vertexes
    if(proximalPos == distalPos) {return;}
    assertOwnership();
    auto key = std::make_pair(proximalPos, distalPos);
    ensureValidLookup();
    auto it = vertex_lookup.find(key);

    igraph_int_t vid = igraph_vcount(&graph);

    std::vector<igraph_int_t> nonAdjVertices;
    if (it != vertex_lookup.end()) {
        // 1. Vertex pair already exists: Merge data into the existing vertex
        vid = it->second;
        
        // The combined IsSplit is true if either is true
        VertexProps props = this->get_vertex_properties(vid);
        props.isSplit = props.isSplit || isSplit;
        
        //Construct a string representing the sorted, unique fragment names associated
        //with this vertex
        std::set<std::string> fragNameSet;
        fragNameSet.insert(props.assocFragments.begin(),props.assocFragments.end());
        fragNameSet.insert(assocFragment);
        std::string fragStr = strjoin(  fragNameSet.begin(),
                                        fragNameSet.end(), DupDelim);

        // Update underlying igraph C attributes
        SETVAB(&graph, "IsSplit", vid, props.isSplit);
        SETVAS(&graph, "assocFragments", vid, fragStr.c_str());

        //Collect all verticies currently not adjacent to this one
        if(bOnline){ //Only if we will be adding edges immediately
            igraph_vs_t vs;
            igraph_vs_nonadj(&vs,vid,IGRAPH_ALL);
            igraph_vit_t vit;
            igraph_vit_create(&graph,vs,&vit);
            while(!IGRAPH_VIT_END(vit)){
                nonAdjVertices.push_back(IGRAPH_VIT_GET(vit));
                IGRAPH_VIT_NEXT(vit);
            }
            igraph_vit_destroy(&vit);
            igraph_vs_destroy(&vs);
        }

        //Remove any edges attached to this vertex which connect
        // to a vertex which now no longer has idependent support
        // A vs B (independent) -> A vs AB (not independent)
        //this->removeSharedFragEdges(assocFragment,vid);
    } 
    else {
        // 2. Vertex pair is unique: Create a brand new vertex
        igraph_add_vertices(&graph, 1, nullptr);
        flag &= ~VALID_CLIQUES;

        vertex_lookup[key] = vid;

        // Set underlying igraph C attributes
        SETVAN(&graph, "ProximalPos", vid, proximalPos);
        SETVAN(&graph, "DistalPos", vid, distalPos);
        SETVAB(&graph, "IsSplit", vid, isSplit);
        SETVAS(&graph, "cliques", vid, "0");
        SETVAS(&graph, "assocFragments", vid, assocFragment.c_str());

        if(bOnline){ //Only if we will be adding edges immediately
            nonAdjVertices.resize(vid);
            std::iota(nonAdjVertices.begin(),nonAdjVertices.end(),0);
        }
    }

    if(!bOnline){
        flag &= ~VALID_EDGES;
    }
    //nonAdjVerticies should be empty if offline
    if(!nonAdjVertices.size()) { return; }

    //Regardless of whether a new vertex was added, new edges may be formed
    //  OR
    // (AB vs A (no-independent support) becomes AB vs AC (independent support)
    // Check against all pre-existing non-adjacent vertices to
    // evaluate edge creation
    for (igraph_integer_t old_vid : nonAdjVertices) {
        checkAndCreateEdge(old_vid, vid);
    }

}

void CBPGraph::BronKerbosh2 (   std::set<igraph_int_t> R,
                                std::set<igraph_int_t> P,
                                std::set<igraph_int_t> X,
                                std::vector<std::set<igraph_int_t>> & res ) const 
{
    //If there are no more candidate nodes to add
    //this clique is maximal
    if(P.size() + X.size() == 0) {
        std::set<igraph_int_t> clique;
        clique.insert(R.begin(),R.end());
        res.push_back(clique);
    }
    if(!P.size()) { return; }
    std::set<igraph_int_t> Q;
    //The pivot index returned is discarded
    /*igraph_int_t pivot =*/
    this->selectPivot(P,Q);
    for( igraph_int_t v : Q) {
        //Identify neighbours (N) of the vertex
        igraph_vs_t vs; // The concept of picking vertices in a graph
        igraph_vit_t vit; // The selection of verteces in this graph
        igraph_vs_adj(&vs,v,IGRAPH_ALL,IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
        igraph_vit_create(&graph, vs, &vit);
        std::set<igraph_int_t> N;
        while(!IGRAPH_VIT_END(vit)) {
            N.insert(IGRAPH_VIT_GET(vit));
            IGRAPH_VIT_NEXT(vit);
        }
        igraph_vit_destroy(&vit);
        igraph_vs_destroy(&vs);
        //Get the updated sets as R + v, Intersect(P,N) and Intersect (X,N)
        std::set<igraph_int_t> Rprime = R;
        Rprime.insert(v);
        std::set<igraph_int_t> Pprime;
        std::set<igraph_int_t> Xprime;
        std::set_intersection(  P.begin(),P.end(),
                                N.begin(),N.end(),
                                std::inserter(Pprime,Pprime.end()));
        std::set_intersection(  X.begin(),X.end(),
                                N.begin(),N.end(),
                                std::inserter(Xprime,Xprime.end()));
        //Make the recursive call
        BronKerbosh2(Rprime,Pprime,Xprime,res);
        //remove v from P
        P.erase(v);
        //add v to x
        X.insert(v);
    }
}


// Evaluates window overlaps and constructs a weighted edge if conditions match
// returns true if an edge was created, false otherwise
bool CBPGraph::checkAndCreateEdge(  igraph_integer_t v1_id,
                                    igraph_integer_t v2_id)
{
    double s1, e1, s2, e2;
    getWindow(VAN(&graph,"ProximalPos",v1_id), VAB(&graph,"IsSplit",v1_id), s1, e1);
    getWindow(VAN(&graph,"ProximalPos",v2_id), VAB(&graph,"IsSplit",v2_id), s2, e2);
    return this->checkAndCreateEdge(v1_id,v2_id,s1,e1,s2,e2);
}

// Used pre-evaluated windows to construct a weighted edge if conditions match
// returns true if an edge was created, false otherwise
bool CBPGraph::checkAndCreateEdge(  igraph_integer_t v1_id,
                                    igraph_integer_t v2_id,
                                    double s1, double e1, double s2, double e2)
{
    assertOwnership();
    //If the sopport for two separate breakpoint is a completely overlapping set
    //  of fragments (alt-mappings of the same fragment),
    //  then they cannot support the same
    //  breakpoint, and therefore no edge should be made
    //if(CBPGraph::fragsets_are_comparable(v1.assocFragments,v2.assocFragments)){
    //    return;
    //}
    //std::set<std::string> fragIntersect;
    //std::set_intersection(  v1.assocFragments.begin(),v1.assocFragments.end(),
    //                        v2.assocFragments.begin(),v2.assocFragments.end(),
    //                        std::inserter(fragIntersect,fragIntersect.end()) );
    //if(fragIntersect.size()){
    //    return;
    //}

    // Compute the overlapping region
    double s_max = std::max(s1, s2);
    double e_min = std::min(e1, e2);

    // Check if Windows do not overlap
    if (s_max > e_min) { return false; }

    igraph_add_edge(&graph, v1_id, v2_id);
    //NOTE: In the current implementation, edge weights are not used
    //weightEdge(igraph_ecount(&graph)-1,s1,e1,s2,e2,e_min,s_max);
    return true;
}


void CBPGraph::constructEdges() {
    assertOwnership();
    //Don't construct if the edges are already valid
    if(flag & VALID_EDGES){ return; }
    //Remove any pre-existing edges
    size_t eCount = this->ecount();
    if(eCount > 0){
        igraph_delete_edges(&graph,igraph_ess_all(IGRAPH_EDGEORDER_ID));
    }
    //There are no edges in an empty or singular graph as loops are forbidden
    if(this->vcount() < 2){
        flag |= VALID_EDGES;
        return;
    }
    
    //Sweep line implementation
    //  pass a line from the start of the first interval to the end of the last interval
    //  whenever the line hits the start of an interval, add that vertex to the
    //  active list.
    //  whenever the line hits the end of an interval, create an edge from the terminated vertex
    //  to all other active vertexes
    struct event_t {
        igraph_int_t vid;
        bool isEnd;
        double pos;
    };
    bool bOpensLeft = GAB(&graph,"OpensLeft");
    //std::cerr << "Pre Event construct\n";
    std::vector<event_t> eventVec;
    eventVec.reserve(this->vcount());
    for(igraph_int_t vid = 0; vid < this->vcount(); vid++){
        std::pair<double,double> window;
        this->getWindow(VAN(&graph,"ProximalPos",vid),VAB(&graph,"IsSplit",vid),
                            window.first,window.second);
        if(bOpensLeft){ //Put the window in the order upstream, downstream
            std::swap(window.first, window.second);
            window.first *= -1;
            window.second *= -1;
        }
        eventVec.push_back({vid,false,window.first});
        eventVec.push_back({vid,true,window.second});
    }
    //Sort the events from most upstream to most downstream
    std::sort(eventVec.begin(),eventVec.end(),
            [](const event_t & a, const event_t & b){
                return a.pos < b.pos;
            } );
    std::unordered_set<igraph_int_t> activeVertexSet;
    //size_t counter = 0;
    //size_t updateAt = 1000;
    //size_t totalEdges = 0;
    std::vector<igraph_int_t> adjVec;
    //Initial guess at the number of edges to be created
    adjVec.reserve(eventVec.size());
    for(auto it = eventVec.begin(); it != eventVec.end(); ) {
        //if(counter > updateAt){
        //    std::cerr << counter << " of " << eventVec.size() << " events processed; " << totalEdges << "edges created so far\r"; 
        //    updateAt = counter + 1000;
        //}
        std::vector<igraph_int_t> closingVertexSet;
        //Find the first event at a position higher than this one
        auto nx = std::next(it);
        while(nx != eventVec.end() && (nx->pos == it->pos)) {
            nx++;
        }
        //Note all vertexes updated by events at this position
        for(; it != nx; it++){
            if(it->isEnd) {
                closingVertexSet.push_back(it->vid);
            } else {
                activeVertexSet.insert(it->vid);
            }
            //counter++;
        }
        //continue if there are no close events
        if(!closingVertexSet.size()) { continue; }
        //Construct an adjacency list for edges to add
        //size_t nEdges = (activeVertexSet.size() - 1) * closingVertexSet.size();
        //totalEdges += nEdges;
        for(igraph_int_t cVid : closingVertexSet){
            for(igraph_int_t aVid : activeVertexSet) {
                //Skip self edges
                if(cVid == aVid) { continue; }
                adjVec.push_back(cVid);
                adjVec.push_back(aVid);
            }
        }
        //Remove closing Vertices from the active set
        for(igraph_int_t vid : closingVertexSet){
            activeVertexSet.erase(vid);
        }
    }
    igraph_vector_int_t adjList;
    igraph_vector_int_init_array(&adjList, adjVec.data(), adjVec.size());
    igraph_add_edges(&graph,&adjList,nullptr);
    igraph_vector_int_destroy(&adjList);

    flag |= VALID_EDGES;
}

//Checks if the unique vertex lookup is valid, and rebuilds it if not
void CBPGraph::ensureValidLookup() {
    if(flag & VALID_LOOKUP) { return; }
    vertex_lookup.clear();
    for(int i = 0; i < this->vcount(); i++){
        int proxPos = std::lround(VAN(&graph,"ProximalPos",i));
        int distPos = std::lround(VAN(&graph,"DistalPos",i));
        vertex_lookup.emplace(std::make_pair(proxPos,distPos),i);
    }
    flag |= VALID_LOOKUP;
}

void CBPGraph::ensureConstructed() {
    if((flag & VALID_EDGES)){ return; }
    //std::cerr << "\t\tStart construction\n";
    this->constructEdges();
    //std::cerr << "\t\tEnd construction\n";
}

//Construct subgraphs of minimum size from each connected component of this graph
//Input - a threshold number of vertexes
//Output a vector of subgraphs, empty if no subgraphs have enough
//vertexes
std::vector<CBPGraph> CBPGraph::decompose(int minVertex) {
    ensureConstructed();
    std::vector<CBPGraph> subGraphVec;
    //Calculate the components
    igraph_graph_list_t components;
    igraph_graph_list_init(&components,0);
    igraph_decompose(&graph, &components, IGRAPH_WEAK, -1, minVertex);
    //Construct the children
    for(igraph_int_t i =0; i < igraph_graph_list_size(&components); i++){
        //The list owns its elements, so if we moved the items out
        //then destroy the list, it might try to free the children
        //so we'll copy instead
        igraph_t child;
        igraph_copy(&child,igraph_graph_list_get_ptr(&components,i));
        subGraphVec.emplace_back(std::move(child),*this);
    }
    igraph_graph_list_destroy(&components);
    return subGraphVec;
}

//Removes vertices which have insufficient support
//  the degree of the vertex is below some threshold, even 
//  when providing a bonus for being a split read
void CBPGraph::filterVertices( double minSupport, double splitBonus) {
    assertOwnership();
    std::vector<igraph_int_t> toFilter;
    for(igraph_int_t i = 0; i < this->vcount(); i++){
        igraph_int_t degree;
        igraph_degree_1(&graph,&degree,i,IGRAPH_ALL,igraph_loops_t(false));
        double support =    1.0 + degree +
                            (VAB(&graph,"IsSplit",i) ? splitBonus : 0);
        if(support < minSupport) {
            toFilter.push_back(i);
        }
    }
    if(toFilter.size()){
        flag &= ~VALID_LOOKUP;
        flag &= ~VALID_CLIQUES;
        igraph_vector_int_t vec;
        igraph_vector_int_init(&vec,toFilter.size());
        for(size_t i = 0; i < toFilter.size();i++) {
            VECTOR(vec)[i] = toFilter[i];
        }
        igraph_delete_vertices(&graph, igraph_vss_vector(&vec));
        igraph_vector_int_destroy(&vec);
    }
}


// Calculates the window boundaries based on direction (OpensLeft) and split status
void CBPGraph::getWindow(int proxPos, bool isSplit, double& start, double& end) const {
    if (!this->opens_left()) {
        // Upstream is numerically smaller (left), Downstream is numerically higher (right)
        start = proxPos - upstreamDist;
        end = proxPos + (isSplit ? readLen : maxInsertSize);
    } else {
        // Upstream is numerically higher (right), Downstream is numerically smaller (left)
        start = proxPos - (isSplit ? readLen : maxInsertSize);
        end = proxPos + upstreamDist;
    }
}

// Initializes the igraph C attribute table (must be called once before using attributes)
void CBPGraph::init_attribute_table() {
    assertOwnership();
    static bool initialized = false;
    if (!initialized) {
        igraph_set_attribute_table(&igraph_cattribute_table);
        initialized = true;
    }
}



//Note: it is assumed that the fragments vectors are sorted
bool CBPGraph::fragsets_are_comparable( std::vector<std::string> fragVec1,
                                        std::vector<std::string> fragVec2)
{ 
    //If both are empty they are comparable
    if(fragVec1.size() == fragVec2.size() && !fragVec1.size()) return true;
    //Determine the smaller vector
    std::vector<std::string> * smaller = &fragVec1;
    std::vector<std::string> * larger = &fragVec2;
    if(fragVec1.size() > fragVec2.size()){
        std::swap(fragVec1,fragVec2);
    }
    return std::includes(   larger->begin(),larger->end(),
                            smaller->begin(),smaller->end());
}

CBPGraph::VertexProps CBPGraph::get_vertex_properties(int id) const {
    if(id >= this->vcount()){
        throw std::invalid_argument("Attempt to get vertex properties for a vertex outside of the graph");
    }
    std::vector<std::string> cliqueStrs = strsplit(VAS(&graph,"cliques",id),DupDelim);
    std::vector<int> cliqueAssignVec;
    for(auto cliqueStr : cliqueStrs){
        try {
        cliqueAssignVec.push_back(std::stoi(cliqueStr));
        } catch (std::invalid_argument &e ) {
            throw std::invalid_argument(std::string(e.what()) + " " + cliqueStr);
        }
    }
    return {    id,
                int(std::lround(VAN(&graph,"ProximalPos",id))),
                int(std::lround(VAN(&graph,"DistalPos",id))),
                VAB(&graph,"IsSplit",id),
                cliqueAssignVec,
                strsplit(VAS(&graph,"assocFragments",id),DupDelim)
    };
}

//Note, edges do not exist between vertexes formed from alternate mappings of
//  the same fragment, therefore any fragment will appear within a clique, just
//  once
//Determines and stores internally all maximal cliques within the current graph
//Output    - true if there are cliques meeting the criteria, false otherwise
bool CBPGraph::maximalCliques( double minVertex, double splitBonus) {
    assertOwnership();
    ensureConstructed();
    std::vector<std::set<igraph_int_t>> cliques;
    std::set<igraph_int_t> nodeIdx; 
    for(igraph_int_t i = 0; i < this->vcount(); i++){
        nodeIdx.insert(i);
    }
    //this->BronKerbosh2({},nodeIdx,{},cliques);
    //The igraph implementation is at least 3x faster...
    igraph_vector_int_list_t cliq;
    igraph_vector_int_list_init(&cliq,0);
    //igraph_set_progress_handler(igraph_progress_handler_stderr);
    igraph_maximal_cliques(&graph,&cliq,int(minVertex - splitBonus),IGRAPH_UNLIMITED,IGRAPH_UNLIMITED);
    //igraph_set_progress_handler(nullptr);
    for(int i = 0; i < igraph_vector_int_list_size(&cliq); i++){
        cliques.push_back(std::set<igraph_int_t>());
        igraph_vector_int_t * ptr = igraph_vector_int_list_get_ptr(&cliq,i);
        for(int j = 0; j < igraph_vector_int_size(ptr); j++){
            cliques[i].insert(VECTOR(*ptr)[j]);
        }

    }
    igraph_vector_int_list_destroy(&cliq);
    //Filter Cliques which are too small
    for(auto it = cliques.begin(); it != cliques.end();){
        bool bSplit = false;
        for(igraph_int_t id : *it){
            if(VAB(&graph,"IsSplit",id)) {
                bSplit = true;
                break;
            }
        }
        if(it->size() + ((bSplit) ? splitBonus : 0.0) < minVertex) {
            it = cliques.erase(it);
        } else {
            it++;
        }
    }
    //Determine all cliques to which a given vertex belongs
    std::multimap<igraph_int_t,size_t> cliqueMap;
    for(size_t i = 0; i < cliques.size(); i++){
        for(igraph_int_t id : cliques[i]) {
            cliqueMap.insert({id,i});
        }
    }
    //Store the clique information
    for(igraph_int_t id = 0; id < this->vcount(); id++){
        //Set the clique to empty for any vertexes not in a valid clique
        if(!cliqueMap.count(id)) { 
            SETVAS(&graph,"cliques",id,"0");
            continue;
        }
        auto range = cliqueMap.equal_range(id);
        std::string str = to_strjoin(
                range.first, range.second, DupDelim,
                [](const std::pair<igraph_int_t,size_t> & item){
                    return std::to_string(item.second);
                    } );

        SETVAS(&graph,"cliques",id,str.c_str());
    }
    flag |= VALID_CLIQUES;
    return bool(cliques.size());
}

void CBPGraph::removeSharedFragEdges(std::string frag, igraph_int_t vid) {
    assertOwnership();
    //Edge Selector for all edges on this vertex
    igraph_es_t es;
    igraph_es_incident(&es,vid,IGRAPH_ALL,IGRAPH_NO_LOOPS);
    igraph_eit_t eit;
    igraph_eit_create(&graph,es, &eit);
    //Iterate over edges and find vertices sharing the fragment
    std::set<igraph_int_t> toRemove;
    while(!IGRAPH_EIT_END(eit)){
        igraph_int_t other_vid = IGRAPH_OTHER(&graph,IGRAPH_EIT_GET(eit),vid);
        if(!this->vertexesHaveIndependentSupport(vid,other_vid)) {
                toRemove.insert(IGRAPH_EIT_GET(eit));
        }
        IGRAPH_EIT_NEXT(eit);
    }
    //Remove noted edges if there are any
    if(toRemove.size()) {
        igraph_vector_int_t vec;
        igraph_vector_int_init(&vec,toRemove.size());
        int i = 0;
        for(igraph_int_t id : toRemove){
            VECTOR(vec)[i++] = id;
        }
        igraph_delete_edges(&graph,igraph_ess_vector(&vec));
        igraph_vector_int_destroy(&vec);
    }
    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);
}

igraph_int_t CBPGraph::selectPivot( const std::set<igraph_int_t> & P,
                                    std::set<igraph_int_t> &symDiff) const
{
    symDiff = P;
    igraph_int_t bestPivot = *P.begin();
    for(igraph_int_t pivot : P){
        //Identify neighbours (N) of the pivot, and remove them from P
        igraph_vs_t vs; // The concept of picking vertices in a graph
        igraph_vit_t vit; // The selection of verteces in this graph
        igraph_vs_adj(&vs,pivot,IGRAPH_ALL,IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
        igraph_vit_create(&graph, vs, &vit);
        std::set<igraph_int_t> PnonPivot = P;
        while(!IGRAPH_VIT_END(vit)) {
            PnonPivot.erase(IGRAPH_VIT_GET(vit));
            IGRAPH_VIT_NEXT(vit);
        }
        igraph_vit_destroy(&vit);
        igraph_vs_destroy(&vs);
        if(PnonPivot.size() < symDiff.size()){
            bestPivot = pivot;
            symDiff = std::move(PnonPivot);
        }
    }
    return bestPivot;
}


void CBPGraph::weightEdge(  igraph_int_t eid, igraph_int_t v1_id, igraph_int_t v2_id,
                            double s1, double e1, double s2,
                            double e2, double eMin, double sMax)
{
    double overlap = eMin - sMax;
    double len1 = e1 - s1;
    double len2 = e2 - s2;

    //Empty intervals are prevented in addOrUpdateVertex
    //  so len1 and len2 are always positive
    double prop1 = overlap / len1;
    double prop2 = overlap / len2;

    double weight = prop1 + prop2;
    
    // Apply SplitFactor multiplier sequentially for each vertex that has isSplit == true
    if (VAB(&graph,"IsSplit",v1_id)) weight *= splitFactor;
    if (VAB(&graph,"IsSplit",v2_id)) weight *= splitFactor;

    // Set edge weight attribute
    SETEAN(&graph, "weight", eid, weight);
}

#endif // BREAKPOINT_GRAPH_H
