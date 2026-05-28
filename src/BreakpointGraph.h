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
#include <stdexcept>
#include <set>
#include <string>
#include "str_utils.h"
#include <vector>


//Note: Current implementation stores all graph and vertex attributes twice
//  in the object, and in the underyling graph
class CBPGraph {
public:
    enum GRAPH_STATES {
        OWNS_GRAPH = 0x1,
        VALID_LOOKUP = 0x2,
    };
// Structure to cache vertex properties internally for quick lookup and manipulation
struct VertexProps {
    igraph_integer_t id;
    int proximalPos;
    int distalPos;
    bool isSplit;
    std::string assocFragments;
};
    //Members
public:
    static const char dupDelim = 29;
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
                            const std::string & assocFragments);
    void filterVertices( double minDegree, double splitBonus);
    std::vector<std::set<int>> maximalCliques(  double minVertex,
                                                double splitBonus);
private:
    void assertOwnership();
    void BronKerbosh2 ( std::set<igraph_int_t> R,
                        std::set<igraph_int_t> P,
                        std::set<igraph_int_t> X,
                        std::vector<std::set<int>> & res, int callDepth);
    void checkAndCreateEdge(igraph_integer_t v1_id, igraph_integer_t v2_id);
    void ensureValidLookup();
    void getWindow(int proxPos, bool isSplit, double& start, double& end) const;
    void init_attribute_table();
};

//DEFINITIONS

//Constructor
CBPGraph::CBPGraph(const std::string& chrom, bool opensLeftVal, 
             int upsDist, int rLen, int maxInsert, double sFactor)
        : flag(OWNS_GRAPH | VALID_LOOKUP),
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
void CBPGraph::assertOwnership() { 
        if(!(flag & OWNS_GRAPH)) {
            throw std::logic_error("Attempt to call non-const function from moved graph");
        }
    }

/**
     * Adds a vertex if the (proximalPos, distalPos) pair is unique.
     * If it already exists, merges attributes with the existing vertex.
     */
void CBPGraph::addOrUpdateVertex(   int proximalPos, int distalPos, bool isSplit,
                                    const std::string& assocFragment)
{
    assertOwnership();
    auto key = std::make_pair(proximalPos, distalPos);
    ensureValidLookup();
    auto it = vertex_lookup.find(key);

    if (it != vertex_lookup.end()) {
        // 1. Vertex pair already exists: Merge data into the existing vertex
        igraph_integer_t vid = it->second;
        
        // The combined IsSplit is true if either is true
        VertexProps props = this->get_vertex_properties(vid);
        props.isSplit = props.isSplit || isSplit;
        
        //Construct a string representing the sorted, unique fragment names associated
        //with this vertex
        std::vector<std::string> fragNames = strsplit(props.assocFragments,dupDelim);
        std::set<std::string> fragNameSet;
        fragNameSet.insert(fragNames.begin(),fragNames.end());
        fragNameSet.insert(assocFragment);
        std::string fragStr = *fragNameSet.begin();
        for(auto it = fragNameSet.begin(); it != fragNameSet.end(); it++){
            if(it == fragNameSet.begin()) { continue; }
            fragStr += dupDelim + *it; 
        }
        props.assocFragments = fragStr;

        // Update underlying igraph C attributes
        SETVAB(&graph, "IsSplit", vid, props.isSplit);
        SETVAS(&graph, "assocFragments", vid, props.assocFragments.c_str());
    } 
    else {
        // 2. Vertex pair is unique: Create a brand new vertex
        igraph_integer_t new_vid = igraph_vcount(&graph);
        igraph_add_vertices(&graph, 1, nullptr);

        vertex_lookup[key] = new_vid;

        // Set underlying igraph C attributes
        SETVAN(&graph, "ProximalPos", new_vid, proximalPos);
        SETVAN(&graph, "DistalPos", new_vid, distalPos);
        SETVAB(&graph, "IsSplit", new_vid, isSplit);
        SETVAS(&graph, "assocFragments", new_vid, assocFragment.c_str());

        // Check against all pre-existing vertices to evaluate edge creations
        for (igraph_integer_t old_vid = 0; old_vid < new_vid; ++old_vid) {
            checkAndCreateEdge(old_vid, new_vid);
        }
    }
}


// Evaluates window overlaps and constructs a weighted edge if conditions match
void CBPGraph::checkAndCreateEdge(  igraph_integer_t v1_id,
                                    igraph_integer_t v2_id)
{
    VertexProps v1 = this->get_vertex_properties(v1_id);
    VertexProps v2 = this->get_vertex_properties(v2_id);

    double s1, e1, s2, e2;
    getWindow(v1.proximalPos, v1.isSplit, s1, e1);
    getWindow(v2.proximalPos, v2.isSplit, s2, e2);

    // Compute the overlapping region
    double s_max = std::max(s1, s2);
    double e_min = std::min(e1, e2);

    if (s_max <= e_min) { // Windows overlap
        double overlap = e_min - s_max;
        double len1 = e1 - s1;
        double len2 = e2 - s2;

        if (len1 > 0 && len2 > 0) {
            double prop1 = overlap / len1;
            double prop2 = overlap / len2;

            double weight = prop1 + prop2;
            
            // Apply SplitFactor multiplier sequentially for each vertex that has isSplit == true
            if (v1.isSplit) weight *= splitFactor;
            if (v2.isSplit) weight *= splitFactor;

            // Add edge to igraph topology
            igraph_add_edge(&graph, v1_id, v2_id);
            igraph_integer_t new_eid = igraph_ecount(&graph) - 1;

            // Set edge weight attribute
            SETEAN(&graph, "weight", new_eid, weight);
        }
    }
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
    static bool initialized = false;
    if (!initialized) {
        igraph_set_attribute_table(&igraph_cattribute_table);
        initialized = true;
    }
}

void CBPGraph::BronKerbosh2 (   std::set<igraph_int_t> R,
                                std::set<igraph_int_t> P,
                                std::set<igraph_int_t> X,
                                std::vector<std::set<int>> & res, int callDepth)
{
    std::cerr << callDepth << ") " << R.size() << " " << P.size() << " " << X.size() << "\t" << res.size() << "\n";
    //If there are no more candidate nodes to add
    //this clique is maximal
    if(P.size() + X.size() == 0) {
        std::set<int> clique;
        clique.insert(R.begin(),R.end());
        res.push_back(clique);
    }
    if(!P.size()) { return; }
    //Select a pivot //TODO Do it smarter
    igraph_int_t pivot = *P.begin();
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
    std::cerr << "\t" << P.size() << "\n";
    while(PnonPivot.size()){
        igraph_int_t v = *P.begin();
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
        BronKerbosh2(Rprime,Pprime,Xprime,res,callDepth+1);
        //remove v from P
        P.erase(v);
        PnonPivot.erase(v);
        //add v to x
        X.insert(v);
    }
}


std::vector<std::set<int>> CBPGraph::maximalCliques(    double minVertex,
                                                        double splitBonus) {
    std::vector<std::set<int>> cliques;
    std::set<igraph_int_t> nodeIdx; 
    for(igraph_int_t i = 0; i < this->vcount(); i++){
        nodeIdx.insert(i);
    }
    this->BronKerbosh2({},nodeIdx,{},cliques,0);
    for(auto it = cliques.begin(); it != cliques.end(); it++){
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
    return cliques;
}


CBPGraph::VertexProps CBPGraph::get_vertex_properties(int id) const {
    return {    id,
                int(std::lround(VAN(&graph,"ProximalPos",id))),
                int(std::lround(VAN(&graph,"DistalPos",id))),
                VAB(&graph,"IsSplit",id),
                VAS(&graph,"assocFragments",id)
    };
}

#endif // BREAKPOINT_GRAPH_H
