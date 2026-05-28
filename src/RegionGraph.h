#ifndef REGION_GRAPH_H
#define REGION_GRAPH_H

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
class CRegionGraph {
public:
    enum GRAPH_STATES {
        OWNS_GRAPH = 0x1,
        VALID_LOOKUP = 0x2,
    };
// Structure to cache vertex properties internally for quick lookup and manipulation
struct VertexProps {
    igraph_integer_t id;
    std::string chromosome;
    bool opensLeft;
    bool fromSplit;
    bool isHost;
    size_t left;
    size_t right;
    std::vector<std::string> assocFragGroups;
};

struct EdgeProps {
    igraph_integer_t id;
    double weight;
    bool fromSplit;
    std::vector<std::string> assocFragGroups;
};
    //Members
public:
    static const char dupDelim = 29;
    static const char fragDelim = 30;
protected:
    igraph_t graph;
    uint8_t flag;
    // Internal caches for O(log N) vertex uniqueness checks and property tracking
    std::multimap<std::string, igraph_integer_t> vertex_by_fragment;
    //std::vector<VertexProps> vertices;
    
//Con-/Destruction
public:
    CRegionGraph();
    ~CRegionGraph() { if(flag & OWNS_GRAPH) {igraph_destroy(&graph); } }
    // Delete copy semantics to prevent double-freeing the underlying igraph_t resource
    CRegionGraph(const CRegionGraph&) = delete;
    CRegionGraph& operator=(const CRegionGraph&) = delete;
    CRegionGraph(CRegionGraph&& other);
    CRegionGraph& operator=(CRegionGraph&& other);
//Accessors
public:
    igraph_t* get_igraph() { assertOwnership(); return &graph; }
    EdgeProps get_edge_properties(int id) const;
    VertexProps get_vertex_properties(int id) const;
    int vcount() const { return igraph_vcount(&graph); }
    int ecount() const { return igraph_ecount(&graph); }
//Methods:
public:
    void addOrUpdateVertex( std::string chromosome, bool opensLeft,
                            bool fromSplit, bool IsHost, size_t left,
                            size_t right, std::vector<std::string>);
    void addOrUpdateVertex(const VertexProps & prop);
    void filterEdges( double minWeight, double splitBonus);
private:
    void assertOwnership();
    void checkAndCreateEdge(igraph_integer_t v1_id, igraph_integer_t v2_id);
    void ensureValidLookup();
    void init_attribute_table();
};

//DEFINITIONS

//Constructor
CRegionGraph::CRegionGraph() : flag(OWNS_GRAPH | VALID_LOOKUP) {
    init_attribute_table();
    // Initialize an undirected graph
    if (igraph_empty(&graph, 0, IGRAPH_UNDIRECTED) != IGRAPH_SUCCESS) {
        throw std::runtime_error("Failed to initialize igraph object.");
    }
    // Set graph-level attributes
    // There aren't any
}

//Move Constructor
CRegionGraph::CRegionGraph(CRegionGraph&& other) :
    graph(std::move(other.graph)),
    flag(other.flag),
    vertex_by_fragment(std::move(other.vertex_by_fragment))
{
    other.flag &= ~OWNS_GRAPH;
}

//Move Assignment Operator
CRegionGraph& CRegionGraph::operator=(CRegionGraph&& other) {
    if(this == &other) { return *this; }
    if(flag & OWNS_GRAPH) { igraph_destroy(&graph); }
    graph = std::move(other.graph);
    flag = other.flag;
    vertex_by_fragment = std::move(other.vertex_by_fragment);
    other.flag &= ~OWNS_GRAPH;
    return *this;
}


//Checks that this object owns its underlying graph and can make changes
void CRegionGraph::assertOwnership() { 
        if(!(flag & OWNS_GRAPH)) {
            throw std::logic_error("Attempt to call non-const function from moved graph");
        }
    }

/**
     * Adds a vertex if the (proximalPos, distalPos) pair is unique.
     * If it already exists, merges attributes with the existing vertex.
     */
void CRegionGraph::addOrUpdateVertex(   std::string chromosome, bool opensLeft,
                                        bool fromSplit, bool isHost,
                                        size_t left, size_t right,
                                        std::vector<std::string> assocFragGroups)
{
    assertOwnership();
    ensureValidLookup();

    igraph_integer_t new_vid = igraph_vcount(&graph);
    igraph_add_vertices(&graph, 1, nullptr);

    //Get the set of old verticies which also have these fragments
    std::set<igraph_int_t> vertexWithFragmentSet;
    for(auto & fragGrp : assocFragGroups) {
        for(std::string frag : strsplit(fragGrp,dupDelim)) {
            auto range = vertex_by_fragment.equal_range(frag);
            for(auto i = range.first; i != range.second; i++){
                vertexWithFragmentSet.insert(i->second);
            }
            //Update the lookup with the new verticies
            vertex_by_fragment.insert({frag,new_vid});
        }
    }

    std::string assocFragStr = assocFragGroups.front();
    for(size_t i = 1; i <= assocFragGroups.size(); i++){
        assocFragStr += "\t" + assocFragGroups[i];
    }

    // Set underlying igraph C attributes
    SETVAS(&graph, "Chromosome", new_vid, chromosome.c_str());
    SETVAB(&graph, "OpensLeft", new_vid, opensLeft);
    SETVAB(&graph, "FromSplit", new_vid, fromSplit);
    SETVAB(&graph, "IsHost", new_vid, isHost);
    SETVAN(&graph, "Left", new_vid, left);
    SETVAN(&graph, "Right", new_vid, right);
    SETVAS(&graph, "assocFragGrps", new_vid, assocFragStr.c_str());

    // Check against relevant pre-existing vertices to evaluate edge creations
    for (igraph_int_t old_vid : vertexWithFragmentSet){
        checkAndCreateEdge(old_vid, new_vid);
    }
}


void CRegionGraph::addOrUpdateVertex(const VertexProps & prop) {
    this->addOrUpdateVertex(prop.chromosome, prop.opensLeft, prop.fromSplit,
                            prop.isHost, prop.left, prop.right,
                            prop.assocFragGroups);
}


// Evaluates window overlaps and constructs a weighted edge if conditions match
void CRegionGraph::checkAndCreateEdge(  igraph_integer_t v1_id,
                                    igraph_integer_t v2_id)
{
    //No self edges
    if(v1_id == v2_id) { return; }
    VertexProps v1 = this->get_vertex_properties(v1_id);
    VertexProps v2 = this->get_vertex_properties(v2_id);

    //Skip edges within parts of the graph
    if(v1.isHost == v2.isHost) { return; }

    //Determine if the edge already exists
    //Adding edges always follows adding a vertex so edges cannot already exist 
    //Also edges are only checked between verticies which share fragments
    //So as long as they are not in the same graph part, we are a go
    igraph_integer_t new_eid = igraph_ecount(&graph);
    igraph_add_edge(&graph, v1_id, v2_id);

    std::set<std::string> theIntersect;
    std::set_intersection(  v1.assocFragGroups.begin(),v2.assocFragGroups.end(),
                            v2.assocFragGroups.begin(),v2.assocFragGroups.end(),
                            std::inserter(theIntersect,theIntersect.end()));
    std::string assocFrag = *theIntersect.begin();
    for(auto it = theIntersect.begin(); it != theIntersect.end(); it++){
        if( it == theIntersect.begin()) { continue; }
        assocFrag += fragDelim + *it;
    }
    SETEAN(&graph, "weight", new_eid, theIntersect.size());
    SETEAN(&graph, "fromSplit", new_eid, v1.fromSplit || v2.fromSplit);
    SETEAS(&graph, "assocFragGrp", new_eid, assocFrag.c_str());
}

//Checks if the unique vertex lookup is valid, and rebuilds it if not
void CRegionGraph::ensureValidLookup() {
    if(flag & VALID_LOOKUP) { return; }
    vertex_by_fragment.clear();
    for(int i = 0; i < this->vcount(); i++){
        for(std::string fragGrp : strsplit(VAS(&graph,"assocFragGrps",i),fragDelim)){
            for( std::string fragName : strsplit(fragGrp,dupDelim)){
                vertex_by_fragment.insert({fragName,i});
            }
        }
    }
    flag |= VALID_LOOKUP;
}

//Removes edges which do not have sufficient support
void CRegionGraph::filterEdges( double minWeight, double splitBonus) {
    std::vector<igraph_int_t> toFilter;
    for(igraph_int_t id = 0; id < this->ecount(); id++){
        double weight = EAN(&graph,"weight",id);
        if(EAB(&graph,"FromSplit",id)){ weight +=1; }
        if(weight < minWeight){
            toFilter.push_back(id);
        }
    }
    if(toFilter.size()){
        flag &= ~VALID_LOOKUP;
        igraph_vector_int_t vec;
        igraph_vector_int_init(&vec,toFilter.size());
        for(size_t i = 0; i < toFilter.size();i++){
            VECTOR(vec)[i] = toFilter[i];
        }
        igraph_delete_edges(&graph, igraph_ess_vector(&vec));
        igraph_vector_int_destroy(&vec);
    }
}

////Removes vertices which have insufficient support
////  the degree of the vertex is below some threshold, even 
////  when providing a bonus for being a split read
//void CRegionGraph::filterVertices( double minSupport, double splitBonus) {
//    assertOwnership();
//    std::vector<igraph_int_t> toFilter;
//    for(igraph_int_t i = 0; i < this->vcount(); i++){
//        igraph_int_t degree;
//        igraph_degree_1(&graph,&degree,i,IGRAPH_ALL,igraph_loops_t(false));
//        double support =    1.0 + degree +
//                            (VAB(&graph,"IsSplit",i) ? splitBonus : 0);
//        if(support < minSupport) {
//
//            toFilter.push_back(i);
//        }
//    }
//    if(toFilter.size()){
//        flag &= ~VALID_LOOKUP;
//        igraph_vector_int_t vec;
//        igraph_vector_int_init(&vec,toFilter.size());
//        for(size_t i = 0; i < toFilter.size();i++) {
//            VECTOR(vec)[i] = toFilter[i];
//        }
//        igraph_delete_vertices(&graph, igraph_vss_vector(&vec));
//        igraph_vector_int_destroy(&vec);
//    }
//}


// Initializes the igraph C attribute table (must be called once before using attributes)
void CRegionGraph::init_attribute_table() {
    static bool initialized = false;
    if (!initialized) {
        igraph_set_attribute_table(&igraph_cattribute_table);
        initialized = true;
    }
}

CRegionGraph::EdgeProps CRegionGraph::get_edge_properties(int id) const {
    return {    id,
                EAN(&graph,"weight",id),
                EAB(&graph,"FromSplit",id),
                strsplit(EAS(&graph,"assocFragGrps",id),fragDelim)
    };
}

CRegionGraph::VertexProps CRegionGraph::get_vertex_properties(int id) const {
    return {    id,
                VAS(&graph,"Chromosome",id),
                VAB(&graph,"OpensLeft",id),
                VAB(&graph,"FromSplit",id),
                VAB(&graph,"IsHost",id),
                size_t(std::lround(VAN(&graph,"Left",id))),
                size_t(std::lround(VAN(&graph,"Right",id))),
                strsplit(VAS(&graph,"assocFragGrps",id),fragDelim)
    };
}



#endif // REGION_GRAPH_H
