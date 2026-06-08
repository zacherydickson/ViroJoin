#ifndef REGION_GRAPH_H
#define REGION_GRAPH_H

#include <algorithm>
#include "igraph/igraph.h"
#include <iostream>
#include <map>
#include <numeric>
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
    std::string to_string() {
        return  std::to_string(id) + ") " + 
                std::string((isHost) ? "H" : "V") + ":" + chromosome + "\t" +
                std::string((opensLeft) ? ((fromSplit) ? "|" : "<") : "") +
                std::to_string(left) + "-" + std::to_string(right) + 
                std::string((!opensLeft) ? ((fromSplit) ? "|" : ">") : "") +
                "\t" + 
                strjoin(assocFragGroups.begin(),assocFragGroups.end(),',');
    }
    bool testOverlap(const VertexProps & other) const {
        if(chromosome != other.chromosome) { return false; }
        if(opensLeft != other.opensLeft) { return false; }
        size_t ml = (left < other.left) ? other.left  : left;
        size_t r = (left < other.left) ? right : other.right;
        return ml < r;
    }
    bool testComparableFragments(const VertexProps & other) const {
        std::set<std::string> thisFragSet;
        std::set<std::string> otherFragSet;
        for(std::string fragGrp : this->assocFragGroups){
            std::vector<std::string> fragList = strsplit(fragGrp,DupDelim);
            thisFragSet.insert(fragList.begin(),fragList.end());
        }
        for(std::string fragGrp : other.assocFragGroups){
            std::vector<std::string> fragList = strsplit(fragGrp,DupDelim);
            otherFragSet.insert(fragList.begin(),fragList.end());
        }
        std::set<std::string> * smaller = &thisFragSet;
        std::set<std::string> * larger = &otherFragSet;
        if(thisFragSet.size() > otherFragSet.size()){
            std::swap(smaller,larger);
        }
        return std::includes(   larger->begin(),larger->end(),
                                smaller->begin(),smaller->end());
    }
    int compare(const VertexProps & other) const { // <, ==, > :::: -1,0,1
        if(isHost != other.isHost){ //Host before virus
            return (isHost) ? -1 : 1;
        }
        if(chromosome != other.chromosome) { //string compare
            return (chromosome < other.chromosome) ? -1 : 1;
        }
        if(opensLeft != other.opensLeft) { //Left before right
            return (opensLeft) ? -1 : 1;
        }
        if(left != other.left){ //Leftmost first
            return (left < other.left) ? -1 : 1;
        }
        if(right != other.right){ //Rightmost first
            return (right < other.right) ? -1 : 1;
        }
        if(fromSplit != other.fromSplit) { //Split before unsplit
            return (fromSplit) ? -1 : 1;
        }
        if(assocFragGroups.size() != other.assocFragGroups.size()){ //Less frags before more
            return (assocFragGroups.size() < other.assocFragGroups.size()) ? -1 : 1;
        }
        std::string thisStr = strjoin(assocFragGroups.begin(),assocFragGroups.end(),FragDelim);
        std::string otherStr = strjoin(other.assocFragGroups.begin(),other.assocFragGroups.end(),FragDelim);
        if(thisStr != otherStr){ //string compare
            return (thisStr < otherStr) ? -1 : 1;
        }
        return 0;
    }
};

struct EdgeProps {
    igraph_integer_t id;
    double weight;
    bool fromSplit;
    std::vector<std::string> assocFragGroups;
    std::pair<igraph_int_t,igraph_int_t> endpoints;
};

    //Members
public:
    static const char DupDelim = 29;
    static const char FragDelim = 30;
protected:
    igraph_t graph;
    uint8_t flag;
    // Internal caches for O(log N) vertex uniqueness checks and property tracking
    std::multimap<std::string, igraph_integer_t> vertex_by_fragment;
    std::vector<EdgeProps> queuedEdges;
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
    std::pair<igraph_int_t,igraph_int_t> get_edge_endpoints(igraph_int_t eid) const;
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
    void mergeUninformitiveOverlap();
    void ensureConstructed();
    static CRegionGraph merge_graphs(const std::vector<CRegionGraph> & graphVec);
    void write_edgelist(FILE * outstream) const;
private:
    void assertOwnership();
    void checkAndQueueEdge(igraph_integer_t v1_id, igraph_integer_t v2_id);
    void ensureValidLookup();
    void init_attribute_table();
    static std::vector<std::string> intersectFragmentGroups(
            std::vector<std::string>, std::vector<std::string>);
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
        for(std::string frag : strsplit(fragGrp,DupDelim)) {
            auto range = vertex_by_fragment.equal_range(frag);
            for(auto i = range.first; i != range.second; i++){
                vertexWithFragmentSet.insert(i->second);
            }
            //Update the lookup with the new verticies
            vertex_by_fragment.insert({frag,new_vid});
        }
    }

    std::string assocFragStr = strjoin( assocFragGroups.begin(),
                                        assocFragGroups.end(), FragDelim);

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
        checkAndQueueEdge(old_vid, new_vid);
    }
}


void CRegionGraph::addOrUpdateVertex(const VertexProps & prop) {
    this->addOrUpdateVertex(prop.chromosome, prop.opensLeft, prop.fromSplit,
                            prop.isHost, prop.left, prop.right,
                            prop.assocFragGroups);
}


// Evaluates window overlaps and constructs a weighted edge if conditions match
void CRegionGraph::checkAndQueueEdge(  igraph_integer_t v1_id,
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
    igraph_integer_t new_eid = igraph_ecount(&graph) + queuedEdges.size();
    queuedEdges.emplace_back();
    EdgeProps & props = queuedEdges.back();
    props.id = new_eid;
    props.endpoints = {v1_id,v2_id};

    std::vector<std::string> fGSVec = CRegionGraph::intersectFragmentGroups(
                                        v1.assocFragGroups,v2.assocFragGroups);
    //No shared fragment groups 
    if(!fGSVec.size()){ return; }


    props.weight = fGSVec.size();
    props.fromSplit = v1.fromSplit || v2.fromSplit;
    props.assocFragGroups = fGSVec;
}

//If there are any queued edges
//the graph is updated with the edges
void CRegionGraph::ensureConstructed() {
    //No edges to update
    if(!queuedEdges.size()){ return; }
    assertOwnership();
    igraph_vector_int_t adjList;
    igraph_vector_int_init(&adjList,queuedEdges.size() * 2);
    size_t counter = 0;
    for(auto & props : queuedEdges){
        VECTOR(adjList)[counter++] = props.endpoints.first;
        VECTOR(adjList)[counter++] = props.endpoints.second;
    }
    igraph_add_edges(&graph,&adjList,nullptr);
    igraph_vector_int_destroy(&adjList);

    for(auto & props : queuedEdges){
        std::string assocFrag = strjoin(props.assocFragGroups.begin(),
                                        props.assocFragGroups.end(), FragDelim);
        SETEAN(&graph, "weight", props.id, props.weight);
        SETEAB(&graph, "FromSplit", props.id, props.fromSplit);
        SETEAS(&graph, "assocFragGrps", props.id, assocFrag.c_str());
    }
    queuedEdges.clear();
}

//Checks if the unique vertex lookup is valid, and rebuilds it if not
void CRegionGraph::ensureValidLookup() {
    if(flag & VALID_LOOKUP) { return; }
    vertex_by_fragment.clear();
    for(int i = 0; i < this->vcount(); i++){
        for(std::string fragGrp : strsplit(VAS(&graph,"assocFragGrps",i),FragDelim)){
            for( std::string fragName : strsplit(fragGrp,DupDelim)){
                vertex_by_fragment.insert({fragName,i});
            }
        }
    }
    flag |= VALID_LOOKUP;
}

//Removes edges which do not have sufficient support
void CRegionGraph::filterEdges( double minWeight, double splitBonus) {
    ensureConstructed();
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



std::vector<std::string> CRegionGraph::intersectFragmentGroups(
        std::vector<std::string> fg1 , std::vector<std::string> fg2)
{
    std::map<std::string,std::pair<std::set<size_t>,std::set<size_t>>> profileByFragMap;
    //Construct the profiles for each fragment
    const std::vector<std::string> * fg_ptr[2] = {&fg1,&fg2};
    for(int i = 0; i < 2; i++){
       for(size_t grIdx = 0; grIdx < fg_ptr[i]->size(); grIdx++){
           for(const std::string & frag :
                   strsplit((*fg_ptr[i])[grIdx],CRegionGraph::DupDelim))
           {
                if(i == 0){
                    profileByFragMap[frag].first.insert(grIdx);
                } else {
                    profileByFragMap[frag].second.insert(grIdx);
                }
           }
       }
    }
    //Group Fragments with the same profile together
    std::map<std::pair<std::set<size_t>,std::set<size_t>>,std::set<std::string>> fragGrpByProfileMap;
    for( const auto & pair : profileByFragMap ){ //frag,profile pair(reg1 profile, reg2 profile)
        //ensure that the fragment appears in at least one fragGroup in both regions
        if(pair.second.first.size() * pair.second.second.size() > 0){
            fragGrpByProfileMap[pair.second].insert(pair.first);
        }
    }
    std::vector<std::string> fragGroupStrVec;
    for( const auto & pair : fragGrpByProfileMap){
        fragGroupStrVec.push_back( strjoin( pair.second.begin(),
                                            pair.second.end(),
                                            CRegionGraph::DupDelim) );
    }
    return fragGroupStrVec;
}

//Returns the vertex id's for the endpoints of an edge
//Always returns them such that the first vertex is the host vertex
std::pair<igraph_int_t,igraph_int_t> CRegionGraph::get_edge_endpoints(
        igraph_int_t eid) const
{
    std::pair<igraph_int_t,igraph_int_t> endpoints;
    igraph_edge(&graph,eid,&endpoints.first,&endpoints.second);
    if(VAB(&graph,"IsHost",endpoints.second)){
        std::swap(endpoints.first,endpoints.second);
    }
    return endpoints;
}

CRegionGraph::EdgeProps CRegionGraph::get_edge_properties(int id) const {
    return {    id,
                EAN(&graph,"weight",id),
                EAB(&graph,"FromSplit",id),
                strsplit(EAS(&graph,"assocFragGrps",id),FragDelim)
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
                strsplit(VAS(&graph,"assocFragGrps",id),FragDelim)
    };
}

//Goes through verticies and merges together those which represent overlapping regions
//for which the set of fragments is comparable (|A + B| = max(|A|,|B|)) 
void CRegionGraph::mergeUninformitiveOverlap() {
    assertOwnership();
    ensureConstructed();
    std::vector<VertexProps> sortedVertexProps;
    for(igraph_int_t i = 0; i < this->vcount(); i++){
        sortedVertexProps.push_back(this->get_vertex_properties(i));
    }
    std::sort(  sortedVertexProps.begin(),sortedVertexProps.end(),
                [](const VertexProps & a, const VertexProps & b) {
                    return (a.compare(b) == -1);
                } );
    std::set<igraph_int_t> toRemoveSet;
    std::vector<VertexProps> toAdd;
    const VertexProps * l = &sortedVertexProps.front();
    VertexProps m;
    m.id=-1; // Use an id of -1 to indicate an unititialized merged vertex
    for(size_t i = 1; i < sortedVertexProps.size(); i++){
        const VertexProps & r = sortedVertexProps[i];
        //If Overlapping and one contains all the fragments in the other
        if(l->testOverlap(r) && l->testComparableFragments(r)){
            //These two can be combined
            if(m.id == -1) {
                //Copy chromosome, opensLeft, isHost
                //Also copy the id which is used to indicate if m contains
                //  a merged vertex or not
                m = *l;
            }
            m.left = std::min(l->left,r.left);
            m.right = std::max(l->right,r.right);
            m.fromSplit = l->fromSplit || r.fromSplit;
            //keep unique fragment groups
            std::set<std::string> fragGrpSet;
            fragGrpSet.insert(l->assocFragGroups.begin(),l->assocFragGroups.end());
            fragGrpSet.insert(r.assocFragGroups.begin(),r.assocFragGroups.end());
            m.assocFragGroups.clear();
            m.assocFragGroups.insert(   m.assocFragGroups.end(),
                                        fragGrpSet.begin(),fragGrpSet.end());
            //Record the old vertexes to remove
            toRemoveSet.insert(l->id);
            toRemoveSet.insert(r.id);
            //Set the merged vertex to be the left for the next comparison
            l = &m;
        } else { // No merge
            if(m.id != -1){ //If there was a a merged vertex, record it
                toAdd.push_back(m);
                m.id = -1;
            }
            //Set the right vertex to be the elft for the enxt comparison
            l = &r;
        }
    }
    if(m.id != -1){ //Check if there is a merged vertex awaiting recording (last compare resulted in a merge)
        toAdd.push_back(m);
    }
    if(!toRemoveSet.size()){
        return; // There were no merges
    }
    //Invalidate the old lookup table
    flag &= ~VALID_LOOKUP;
    //Delete the old verticies
    igraph_vector_int_t toRemove;
    igraph_vector_int_init(&toRemove,toRemoveSet.size());
    size_t counter = 0;
    for(igraph_int_t id : toRemoveSet){
        VECTOR(toRemove)[counter++] = id;
    }
    igraph_delete_vertices(&graph,igraph_vss_vector(&toRemove));
    igraph_vector_int_destroy(&toRemove);
    //Add the merged vertices back in
    for(const auto & prop : toAdd){
        this->addOrUpdateVertex(prop);
    }
}


CRegionGraph CRegionGraph::merge_graphs(const std::vector<CRegionGraph> & graphVec) {
    //Collate information from the graphs
    size_t totalVertices = 0;
    std::vector<igraph_int_t> adjVec;
    std::vector<VertexProps> vPropVec;
    std::vector<EdgeProps> ePropVec;
    //Iterate over each graph and store the adjacency information
    //  as well as vertex and edge properties
    for(const CRegionGraph & graphObj : graphVec) {
        igraph_adjlist_t adjList;
        //documentation says performance hit from not using LOOPS, or MULTIPLE
        //Construction of graphs should prevent such edges anyway
        igraph_adjlist_init(&graphObj.graph,&adjList,IGRAPH_ALL,IGRAPH_LOOPS_TWICE,IGRAPH_MULTIPLE);
        igraph_int_t adjListSize = igraph_adjlist_size(&adjList);
        //Copy the adjacency list, with vertices offset
        //Also retain the vertex and edge attributes
        size_t futureTotalV = vPropVec.size() + adjListSize;
        vPropVec.resize(futureTotalV);
        //Iterate over node ids in the graph
        for(igraph_int_t i = 0; i < adjListSize; i++){
            igraph_int_t vid = i+totalVertices;
            //Store vertex properties
            vPropVec[vid] = graphObj.get_vertex_properties(i);
            vPropVec[vid].id = vid;
            //Resize vectors
            igraph_vector_int_t * vec_ptr = igraph_adjlist_get(&adjList,i);
            igraph_int_t vecSize = igraph_vector_int_size(vec_ptr);
            //Iterate over neighbours
            for(igraph_int_t j = 0; j < vecSize; j++){
                //Get the vid of the neighbour
                igraph_int_t n = VECTOR(*vec_ptr)[j];
                //Only count undirected edges once (when the lower index vertex is checked)
                if(n <= i){ continue; }
                adjVec.push_back(i + totalVertices);
                adjVec.push_back(n + totalVertices);
                //Get the edge info for the edge between node i and n
                igraph_vector_int_t eids;
                igraph_get_all_eids_between(&(graphObj.graph),&eids,i,n,
                                            IGRAPH_UNDIRECTED);
                //Assumes there is only one edge between any given pair of nodes
                igraph_int_t eid;
                igraph_get_eid(&(graphObj.graph),&eid,i,n,IGRAPH_UNDIRECTED,false);
                if(eid == -1) {continue; }
                igraph_int_t mergedEid = ePropVec.size();
                ePropVec.push_back(graphObj.get_edge_properties(eid));
                ePropVec.back().id = mergedEid;
            }
        }
        totalVertices = futureTotalV;
        //Cleanup
        igraph_adjlist_destroy(&adjList);
    }

    //Construct a combined igraph_adjlist_t object
    //igraph_adjlist_t merged_adjlist;
    //igraph_adjlist_init_empty(&merged_adjlist,totalVertices);
    //for(size_t vid1 = 0; vid1 < mergedAdjList.size(); vid1++){
    //    igraph_vector_int_t * vec_ptr = igraph_adjlist_get(merged_adjlist,vid1);
    //    for(igraph_int_t vid2 : mergedAdjList[vid1]){
    //        igraph_vector_int_push_back(vec_ptr,vid2);
    //    }
    //}
    //Create an igraph_vector_int_t object
    igraph_vector_int_t edges = igraph_vector_int_view(adjVec.data(),adjVec.size());
    //Construct the merged graph
    CRegionGraph mergedRegGraph;
    mergedRegGraph.flag &= ~VALID_LOOKUP;
    igraph_t & mergedGraph = mergedRegGraph.graph;
    igraph_empty(&mergedGraph,totalVertices,IGRAPH_UNDIRECTED);
    igraph_add_edges(&mergedGraph, &edges, nullptr);
    //Update the vertex information of the merged graph
    for(const VertexProps & vProp : vPropVec){
        std::string assocFragStr = strjoin(vProp.assocFragGroups.begin(),
                                           vProp.assocFragGroups.end(),
                                           FragDelim);
        SETVAS(&mergedGraph, "Chromosome", vProp.id, vProp.chromosome.c_str());
        SETVAB(&mergedGraph, "OpensLeft", vProp.id, vProp.opensLeft);
        SETVAB(&mergedGraph, "FromSplit", vProp.id, vProp.fromSplit);
        SETVAB(&mergedGraph, "IsHost", vProp.id, vProp.isHost);
        SETVAN(&mergedGraph, "Left", vProp.id, vProp.left);
        SETVAN(&mergedGraph, "Right", vProp.id, vProp.right);
        SETVAS(&mergedGraph, "assocFragGrps", vProp.id, assocFragStr.c_str());
    }

    //Update the edge information of the merged graph
    for(const EdgeProps & eProp  : ePropVec){
        std::string assocFrag = strjoin(eProp.assocFragGroups.begin(),
                                        eProp.assocFragGroups.end(), FragDelim);
        SETEAN(&mergedGraph, "weight", eProp.id, eProp.weight);
        SETEAB(&mergedGraph, "FromSplit", eProp.id, eProp.fromSplit);
        SETEAS(&mergedGraph, "assocFragGrps", eProp.id, assocFrag.c_str());
    }
    //Return the merged graph
    return mergedRegGraph;
}

void CRegionGraph::write_edgelist(FILE * outstream) const {
    if(queuedEdges.size()){
        throw std::logic_error("Attempt to write edge list before construction\n");
    }
    igraph_write_graph_edgelist(&graph,outstream);
}

#endif // REGION_GRAPH_H
