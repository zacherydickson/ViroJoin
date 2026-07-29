#include <algorithm>
#include <array>
#include <iostream>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <filesystem>
#include <forward_list>
#include <memory>
#include <regex>
#include <unordered_map>
#include <unordered_set>

#include "config.h"
#include "ChimericFragment.h"
#include "edge_utils.h"
#include "utils.h"
#include "sam_utils.h"
#include "str_utils.h"
#include <ssw.h>
#include <ssw_cpp.h>
#include <cptl_stl.h>
#include "BranchedQueue.h"

//==== TYPE DECLARATIONS

//Read_t, Edge_t, Region_t, SQPair_t, ReadPair_t, ReadPairAlnSummary_t
//As well as vectors, pointers, maps, and sets of these types are defined in 
//edge_utils.h

struct RRLabelAssoc_t {
    std::string readName;
    long int regionId;
    uint16_t flag;
}; 

struct AlignmentTableRow_t {
    ReadPair_pt label;
    std::string seq;
    size_t nFill;
};

typedef std::vector<AlignmentTableRow_t> AlignmentTable_t;

//CBranchedQueue defined in BranchedQueue.hpp

using CBranchedEdgeQueue = CBranchedQueue<Edge_t,ReadPair_pt,ReadPairSet_t,AlignmentMap_t>;

//==== GLOBAL VARIABLE DECLARATIONS

//std::unordered_set<std::string> VirusNameSet;
config_t Config;
stats_t Stats;
bam_hdr_t* JointHeader;
//For an alignment to pass it must have a score of at least 30
//NOTE:The Built in Filter is not used as it leads to alignments without cigar strings
//which explode when we go to test them, by using a default filter, everything
//has a cigar 
StripedSmithWaterman::Filter AlnFilter;//(true,true,30,32767);
StripedSmithWaterman::Aligner Aligner(1,4,6,1,false);
int32_t AlnMaskLen;
bool ExploratoryDeduplication = false;
static size_t MinimumReads = 4;
static size_t SplitBonus = 1;
static double MaxDiffRate = 0.06;

std::mutex Mtx;

//==== FUNCTION DECLARATIONS

void AlignRead( int id, const Read_pt & read, const RegionSet_t & regSet,
                AlignmentMap_t & alnMap);
void AlignReads(const Read2RegionsMap_t &regMap,
                AlignmentMap_t & alnMap);
bool AreConsistentCigars(   std::vector<uint32_t> vec1,
                            std::vector<uint32_t> vec2,
                            bool bFromBack);
AlignmentTable_t BuildAlignmentTable(   const Edge_t & edge,
                                        const AlignmentMap_t & alnMap);
EdgeVec_t ConsensusSplitEdge( Edge_t & edge, const AlignmentMap_t & alnMap);
bool ConstructBamEntry( const Read_pt & query, const Region_pt & subject,
                        bool isSupplemental,
                        const Read_pt & mate, const Region_pt & mateSubject, 
                        const AlignmentMap_t & alnMap,
                        bam1_t* entry);
breakpoint_t ConstructBreakpoint(const Region_pt & reg,size_t offset);
call_t ConstructCall(   int id, const Edge_t & edge,
                        const AlignmentMap_t & alnMap,
                        const ReadPairSet_t & used);
void ConstructEdgeQueue(const EdgeVec_t & edgeVec,
                        CBranchedEdgeQueue & edgeQueue);
//JunctionInterval_t ConstructJIV(char strand, bool isVirus,
//                                const StripedSmithWaterman::Alignment & aln);
//std::vector<uint32_t> ConstructJointModCigar(
//    const StripedSmithWaterman::Alignment & hAln,
//    const StripedSmithWaterman::Alignment & vAln,
//    bool bHRev, bool bVRev);
void DeduplicateEdge(Edge_t & edge,const AlignmentMap_t & alnMap);
size_t FillStringFromAlignment( std::string & outseq,
                                const std::string & inseq,
                                size_t offset, size_t maxLen,
                                const std::vector<uint32_t> & cigarVec);
void FilterEdgeVec(EdgeVec_t & edgeVec, const ReadPairSet_t * used = nullptr);
void FilterHighInsertReads(Edge_t & edge, const AlignmentMap_t & alnMap);
void FilterSuspiciousReads(Edge_t & edge, const AlignmentMap_t & alnMap);
template<class T>
void FilterVector(  std::vector<T> & vec,
                    const std::unordered_set<size_t> & idxSet);
//std::string GenerateConsensus(const std::vector<std::string> & rowVec,
//                                std::vector<size_t> * diffVec = nullptr);
std::string GenerateConsensus(const AlignmentTable_t & alnTable,
                                std::vector<size_t> * diffVec = nullptr);
std::string GetAlignedSequence( const Edge_t & edge, const ReadPair_pt & rp,
                                const AlignmentMap_t & alnMap, size_t & nFill);
EdgeVec_t LoadEdges(std::string edgeFName, std::string feFName, 
                    const Name2ReadPairMap_t & rpMap,
                    const RegID2RegionMap_t & regMap,
                    const AlignmentMap_t & alnMap);
std::vector<RRLabelAssoc_t> LoadReadRegionAssoc(const std::string & rrFName);
ReadSet_t LoadReads(const std::string & bamFName,
                    const std::unordered_set<std::string> * readNames = nullptr);
ReadVec_t LoadReadsInRef(
        const std::string & alnFileName, int tid,hts_pos_t beg, hts_pos_t end,
        const std::unordered_set<std::string> * rNames = nullptr);
RegID2RegionMap_t LoadRegions(const std::string jointRefFName,
                        const std::string regCandFName);
void LoadRegionSeq( const std::string & regionsFName,
                    Region2ReadsMap_t & reg2readSetMap,
                    Name2RegionMap_t & nameMap);
void OrderEdges(EdgeVec_t & edgeVec,const AlignmentMap_t & alnMap);
void OutputEdgeCall(int id, const Edge_t & edge, const AlignmentMap_t & alnMap,
                std::ofstream & out, const ReadPairSet_t & used);
void OuputEdgeReads(int id, const Edge_t & edge, const AlignmentMap_t & alnMap,
                const std::string & readDir, const ReadPairSet_t & used);
void OutputEdgeBP(  int id, std::ofstream & hostOut, std::ofstream & virusOut,
                    const Edge_t & edge, const AlignmentMap_t & alnMap,
                    const ReadPairSet_t & used);
void OutputEdgesByQ(EdgeVec_t & edgeVec,const AlignmentMap_t & alnMap,
                    const Name2ReadMap_t & readNameMap,
                    const std::string & resFName, const std::string & readDir,
                    const std::string & hostbpFName,
                    const std::string & virusbpFName);
bool PassesEffectiveReadCount(  const Edge_t & edge,
                                const ReadPairSet_t * used = nullptr);
void ProcessEdge(int id,Edge_t & edge, const AlignmentMap_t & alnMap);
void ProcessEdges(EdgeVec_t & edgeVec, const AlignmentMap_t & alnMap);
//EdgeVec_t RecursiveSplitEdge(Edge_t & edge, std::vector<ReadPair_pt> rowLabelVec,
//                        std::vector<std::string> rowSeqVec,
//                        std::vector<size_t> nFillVec);
template<class T>
EdgeVec_t RecursiveSplitEdge(   Edge_t & edge,
                                std::vector<T> vec,
                                std::vector<bool> (*globalTest)(const std::vector<T> &),
                                bool (*pairwiseTest)(const T &, const T &)
                                );
void RemoveUnalignedReads(EdgeVec_t & edgeVec,const AlignmentMap_t & alnMap);
void SortEdgeVec(   EdgeVec_t & edgeVec, const AlignmentMap_t & alnMap,
                    const ReadPairSet_t & used);
EdgeVec_t SplitEdges(   EdgeVec_t & edgeVec, const AlignmentMap_t & alnMap,
                        EdgeVec_t (*edgeSplitter)(Edge_t &, const AlignmentMap_t &));
std::vector<bool> TestConsistencyGlobal( const AlignmentTable_t & alnTable);
bool TestConsistencyPairwise( const AlignmentTableRow_t & a,
                              const AlignmentTableRow_t & b);

//==== MAIN

//Parse Fasta Files for regions and reads as well as a file defining
//edges
//Regions sequences are expected to be the same strand as labelled
//Region names should be in the format CONTIG,OFF,END,STRAND(:+) 
//Reads are always on the reference strand regardless of fwd/rev
//Read names should be in the format READID_[12]
//
//After loading each read is mapped to all associated regions with
//poor alignments filtered away
//
//Edges are then processed, before each step checking if enough reads
//remain
//  Deduplication (based on alignment positions and cigars)
//  Split Edges By Consensus Sequence (Ex. R1->S1 R2->S2, R1->S1) -> 2 edges
//  Sum the score for each edge (conditional(from shared reads), and
//        unconditional (score from edge specific reads)
//  Sort edges on unconditional score, then by conditional score
//  Accept from best to worst until none pass anymore
int main(int argc, const char* argv[]) {

    //## Parse Inputs
    std::string joint_ref_file_name = argv[1];
    std::string workdir = argv[2];
    std::string workspace = argv[3];
    //## Files to be used from the workdir
    std::string stats_file_name = workspace + "/stats.txt";
    std::string config_file_name = workdir + "/config.txt";
    //std::string region_fasta_file_name = workdir + "/regions.fna";
    //std::string read_fasta_file_name = workdir + "/edge_reads.fna";
    //std::string edge_file_name = workdir + "/edges.tab";


    std::string bam_file_name = workspace + "/retained-pairs.remapped.cs.bam";
    std::string edge_file_name = workdir + "/edge-candidates.tab";
    std::string fragment_edge_file_name = workdir + "/fragment-edge-associations.tab";
    std::string read_region_file_name = workdir + "/read-region-associations.tab";
    std::string region_bed_file_name = workdir + "/region-candidates.bed";
    //## Output Files
    std::string res_file_name = workdir + "/results.txt";
    std::string reads_dir = workdir + "/readsx";
    std::string hostbp_file_name = workdir + "/host_bp_seqs.fa";
    std::string virusbp_file_name = workdir + "/virus_bp_seqs.fa";
    //Quick test for required files
    for(const std::string & path : {
            joint_ref_file_name, bam_file_name, edge_file_name,
            fragment_edge_file_name, read_region_file_name,
            region_bed_file_name, reads_dir} )
    {
        if(!std::filesystem::exists(path)){
            fprintf(stderr,"[ERROR] %s does not exist",path.c_str());
            return 1;
        }
    }
    //## Configuration
    //LoadVirusNames(virus_ref_file_name,VirusNameSet);
    Config = parse_config(config_file_name);
    AlnMaskLen =  Config.read_len/2;
    //Set global variable across edges
    Edge_t::MinimumClipLen = Config.min_sc_size;
    ExploratoryDeduplication = Config.explore;
    Stats = parse_stats(stats_file_name);
    JointHeader = sam_hdr_read(sam_open(bam_file_name.c_str(),"r"));
    //## Raw Data
    Read2RegionsMap_t  read2regSetMap;
    Region2ReadsMap_t  reg2readSetMap;
    Name2ReadMap_t readNameMap;
    
    std::vector<RRLabelAssoc_t> rrLabelAssocVec =
        LoadReadRegionAssoc(read_region_file_name);
    std::unordered_set<std::string> readNameSet;
    //Get Set of unique read names
    for(const auto & assoc : rrLabelAssocVec){
        readNameSet.insert(assoc.readName);
    }
    ReadSet_t readSet = LoadReads(bam_file_name,&readNameSet);
    RegID2RegionMap_t regIdtoRegionMap = LoadRegions(joint_ref_file_name,
                                        region_bed_file_name);
    fprintf(stderr,"Indexing Read Pairs ...\n");
    //Index reads by read name - These are also the fragment objects
    Name2ReadPairMap_t rNametoReadPairMap;
    for(const auto & read : readSet){
        if(!rNametoReadPairMap.count(read->name)) {
            rNametoReadPairMap.insert({read->name,ReadPair_pt( new ReadPair_t)});
        }
        rNametoReadPairMap[read->name]->getRead(read->isR1) = read;
    }
    fprintf(stderr,"Indexed %lu Read Pairs\n",rNametoReadPairMap.size());
    fprintf(stderr,"Associating read pairs with regions...\n");
    for(RRLabelAssoc_t & labelAssoc : rrLabelAssocVec) {
        auto & readPair = rNametoReadPairMap.at(labelAssoc.readName);
        auto & region = regIdtoRegionMap.at(labelAssoc.regionId);
        for(uint32_t pairBit : {BAM_FREAD1, BAM_FREAD2} ){
            //Skip if the particular read in the pair isn't associated with the region
            if(!(labelAssoc.flag & pairBit)) { continue; }
            bool checkingR1 = (pairBit == BAM_FREAD1);
            read2regSetMap[readPair->getRead(checkingR1)].insert(region);
            reg2readSetMap[region].insert(readPair->getRead(checkingR1));
        }
    }
    fprintf(stderr,"Bi-directionally Associated %lu ReadPairs and %lu Regions...\n", read2regSetMap.size(), reg2readSetMap.size());
    
    //Perform alignments
    AlignmentMap_t alnMap;
    AlignReads(read2regSetMap,alnMap);

    //Load Edges
    EdgeVec_t edgeVec = LoadEdges(  edge_file_name, fragment_edge_file_name,
                                    rNametoReadPairMap, regIdtoRegionMap,
                                    alnMap);
    //Note: Some of this effort could be performed during the loading step
    ////## Alignments 
    RemoveUnalignedReads(edgeVec,alnMap);
    ////## Edge Processing
    OrderEdges(edgeVec,alnMap);
    ////## Output
    //TODO: FIXME: The same readpair has been found supporting multiple breakpoints
    //      This is unacceptable
    OutputEdgesByQ(edgeVec,alnMap,readNameMap,res_file_name,reads_dir,
                hostbp_file_name,virusbp_file_name);
    ////## Cleanup
    bam_hdr_destroy(JointHeader);
    fprintf(stderr,"edge_mapper Done\n");
    return 0;
}

//==== FUNCTION DEFINITIONS

//Aligns a read to all regions with which it is associated
//  filters the alignments based on a minimum score, and a minimum score
//  relative to the best alignment
//Inputs - an id, used by thead pool
//         - a read object
//         - a set of regions
//         - a constant reference to a mapping of subject query pairs generated
//                alignment objects
//Output - None, modifies the alignment mapping
void AlignRead( const Read_pt & read, const RegionSet_t & regSet,
                AlignmentMap_t & alnMap)
{
    //Need to track host and virus scores separately
    //0th element is host, 1st element is virus
    std::array<uint16_t,2> bestScore = {0,0};
    //Iterate over regions and do the alignments
    std::vector<const SQPair_t *> sqPairVec;
    for( const Region_pt & reg : regSet){
        StripedSmithWaterman::Alignment aln;
        aln.sw_score=0;
        for(char strand : {'-' , '+'}){ //Align in both the fwd and reverse orientations
            const std::string & query = read->seq.get(strand == '+');
            StripedSmithWaterman::Alignment curAln;
            int step = 0;
            bool bPass = true;
            uint32_t opdata;
            bool bClipped = false;
            do {
                switch (step) {
                    case 0: //Does the read align?
                        bPass = Aligner.Align(  query.c_str(),
                                                reg->sequence.c_str(),
                                                reg->sequence.length(),
                                                AlnFilter, &(curAln),AlnMaskLen);
                        //Co-opting unused variable in structure to store the strand of the read
                        curAln.sw_score_next_best = (uint16_t) strand;
                        break;
                    case 1: //Is alignment better than the other strand? (or first)
                        bPass = curAln.sw_score > aln.sw_score;
                        break;
                    case 2: // alignment length, score, clippedness
                        //curAln is better than previous best
                        bPass = accept_alignment(curAln, Config.min_sc_size);
                        break;
                    case 3: // read low complexity filter
                        bPass = !is_low_complexity(query.c_str(),
                                            curAln.query_begin,curAln.query_end);
                        break;
                    case 4: // region low complexity filter
                        bPass = !is_low_complexity(reg->sequence.c_str(),
                                            curAln.ref_begin,curAln.ref_end);
                        break;
                    case 5: //Split side filter
                        //Get the cigar op on the distal side of the read
                        opdata = curAln.cigar[  (reg->opensLeft()) ?  
                                                curAln.cigar.size()-1 : 0 ];
                        bClipped =  (cigar_int_to_op(opdata) == 'S') && 
                                    ( cigar_int_to_len(opdata) >=
                                        uint32_t(Config.min_sc_size) );
                        //Fail if an alignment is split on the distal side
                        bPass = !bClipped;
                }
            } while(++step < 6 && bPass);
            if(bPass) { aln = curAln; }
        }
        // Do not store failed alignments
        if(aln.sw_score <= 0) { continue; }
        // Store the alignment
        Mtx.lock();
        auto res = alnMap.insert(std::make_pair(  SQPair_t(reg,read), aln));
        Mtx.unlock();
        //Track all passing alignments for future filtering against
        //highest observed scores across all regions in each genome
        sqPairVec.push_back(&(res.first->first));
        if(aln.sw_score > bestScore[reg->isViral()]){
           bestScore[reg->isViral()] = aln.sw_score;
        }
    }
    std::array<double,2> minScore = {        0.75 * double(bestScore[0]),
                                        0.75 * double(bestScore[1])};
    //Iterate over sq pairs and erase any which are below threshold
    for( const SQPair_t * & pPair : sqPairVec){
        Mtx.lock();
        auto it = alnMap.find(*pPair);
        if(it->second.sw_score < minScore[pPair->subject->isViral()]){
            alnMap.erase(it);
        }
        Mtx.unlock();
    }
}

//Iterate over each read and align it to all associated regions
//  Filters out low scoring alignments
//  Paralleizes on Reads
//Inputs- a reference to a mapping from reads to regions
//	- a reference to a mapping from read-region pairs to alignments
//Output - none, Modifies the alignment map
void AlignReads(const Read2RegionsMap_t &regMap,
                AlignmentMap_t & alnMap) {
    fprintf(stderr,"Aligning reads ...\n");
    ctpl::thread_pool threadPool (Config.threads);
    std::vector<std::future<void>> futureVec;
    for(const auto & pair : regMap){
        std::future<void> future = threadPool.push( 
                [&pair,&alnMap](int id) {
                    AlignRead(pair.first,pair.second,alnMap);
                } );
        futureVec.push_back(std::move(future));
    }
    int pert = 0;
    size_t complete =0;
    for( auto & future : futureVec){
        future.get();
        complete++;
        double progress = complete / double(futureVec.size());
        if(1000.0 * progress > pert+1){
            pert = 1000 * progress;
            fprintf(stderr,"Progress: %0.1f%%\r",progress*100.0);
        }
    }
    //Single threaded
    //for(const auto & pair : regMap){
    //    AlignRead(0,pair.first,pair.second,alnMap);
    //}
    fprintf(stderr,"Passing Alignments: %zu\n",alnMap.size());
}

bool AreConsistentCigars(   std::vector<uint32_t> vec1,
                            std::vector<uint32_t> vec2,
                            bool bFromBack)
{
    size_t max = (vec1.size() < vec2.size()) ? vec2.size() : vec1.size(); 
    if(bFromBack){
        std::reverse(vec1.begin(),vec1.end());
        std::reverse(vec2.begin(),vec2.end());
    }
    //Iterate over both vectors
    for(size_t i = 0; i < max; i++) {
        uint32_t op1 = vec1[i], op2=vec2[i];
        //If equal pass
        if(op1 == op2) continue;
        //If last comparison and same op, Full Pass
        if( i+1 == max && bam_cigar_opchr(op1) == bam_cigar_opchr(op2))
            return true;
        return false;
    }
    return true;
}


AlignmentTable_t BuildAlignmentTable(   const Edge_t & edge,
                                        const AlignmentMap_t & alnMap)
{
    //Build the Table of aligned sequences
    AlignmentTable_t alignmentTable;
    //To track which rows are still be processed
    for(const ReadPair_pt & rp : edge.getSupport()){
        AlignmentTableRow_t tableRow;
        tableRow.label = rp;
        tableRow.nFill = 0;
        tableRow.seq = GetAlignedSequence(edge,rp,alnMap,tableRow.nFill);
        alignmentTable.push_back(tableRow);
    }
    return alignmentTable;
}

//Splits an edge into a number of edges for each unique consensus
//sequence of reads observed
//Inputs - an edge to process
//         - an alignment map
//         - a vector in which to store newly created edges
//Output - None, modifies the newEdges vector and edge object
EdgeVec_t ConsensusSplitEdge( Edge_t & edge, const AlignmentMap_t & alnMap) {
    AlignmentTable_t alnTable = BuildAlignmentTable(edge,alnMap);
    return RecursiveSplitEdge(edge,alnTable,&TestConsistencyGlobal,&TestConsistencyPairwise);//,rowLabelVec,rowSeqVec,nFillVec); 
}



//Given a read region pair, and alignment info, construct a bam entry for
//the pair's alignment
//Inputs - a read-region pair
//         - an alignment Map
//         - a bam1_t pointer to store the result in
//Output - Boolean whether the construction was successful or not
//         - Also modifes the entry object
bool ConstructBamEntry( const Read_pt & query, const Region_pt & subject,
                        bool isSupplemental,
                        const Read_pt & mate, const Region_pt & mateSubject,
                        const AlignmentMap_t & alnMap,
                        bam1_t* entry)
{
    const StripedSmithWaterman::Alignment & aln =
        alnMap.at(SQPair_t(subject,query));
    //sw_score_next_best has been co-opted to store whether the cannonical or non-cannonical representation of
    //the read aligned against the subject
    bool isCannonicalQuery = ((char) aln.sw_score_next_best == '+');
    //bool bRev = subject.opensLeft();
    //Information for all alignments
    uint16_t flag = query->isR1 ? BAM_FREAD1 : BAM_FREAD2;
    entry->core.qual = 255;
    entry->core.l_extranul = (4 - (query->name.length() % 4)) % 4;
    entry->core.l_qname = query->name.length() + entry->core.l_extranul;
    std::string qSeq = query->seq.get(isCannonicalQuery);
    int l_qseq = qSeq.length();
    std::vector<char> qual(l_qseq,'<');
    int l_aux = 0;
    auto ql = bam_cigar2qlen(aln.cigar.size(),aln.cigar.data());
    if(l_qseq != ql) return false;
    std::vector<uint32_t> cigar = aln.cigar;
    hts_pos_t pos = subject->offset + aln.ref_begin;
    int32_t mtid = -1;
    hts_pos_t mpos = -1;
    hts_pos_t isize = 0;
    //Extra info
    if(isSupplemental){ // Basic info only
        flag |= BAM_FSUPPLEMENTARY;
        //Supplementary alignments point away from the junction
        if(!subject->opensLeft()) {flag |= BAM_FREVERSE;}
    } else { // Add Mate Information
        //Primary Alignments point towards the junction
        if(subject->opensLeft()) {flag |= BAM_FREVERSE;}
        SQPair_t msqp(mateSubject,mate);
        auto it = alnMap.find(msqp);
        if(it != alnMap.end()){
            flag |= BAM_FPAIRED;
            if(mateSubject->opensLeft()) {flag |= BAM_FMREVERSE;}
            mtid = sam_hdr_name2tid(JointHeader,mateSubject->chromosome.c_str());
            const StripedSmithWaterman::Alignment & mateAln = it->second;
            mpos = mateSubject->offset + mateAln.ref_begin;
        } else {
            //A primary alignment's mate SHOULD have a primary alignment
            throw std::runtime_error("Missing primary mate aln");
        }
    }
    bam_set1(   entry,
                query->name.length(),query->name.c_str(),
                flag, sam_hdr_name2tid(JointHeader,subject->chromosome.c_str()), pos, 255,
                cigar.size(), cigar.data(),
                mtid, mpos, isize,
                l_qseq, qSeq.c_str(), qual.data(),
                l_aux);
    return true;
}

//Constructs a breakpoint from a known region
//Inputs - an offset describing where in the region the offset is
//         - a region
//Output - a breakpoint_t object (see util.h)
breakpoint_t ConstructBreakpoint(const Region_pt & reg,size_t offset){
    //HL -, HR +, VL +, VR -
    bool bRev = (reg->opensLeft() != reg->isViral());
    int pos = reg->offset + offset + 1;
    return breakpoint_t(reg->chromosome,pos,pos,bRev);
}

//Given an edge and alignment information, calculates summary stats
//  hostPBS - the average score/aligned base on the host side
//  coverage - the half the total length of the host and virus sides as a
//                proportion of the max insert size
//Then constructs a call_t objects which can be output
//Inputs - an identifier for the output junction
//         - an edge object
//         - an alignment map
//Output - a call_t object (see utils.h)
call_t ConstructCall(int id, const Edge_t & edge, const AlignmentMap_t & alnMap,
                        const ReadPairSet_t & used)
{
    size_t nReads = 0;
    size_t nSplit = 0;
    std::pair<size_t,size_t> offsets = edge.getOffsets();
    breakpoint_t hostBP = ConstructBreakpoint(edge.hostRegion,offsets.first);
    breakpoint_t virusBP = ConstructBreakpoint(edge.virusRegion,offsets.second);
    double hostPBS = 0, virusPBS = 0;
    double hostCov = 0, virusCov = 0;
    int32_t hostLeft = edge.hostRegion->sequence.length(), hostRight = 0;
    int32_t virusLeft = edge.virusRegion->sequence.length(), virusRight = 0;
    int score = 0;
    for(const auto & pair : edge.getSupportSummaryMap()){
        const ReadPair_pt & frag = pair.first;
        const ReadPairAlnSummary_t & summary = pair.second;
        if(used.count(frag)) { continue; }
        nReads++;
        if(summary.isSplit) { nSplit++; }
        if(summary.hostLeft < hostLeft) { hostLeft = summary.hostLeft; }
        if(summary.hostRight > hostRight) { hostRight = summary.hostRight; }
        if(summary.virusLeft < virusLeft) { virusLeft = summary.virusLeft; }
        if(summary.virusRight > virusRight) { virusRight = summary.virusRight; }
        double hLen = summary.hostQAlnBases;
        double vLen = summary.virusQAlnBases;
        hostPBS += double(summary.hostScore) / hLen;
        virusPBS += double(summary.virusScore) / vLen;
        score += summary.score();
    }
    hostPBS /= double(nReads);
    virusPBS /= double(nReads);
    if(hostLeft <= hostRight) {
        hostCov = double(hostRight - hostLeft) / (Stats.max_is - MinimumAlignmentLength);
    }
    if(virusLeft <= virusRight) {
        virusCov = double(virusRight - virusLeft) / (Stats.max_is - MinimumAlignmentLength);
    }
    return call_t(id,hostBP,virusBP,nReads,nReads,nSplit,0,0,
            score,hostPBS,virusPBS,hostCov,virusCov);
}

//Takes a vector of edges and copies them into a branched edge queue
//Which stores the edges and their reads in such a way that the top of
//the queue is always the next best edge
//Input - a const reference to a vector of edges to put into the queue
//          It isn't required, but ideally this edge vector is sorted in
//          decreasing order
//      - a 
void ConstructEdgeQueue(const EdgeVec_t & edgeVec,
                        CBranchedEdgeQueue & edgeQueue)
{
    fprintf(stderr,"Constructing edgeQueue ...\n");
    int pert = 0;
    size_t processed = 0;
    //Iterate over the edge vector from back to front
    for(auto it = edgeVec.rbegin(); it != edgeVec.rend(); it++){
        edgeQueue.addEdge(*it);
        double progress = (++processed) / double(edgeVec.size());
        if(1000.0 * progress > pert+1){
            pert = 1000 * progress;
            fprintf(stderr,"Progress: %0.1f%%\r",progress*100.0);
        }
    }
    fprintf(stderr,"Edge Queue has %lu edges across %lu branches\n",
            edgeQueue.queueSize(),edgeQueue.size());
}

////Given the an alignment from either/host or virus 
////Assigns the  reference begin/end offsets into juncton proximal/distal offsets
////  depending on the strand of the reference
//// Host alignments are distal-proximal, unless reversed
//// Viral alignments are proximal-distal, unless reversed
//// Boils down to an XOR operation on the viralness and reversedness
////  if one is true proximal-distal, else distal-proximal
//// At the end the reversed status may flip dist/proximal
////Inputs - the strand of the reference
////         - whether this alignment was against a viral reference
////         - the alignment in question
////Output - a structure containing the junction distal and proximal positions
//JunctionInterval_t ConstructJIV(char strand, bool isVirus,
//                                const StripedSmithWaterman::Alignment & aln)
//{
//
//    //sw_score_next_best has been co-opted to store the strand of the read's alignment
//    //against the subject
//    char queryStrand = (char) aln.sw_score_next_best;
//    bool bRev = (strand != queryStrand);
//    JunctionInterval_t jIV;
//    //!= is an XOR operation on boolean values
//    if(bRev != isVirus){
//        jIV.proximal = aln.ref_begin;
//        jIV.distal = aln.ref_end;
//    } else {
//        jIV.distal = aln.ref_begin;
//        jIV.proximal = aln.ref_end;
//    }
//    if(bRev) {
//        std::swap(jIV.distal,jIV.proximal);
//    }
//    return jIV;
//}

////Construct the joint cigar vector of the host and virus side in
////the order of human then virus, the proximal soft clip is
////truncated to a length of 1
////Input - an alignment to the host side
////      - an alignment to the virus side
////      - a boolean of whether the host alignment is arranged distal-proximal
////      - a boolean of whether the virus alignment is arranged proximal-distal
////Output - a vecotr of cigar ops
//std::vector<uint32_t> ConstructJointModCigar(
//    const StripedSmithWaterman::Alignment & hAln,
//    const StripedSmithWaterman::Alignment & vAln,
//    bool bHRev, bool bVRev)
//{
//    std::vector<uint32_t> opVec = hAln.cigar;
//    std::vector<uint32_t> vOpVec = vAln.cigar;
//    if(bHRev){
//        std::reverse(opVec.begin(),opVec.end());
//    }
//    if(bam_cigar_opchr(opVec.back()) == 'S'){
//        opVec.back() = bam_cigar_gen(1,'S');
//    }
//    if(bVRev){
//        std::reverse(vOpVec.begin(),vOpVec.end());
//    }
//    if(bam_cigar_opchr(vOpVec.front()) == 'S'){
//        vOpVec.front() = bam_cigar_gen(1,'S');
//    }
//    //Append the virus cigar vector
//    opVec.insert(opVec.end(),vOpVec.begin(),vOpVec.end());
//    //Convert match/mismatch (X/=) to aligned(M)
//    for(auto it = opVec.begin(); it != opVec.end();it++){
//        char opChr = bam_cigar_opchr(*it);
//        if(opChr == '=' || opChr == 'X'){
//            *it = bam_cigar_gen(bam_cigar_oplen(*it),'M');
//        }
//    }
//    //Collapse together consecutive align ops
//    for(auto it = opVec.begin(),nx=std::next(it); nx != opVec.end();){
//        char opChr = bam_cigar_opchr(*it);
//        char nextOpChr = bam_cigar_opchr(*nx);
//        if(opChr == nextOpChr && opChr == 'M'){
//            *it = bam_cigar_gen(bam_cigar_oplen(*it)+bam_cigar_oplen(*(it+1)),
//                                'M');
//            nx = opVec.erase(nx);
//        } else {
//            it = nx;
//            nx++;
//        }
//    }
//    return opVec;
//}

//Identifies Duplicate reads at an edge, duplicates are defined as having
//the same left map in the host and right map in the virus
//OR
//  One end matches and a joint cigar string of the host and virus
//  mapping segments matches
//Inputs - an edge to process
//         - an alignment map to inform the deduplication
//Output - None, modifies the given object
void DeduplicateEdge(Edge_t & edge ,const AlignmentMap_t & alnMap) {
    std::unordered_set<int32_t> hostPositionSet;
    std::unordered_set<int32_t> virusPositionSet;
    ReadPairSet_t toRemoveSet;
    typedef const std::pair<const ReadPair_pt,ReadPairAlnSummary_t>* rpp_pt;
    std::vector<rpp_pt> supportVec;
    for( const auto & pair : edge.getSupportSummaryMap() ) {
        supportVec.push_back(&pair);
    }
    //Sort such that the highest scoring read pair comes first
    std::sort(supportVec.begin(),supportVec.end(),
            [](rpp_pt a, rpp_pt b) {
                return a->second.score() > b->second.score();
            }
        );
    //Retain the first (highest scoring) read pair for any given pair of distal coordinates
    for( rpp_pt pair_ptr : supportVec ) {
        const ReadPair_pt & frag = pair_ptr->first;
        const ReadPairAlnSummary_t & summary = pair_ptr->second;
        bool bSeenHuman = hostPositionSet.count(summary.hostDistal);
        bool bSeenVirus = virusPositionSet.count(summary.virusDistal);
        //If Exploratory, then one end may be non-unique
        //otherwise both must be unique
        if( (!ExploratoryDeduplication && (bSeenHuman || bSeenVirus)) ||
            (bSeenHuman && bSeenVirus) ) 
        {
            toRemoveSet.insert(frag);
            continue;
        } 
        hostPositionSet.insert(summary.hostDistal);
        virusPositionSet.insert(summary.virusDistal);
    }
    for(const auto & frag : toRemoveSet){
        edge.removeSupport(frag);
    }
}

//Sets the characters of an output string to the appropriate characters
//from an input string according to a set of cigar operations
//Only up to maxpos in the output string will be altered
//The output is always in reference coordinates, as a result insertions
//are treated as alignments, that is they consume both the reference and
//query
//This may lead to the sequence being truncated at maxpos, which is fine
//The only situation this is a problem in is if a sequencing error exists
//at the end of a sequence with a long insertion *shrug*
//Inputs - a sequence to fill in
//         - a sequence to fill from
//         - where in the outseq to begin filling
//         - where in the outsew to stop filling
//         - a vector of bam formated cigar information
//Output - The number of defined characters filled in
size_t FillStringFromAlignment( std::string & outseq,
                                const std::string & inseq,
                                size_t offset, size_t maxpos,
                                const std::vector<uint32_t> & cigarVec) {
    size_t pos = offset;
    size_t qpos = 0;
    size_t nFill = 0;
    //Iterate over cigar operations and set characters in the outseq
    for(size_t i = 0; i < cigarVec.size(); i++){
        size_t opLen = bam_cigar_oplen(cigarVec[i]);
        char opchr = bam_cigar_opchr(cigarVec[i]);
        switch (opchr) {
            case 'H':
            case 'S':
                qpos += opLen;
                break;
            case 'M':
            case '=':
            case 'X':
                for(size_t j = 0;
                    j < opLen && pos < maxpos && qpos < inseq.size();
                    j++)
                {
                    if(outseq[pos] == 'N'){
                        nFill++;
                    }
                    outseq[pos++] = inseq[qpos++];
                }
                break;
            case 'I': //Insertions are treated as internal soft clips
                //TODO: This should eventually be addressed, for now it is the behaviour
                qpos += opLen;
                break;
            case 'D':
                for(size_t j = 0;
                    j < opLen && pos < maxpos && qpos < inseq.size();
                    j++)
                {
                    if(outseq[pos] == 'N'){
                        nFill++;
                    }
                    outseq[pos++] = '-';
                }
                break;
        }
    }
    return nFill;
}

//Filters an edge Vector to only contain edges which have enough effective
//reads
//IMPORTANT: does not maintain element order, follow up with a sort!
//Inputs - an edge vector
//         - a pointer to a set of used reads, may be null
//Output - none, modifies the given edge vector
void FilterEdgeVec(EdgeVec_t & edgeVec, const ReadPairSet_t * used){
    //size_t filtered = 0;
    for(auto it = edgeVec.begin(); it != edgeVec.end(); ){
        if(PassesEffectiveReadCount(*it,used)){
            it++;
        } else {
            //Delete the element by overwriting with last element,
            //        then delete the last element
            //Does not maintain element order, but is fast
            *it = std::move(edgeVec.back());
            edgeVec.pop_back();
            //filtered++;
        }
    }
}


//Identifies read pairs with an apparent insert size which is too large
//and removes them
//Inputs - an edge to process
//         - an alignment map 
//Output - None, modifies the edge object
void FilterHighInsertReads(Edge_t & edge, const AlignmentMap_t & alnMap){
    std::vector<ReadPair_pt> toRemoveVec;
    for(auto & pair : edge.getSupportSummaryMap()){
        const ReadPair_pt & frag = pair.first;
        const ReadPairAlnSummary_t & summary = pair.second;
        if(summary.calcIS() > Stats.max_is){
            toRemoveVec.push_back(frag);
        }
    }
    for(const ReadPair_pt & frag : toRemoveVec){
        edge.removeSupport(frag);
    }
}

//Remove reads that have suspicious alignments, alignments are considered
//suspicious if:
//  a split read's aligned position is too far from the breakpoint
//  one mate's alignment is strictly inside another ( it doesn't have a junction
//   distal position)
//Inputs - an edge to process
//         - an alignment map
//Output - None, modifies the edge object
void FilterSuspiciousReads(Edge_t & edge, const AlignmentMap_t & alnMap) {
    ReadPairSet_t toRemoveSet;
    std::pair<size_t,size_t> offsets = edge.getOffsets();
    //First pass to id strictly nested alignments
    for(const auto & pair : edge.getSupportSummaryMap()) {
        const ReadPair_pt & frag = pair.first;
        const ReadPairAlnSummary_t & summary = pair.second;
        if (!(summary.distalContribFlag & ReadPairAlnSummary_t::HAS_R1) || // R1 is inside R2
            !(summary.distalContribFlag & ReadPairAlnSummary_t::HAS_R2) ) // R2 is inside R1
        {
            toRemoveSet.insert(frag);
        }
    }
    //Second pass to id reads with split position too far
    for(const Region_pt & reg : {edge.hostRegion,edge.virusRegion}){
        long int curOffset =  (reg == edge.hostRegion) ?
                              offsets.first : offsets.second;
        //Set the most distal acceptable clip position
        for(const auto & pair : edge.getSupportSummaryMap()) {
            const ReadPair_pt & frag = pair.first;
            //Skip already excluded  fragments
            if( toRemoveSet.count(frag) ) { continue; }
            //Check Each Read
            for( bool checkR1 : {true, false}) {
                SQPair_t sqp (reg,frag->getRead(checkR1));
                //Skip unaligned reads
                auto it = alnMap.find(sqp);
                if(it == alnMap.end()) { continue; }
                auto & aln = it->second;
                //Skip unsplit alignments
                if(!Edge_t::AlignmentIsSplit(reg->opensLeft(),aln)){ continue; }
                long int proxPos =  reg->opensLeft() ?
                                    aln.ref_begin : aln.ref_end;
                //As the offset is defined as the most extreme position across frags,
                // any individual fragment's proximal position will be at least as
                // distal (no chance of underflow)
                if(std::abs(proxPos - curOffset) > Config.max_sc_dist) {
                    toRemoveSet.insert(frag);
                }
            }
        }
    }
    for(const ReadPair_pt & frag : toRemoveSet){
        edge.removeSupport(frag);
    }
}

//Generic function for filtering a vector to only a given set of indexes
//Inputs - a vector to process
//         - a set of indexes
//Output - none, modifes the given vector
template<class T>
void FilterVector(  std::vector<T> & vec,
                    const std::unordered_set<size_t> & idxSet)
{
    size_t idx = 0;
    for( auto it = vec.begin(); it != vec.end(); idx++){
        if(!idxSet.count(idx)){
            it = vec.erase(it);
        } else {
            it++;
        }
    }
}

////Given a vecotr of strings (all assumed to be the same length) construct
////a majority rule consensus string
////Also record the number of differences from the generated consensus for
////each string
////Inputs - a vector of strings
////         - a reference to a vector of size_t, will be scaled to rowVec
////         Size, and will contain the # of sites which differ from the
////         conensus for each row
////Output - a majority rule consensus sequence
//std::string GenerateConsensus(const std::vector<std::string> & rowVec,
//                                std::vector<size_t> * diffVec)
//{
//    
//    if(!rowVec.size()) return std::string();
//    if(diffVec) diffVec->assign(rowVec.size(),0);
//    std::string cons(rowVec.front().length(),'N');
//    size_t nCol = rowVec[0].length();
//    for(size_t col = 0; col < nCol; col++){
//        std::unordered_map<char,size_t> nucCount;
//        size_t max = 0;
//        char best = 'N';
//        //Iterate to determine consensus residue
//        for(const std::string & seq : rowVec){
//            char nuc = seq.at(col);
//            size_t count = ++nucCount[nuc];
//            if(count > max && nuc != 'N'){
//                best = nuc;
//                max = count;
//            }
//        }
//        cons[col] = best;
//        if(diffVec) {
//            //Iterate to count differences from consensus
//            for(size_t row = 0; row < rowVec.size(); row++){
//                char nuc = rowVec.at(row).at(col);
//                if(nuc != best && nuc != 'N'){
//                    (*diffVec)[row]++;
//                }
//            }
//        }
//    }
//    return cons;
//}

//Given a vecotr of strings (all assumed to be the same length) construct
//a majority rule consensus string
//Also record the number of differences from the generated consensus for
//each string
//Inputs - a vector of strings
//         - a reference to a vector of size_t, will be scaled to rowVec
//         Size, and will contain the # of sites which differ from the
//         conensus for each row
//Output - a majority rule consensus sequence
std::string GenerateConsensus(const AlignmentTable_t & alnTable,
                                std::vector<size_t> * diffVec)
{
    if(!alnTable.size()) return std::string();
    if(diffVec) diffVec->assign(alnTable.size(),0);
    size_t nCol = alnTable.front().seq.length();
    std::string cons(nCol,'N');
    for(size_t col = 0; col < nCol; col++){
        std::unordered_map<char,size_t> nucCount;
        size_t max = 0;
        char best = 'N';
        //Iterate to determine consensus residue
        for(const AlignmentTableRow_t & row : alnTable){
            const std::string & seq = row.seq;
            char nuc = seq.at(col);
            size_t count = ++nucCount[nuc];
            if(count > max && nuc != 'N'){
                best = nuc;
                max = count;
            }
        }
        cons[col] = best;
        if(diffVec) {
            //Iterate to count differences from consensus
            for(size_t row = 0; row < alnTable.size(); row++){
                char nuc = alnTable.at(row).seq.at(col);
                if(nuc != best && nuc != 'N'){
                    (*diffVec)[row]++;
                }
            }
        }
    }
    return cons;
}

//Uses the region information from an edge, and alignment information to
//generate the aligned sequence against a concatenation of both regions
//The Output will have N's in positions with N's or no coverage
//Inputs - an edge, containing the host and virus region information
//         - a read for which to contruct the sequence
//         - an alignment map containing alignment information for the read
//         - a reference to a size_t to store the number of filled characters
//Output - a string which represents the alignment to the host region
//            concatenated to the alignmnet for the virus region
std::string GetAlignedSequence( const Edge_t & edge, const ReadPair_pt & rp,
                                const AlignmentMap_t & alnMap, size_t & nFill)
{
    nFill = 0;
    std::string outseq("");
    for( const Region_pt & reg : { edge.hostRegion, edge.virusRegion}) {
        size_t regLen = reg->sequence.length();
        size_t offset = outseq.length();
        outseq += std::string(regLen,'N');
        //Fill in R2 then R1, which gives priority to R1 bases
        for(bool checkR1 : {false, true} ){
            const Read_pt & read = rp->getRead(checkR1);
            SQPair_t sqp(reg,read);
            auto it = alnMap.find(sqp);
            if(it == alnMap.end()) {continue; }
            const StripedSmithWaterman::Alignment & aln = it->second;
            //Determine which strand of was aligned to the region
            //sw_score_next_best has been co-opted to store the strand of the read's alignment
            //against the subject
            bool isCanonical = (char(aln.sw_score_next_best) == '+');
            const std::string & seq = read->seq.get(isCanonical);
            //Fill in the appropriate side of the alignment
            nFill += FillStringFromAlignment(   outseq,seq,aln.ref_begin+offset,
                                                outseq.length(),aln.cigar);
        }
    }
    return outseq;
}


//Given paths to the edge and fragment edge associations as well as a mapping from fragment names to 
//  ReadPair objects, construct a vector of Edges
//Inputs - a path to a tab delim file with edge ids, and region ids
//       - a path to a tab delim file with fragment names, and edge ids
//       - a const reference to a mapping between fragment names and read pairs
//Output - A vector of edge objects
EdgeVec_t LoadEdges(std::string edgeFName, std::string feFName, 
                    const Name2ReadPairMap_t & rpMap,
                    const RegID2RegionMap_t & regMap,
                    const AlignmentMap_t & alnMap)
{
    fprintf(stderr,"Loading Edges from %s...\n",edgeFName.c_str());
    EdgeVec_t edgeVec;
    std::ifstream edgeFile(edgeFName);
    long int edgeID;
    long int hostRegId, virusRegId;
    size_t support;
    std::unordered_map<long int, size_t> edgeID2VecIdxMap;
    //Construct the edges between the host and viral regions
    while(edgeFile >> edgeID >> hostRegId >> virusRegId >> support) {
        edgeID2VecIdxMap[edgeID] = edgeVec.size();
        edgeVec.emplace_back(regMap.at(hostRegId),regMap.at(virusRegId));
        edgeVec.back().id = edgeID;
    }
    fprintf(stderr,"Loaded %lu Edges\n",edgeVec.size());
    fprintf(stderr,"Loading Fragment-edge associations from %s ...\n",
            feFName.c_str());
    std::ifstream feFile(feFName);
    std::string name;
    long int groupID;
    size_t counter = 0;
    std::unordered_map<size_t,ReadPairSet_t> vecIdx2ReadPairSetMap;
    while(feFile >> name >> edgeID >> groupID){
        size_t vecIdx = edgeID2VecIdxMap.at(edgeID);
        vecIdx2ReadPairSetMap[vecIdx].insert(rpMap.at(name));
        counter++;
    }
    fprintf(stderr,"Loaded %lu fragment-edge associations\n",counter);
    fprintf(stderr,"Assigning fragments to edges ...\n");
    ctpl::thread_pool threadPool (Config.threads);
    std::vector<std::future<void>> futureVec;
    for(const auto & pair : vecIdx2ReadPairSetMap){
        auto future = threadPool.push(
                [&edgeVec,&pair,&alnMap](int id) {
                    for(const auto & frag : pair.second){
                        edgeVec[pair.first].addSupport(frag,alnMap);
                    }
                } );
        futureVec.push_back(std::move(future));
    }
    int pert = 0;
    size_t complete = 0;
    for (auto & future : futureVec) {
        future.get();
        complete++;
        double progress = complete / double(futureVec.size());
        if(1000.0 * progress > pert+1){
            pert = 1000 * progress;
            fprintf(stderr,"Progress: %0.1f%%\r",progress*100.0);
        }
    }
    fprintf(stderr,"Assigned %lu fragment-edge associations and %lu edges\n",
            counter,edgeVec.size() );
    return edgeVec;
}

//Parses a Read Region Association File and constructs a vector of associations
//  this only associates labels together
std::vector<RRLabelAssoc_t> LoadReadRegionAssoc(const std::string & rrFName) {
    fprintf(stderr,
            "Loading read region associations from %s ...\n", rrFName.c_str() );
    std::vector<RRLabelAssoc_t> rrLabelAssocVec;
    std::ifstream in(rrFName);
    RRLabelAssoc_t obj;
    while(in >> obj.readName >> obj.regionId >> obj.flag){
        rrLabelAssocVec.push_back(obj);
    }
    fprintf(stderr,
            "Loaded %lu read-region associations\n", rrLabelAssocVec.size() );
    return rrLabelAssocVec;
}

//Given a bam file, loads r1 and r2 reads,
//  a set of read names to restrict to may optionally be provided
//Inputs - a path to a bam file
//       - [optional] an allow filter set of strings - Default allow all
//Output - a set of shared read pointers
ReadSet_t LoadReads(const std::string & bamFName,
                    const std::unordered_set<std::string> * allowList_ptr)
{
    fprintf(stderr,"Loading candidate reads from %s ...\n", bamFName.c_str());
    ReadSet_t readSet;
    ctpl::thread_pool threadPool(Config.threads);
    //Get the number of references to process
    open_samFile_t* alnFile = open_samFile(bamFName.c_str(),false,true);
    int nref = sam_hdr_nref(alnFile->header);
    //First pass determine the average number of mapped reads per reference
    // as well as the number mapped per region
    uint64_t unmapped = 0;
    uint64_t total = 0;
    std::vector<uint64_t> mappedVec(nref,0);
    for(int tid = 0; tid < nref; tid++){
        hts_idx_get_stat(alnFile->idx,tid,&mappedVec[tid],&unmapped);
        total += mappedVec[tid];
    }
    double target = double(total) / double(nref);
    std::vector<std::future<ReadVec_t>> futureVec;
    for(int tid = 0; tid < nref; tid++){
        hts_pos_t nPart = std::ceil(double(mappedVec[tid]) / target);
        hts_pos_t refLen = sam_hdr_tid2len(alnFile->header,tid);
        hts_pos_t nBases = std::ceil(double(refLen) / double(nPart));
        for(int part = 0; part < nPart; part++){
            hts_pos_t beg = part * nBases;
            hts_pos_t end = (part+1) * nBases;
            if(end > refLen) {end = refLen;}
            auto future = threadPool.push(
                    [&allowList_ptr,tid,&bamFName,beg,end](int id ) {
                        return LoadReadsInRef(bamFName,tid,beg,end,allowList_ptr);
                    } );
            futureVec.push_back(std::move(future));
        }
    }
    close_samFile(alnFile);
    double counter = 0;
    int perc = 0;
    for(auto & future : futureVec){
        auto readVec = future.get();
        readSet.insert(readVec.begin(),readVec.end());
        double progress = (counter++)/double(futureVec.size());
        if(100.0* progress > perc+1) {
            perc = 100.0 * progress;
            fprintf(stderr,"Progress: %0.1f%%\r",progress*100.0);
        }
    }
    fprintf(stderr,"Loaded %lu reads\n",readSet.size());
    return readSet;
}

ReadVec_t LoadReadsInRef(const std::string & alnFileName,
                            int tid, hts_pos_t beg, hts_pos_t end,
                            const std::unordered_set<std::string> * allowList_ptr)
{
    ReadVec_t readVec;
    open_samFile_t* alnFile = open_samFile(alnFileName.c_str(),false,true);
    //hts_pos_t refLen = sam_hdr_tid2len(alnFile->header,tid);
    hts_itr_t * itr = sam_itr_queryi(alnFile->idx,tid,beg,end);
    bam1_t* entry = bam_init1();
    while(sam_itr_next(alnFile->file,itr,entry) >= 0){
        //If an allow list was provided, only keep reads in that allow list
        if(!allowList_ptr || allowList_ptr->count(bam_get_qname(entry))){
            readVec.push_back(std::make_shared<Read_t>(entry));
        }
    }
    bam_destroy1(entry);
    hts_itr_destroy(itr);
    close_samFile(alnFile);
    return readVec;
}



//Given a pair of files which contain the regions of interest, and the raw
//  fasta sequences from which those sequences were drawn, extracts the info
//  and constructs a set of region objects
//
//Inputs - a path to a joint fasta reference
//       - a path to a region candidate bed file
RegID2RegionMap_t LoadRegions(const std::string jointRefFName,
                        const std::string regCandFName) 
{
    fprintf(stderr,
            "Loading candidate regions from %s and %s ...\n", 
            jointRefFName.c_str(), regCandFName.c_str());
    RegionSet_t regSet;
    RegID2RegionMap_t regIdtoRegionMap;
    std::ifstream regCandFile(regCandFName);
    faidx_t * jointRefFai = fai_load(jointRefFName.c_str());
    if(!jointRefFai) {
        throw std::runtime_error(   "Failure to open index joint reference: " + 
                                    jointRefFName);
    }
    std::string chr;
    size_t off, end, id;
    int flag;
    char strand;
    while(regCandFile >> chr >> off >> end >> id >> flag >> strand){
        bool bViral = !(flag & (1 << ChimericFragment_t::FLAG_BITS));
        bool opensLeft = (flag & ChimericFragment_t::OPENS_LEFT);
        std::string reg =   "{" + chr + "}:" + std::to_string(off +1) + "-" +
                            std::to_string(end);
        hts_pos_t regLen;
        char* regSeq = fai_fetch64(jointRefFai,reg.c_str(), &regLen);
        if(!regSeq || regLen < 0) {
            throw std::runtime_error(   "Failure to extract region (" + reg + 
                                        ") from joint reference: " +
                                        jointRefFName);
        }
        std::string regSeqStr(regSeq,regLen);
        auto pair = regSet.insert(std::make_shared<Region_t>( 
                        id, chr, regSeqStr, off, end, bViral, opensLeft) );
        //the region might not be inserted if it is the same, but with a different id
        //the first elem of pair is an iterator to the inserted (or blocking) region
         regIdtoRegionMap[id] = *(pair.first);
    }
    fai_destroy(jointRefFai);
    fprintf(stderr,"Loaded %lu candidate regions with %lu ids\n",regSet.size(),regIdtoRegionMap.size());
    return regIdtoRegionMap;
}

//Process Edges and puts them into an order from most likely to be real to
//least
//Processing includes:
//  Spliting Edges which appear to have multiple
//  consensus sequences
//  Finding the actual breakpoints
//  Deduplicating
//  Removing High insert size reads
//  Identifying which if any reads are unique to the edge
//Inputs - a reference to an edge vector
//         - a const reference to an alignment map
//Output - None, modifies the edge vector
void OrderEdges(EdgeVec_t & edgeVec,const AlignmentMap_t & alnMap) {
    fprintf(stderr,"Ordering Edges ...\n");
    EdgeVec_t newEdges = SplitEdges(edgeVec,alnMap,&ConsensusSplitEdge); 
    //Eliminate Edges with low read counts
    FilterEdgeVec(edgeVec);
    //Add any new edges back in (these are already filtered)
    edgeVec.insert(edgeVec.end(),newEdges.begin(),newEdges.end());
    ////Process all edges
    ProcessEdges(edgeVec,alnMap);
    //Remove insufficiently supported Edges
    //Find the edges which are unique to a particular edge
    ReadPairSet_t used; //Sort Edge Needs an object to work with
    SortEdgeVec(edgeVec,alnMap,used);
    fprintf(stderr,"Ordered %zu edges\n",edgeVec.size());
}

void OutputEdgeCall(int id, const Edge_t & edge, const AlignmentMap_t & alnMap,
                std::ofstream & out, const ReadPairSet_t & used)
{
    //Output the call
    call_t call = ConstructCall(id, edge,alnMap, used);
    out << call.to_string() << "\n"; 
}


void OutputEdgeReads(int id, const Edge_t & edge, const AlignmentMap_t & alnMap,
                const std::string & readDir, const ReadPairSet_t & used)
{
    //Output the reads
    std::string fName = std::to_string(id)+".bam";
    samFile* writer = open_bam_writer(readDir, fName, JointHeader);
    std::vector<bam1_t*> entryVec;
    for(const auto & pair : edge.getSupportSummaryMap()){
        const ReadPair_pt & frag = pair.first;
        if(used.count(frag)) continue;
        const ReadPairAlnSummary_t & summary = pair.second;
        //Order the regions such that R1 will be primarily be against the first
        // and R2 against the second
        std::vector<Region_pt> regionOrder;
        if( (summary.distalContribFlag & 0b0110) == 0b0110) { //R1 V primary, R2 H primary
            regionOrder = {edge.virusRegion,edge.hostRegion};
        } else if ( (summary.distalContribFlag & 0b1001) == 0b1001) { //R1 H primary, R2 V primary
            regionOrder = {edge.hostRegion,edge.virusRegion};
        } else { // Impossible cases of missing info or subsumed reads
            throw std::runtime_error(   "Summary - Aln mismatch; case: " +
                                        std::to_string(summary.distalContribFlag));
        }
        //Construct the Bam Entries
        for(bool checkR1 : {true, false} ){
            const Read_pt read = frag->getRead(checkR1);
            //The first region processed is the primary, all supsequent are supp
            for(bool bSupplementary : {true,false} ){
                //R1 takes the 0th for primary, and 1st for secondary
                //R2 takes the 1st for primary, and 0th for secondary
                bool regIdx = bSupplementary ? checkR1 : !checkR1;
                const Region_pt & reg = regionOrder[regIdx];
                SQPair_t sqp(reg,read);
                //Skip unaligned reads
                if(!alnMap.count(sqp)){ continue; }
                const Region_pt & mateReg = regionOrder[!regIdx];
                entryVec.push_back(bam_init1());
                bool bCons = ConstructBamEntry(read,reg,bSupplementary,
                                               frag->getRead(!checkR1), mateReg,
                                               alnMap, entryVec.back() );
                if(!bCons) throw std::runtime_error("Cigar failure");
            }
        }
    }
    //Sort the reads by position
    std::sort(  entryVec.begin(),entryVec.end(),
                [] (bam1_t* & a, bam1_t* & b) {
                    return compareBamByPos(a,b) == -1;
                });
    //Iterate over sorted entries and write them
    for(bam1_t * & entry : entryVec){
        int ok = sam_write1(writer,JointHeader,entry);
        if(ok < 0) throw std::runtime_error("Failed to write to " +
                                            std::string(writer->fn));
        bam_destroy1(entry);
    }
    sam_close(writer);
    //Construct an index for the bam file
    std::string fullFName = readDir + '/' + fName;
    int code = sam_index_build(fullFName.c_str(),0);
    if( code != 0 ){
        throw std::runtime_error("Failed to index " + fullFName);
    }
}

//Constructs a consensus sequence for the edge and outputs the host and
//viral sides
//Inputs - an id for the junction
//         - output file streams fro the host and virus
//         - an edge to process
//         - an alignment map
//Output - None, writes sequences to the file streams
void OutputEdgeBP(  int id, std::ofstream & hostOut, std::ofstream & virusOut,
                    const Edge_t & edge, const AlignmentMap_t & alnMap,
                    const ReadPairSet_t & used)
{
    //Build the Table of aligned sequences
    AlignmentTable_t alnTable = BuildAlignmentTable(edge,alnMap);
    std::string consensus = GenerateConsensus(alnTable);
    //Split the consensus and strip off leading and trailing N's
    std::regex rgx("^N+|N+$");
    std::string hostSeq = std::regex_replace(
                            consensus.substr(0,edge.hostRegion->sequence.length()),
                            rgx,"");
    std::string virusSeq = std::regex_replace(
                            consensus.substr(edge.hostRegion->sequence.length()),
                            rgx,"");
    char hostSuffix = (edge.hostRegion->opensLeft()) ? 'R' : 'L';
    char virusSuffix = (edge.virusRegion->opensLeft()) ? 'R' : 'L';
    hostOut << '>' << id << '_' << hostSuffix << '\n' <<
            hostSeq << '\n';
    virusOut << '>' << id << '_' << virusSuffix << '\n' <<
            virusSeq << '\n';
}

//Proceeds from high confidence edges to low confidence edges, ensuring
//each read is used exactly once
//  The edges may be reordered as 
//The output order is determined by constructing a queue from the reads
//  using a branched queue structure
//Inputs - a reference to a vector of edges sorted from most to
//            least confident
//         - a const reference to an alignment map
//         - a string representing the results file
//         - a string representing the reads directory 
//Output - None, prints to outfile
void OutputEdgesByQ(   EdgeVec_t & edgeVec,const AlignmentMap_t & alnMap,
                    const Name2ReadMap_t & readNameMap,
                    const std::string & resFName, const std::string & readDir,
                    const std::string & hostbpFName,
                    const std::string & virusbpFName)
{
    fprintf(stderr,"Outputting Edges ensuring reads support only one edge...\n");
    ReadPairSet_t used;
    std::ofstream out(resFName);
    std::ofstream hbpOut(hostbpFName);
    std::ofstream vbpOut(virusbpFName);
    //Initialize the Edge Queue
    CBranchedEdgeQueue edgeQueue(&alnMap,&used);
    ConstructEdgeQueue(edgeVec,edgeQueue);
    //Pull Elements off the edge queue
    fprintf(stderr,"Processing Edge Queue...\n");
    int nextJunctionID = 0;
    size_t start = edgeQueue.queueSize();
    int pert = 0;
    while(!edgeQueue.empty()){
        const Edge_t & edge = edgeQueue.top();
        if(PassesEffectiveReadCount(edge,&used)){
            //Output to the res file, bam files, and fasta files
            OutputEdgeCall(nextJunctionID,edge,alnMap,out,used);
            OutputEdgeReads(nextJunctionID,edge,alnMap,readDir,used);
            OutputEdgeBP(nextJunctionID,hbpOut,vbpOut,edge,alnMap,used);
            //Update the used Set
            used.insert(edge.getSupport().begin(),edge.getSupport().end());
            nextJunctionID++;
            edgeQueue.pop();
        } else {
            //Eat the top edge and feed its reads to edges beneath it in
            //the queue
            edgeQueue.cannabalize();
        }
        double progress = (start - edgeQueue.queueSize()) / double(start);
        if(1000.0 * progress > pert+1){
            pert = 1000 * progress;
            fprintf(stderr,"Progress: %0.1f%%\r",progress*100.0);
        }
    }
    fprintf(stderr,"Output %d Edges...\n",nextJunctionID);
}

//Given an edge reports wheteher it has enough effective reads
//Inputs - an Edge
//Output - boolean wether it has enough reads
bool PassesEffectiveReadCount(  const Edge_t & edge,
                                const ReadPairSet_t * used)
{
    //size_t count = edge.supportSet.size() + ((edge.nSplit) ? SplitBonus : 0);
    double count = 0;
    bool bSplit = false;
    //null or empty used means we can take the edge at its word
    if(!used || used->empty()){
        count = edge.getSupport().size();
        if(edge.isSplit()) { count += SplitBonus; }
        return (count >= MinimumReads);
    }
    //Used Exists and is non-empty, we have to check the fragments
    for(const auto & pair : edge.getSupportSummaryMap()){
        if(used->count(pair.first)) { continue; }
        count += 1.0;
        if(pair.second.isSplit) {bSplit = true; }
    }
    if(bSplit) count += SplitBonus;
    return (count >= MinimumReads);
}

///Performs all filtering steps on the edge
//  Identifying break point locations
//  Deduplicating reads
//  Removing High insert size reads
//Inputs - an id, used by thread_pool
//         - an edge
//         - an alignment map
void ProcessEdge(int id,Edge_t & edge, const AlignmentMap_t & alnMap){
    FilterHighInsertReads(edge,alnMap);
    FilterSuspiciousReads(edge,alnMap);
    DeduplicateEdge(edge,alnMap);
    //Explicit call to ensure the offsets are ready when needed
    edge.getOffsets();
}

void ProcessEdges(EdgeVec_t & edgeVec, const AlignmentMap_t & alnMap){
    fprintf(stderr,"Processing %zu Edges ...\n",edgeVec.size());
    ctpl::thread_pool threadPool (Config.threads);
    std::vector<std::future<void>> futureVec;
    for( Edge_t & edge : edgeVec){
        auto future = threadPool.push(        ProcessEdge,std::ref(edge),
                                        std::cref(alnMap));
        futureVec.push_back(std::move(future));
    }
    int pert = 0;
    size_t complete = 0;
    for (auto & future : futureVec){
        future.get();
        complete++;
        double progress = complete / double(futureVec.size());
        if(1000.0 * progress > pert+1){
            pert = 1000 * progress;
            fprintf(stderr,"Progress: %0.1f%%\r",progress*100.0);
        }
    }
    FilterEdgeVec(edgeVec);
    fprintf(stderr,"Processed and retained %zu Edges\n",edgeVec.size());
}

////Recursivly processes prepared data describing the sequences of an edge
////First a consensus sequence is generated for the edge
////then all reads which are too different from the consensus (adjusted for
////the length of the read) are identified
////These reads are removed fromt he parent edge and moved to a child edge
////Additionally any reads from the parent edge which are completely
////consistent with any of the discarded reads are also included
////Overall this allows the potential for multiple alleles of junctions
////Any child edges with too fiew reads are ignored
////Continues until no valid child edge is made
////Inputs - an edge to process
////         - a vector of the reads in the edge
////         - a vector of the aligned sequences of the reads
////         - a vector of the aligned length of the reads
////         - a reference to a vector of edges in which to store new edges
////Output - None, modifies all inputs
//EdgeVec_t RecursiveSplitEdge(Edge_t & edge, std::vector<ReadPair_pt> rowLabelVec,
//                        std::vector<std::string> rowSeqVec,
//                        std::vector<size_t> nFillVec)
//{
//    size_t nRowIn = rowLabelVec.size();
//    std::vector<size_t> diffCount;
//    GenerateConsensus(rowSeqVec,&diffCount);
//    std::unique_ptr<Edge_t> newEdge_p(nullptr);
//    std::unordered_set<size_t> roiSet;
//    for(size_t a = 0; a < rowSeqVec.size(); a++){
//        //Calculate the # of diffs per defined site
//        double diffRate = double(diffCount[a]) / double(nFillVec[a]);
//        if(diffRate < MaxDiffRate) continue;
//        if(!newEdge_p){
//            newEdge_p = std::make_unique<Edge_t>(edge.hostRegion,edge.virusRegion);
//            newEdge_p->id = edge.id;
//        }
//        //Move the fragment to the new edge
//        edge.transferSupport(rowLabelVec[a],*newEdge_p);
//        roiSet.insert(a);
//        //Any reads consistent with this read will be included
//        for(size_t b = a + 1; b < rowSeqVec.size(); b++){
//            if(!IsConsistent(rowSeqVec[a],rowSeqVec[b])) continue;
//            edge.transferSupport(rowLabelVec[b],*newEdge_p,false);
//            roiSet.insert(b);
//        }
//    }
//    //We are done if no new edge was created
//    if(!newEdge_p) return EdgeVec_t();
//    //We are also done if the new edge is too small
//    if(!PassesEffectiveReadCount(*newEdge_p)){
//        return EdgeVec_t();
//    }
//    //Reduce the vectors to only the rows of interest for the new edge
//    FilterVector(rowLabelVec,roiSet); 
//    FilterVector(rowSeqVec,roiSet); 
//    FilterVector(nFillVec,roiSet); 
//    //Prevent infinite recursion by requiring that the recursion stops if
//    //the next round isn't smaller
//    if(rowLabelVec.size() >= nRowIn) return EdgeVec_t(1,*newEdge_p);
//    EdgeVec_t res = RecursiveSplitEdge(*newEdge_p,rowLabelVec,rowSeqVec,nFillVec);
//    //Check if the splitting process left the created edge large enough
//    if(PassesEffectiveReadCount(*newEdge_p)){
//        res.insert(res.begin(),*newEdge_p);
//    }
//    return res;
//}

//Recursivly processes prepared data describing the sequences of an edge
//First a consensus sequence is generated for the edge
//then all reads which are too different from the consensus (adjusted for
//the length of the read) are identified
//These reads are removed fromt he parent edge and moved to a child edge
//Additionally any reads from the parent edge which are completely
//consistent with any of the discarded reads are also included
//Overall this allows the potential for multiple alleles of junctions
//Any child edges with too fiew reads are ignored
//Continues until no valid child edge is made
//Inputs - an edge to process
//         - a vector of the reads in the edge
//         - a vector of the aligned sequences of the reads
//         - a vector of the aligned length of the reads
//         - a reference to a vector of edges in which to store new edges
//Output - None, modifies all inputs
template<class T>
EdgeVec_t RecursiveSplitEdge(   Edge_t & edge,
                                std::vector<T> vec,
                                std::vector<bool> (*globalTest)(const std::vector<T> &),
                                bool (*pairwiseTest)(const T &, const T &)
                                )

{
    size_t nRowIn = vec.size();
    std::vector<bool> rowConsisVec = globalTest(vec);
    std::unique_ptr<Edge_t> newEdge_p(nullptr);
    std::unordered_set<size_t> roiSet;
    for(size_t a = 0; a < vec.size(); a++){
        const T & rowA = vec[a];
        bool bConsistent = rowConsisVec[a];
        if(bConsistent) continue;
        if(!newEdge_p){
            newEdge_p = std::make_unique<Edge_t>(edge.hostRegion,edge.virusRegion);
            newEdge_p->id = edge.id;
        }
        //Move the fragment to the new edge
        edge.transferSupport(rowA.label,*newEdge_p);
        roiSet.insert(a);
        //Any reads consistent with this read will be included
        for(size_t b = a + 1; b < vec.size(); b++){
            const T & rowB = vec[b];
            if(!pairwiseTest(rowA,rowB)) continue;
            edge.shareSupport(rowB.label,*newEdge_p);
            roiSet.insert(b);
        }
    }
    //We are done if no new edge was created
    if(!newEdge_p) return EdgeVec_t();
    //We are also done if the new edge is too small
    if(!PassesEffectiveReadCount(*newEdge_p)){
        return EdgeVec_t();
    }
    //Reduce the vectors to only the rows of interest for the new edge
    FilterVector(vec,roiSet); 
    //Prevent infinite recursion by requiring that the recursion stops if
    //the next round isn't smaller
    if(vec.size() >= nRowIn) return EdgeVec_t(1,*newEdge_p);
    EdgeVec_t res = RecursiveSplitEdge(*newEdge_p,vec,globalTest,pairwiseTest);
    //Check if the splitting process left the created edge large enough
    if(PassesEffectiveReadCount(*newEdge_p)){
        res.insert(res.begin(),*newEdge_p);
    }
    return res;
}

//Some reads may have failed during alignment, remove them from the edges
//Inputs - a vector of edges to modify
//         - an alignment map to check
//Output - None, modifes the edge vector
void RemoveUnalignedReads(EdgeVec_t & edgeVec,const AlignmentMap_t & alnMap){
    fprintf(stderr,"Removing Unaligned reads from %zu edges...\n",edgeVec.size());
    for( Edge_t & edge : edgeVec){
        std::vector<ReadPair_pt> toRemoveVec;
        for( const auto & pair : edge.getSupportSummaryMap()){
            const ReadPair_pt & frag = pair.first;
            const ReadPairAlnSummary_t & summary = pair.second;
            //A paired alignment passes if at least one alignment exists for each of the following
            //  R1, R2, Host, Virus
            bool bPass =
                (summary.hasAlnFlag & ReadPairAlnSummary_t::HAS_R1) &&
                (summary.hasAlnFlag & ReadPairAlnSummary_t::HAS_R2) &&
                (summary.hasAlnFlag & ReadPairAlnSummary_t::HAS_HOST) &&
                (summary.hasAlnFlag & ReadPairAlnSummary_t::HAS_VIRUS);
            if(!bPass){
                toRemoveVec.push_back(frag);
            }
        }
        for( const ReadPair_pt & frag : toRemoveVec){
            edge.removeSupport(frag);
        }
    }
    FilterEdgeVec(edgeVec);
    fprintf(stderr,"Edges remaining: %zu\n",edgeVec.size());
}


void SortEdgeVec(   EdgeVec_t & edgeVec, const AlignmentMap_t & alnMap,
                    const ReadPairSet_t & used) {
    std::sort(  edgeVec.begin(), edgeVec.end(),
                [&alnMap,&used](Edge_t & a, Edge_t & b){
                    //Sort in descending order (we'll process from the back)
                    double aScore = a.score(alnMap,used);
                    double bScore = b.score(alnMap,used);
                    if(aScore != bScore){
                        return aScore < bScore;
                    }
                    return false;
                    //return a.uniqueReadSet.size() < b.uniqueReadSet.size();
                });
}

//Given a vector of edges, in parallel splits each into edges for
//each consensus sequence present
//Inputs - a vector of edges
//         - an alignment map
//Output - a vector of new edges, also modifies the edges in the input vecto
//
EdgeVec_t SplitEdges(   EdgeVec_t & edgeVec, const AlignmentMap_t & alnMap,
                        EdgeVec_t (*edgeSplitter)(Edge_t &, const AlignmentMap_t &))
{
    fprintf(stderr,"Splitting Edges based on consensus sequences ...\n");
    //Multithreaded
    ctpl::thread_pool threadPool (Config.threads);
    std::vector<std::future<EdgeVec_t>> futureVec;
    for( Edge_t & edge : edgeVec){
        auto future = threadPool.push(
                [edgeSplitter,&edge,&alnMap] (int id) {
                    return edgeSplitter(edge,alnMap);
                } );
        futureVec.push_back(std::move(future));
    }
    EdgeVec_t newEdges;
    for(auto & future : futureVec){
        EdgeVec_t localNew = future.get();
        newEdges.insert(newEdges.end(),localNew.begin(),localNew.end());
    }
    ////Single Threaded version
    //EdgeVec_t newEdges;
    //for( Edge_t & edge : edgeVec){
    //    EdgeVec_t localNew = ConsensusSplitEdge(1,edge,alnMap);
    //    newEdges.insert(newEdges.end(),localNew.begin(),localNew.end());
    //}
    fprintf(stderr,"Identified %zu new Edges ...\n", newEdges.size());
    return newEdges;
}

//Specialized function for testing whether all rows in an aignment table are consistent with
// the consensus of that alignment table
//Inputs - a vector of alignment table rows (cref)
//Output - a vector of boolean results (one for each row), true if consistent 
std::vector<bool> TestConsistencyGlobal (
                    const AlignmentTable_t & alnTable)
{
    std::vector<bool> resVec;
    std::vector<size_t> diffCount;
    GenerateConsensus(alnTable,&diffCount);
    for(size_t a = 0; a < alnTable.size(); a++){
        const AlignmentTableRow_t & rowA = alnTable[a];
        //Calculate the # of diffs per defined site
        double diffRate = double(diffCount[a]) / double(rowA.nFill);
        bool bRes = (diffRate < MaxDiffRate) ? true : false;
        resVec.push_back(bRes); 
    }
    return resVec;
}


//Given two partial DNA sequences tests if the two have consistent
//sequences: that is they match at all non-N positions
//Inputs - two strings representing the two sequences
//Output - a boolean of whether they are consistent or not
bool TestConsistencyPairwise(const AlignmentTableRow_t & a, const AlignmentTableRow_t & b){
    const std::string & seq1 = a.seq;
    const std::string & seq2 = b.seq;
    if(seq1.length() != seq2.length()) return false; 
    for(size_t i = 0; i < seq1.length(); i++){
        char c1 = seq1.at(i);
        char c2 = seq2.at(i);
        if(c1 != c2 && c1 != 'N' && c1 != 'N') return false;
    }
    return true;
}






