#ifndef SURVERYOR_EDGE_UTILS_H
#define SURVERYOR_EDGE_UTILS_H

#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <htslib/kseq.h>
#include <list>
#include <vector>

#include "utils.h"
#include "sam_utils.h"
#include <ssw.h>
#include <ssw_cpp.h>


struct CannonicalSeq_t {
    protected:
    std::string cannon;
    std::string noncannon;
    public:
    CannonicalSeq_t(std::string s) : 
        cannon(s), noncannon(get_seqrc(s))
    {
        if(!is_canonical(s)){
            std::swap(cannon,noncannon);
        }
    }
    const std::string & get(bool cannonical = true) const {
        return (cannonical) ? cannon : noncannon;
    }
    //Determines if a DNA sequence is in canonical form (the lexicographically earlier sequence)
    static bool is_canonical(const std::string & seq) {
        int cmp = 0;
        auto it = seq.begin();
        auto rit = seq.rbegin();
        while(&(*it) != &(*rit) && cmp == 0){
            char r = 'N';
            switch(*rit) {
                case 'A': r = 'T'; break;
                case 'C': r = 'G'; break;
                case 'T': r = 'A'; break;
                case 'G': r = 'C'; break;
            }
            if(*it != r) { cmp = (*it < r) ? -1 : 1; }
            it++;
            rit++;
        }
        return cmp <= 0;
    }

};

struct Read_t;

typedef std::shared_ptr<Read_t> Read_pt;

struct Read_t {
    Read_t(std::string nm, std::string s, bool isr1) :
        name(nm), seq(to_upper(s)), isR1(isr1) {}
    Read_t(bam1_t* aln) : Read_t(   bam_get_qname(aln),
                                    get_sequence(aln),
                                    aln->core.flag & BAM_FREAD1) {}
    const std::string name;
    const CannonicalSeq_t seq;
    const bool isR1;
    int compare(const Read_t & other) const {
        if(this->name != other.name) {
            return (this->name < other.name) ? -1 : 1;
        }
        if(this->isR1 != other.isR1) {
            return (this->isR1) ? -1 : 1;
        }
        if(this->seq.get() != other.seq.get()){
            return (this->seq.get() < other.seq.get()) ? -1 : 1;
        }
        return 0;
    }
    std::string to_string(bool seq = false) const {
        std::string str;
            str =   '-' + name + " R" + ((isR1) ? '1' : '2') +
                    this->seq.get(); 
        return str;
    }
};

struct Read_pt_EqFunctor {
    bool operator()(const Read_pt & a, const Read_pt & b) const {
        return a->compare(*b) == 0;
    }
};

struct Read_pt_HashFunctor {
    size_t operator()(const Read_pt & a) const {
        return std::hash<std::string>{}(a->to_string());
    }
};


//Object describing a read
//struct Read_t {
//    Read_t(std::string nm, bool bSplit, bool bViralR1) :
//        name(nm), isSplit(bSplit), viralR1(bViralR1) {}
//    std::string name;
//    std::string hostSegment;
//    std::string hostRC;
//    std::string virusSegment;
//    std::string virusRC;
//    bool isSplit = false;
//    bool viralR1 = false;
//    Read_pt mate = nullptr;
//    int compare(const Read_t & other) const {
//        if(this->name != other.name) {
//            return (this->name < other.name) ? -1 : 1;
//        }
//        if(this->hostSegment.size() != other.hostSegment.size()){
//            return (this->hostSegment.size() < other.hostSegment.size())
//                    ? -1 : 1;
//        }
//        if(this->virusSegment.size() != other.virusSegment.size()){
//            return (this->virusSegment.size() < other.virusSegment.size())
//                    ? -1 : 1;
//        }
//        if(this->isSplit != other.isSplit){
//            return (other.isSplit) ? -1 : 1;
//        }
//        if(this->viralR1 != other.viralR1){
//            return (this->viralR1) ? -1 : 1;
//        }
//        return 0;
//    }
//    std::string to_string(bool seq = false) const {
//        std::string str;
//        if(!seq){
//            str =   name + '-' + std::to_string(hostSegment.size()) + '-' + 
//                    std::to_string(virusSegment.size()) + '-' +
//                    std::to_string(isSplit) + std::to_string(viralR1);
//        } else {
//            if(hostSegment.size()){
//                str += std::to_string('>') + name + '_' + std::to_string(1);
//                str += std::to_string('\n') + hostSegment;
//            }
//            if(virusSegment.size()){
//                str += std::to_string('>') + name + '_' + std::to_string(2);
//                str += std::to_string('\n') + virusSegment;
//            }
//        }
//        return str;
//    }
//    const std::string & getSegment(bool bVirus, bool bRC) const {
//        if(bVirus){
//           return (bRC) ? this->virusRC : this->virusSegment;
//        } else {
//           return (bRC) ? this->hostRC : this->hostSegment;
//        }
//    }
//};
//
//struct Read_pt_EqFunctor {
//    bool operator()(const Read_pt & a, const Read_pt & b) const {
//        return a->compare(*b) == 0;
//    }
//};
//
//struct Read_pt_HashFunctor {
//    size_t operator()(const Read_pt & a) const {
//        return std::hash<std::string>{}(a->to_string());
//    }
//};

struct ReadPair_t {
    ReadPair_t() : R1(nullptr), R2(nullptr) {}
    Read_pt R1;
    Read_pt R2;
    Read_pt & operator[](bool bR1){ return (bR1) ? R1 : R2; }
    Read_pt & getRead(bool bR1){ return (bR1) ? R1 : R2; }
};


typedef std::shared_ptr<ReadPair_t> ReadPair_pt;
typedef std::unordered_map<std::string,ReadPair_pt> Name2ReadPairMap_t;
typedef std::unordered_set<ReadPair_pt> ReadPairSet_t;
typedef std::unordered_map<std::string,Read_pt> Name2ReadMap_t;

//Object describing a genomic location containing a
//candidate junction
struct Region_t {
    enum REGION_INFO_BITS {
        IS_VIRAL = 0x1,
        OPENS_LEFT = 0x2
    };
    const long int id;
    const std::string chromosome;
    const std::string sequence;
    const size_t offset;
    const size_t end;
    const uint8_t flag;
    Region_t(   long int id, const std::string & chr, const std::string &seq,
                size_t o, size_t e, bool bV, bool oL) :
        id(id), chromosome(chr), sequence(to_upper(seq)),
        offset(o), end(e), flag(flag_from_bool(bV,oL))
    {}
    //{
    //    std::vector<std::string> fields = strsplit(name,',');
    //    std::vector<std::string> coords = strsplit(coordStr,'-');
    //    this->chr = fields[0];
    //    this->left = std::stoul(fields[1]);
    //    this->right = std::stoul(fields[2]);
    //    this->seqLeft = std::stoul(coords[0]);
    //    this->seqRight = std::stoul(coords[1]);
    //    this->strand = fields[3][0];
    //    for (auto & c: this->sequence) c = (char)toupper(c);
    //}
    //Region_t(kseq_t *seq, bool bVirus,const std::string & coordStr) : 
    //    Region_t(std::string(seq->name.s),std::string(seq->seq.s),bVirus,coordStr)
    //{}
    //Region_t(const Region_t & other) :
    //    left(other.left), right(other.right),
    //    seqLeft(other.seqLeft), seqRight(other.seqRight),
    //    chr(other.chr), strand(other.strand), sequence(other.sequence) {}
    static uint8_t flag_from_bool(bool bVirus, bool opensLeft) {
        uint8_t flag = 0;
        if(bVirus) { flag |= IS_VIRAL; }
        if(opensLeft) { flag |= OPENS_LEFT; }
        return flag;
    }
    //size_t left, right; //These are labels defining the region
    ////These are the actual genomic corrdinates of the sequence associated with this region
    //size_t seqLeft, seqRight; 
    //std::string chr;
    //char strand;
    //std::string sequence;
    //bool isVirus;
    int compare(const Region_t other) const {
        if(this->chromosome != other.chromosome){ 
            return (this->chromosome < other.chromosome) ? -1 : 1;
        }
        //host less than viral, opens left less than opens right
        if(this->flag != other.flag){
            return (this->flag < other.flag) ? -1 : 1;
        }
        if(this->offset != other.offset){
            return (this->offset < other.offset) ? -1 : 1;
        }
        if(this->end != other.end) {
            return (this->end < other.end) ? -1 : 1;
        }
        return 0;
    }
    bool opensLeft() const { return flag & OPENS_LEFT; }
    bool isViral() const { return flag & IS_VIRAL; }
    char strand() const {
        return (bool(flag & OPENS_LEFT) == bool(flag & IS_VIRAL)) ? '+' : '-';
    }
    std::string to_string(bool seq=false) const {
        std::string str=    '>' + chromosome + ":(" + 
                            std::to_string(int(flag)) + ")" +
                            std::to_string(offset+1) + "-" +
                            std::to_string(end);
        if(seq) str += '\n' + sequence;
        return str;
    }
};

////Object describing a genomic location containing a
////candidate junction
//struct Region_t {
//    Region_t(   const std::string & name, const std::string &sequence,
//                bool bVirus, const std::string & coordStr) :
//        sequence(sequence), isVirus(bVirus)
//    {
//        std::vector<std::string> fields = strsplit(name,',');
//        std::vector<std::string> coords = strsplit(coordStr,'-');
//        this->chr = fields[0];
//        this->left = std::stoul(fields[1]);
//        this->right = std::stoul(fields[2]);
//        this->seqLeft = std::stoul(coords[0]);
//        this->seqRight = std::stoul(coords[1]);
//        this->strand = fields[3][0];
//        for (auto & c: this->sequence) c = (char)toupper(c);
//    }
//    Region_t(kseq_t *seq, bool bVirus,const std::string & coordStr) : 
//        Region_t(std::string(seq->name.s),std::string(seq->seq.s),bVirus,coordStr)
//    {}
//    Region_t(const Region_t & other) :
//        left(other.left), right(other.right),
//        seqLeft(other.seqLeft), seqRight(other.seqRight),
//        chr(other.chr), strand(other.strand), sequence(other.sequence) {}
//    size_t left, right; //These are labels defining the region
//    //These are the actual genomic corrdinates of the sequence associated with this region
//    size_t seqLeft, seqRight; 
//    std::string chr;
//    char strand;
//    std::string sequence;
//    bool isVirus;
//    int compare(const Region_t other){
//        if(this->chr != other.chr){ 
//            return (this->chr < other.chr) ? -1 : 1;
//        }
//        if(this->strand != other.strand){
//            return (this->strand == '-') ? -1 : 1;
//        }
//        if(this->left != other.left){
//            return (this->left < other.left) ? -1 : 1;
//        }
//        if(this->right != other.right) {
//            return (this->right < other.right) ? -1 : 1;
//        }
//        return 0;
//    }
//    std::string to_string(bool seq=false) const {
//        std::string str=    '>' + chr + ":" + strand +
//                            std::to_string(left) + "-" +
//                            std::to_string(right);
//        if(seq) str += sequence;
//        return str;
//    }
//};

typedef std::shared_ptr<Region_t> Region_pt;
typedef std::unordered_map<std::string,Region_pt> Name2RegionMap_t;
typedef std::unordered_map<long int,Region_pt> RegID2RegionMap_t;

struct Region_pt_EqFunctor {
    bool operator()(const Region_pt & a, const Region_pt & b) const{
        return a->compare(*b) == 0;
    }
};

struct Region_pt_HashFunctor {
    size_t operator()(const Region_pt & a) const {
        return std::hash<std::string>{}(a->to_string());
    }
};



//Sets of shared pointers to reads and regions
typedef std::unordered_set<Read_pt,Read_pt_HashFunctor,Read_pt_EqFunctor>
            ReadSet_t;
typedef std::unordered_set<Region_pt,Region_pt_HashFunctor,Region_pt_EqFunctor>
            RegionSet_t;


typedef std::unordered_map< Read_pt,RegionSet_t,
                            Read_pt_HashFunctor,Read_pt_EqFunctor>
            Read2RegionsMap_t;
typedef std::unordered_map< Region_pt, ReadSet_t,
                            Region_pt_HashFunctor,Region_pt_EqFunctor>
            Region2ReadsMap_t;
typedef std::unordered_map< Read_pt, ReadSet_t,
                            Read_pt_HashFunctor,Read_pt_EqFunctor>
            Read2ReadsMap_t;

struct Edge_t; 

typedef std::vector<Edge_t> EdgeVec_t;
typedef std::list<Edge_t> EdgeList_t;

//Object associating a query read with a subject region
struct SQPair_t {
    SQPair_t(const Region_pt & s, const Read_pt & q) : subject(s), query(q) {}
    SQPair_t(const SQPair_t & other) :
        subject(other.subject), query(other.query) {}
    Region_pt subject;
    Read_pt query;
    int compare(const SQPair_t & other) const {
        int res = this->subject->compare(*other.subject);
        if(res != 0) return res; 
        res = this->query->compare(*other.query);
        if(res != 0) return res;
        return 0;
    }
    std::string to_string() const {
        return        this->subject->to_string() + " vs " +
                this->query->to_string();
    }
};

struct SQPair_EqFunctor {
    bool operator()(const SQPair_t & a, const SQPair_t & b) const {
        return a.compare(b) == 0;
    }
};

struct SQPair_HashFunctor {
    size_t operator()(const SQPair_t &a) const {
        return std::hash<std::string>{}(a.to_string());
    }
};

//Mapping from a subject query pair to the resulting alignment
typedef std::unordered_map< SQPair_t,StripedSmithWaterman::Alignment,
                            SQPair_HashFunctor,SQPair_EqFunctor>
            AlignmentMap_t;


struct ReadPairAlnSummary_t {
    ReadPairAlnSummary_t() :
        isSplit(false) , hostScore(0.0), virusScore(0.0),
        hostLeft(-1), hostRight(-1), virusLeft(-1), virusRight(-1),
        hostQAlnBases(-1), virusQAlnBases(-1)
    {}
    bool isSplit;
    double hostScore, virusScore;
    int32_t hostLeft, hostRight, virusLeft, virusRight;
    int32_t hostQAlnBases, virusQAlnBases;
    double score() const { return hostScore + virusScore; }
    int32_t calcIS() const {
        return (hostRight - hostLeft) + (virusRight - virusLeft);
    }
};
typedef std::unordered_map<ReadPair_pt,ReadPairAlnSummary_t> ReadPairAlnSummaryMap_t;

//Object associating a pair of regions and the reads spanning the pair
struct Edge_t {
    //#TODO: Have support have a way to group duplicates together
    static size_t MinimumClipLen;
    public:
    long int id;
    Region_pt hostRegion;
    Region_pt virusRegion;
    protected:
    ReadPairSet_t supportSet;
    ReadPairAlnSummaryMap_t supportAlnSummaryMap;
    bool validOffsets = false;
    size_t hostOffset;
    size_t virusOffset;
    size_t nSplit = 0;
    //double lastScore = -1;
    double lastScore = 0;
    double m_cachedScore = -1;
    public:
    Edge_t() : Edge_t(nullptr,nullptr) {}
    //Edge_t(const std::string & regStr, const std::string & readStr);
    Edge_t(Region_pt hostReg, Region_pt virusReg) :
            id(-1),
            hostRegion(hostReg), virusRegion(virusReg), supportSet(),
            supportAlnSummaryMap(), validOffsets(false),
            hostOffset(0), virusOffset(0), nSplit(0), lastScore(0),
            m_cachedScore(-1)
    {}
    Edge_t(const Edge_t & other) :
        id(other.id),
        hostRegion(other.hostRegion), virusRegion(other.virusRegion),
        supportSet(other.supportSet),
        supportAlnSummaryMap(other.supportAlnSummaryMap),
        validOffsets(other.validOffsets),
        hostOffset(other.hostOffset), virusOffset(other.virusOffset), 
        nSplit(other.nSplit), lastScore(other.lastScore),
        m_cachedScore(other.m_cachedScore)
    {}
    public:
    const ReadPairSet_t & getSupport() const { return this->supportSet; }
    std::pair<size_t,size_t> getOffsets() {
        return (validOffsets) ? std::make_pair(hostOffset,virusOffset) :
                                determineOffsets();
    }
    std::pair<size_t,size_t> getOffsets() const {
        if(!validOffsets) {
            throw std::logic_error("Call for const getOffsets prior to determineOffsets");
        }
        return std::make_pair(hostOffset,virusOffset);
    }
    const ReadPairAlnSummaryMap_t & getSupportSummaryMap() const {
        return this->supportAlnSummaryMap;
    }
    bool addSupport(const ReadPair_pt & frag, ReadPairAlnSummary_t summary) {
        auto res = this->supportSet.insert(frag);
        if(!res.second){ return false; }
        validOffsets = false;
        if(summary.isSplit) nSplit++;
        this->lastScore += summary.score();
        this->supportAlnSummaryMap[frag] = summary;
        m_cachedScore = -1;
        return true;
    }
    bool addSupport(const ReadPair_pt & frag, const AlignmentMap_t alnMap){
        return addSupport(frag,this->getRPAlnSummary(frag,alnMap));
    }
    bool transferSupport(const ReadPair_pt & frag, Edge_t & other,bool bRemove = true){
        //Cannot transfer support an edge does not have
        if(!this->supportSet.count(frag)) { return false; }
        //Cannot transfer support to an edge with different regions
        if( (this->hostRegion->compare(*other.hostRegion) != 0) ||
            (this->virusRegion->compare(*other.virusRegion) != 0) )
        {
            return false;
        }
        //
        //Add the support from this fragment to the new edge
        bool retVal = other.addSupport(frag,this->supportAlnSummaryMap[frag]);
        if(bRemove) {
            retVal = retVal && (this->removeSupport(frag));
        }
        return retVal;
    }
    protected:
    //Iterates over supporting fragments and identifies the most junction proximal
    //  position observed; This is cached for the future;
    std::pair<size_t,size_t> determineOffsets() {
        //Iterate over reads to find the extremes
        for ( const Region_pt & reg : {hostRegion, virusRegion} ) { 
            int32_t minLeft = reg->sequence.length();
            int32_t maxRight = 0;
            for( const auto & pair :  supportAlnSummaryMap) {
                const ReadPairAlnSummary_t & summary = pair.second;
                //Get the left and right positions of the alignment,
                // and update the overall positions for the pair
                const int32_t * left_ptr = (reg == this->hostRegion) ?
                                        &(summary.hostLeft) :
                                        &(summary.virusLeft);
                const int32_t * right_ptr = (reg == this->hostRegion) ?
                                        &(summary.hostRight) :
                                        &(summary.virusRight);
                if(*left_ptr < minLeft) { minLeft = *left_ptr; } 
                if(*right_ptr > maxRight) { maxRight = *right_ptr; } 
            }
            //Set the appropriate offset
            auto * offset_ptr = (reg == hostRegion) ?   &(hostOffset) :
                                                        &(virusOffset);
            *offset_ptr = (reg->opensLeft()) ? minLeft : maxRight;
        }
        validOffsets=true;
        return std::make_pair(hostOffset,virusOffset);
    }
    //Given alignment information across all read-region combis,
    //a summary is extracted for the pair of reads mapped to this edges' regions
    //This is cached for future reference
    ReadPairAlnSummary_t getRPAlnSummary(   const ReadPair_pt & frag,
                                            const AlignmentMap_t alnMap)
    {
        ReadPairAlnSummary_t summary;
        if(!this->hostRegion || !this->virusRegion) { return summary; }
        summary.hostLeft = this->hostRegion->sequence.length();
        summary.hostRight = 0;
        summary.virusLeft = this->virusRegion->sequence.length();
        summary.virusRight = 0;
        summary.hostQAlnBases = 0;
        summary.virusQAlnBases = 0;
        //Check each combination of Read vs Region and calculate the total score
        //  as well as whether the alignment is split
        for ( bool checkR1 : {true, false} ){
            int splitCount = 0;
            for ( const Region_pt & curReg :
                    {this->hostRegion, this->virusRegion})
            {
                SQPair_t pair(curReg,frag->getRead(checkR1));
                //Skip read-region pairs with no alignment
                if(!alnMap.count(pair)) { continue; }
                const StripedSmithWaterman::Alignment & aln =
                    alnMap.at(pair);
                int opIdx = curReg->opensLeft() ? 0 : aln.cigar.size()-1;
                uint32_t c = aln.cigar[opIdx];
                if( cigar_int_to_op(c) == 'S' &&
                    cigar_int_to_len(c) >= Edge_t::MinimumClipLen)
                {
                    splitCount++;
                }
                double * score_ptr =    (curReg->isViral()) ?
                                        &(summary.virusScore) :
                                        &(summary.hostScore);
                *score_ptr += aln.sw_score;
                //Get the left and right positions of the alignment,
                // and update the overall positions for the pair
                int32_t * left_ptr = (curReg == this->hostRegion) ?
                                        &(summary.hostLeft) :
                                        &(summary.virusLeft);
                int32_t * right_ptr = (curReg == this->virusRegion) ?
                                        &(summary.hostRight) :
                                        &(summary.virusRight);
                if(aln.ref_begin < *left_ptr) { *left_ptr = aln.ref_begin; } 
                if(aln.ref_end > *right_ptr) { *right_ptr = aln.ref_end; } 
                //Total the number of query bases mapped to each region
                int32_t * qAlnBasesPtr =    curReg->isViral() ? 
                                            &(summary.virusQAlnBases) :
                                            &(summary.hostQAlnBases);
                *qAlnBasesPtr += aln.query_end - aln.query_begin + 1;
            }
            //For a read to be split it must have a split alignment to both the
            // host and viral regions
            if(splitCount == 2) {
                summary.isSplit = true;
            }
        }
        return summary;
    }
    public:
    size_t splitCount() const { return this->nSplit; }
    bool isSplit() const { return this->nSplit > 0; }
    //Retained for backwards compatibility
    double cachedScore(   const AlignmentMap_t & alnMap,
                    const ReadPairSet_t & used)
    {
        if(this->m_cachedScore == -1) {
            this->m_cachedScore = this->score(alnMap,used);
        }
        return this->m_cachedScore;
    }
    bool removeSupport(const ReadPair_pt & frag){
        //if(this->readSet.empty()) return false;
        if(!this->supportSet.erase(frag)) {
            return false;
        }
        ReadPairAlnSummary_t summary = this->supportAlnSummaryMap.at(frag);
        this->supportAlnSummaryMap.erase(frag);
        if(summary.isSplit && nSplit) nSplit--;
        this->lastScore -= summary.score();
        validOffsets = false;
        m_cachedScore = -1;
        return true;
    }
    //alnMap Retained for backwards compatibility
    double score(   const AlignmentMap_t & alnMap,
                    const ReadPairSet_t & used) const
    {
        if(used.empty()) { return lastScore; }
        double my_score = 0;
        for(const auto & pair: this->supportAlnSummaryMap ) {
            const ReadPair_pt & frag = pair.first;
            if(used.count(frag)) { continue; }
            my_score += pair.second.score();
            //for ( bool checkR1 : {true, false} ){
            //    for ( const Region_pt & curReg :
            //            {this->hostRegion, this->virusRegion})
            //    {
            //        SQPair_t pair(curReg,frag->getRead(checkR1));
            //        const StripedSmithWaterman::Alignment & aln =
            //            alnMap.at(pair);
            //        my_score += aln.sw_score;
            //    }
            //}
        }
        return my_score;
    }
    //private:
    //void parseRegString(const std::string & regStr);
    //void parseReadString(const std::string & readStr);
};


size_t Edge_t::MinimumClipLen = 20;

#endif //SURVERYOR_EDGE_UTILS_H



