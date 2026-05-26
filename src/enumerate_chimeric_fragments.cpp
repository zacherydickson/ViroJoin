#include <iostream>
#include <list>
#include <memory>
#include <vector>
#include <array>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <stdexcept>

#include "sam_utils.h"
#include "str_utils.h"
#include "ChimericFragment.h"
#include "config.h"
#include <cptl_stl.h>
#include "utils.h"
#include <mutex>

std::mutex mu;

class ClippedCXA : public CXA {
    public:
        enum FLAG_BITS {
            IS_CLIPPED = 0x1,
            IS_LEFT = 0x2,
            IS_R1 = 0x4,
            PRIMARY_IS_REVERSED=0x8
        };
        uint8_t flag;
        ClippedCXA(std::string xaStr,uint8_t flag) : CXA(xaStr), flag(flag) {}
        ClippedCXA(const CXA & other,uint8_t flag) : CXA(other), flag(flag) {}
        bool is_clipped() const { return flag & IS_CLIPPED; }
        bool is_left() const { return flag & IS_LEFT; }
        bool is_r1() const { return flag & IS_R1; }
        bool primary_is_rev() const { return flag & PRIMARY_IS_REVERSED; }
        std::string to_string() const {
            return CXA::to_string() + " " + std::to_string(flag); 
        }
};

typedef std::shared_ptr<ClippedCXA> ClippedCXA_spt;

struct MateAlignmentInfo_t {
    enum ALN_PARTS {
        UNSPLIT, ANCHOR_LEFT, ANCHOR_RIGHT, CLIP_LEFT, CLIP_RIGHT
    };
    static const int NPartTypes = 5;
    protected:
    ClippedCXA_spt parts[NPartTypes];
    public:
    MateAlignmentInfo_t() {
        for(int i = 0; i < NPartTypes; i++){
            parts[i] = nullptr;
        }
    }
    MateAlignmentInfo_t(const MateAlignmentInfo_t & other) {
        for(int i = 0; i < NPartTypes; i++){
            this->parts[i] = other.parts[i];
        }
    }
    static bool Part_Is_Anchor(ALN_PARTS partType) {
        return partType == ANCHOR_LEFT || partType == ANCHOR_RIGHT ;
    }
    static bool Part_Is_Clip(ALN_PARTS partType) {
        return partType == CLIP_LEFT || partType == CLIP_RIGHT;
    }
    static bool Part_Is_Left(ALN_PARTS partType) {
        return partType == ANCHOR_LEFT || partType == CLIP_LEFT;
    }
    static bool Part_Is_Right(ALN_PARTS partType) {
        return partType == ANCHOR_RIGHT || partType == CLIP_RIGHT;
    }
    
    bool add(ClippedCXA_spt part, ALN_PARTS partType) {
        //Cannot add part if it already exists
        if(this->parts[partType]) { return false; }
        if(partType == UNSPLIT) { return do_add(part,partType); }
        //Part type is split
        //Cannot add split parts if an unsplit part is already present
        if(this->parts[UNSPLIT]) { return false; }
        //Cannot have more than one left part
        if(this->has_left() && Part_Is_Left(partType)) { return false;}
        //Cannot have more than one right part
        if(this->has_right() && Part_Is_Right(partType)) { return false;}
        if(Part_Is_Anchor(partType)) {
            //Cannot have more than one anchor
            if(this->has_anchor()) { return false; }
            return do_add(part,partType);
        }
        //Part is a clip
        //Cannot have more than one clip
        if(this->has_clip()) { return false; }
        return do_add(part,partType);
    }
    protected:
    bool do_add(ClippedCXA_spt part, ALN_PARTS partType) {
        parts[partType] = part;
        return true;
    }
    public:
    bool has_anchor () const {
        return parts[ANCHOR_LEFT] || parts[ANCHOR_RIGHT];
    }
    bool has_clip () const {
        return parts[CLIP_LEFT] || parts[CLIP_RIGHT];
    }
    bool has_left() const {
        return parts[ANCHOR_LEFT] || parts[CLIP_LEFT];
    }
    bool has_right() const {
        return parts[ANCHOR_RIGHT] || parts[CLIP_RIGHT];
    }
    bool is_right_anchored() const { return bool(parts[ANCHOR_RIGHT]); }
    bool is_split () const { return !(parts[UNSPLIT]); }
    bool is_complete() const {
        return parts[UNSPLIT] || (this->has_anchor() && this->has_clip());
    }
    ClippedCXA_spt get_alignment(ALN_PARTS partType) const {
        return parts[partType];
    }
    ClippedCXA_spt get_anchor() const {
        return (parts[ANCHOR_LEFT]) ? parts[ANCHOR_LEFT] : parts[ANCHOR_RIGHT];
    }
    ClippedCXA_spt get_clip() const {
        return (parts[CLIP_LEFT]) ? parts[CLIP_LEFT] : parts[CLIP_RIGHT];
    }
    std::vector<ClippedCXA_spt> get_alignments() const {
        std::vector<ClippedCXA_spt> alnVec;
        if(this->is_split()){
            alnVec.push_back(
                (parts[ANCHOR_LEFT]) ? parts[ANCHOR_LEFT] : parts[ANCHOR_RIGHT]
            );
            alnVec.push_back(
                (parts[CLIP_LEFT]) ? parts[CLIP_LEFT] : parts[CLIP_RIGHT]
            );
        } else {
            alnVec.push_back(parts[UNSPLIT]);
        }
        return alnVec;
    }
    std::string to_string() const {
        std::string str;
        for(int i = 0; i < NPartTypes; i++){
            if(i > 0){ str += "\t"; }
            str += std::to_string(i) + ") ";
            str += (parts[i]) ? parts[i]->to_string() : ".";
        }
        return str;
    }
};


typedef std::vector<bam1_t*> AlnVector_t;
typedef std::unique_ptr<AlnVector_t> AlnVector_pt;
typedef std::list<MateAlignmentInfo_t> MateAlnInfoList_t;
typedef std::pair<MateAlnInfoList_t,MateAlnInfoList_t> FragmentInfo_t;

void LogMateAlnInfoList(const MateAlnInfoList_t  & list) {
    for(const auto & aln : list){
        fprintf(stderr,"%s\n",aln.to_string().c_str());
    }
}

std::unordered_set<std::string> VirusNameSet;

//===== Function Declarations

void AddClipsToFragmentInfo (   MateAlnInfoList_t & AlnInfoList,
                                const std::vector<ClippedCXA_spt> & clips );
bool AddMateAlignmentInfoToFragment(const MateAlignmentInfo_t obj,
                                    ChimericFragment_t & frag);
FragmentInfo_t ConstructFragmentInfo( bam_hdr_t* header,
                                        const AlnVector_t & alnVec);
void DestroyAlnVector(AlnVector_pt & alnVec);
void FilterIncompleteFragmentInfo(MateAlnInfoList_t & alnInfoList);
int ParseAlnID(bam1_t* aln, std::string & qname, uint8_t & flag);
void ParseReadXA(bam1_t *read, std::string primaryContig,std::vector<CXA> & out);
void ProcessAlnVec(int id, std::ofstream & outbed, bam_hdr_t* header, AlnVector_pt alnVecPtr);
AlnVector_pt ReadAlnSet(open_samFile_t* alnFile, bam1_t* & read_buf);

// ===== MAIN

//Process BAM files from the working environment to extract all
//candidate fragments from the the BWA alignment
//Output the results in BEDPE format
//  ([chr,off,end] up and down, name, bSplit, [strand] up and down,
//Inputs - A path to the viral reference in fasta format
//	 - A path to the working directory
//	 - A path to the bam workspace
//Outputs - A BEDPE formated file 
//Each entry describes a fragment, the breakpoint is downstream of the
//  first interval, and upstream of the second interval
//  I.e the first interval's offset and end are junction distal and proximal
//      and the second interval's offset and end are junction proximal and distal
//The name column is id of the read-pair supporting the candidate
//A read-pair may infer multiple fragments, each of which may support a distinct
//  breakpoint
//The score field is used to store a flag with lots of information
//  The flag corresponds to a bit vector
//  It is big endian IV1 bits, a  spacer bit, then IV2 bits
//  Within interval bits it is little endian See ChimericFragment_t INFOFLAGBIT
//  Whether the interval directly abutts the breakpoint is in here
//The strand columns inform the breakpoint configuration
//      + = off is distal, end is proximal
//      - = off is proximal, end is distal
//      These are reversed for IV2
int main(int argc, char* argv[]) {
    //##PARSE INPUTS
    std::string virus_names_file = argv[1];
    std::string workdir = argv[2];
    std::string workspace = argv[3];

    //##Files to be used from the workspace
    std::string bam_fname = workspace + "/all_alignments.ns.bam";

    //##Output file
    std::string bed_fname = workdir + "/junction-candidates.bedpe";
    std::ofstream outbed(bed_fname);

    //TODO: Multithreading
    //Set up the thread pool
    int nThread = parse_config_threads(workdir + "/config.txt");
    ctpl::thread_pool thread_pool(nThread);

    //##LOAD DATA INTO GLOBAL VARIABLES
    //Load names of viral contigs
    LoadVirusNames(virus_names_file,VirusNameSet);


    bam1_t* read_buf = nullptr;
    open_samFile_t* alnFile = open_samFile(bam_fname.c_str(), false, false);
    std::vector<std::future<void>> futureVec;
    auto header = alnFile->header;
    for(AlnVector_pt alnVecPtr; (alnVecPtr = ReadAlnSet(alnFile,read_buf)) != nullptr; ){
        std::future<void> future = thread_pool.push(
                [ptr = std::move(alnVecPtr), file=&outbed, h=&header ](int id) mutable {
                    ProcessAlnVec(id, *file,*h,std::move(ptr));
                }
        );
        futureVec.push_back(std::move(future));
        //TODO: Maybe make it sort its output on its own?
    }
    for(auto & future : futureVec){
        future.get();
    }
    close_samFile(alnFile);
    bam_destroy1(read_buf);
}

//===== Function Defintions

bool AddMateAlignmentInfoToFragment(const MateAlignmentInfo_t obj,
                                    ChimericFragment_t & frag)
{
    if(!obj.is_complete()) {return false; }
    std::vector<ClippedCXA_spt> alnVec = obj.get_alignments();
    ChimericFragment_t::IV_IDX ivIdx;
    //Check if there is just one alignment
    if(!obj.is_split()){
        ClippedCXA_spt aln = obj.get_alignment(MateAlignmentInfo_t::UNSPLIT);
        ivIdx = (VirusNameSet.count(aln->chr)) ?    ChimericFragment_t::IV2 :
                                                    ChimericFragment_t::IV1;
        return frag.add_alignment( *aln,false,false,false,aln->is_r1(),
                                    aln->clipSide(),ivIdx);
    }
    //Alignment is Split
    //Try to add the anchor
    ClippedCXA_spt anchor = obj.get_anchor();
    ivIdx = (VirusNameSet.count(anchor->chr)) ? ChimericFragment_t::IV2 :
                                                ChimericFragment_t::IV1;
    bool bAnchorIsRight = obj.is_right_anchored();
    if(!frag.add_alignment( *anchor,false,true,!bAnchorIsRight,
                            anchor->is_r1(),anchor->clipSide(),ivIdx))
    {
        return false;
    }
    //Try to add the clip
    ClippedCXA_spt clip = obj.get_clip();
    ivIdx = (VirusNameSet.count(clip->chr)) ? ChimericFragment_t::IV2 :
                                              ChimericFragment_t::IV1;
    uint8_t clipSide = clip->clipSide();
    //0 LA 0 c+ LeftClipped
    //0 LA 1 c- RightClipped
    //1 RA 0 c+ RightClipped
    //1 RA 1 c- LeftCipped
    clipSide |= (bAnchorIsRight == clip->bRev) ?    CXA::LEFT_CLIPPED :
                                                    CXA::RIGHT_CLIPPED;
    return frag.add_alignment(  *clip,true,false,bAnchorIsRight,clip->is_r1(),
                                clipSide,ivIdx);
}



//Given a set of Alignment Information objects, attempt to add clip information
//to the anchors
//Process: iterate over alnInfoList, skipping over non-anchors and leaving them
// as is. Pull the anchor object from the list as a parent object. Try to create 
// a new object from the parent and each clip. for each on that works insert it
// back into the list before the item previously following the parent
//Note: If an anchor does not have any clipps added to it, it is removed
//          (this probably removes the need for a subsequent filtering step ...)
//Inputs    - a list of MateAlignmentInfo_t objects
//          - a vector of shared Clipped CXA pointers to add
//Output    - None, modifies the given MateAlnInfoList_t Object
void AddClipsToFragmentInfo (   MateAlnInfoList_t & alnInfoList,
                                const std::vector<ClippedCXA_spt> & clips )
{
    for(auto it = alnInfoList.begin(); it != alnInfoList.end(); ){
        //Cannot add a clip to an alignment without an anchor
        if(!it->has_anchor()) { it++; continue; }
        MateAlignmentInfo_t parentObj(*it);
        ClippedCXA_spt anchor = parentObj.get_alignments()[0];
        //Pull the anchor out of the list to build new objects
        it=alnInfoList.erase(it);
        //If the primary alignment was reversed, but this alternative alignment 
        //was not, then the clip sequence corresponds to the other side of the alignment
        bool bAnchorSwap = (anchor->bRev != anchor->primary_is_rev());
        for(const ClippedCXA_spt & clip : clips){
            bool bClipLeft = clip->is_left();
            MateAlignmentInfo_t obj(parentObj);
            MateAlignmentInfo_t::ALN_PARTS partType = 
                (bClipLeft == bAnchorSwap) ?
                    MateAlignmentInfo_t::CLIP_RIGHT :
                    MateAlignmentInfo_t::CLIP_LEFT;
            if(obj.add(clip,partType)){
                alnInfoList.insert(it,obj);
            }
        }
    }
}

FragmentInfo_t ConstructFragmentInfo( bam_hdr_t* header,
                                        const AlnVector_t & alnVec) {
    MateAlnInfoList_t R1AlnInfoList;
    MateAlnInfoList_t R2AlnInfoList;

    std::vector<ClippedCXA_spt> alnMappings;
    std::vector<ClippedCXA_spt> R1Clips;
    std::vector<ClippedCXA_spt> R2Clips;
    std::string qName("");
    //Load up all alternative alignments
    for( bam1_t* aln : alnVec){
        uint8_t flag;
	std::string cname = sam_hdr_tid2name(header,aln->core.tid);
        ParseAlnID(aln,qName,flag);
        if(!(flag & ClippedCXA::IS_CLIPPED)){
            if(aln->core.flag & BAM_FREAD1) { flag |= ClippedCXA::IS_R1; }
            if(aln->core.flag & BAM_FREVERSE) { flag |= ClippedCXA::PRIMARY_IS_REVERSED; }
        }
        MateAlnInfoList_t & mateAlnInfoVec =
            (flag & ClippedCXA::IS_R1) ? R1AlnInfoList : R2AlnInfoList;
        std::vector<ClippedCXA_spt> & mateClips =
            (flag & ClippedCXA::IS_R1) ? R1Clips : R2Clips;

        std::vector<CXA> theseAln;
        ParseReadXA(aln,cname,theseAln);

        for(auto aln : theseAln){
            ClippedCXA_spt ccxa(new ClippedCXA(aln,flag));
            alnMappings.push_back(ccxa);
            //If this alignment is from a clip
            if(ccxa->is_clipped()) {
                //Record it for processing after non-clips are processed
                mateClips.push_back(ccxa);
                continue;
            }
            uint8_t clipSide = ccxa->clipSide();
            for( CXA::CLIP_SIDE side : 
                    {CXA::UNCLIPPED, CXA::LEFT_CLIPPED, CXA::RIGHT_CLIPPED} )
            {
                //Figure out which part of an alignment this CXA represents
                MateAlignmentInfo_t::ALN_PARTS part; 
                switch (side) {
                    case CXA::UNCLIPPED:
                        //Actual clips never hit this category, so no need to skip
                        part = MateAlignmentInfo_t::UNSPLIT;
                        break;
                    case CXA::LEFT_CLIPPED:
                        //Skip left clip if the alignment isn't left clipped
                        if(!(clipSide & CXA::LEFT_CLIPPED)) { continue; }
                        part = MateAlignmentInfo_t::ANCHOR_RIGHT;
                        break;
                    case CXA::RIGHT_CLIPPED:
                        //Skip right clip if the alignment isn't right clipped
                        if(!(clipSide & CXA::RIGHT_CLIPPED)) { continue; }
                        part = MateAlignmentInfo_t::ANCHOR_LEFT;
                        break;
                    case CXA::DOUBLE_CLIPPED:
                        //Side can't be in this state...
                        throw std::logic_error("Attempt to interpret alignment as double clipped");
                        break;
                }
                //Create the alignment info object and add this alignment
                mateAlnInfoVec.emplace_back();
                if(!mateAlnInfoVec.back().add(ccxa,part)){
                    throw std::logic_error("Failure to add ccxa to empty mate info");
                }
            }
        }
    }
    AddClipsToFragmentInfo(R1AlnInfoList,R1Clips);
    AddClipsToFragmentInfo(R2AlnInfoList,R2Clips);
    return std::make_pair(R1AlnInfoList,R2AlnInfoList);
}

//Free the bam1_t objects in an alignment vector
void DestroyAlnVector(AlnVector_pt & alnVec) {
    for (bam1_t * aln : *alnVec){
        bam_destroy1(aln);
    }
}

void FilterIncompleteFragmentInfo(MateAlnInfoList_t & alnInfoList){
    for(auto it = alnInfoList.begin(); it != alnInfoList.end();){
        if(it->is_complete()){
            it++;
        } else {
            it = alnInfoList.erase(it);
        }
    }
}

//parses the query name from an alignment object
//  query names are assumed in the form ([^_]+)(_([LR])_([12]))?
//  where $1 is the raw qname to be extracted
//  if $2 is present then the 0x1 bit is set indicating the alignment is a clip
//  the 0x2 bit is set if clipped and $3 is L
//  the 0x4 bit is set if clipped and $4 is 1
//  IF The input read names end in _[LR]_[12], they will be interpretted as
//  clipped
//Inputs - a reference to a string in which to place the raw query name
//       - a reference to a byte in which to store flags
//Output - error code, 0 for success
int ParseAlnID(bam1_t* aln, std::string & qname, uint8_t & flag){
    qname = "";
    flag = 0;
    if(!aln) { return 1; } //Undefined aln object
    std::string alnName = bam_get_qname(aln);
    std::vector<std::string> nameParts = strsplit(alnName,'_');
    qname = nameParts[0];
    if(nameParts.size() > 1) {
        for(size_t i = 1; i < nameParts.size() - 2; i++){
            qname += "_" + nameParts[i];
        }
        char side = nameParts[nameParts.size() - 2][0];
        char read = nameParts[nameParts.size() - 1][0];
        if( (side == 'L' || side == 'R') &&
            (read == '1' || read == '2')) 
        { // Is clipped
            flag |= ClippedCXA::IS_CLIPPED;
            if(side == 'L'){ flag |= ClippedCXA::IS_LEFT; } //Set the left bit
            if(read == '1'){ flag |= ClippedCXA::IS_R1; } //Set the R1 bit
        } else {
            qname += "_" + nameParts[nameParts.size() - 2];
            qname += "_" + nameParts[nameParts.size() - 1];
        }
    }
    return 0;
}

void ParseReadXA (  bam1_t *read, std::string primaryContig,
		    std::vector<CXA> & out){
    uint8_t * nm = bam_aux_get(read,"NM"); 
    int nmVal = (nm) ? bam_aux2i(nm) : 0;

    std::string xaStr = primaryContig + "," + 
			((read->core.flag & BAM_FREVERSE) ? "-" : "+")  +
			std::to_string(read->core.pos+1) + "," +
			"1M" + "," + std::to_string(nmVal);
    out.emplace_back(xaStr);
    out.front().nCigar = read->core.n_cigar;
    out.front().cigar = (uint32_t*) std::realloc(out.front().cigar,
					    sizeof(uint32_t) * (out.front().nCigar));
    memcpy( out.front().cigar,bam_get_cigar(read),
	    sizeof(uint32_t) * read->core.n_cigar);


    uint8_t * xa = bam_aux_get(read,"XA");
    if(!xa) return;
    std::string xaListStr = bam_aux2Z(xa);
    size_t pos = 0, prev = 0;
    while((pos = xaListStr.find(';',prev+1)) != std::string::npos){
	xaStr = xaListStr.substr(prev,pos-prev);
	prev=pos+1;
	out.emplace_back(xaStr);
    }
}

//Given a vector of alignments, stitches combinations of alignments into 
// consistent fragments and outputs them to the provided ofstream
//Each alignment may have alternative alignments, any of which is considered 
//  Equally valid
//Also considers all valid clipping arrangements for an alignment,
//  only the consistent combinations of alt alignments and clip configurations
//  which have both distal positions corresponding to read termini, and infer
//  a chimeric fragment will be carried forward
//  This isn't the most efficient method, but it does allow considering
//  everything
//Responsible for destroying alignments in tha alnVector
//Inputs - an open output file stream ofstream object
//       - an AlnVector_pt object to process
//Output - None, writes to outbed
void ProcessAlnVec(int id, std::ofstream & outbed, bam_hdr_t* header, AlnVector_pt alnVecPtr) {
    FragmentInfo_t fragInfo = ConstructFragmentInfo(header, *alnVecPtr);
    std::string qName("");
    uint8_t dummy;
    ParseAlnID((*alnVecPtr)[0],qName,dummy);
    //Skip fragments which do not have alignments from both the R1 and the R2
    if(!fragInfo.first.size() || !fragInfo.second.size()) { return; }
    //Construct ChimericFragments
    std::vector<ChimericFragment_t> fragmentVec;
    //Iterate over combinations of R1 and R2 alignment information objects
    for(const MateAlignmentInfo_t & r1 : fragInfo.first){
        ChimericFragment_t parentFrag(qName);
        //Skip fragments where the R1 alignments are inconsistent
        if(!AddMateAlignmentInfoToFragment(r1,parentFrag) ){ continue; }
        for(const MateAlignmentInfo_t & r2 : fragInfo.second){
            //Make a copy of the parent frag to work with
            ChimericFragment_t frag(parentFrag);
            //Skip fragments where R2 alignments are inconsistent
            if(!AddMateAlignmentInfoToFragment(r2,frag)) { continue; }
            fragmentVec.push_back(frag);
        }
    }
    bool bValid = true;
    //Perform check that
    //all fragments are chimeric
    for(const ChimericFragment_t & frag : fragmentVec ){
        if(frag.not_chimeric()){
            bValid = false;
            break;
        }
    }
    if(bValid){
        //Collect Complete fragments, and remove fragments which are
        //dominated by another (same fragment with better coverage/ support)
        //Process fragments until there are no more, and retain the dominant fragments
        std::list<ChimericFragment_t> domFragList;
        for(; fragmentVec.size() > 0; fragmentVec.pop_back()){
            const ChimericFragment_t & frag = fragmentVec.back();
            //Skip incomplete fragments
            if(!frag.is_complete()) { continue; }
            bool bDominated = false;
            //Check if this complete fragment is dominated by any previous
            for(auto it = domFragList.begin();
                    !bDominated && it != domFragList.end(); )
            {
                int cmp = it->dominant_comparison(frag);
                if(cmp == 1){ //The previous fragment dominates this fragment
                    bDominated = true;
                } else if(cmp == -1) {//This fragment dominates the previous,
                    //Remove the dominated fragment
                    it = domFragList.erase(it);
                } else { //No domination
                    it++;
                }
            }
            //Skip dominated fragments
            if(bDominated) { continue; }
            domFragList.push_back(frag);
        }
        for( const ChimericFragment_t & frag : domFragList ){
            std::string bedpeStr = frag.to_bedpe();
            mu.lock();
            outbed << bedpeStr << "\n";
            mu.unlock();
        }
    }
    DestroyAlnVector(alnVecPtr);
}

//Given an open bam file, and a read object to act as lookahead buffer, loads
// reads from the bam file into a vector until a read with a different id is
// loaded, this is retained in the lookahead buffer
// a null ptr is returned if no read can be read from the bam file, and the 
// lookahed buffer is empty
//The caller is responsible for detroying and freeing all bam1_t objects in the
//  vector
//Inputs - an open_samFile_t pointer to a valid open bam file
//       - a valid bam1_t object to act as a lookahead buffer
//Output - an AlnVector_pt object, null in the case of an empty vector
//Exceptions - 
AlnVector_pt ReadAlnSet(open_samFile_t* alnFile, bam1_t* & read_buf) {
    AlnVector_pt alnVector( new AlnVector_t());
    std::string qName = "";
    uint8_t flag = 0;
    int parseRes = 0;
    int readRes = 0;
    if(read_buf){
        alnVector->push_back(bam_dup1(read_buf));
        parseRes = ParseAlnID(read_buf,qName,flag);
        if(parseRes){ 
            char buf[100];
            sprintf(buf,"ID parse error: %d",parseRes);
            throw std::invalid_argument(buf);
        }
    } else {
        read_buf = bam_init1();
    }
    while ((readRes = sam_read1(alnFile->file, alnFile->header, read_buf)) >= 0) {
        std::string curQName("");
        parseRes = ParseAlnID(read_buf,curQName,flag);
        if(parseRes){ 
            char buf[100];
            sprintf(buf,"ID parse error: %d",parseRes);
            throw std::invalid_argument(buf);
        }
        if(qName == ""){
            qName = curQName;
        }
        if (qName == curQName){
            alnVector->push_back(bam_dup1(read_buf));
        } else {
            break;
        }
    }
    if(readRes == -1){
        bam_destroy1(read_buf);
        read_buf = nullptr;
    } else if(readRes < -1){
        char buf[100];
        sprintf(buf,"BAM Read Error: %d",readRes);
        throw std::invalid_argument(buf);
    }
    if(!alnVector->size()){
        return nullptr;
    }
    return alnVector;
}
