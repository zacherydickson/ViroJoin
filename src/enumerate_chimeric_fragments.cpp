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

//enum JunctionSide_t {JS_HOST, JS_VIRUS};
//enum GoodClipType_t {
//    GCT_R1L = 0, //good left clip from R1
//    GCT_R1R = 1, //good right clip from R1
//    GCT_R2L = 2, //good left clip from R2
//    GCT_R2R = 3, //good right clip from R2
//};
//
////Note Junctions fall into one of 8 situations:
////  0		1	2		3
////  (A) H+V+, (B) H+V-, (C) H-V+, and (D) H-V-
////  (a)	V-H-, (b) V+H-, (c) V-H+, and (d) V+H+ 
////  Note that the latter 4 are equivalent to the former 4
////  The reported junction strand will be from the former 4
////00 -> ++, 01 -> +-, 10 -> -+, 11 -> --
//uint8_t JunctionOrientation[32] = {
//    //	    Seq		Cli	Anc	aVir	Left
//    0b00, //   1	+	+	V	L
//    0b11, //   1 	+	+	V	R
//    0b11, //   1 	+	+	H	L
//    0b00, //   1 	+	+	H	R
//    0b11, //   1 	+	-	V	L
//    0b00, //   1 	+	-	V	R
//    0b00, //   1 	+	-	H	L
//    0b11, //   1 	+	-	H	R
//    0b10, //   1 	-	+	V	L
//    0b01, //   1 	-	+	V	R
//    0b10, //   1 	-	+	H	L
//    0b01, //   1 	-	+	H	R
//    0b01, //   1 	-	-	V	L
//    0b10, //   1 	-	-	V	R
//    0b10, //   1 	-	-	H	L
//    0b01, //   1 	-	-	H	R
//    0b11, //   2 	+	+	V	L
//    0b00, //   2 	+	+	V	R
//    0b00, //   2 	+	+	H	L
//    0b11, //   2 	+	+	H	R
//    0b00, //   2 	+	-	V	L
//    0b11, //   2 	+	-	V	R
//    0b11, //   2 	+	-	H	L
//    0b00, //   2 	+	-	H	R
//    0b01, //   2 	-	+	V	L
//    0b10, //   2 	-	+	V	R
//    0b10, //   2 	-	+	H	L
//    0b01, //   2 	-	+	H	R
//    0b10, //   2 	-	-	V	L
//    0b01, //   2 	-	-	V	R
//    0b01, //   2 	-	-	H	L
//    0b10  //   2 	-	-	H	R
//};
//
////Same 4(8) cases as previously but the information is more constrained
////only 3 tests needed
//uint8_t PairedJunctionOrientation[8] = {
//    //R1 Host
//    //	R1 Fwd
//    0b01, //R2 Fwd
//    0b00, //R2 Rev
//    // R1 Rev
//    0b11, //R2 Fwd
//    0b10, //R2 Rev
//    //R1 Virus
//    //	R1 Fwd
//    0b01, //R2 Fwd
//    0b11, //R2 Rev
//    //	R1 Rev
//    0b00, //R2 Fwd
//    0b10, //R2 Rev
//};
//
//typedef std::array<bam1_t*,4> ClipArray_t;
//typedef std::unordered_map<std::string, ClipArray_t> GoodClipMap_t;
//typedef std::unordered_set<std::string> QNameSet_t;


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
        std::cerr << aln.to_string() << "\n";
    }
}

std::unordered_set<std::string> VirusNameSet;
//GoodClipMap_t GoodClipMap;
//QNameSet_t GoodClipSet;
//std::mutex mtx;

//===== Function Declarations

void AddClipsToFragmentInfo (   MateAlnInfoList_t & AlnInfoList,
                                const std::vector<ClippedCXA_spt> & clips );
bool AddMateAlignmentInfoToFragment(const MateAlignmentInfo_t obj,
                                    ChimericFragment_t & frag);
FragmentInfo_t ConstructFragmentInfo( bam_hdr_t* header,
                                        const AlnVector_t & alnVec);
//std::string ConstructCandidateString();
//void DestroyGoodClips();
void DestroyAlnVector(AlnVector_pt & alnVec);
void FilterIncompleteFragmentInfo(MateAlnInfoList_t & alnInfoList);
//char DetermineClipJunctionStrand (uint8_t flag);
//std::array<char,2> DetermineJunctionOrientation (   bool bViralAnchor,
//		    bool isLeftClip, bool bAnchorRev, bool bClipRev, bool isR1);
//std::array<char,2> DeterminePairedJunctionOrientation(bool r1Virus, bool r1Rev,
//		    bool r2Rev);
//std::array<hts_pos_t,2> DeterminePairedJuncRelPos(bool r1Virus, 
//                    bool r1Rev, bool r2Rev, hts_pos_t r1L, hts_pos_t r1R,
//                    hts_pos_t r2L, hts_pos_t r2R);
//bool IsLeftOfJunction(uint8_t flag);
//void LoadAnchorOrientation(std::string fname);
//void LoadGoodClips(std::string fname);
//void OutputBEDEntries(	std::ofstream & outbed, const bam1_t* read,
//			std::string cname, uint8_t clipflag = 0x0,
//			std::string qname = std::string());
int ParseAlnID(bam1_t* aln, std::string & qname, uint8_t & flag);
void ParseReadXA(bam1_t *read, std::string primaryContig,std::vector<CXA> & out);
void ProcessAlnVec(std::ofstream & outbed, bam_hdr_t* header, AlnVector_pt alnVecPtr);
//void ProcessAnchor(bam1_t *read, std::string cname, std::ofstream & outbed);
//void ProcessClip(bam1_t *read, std::string cname, std::ofstream & outbed);
//void ProcessPair(   bam1_t *r1, bam1_t *r2, std::string cname1,
//		    std::string cname2, std::ofstream &outbed);
//void ProcessPairs(std::string fname, std::ofstream & outbed);
//void ProcessSplitRead(	bam1_t *anchor, bam1_t clip, int jSide,
//			std::string primaryContig, std::string clipCName,
//			std::ofstream & outbed);
//void ProcessSplitReads(	std::string anchor_fname, std::string clip_fname, 
//			int jSide, std::ofstream & outbed);
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
    // TODO: Implement the design described above

    /*//Set up the thread pool
    int nThread = parse_config_threads(workdir + "/config.txt");
    ctpl::thread_pool thread_pool(nThread); */

    //##LOAD DATA INTO GLOBAL VARIABLES
    //Load names of viral contigs
    LoadVirusNames(virus_names_file,VirusNameSet);

    bam1_t* read_buf = nullptr;
    open_samFile_t* alnFile = open_samFile(bam_fname.c_str(), false, false);
    //int counter =0;
    for(AlnVector_pt alnVecPtr; (alnVecPtr = ReadAlnSet(alnFile,read_buf)) != nullptr; ){
        //if(counter++ == 1)
        //std::cout << "BEGIN BLOCK\t" << alnVecPtr->size() << "\n";
        ProcessAlnVec(outbed,alnFile->header,std::move(alnVecPtr));
        ////TODO: Process Aln Vec
        //std::string qname;
        //uint8_t flag;
        //for( bam1_t* aln : *alnVecPtr ){
        //    ParseAlnID(aln,qname,flag);
        //    std::cout << qname << "\t" << int(flag) << "\n";
        //}
        //DestroyAlnVector(alnVecPtr);
        //std::cout << "END BLOCK\n";
    }
    close_samFile(alnFile);
    bam_destroy1(read_buf);


    
    ////Load the ids and directions of clips which map properly
    //for (int side = JS_HOST; side <= JS_VIRUS; side++){
    //    LoadGoodClips(clip_bam_fnames[side]);
    //    ProcessSplitReads(  anchor_bam_fnames[side],clip_bam_fnames[side],
    //    		    side,outbed);
    //    DestroyGoodClips();
    //}

    ////Pass over the paired reads to find valid chimeras
    //ProcessPairs(bam_fname,outbed);
}

//===== Function Defintions

////Given information to be printed, constructs a string describing the candidate breakpoint
////Inputs -
////Output - a string
//std::string ConstructCandidateString(   std::string chr, size_t pos,
//                                        std::string qname, char strand)
//{
//    std::string str =  chr + '\t' + std::to_string(pos) + '\t' +
//	                std::to_string(pos+ 1) + '\t' +
//                        qname + "\t.\t" + strand;
//    return str;
//}
//
////Looks up the case in a precalculated table based on 5 boolean values
////This table gives information on if the human and viral sides are in the
////+ or - orientation
////The result can then be matched to the anchor and clip as needed
////Inputs - 5 boolean values defining the orientation
////Output - a string 
//std::array<char,2> DetermineJunctionOrientation (  bool bViralAnchor, bool isLeftClip,
//					    bool bAnchorRev, bool bClipRev,
//					    bool isR1) {
//    char anchorStrand, clipStrand;
//    ////Result is two bits where the one's bit is one if the virus is rev
//    //// and the twos bit is one if the host is rev
//    uint8_t result = 0;
//    bool bRightClip = !isLeftClip;
//    bool bIs5Prime = (bViralAnchor != bRightClip);
//    if(!bIs5Prime) result |= 0b10; //Host is Inverted at the 3' junction
//    //virus matches host in fwd insertions, and doesn't in reverse insertions
//    result |= ((result >> 1) ^ bClipRev); // XOR can act as a negate
//    //If the anchor and clip strand do not match ( 0b01=1 or 0b10=2 )
//    //And the virus is the anchor, the strands need to be flipped
//    if(bClipRev && bViralAnchor){
//	result = (~result) & 0b11;
//    }
//    if(bViralAnchor){ // Anchor takes viral result
//        anchorStrand = (result & 0b01) ? '-' : '+';
//        clipStrand = (result & 0b10) ? '-' : '+';
//    } else {
//        anchorStrand = (result & 0b10) ? '-' : '+';
//        clipStrand = (result & 0b01) ? '-' : '+';
//    }
//    return {anchorStrand,clipStrand};
//}
//
//
//std::array<char,2> DeterminePairedJunctionOrientation(bool r1Virus, bool r1Rev,
//		    bool r2Rev) {
//    char hostStrand, virStrand;
//    uint8_t flag = 0;
//    if(r1Virus) flag |= 0x4;
//    if(r1Rev) flag |= 0x2;
//    if(r2Rev) flag |= 0x1;
//    uint8_t result = PairedJunctionOrientation[flag];
//    hostStrand = (result & 0x2) ? '-' : '+';
//    virStrand = (result & 0x1) ? '-' : '+';
//    return {hostStrand,virStrand};
//}
//
////Procedure have an array of 4 positions, and manipulate it according to the
////specific read characteristics
////The initial assumption is that the pair does support a junction and therefore
////the 5` end of a read is the junction distal position
//std::array<hts_pos_t,2> DeterminePairedJuncRelPos(bool r1Virus,
//                    bool r1Rev, bool r2Rev, hts_pos_t r1L, hts_pos_t r1R,
//                    hts_pos_t r2L, hts_pos_t r2R)
//{
//    hts_pos_t r1prox =  r1Rev ? r1L : r1R;
//    hts_pos_t r1dist = !r1Rev ? r1L : r1R;
//    hts_pos_t r2prox =  r2Rev ? r2L : r2R;
//    hts_pos_t r2dist = !r2Rev ? r2L : r2R;
//    //Start with the assumption that host is r1
//    //so that it goes HD,HP | VP,VD
//    std::array<hts_pos_t,4> p = {r1dist,r1prox,r2prox,r2dist};
//    //If host is r2, everything 
//    if(r1Virus){
//        std::reverse(std::begin(p),std::end(p));
//    }
//    //Proximal positions are the middle two
//    return {p[1],p[2]};
//}
//
//
////Opens a bam file containing mapped clips.
////It is assumed that the bam file only contains primary mappped clips
////(no secondary/supplementary/unmapped)
////All reads in these files are assumed to define good clips
////Each alignment object is stored for later use (the entire clips file is
////loaded into memory)
////Input  - a string reperesenting a file name
////Output - none, modifies the global GoodClipMap
//void LoadGoodClips(std::string fname){
//    open_samFile_t* clips_file = open_samFile(fname.c_str(), false, false);
//    bam1_t* read = bam_init1();
//
//    while (sam_read1(clips_file->file, clips_file->header, read) >= 0) {
//        std::string clip_name = bam_get_qname(read);
//        std::string qname = clip_name.substr(0, clip_name.length()-4);
//	GoodClipSet.insert(qname);
//	if(!GoodClipMap.count(qname)){
//	    GoodClipMap[qname] = {nullptr,nullptr,nullptr,nullptr};
//	}
//	bool isLeftClip(clip_name[clip_name.length()-3] == 'L');
//	bool isR1(clip_name[clip_name.length()-1] == '1');
//	int idx = 0;
//	if(!isLeftClip) idx += 1;
//	if(!isR1) idx += 2;
//	GoodClipMap[qname][idx]	= bam_dup1(read);
//    }
//    close_samFile(clips_file);
//    bam_destroy1(read);
//}
//
//void DestroyGoodClips(){
//    for( auto & pair : GoodClipMap ){
//	for( int i = 0; i < 4; i++){
//	    if(pair.second[i]){
//		bam_destroy1(pair.second[i]);
//		pair.second[i] = nullptr;
//	    }
//	}
//    }
//    GoodClipMap.clear();
//}


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
            //std::cerr << "\tAttempt Add As " << partType << "\n" << clip->to_string() << "\n" << "\tto\n" << obj.to_string() << "\n";
            if(obj.add(clip,partType)){
                //std::cerr << "\tSuccess\n";
                alnInfoList.insert(it,obj);
            }
            //else { std::cerr << "\tFailure\n"; }
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
        //std::cerr << "SAM | " << cname << "\t" << qName << "\t" << aln->core.flag << "\t" << aln->core.pos << "\t" << std::to_string(int(flag)) << "\n";
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

        //std::cerr << "=== Build CXA\n";
        for(auto aln : theseAln){
            ClippedCXA_spt ccxa(new ClippedCXA(aln,flag));
            //std::cerr << "CCXA | " << ccxa->to_string() << "\n";
            alnMappings.push_back(ccxa);
            //If this alignment is from a clip
            if(ccxa->is_clipped()) {
                //std::cerr << "\tClipped\n";
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
                        //Skip attempt UNCLIPPED if the alignment is a clip!
                        //if(flag & ClippedCXA::IS_CLIPPED) { continue; }
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
                //std::cerr << "\tAttempt Add As " << side << "\n";;
                //Create the alignment info object and add this alignment
                mateAlnInfoVec.emplace_back();
                if(!mateAlnInfoVec.back().add(ccxa,part)){
                    throw std::logic_error("Failure to add ccxa to empty mate info");
                }
                //std::cerr << "MAIV | " << mateAlnInfoVec.back().to_string() << "\n";
            }
        }
        //std::cerr << alnMappings.back().size() << "\t";
    }
    //std::cerr << "\n=== Pre-Clip R1 MAIV\n";
    //LogMateAlnInfoList(R1AlnInfoList);
    //std::cerr << "\n=== Pre-Clip R2 MAIV\n";
    //LogMateAlnInfoList(R2AlnInfoList);
    //Add clips to their respective anchors
    //std::cerr << "\n=== Add Clips R1\n";
    AddClipsToFragmentInfo(R1AlnInfoList,R1Clips);
    //std::cerr << "\n=== Add Clips R2\n";
    AddClipsToFragmentInfo(R2AlnInfoList,R2Clips);
    //std::cerr << "\n=== Post-Clip R1 MAIV\n";
    //LogMateAlnInfoList(R1AlnInfoList);
    //std::cerr << "\n=== Post-Clip R2 MAIV\n";
    //LogMateAlnInfoList(R2AlnInfoList);
    ////Filter out incomplete fragments
    //std::cerr << "\n=== Filter Fragments R1 \n";
    //FilterIncompleteFragmentInfo(R1AlnInfoList);
    //std::cerr << "\n=== Filter Fragments R2 \n";
    //FilterIncompleteFragmentInfo(R2AlnInfoList);
    //std::cerr << "\n=== Post-Filter R1 MAIV\n";
    //LogMateAlnInfoList(R1AlnInfoList);
    //std::cerr << "\n=== Post-Filter R2 MAIV\n";
    //LogMateAlnInfoList(R2AlnInfoList);
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
        //std::cerr << "\tCheck for completeness\n" << it->to_string() << "\n";
        if(it->is_complete()){
            //std::cerr << "\tComplete\n";
            it++;
        } else {
            //std::cerr << "\tIncomplete\n";
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
void ProcessAlnVec(std::ofstream & outbed, bam_hdr_t* header, AlnVector_pt alnVecPtr) {
    //std::cerr << "Init Proc Aln\n";
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
        //std::cerr << "Pre-R1\t" << parentFrag.to_bedpe(true) << "\n";
        //Skip fragments where the R1 alignments are inconsistent
        if(!AddMateAlignmentInfoToFragment(r1,parentFrag) ){ /*std::cerr << "failure to add R1\n";*/ continue; }
        //std::cerr << "Post-R1\t"<< parentFrag.to_bedpe(true) << "\n";
        for(const MateAlignmentInfo_t & r2 : fragInfo.second){
            //Make a copy of the parent frag to work with
            ChimericFragment_t frag(parentFrag);
            //std::cerr << "Pre-R2\t"<< frag.to_bedpe(true) << "\n";
            //Skip fragments where R2 alignments are inconsistent
            if(!AddMateAlignmentInfoToFragment(r2,frag)) {  /*std::cerr << "failure to add R2\n";*/ continue; }
            //std::cerr << "Post-R2\t"<< frag.to_bedpe(true) << "\n";
            fragmentVec.push_back(frag);
        }
    }
    //std::vector<std::vector<CXA>> alnMappings;
    //std::string qName("");
    //std::vector<uint8_t> flagVec;
    ////Load up all alternative alignments
    //for( bam1_t* & aln : *alnVecPtr){
    //    alnMappings.push_back(std::vector<CXA>());
    //    flagVec.push_back(0);
    //    std::string cname = sam_hdr_tid2name(header,aln->core.tid);
    //    ParseAlnID(aln,qName,flagVec.back());
    //    //if(qName == "Fake:H+V-:HChr8:144Mbp:chr11:86Mbp:V:96kbp1") {
    //    //    std::cerr << cname << "\t" << aln->core.flag << "\t" << aln->core.pos << "\n";
    //    //}
    //    if(!(flagVec.back() & 0x1) && aln->core.flag & BAM_FREAD1) {
    //        flagVec.back() |= 0x4;
    //    }
    //    ParseReadXA(aln,cname,alnMappings.back());
    //    //std::cerr << alnMappings.back().size() << "\t";
    //}
    //bool bLog = true;//(qName == "Fake:H+V-:HChr8:144Mbp:chr11:86Mbp:V:96kbp1");
    //if(bLog) std::cerr << "\n" << qName << "\t" << alnMappings.size() << "\n";
    ////Construct all fragments which are consistent with the alignments
    //std::vector<ChimericFragment_t> fragmentVec = {ChimericFragment_t(qName)};
    //return;
    ////TODO: Construct Fragments
//  //  fragmentVec[0].name = qName;
    //std::vector<ChimericFragment_t> fragmentVecTmp;
    //if(bLog) std::cerr << "Pre Build Fragments\n";
    //for(size_t i = 0; i < alnMappings.size(); i++){
    //    const std::vector<CXA> & partMappings = alnMappings[i];
    //    const uint8_t & flag = flagVec[i];
    //    while(!fragmentVec.empty()){
    //        ChimericFragment_t & parentFrag = fragmentVec.back();
    //        //bool bMod = false;
    //        int part = 0;
    //        for( const CXA & cxa : partMappings){
    //            part++;
    //            bool isViral = VirusNameSet.count(cxa.chr);
    //            uint8_t clipSide = cxa.clipSide();
    //            if(flag & 0x1){ // If the alignment is from a clip
    //                //If left clip, then the right is clipped away
    //                //Otherwise the left is clipped away
    //                //This flips if the clip maps to the reverse strand
    //                //1Left 0+ = Right 
    //                //0Right 0+ = Left
    //                //1Left 1- = Left
    //                //0Right 1- = Right
    //                //Boils down to whehter they match or not
    //                clipSide |= (bool(flag & 0x2) != cxa.bRev) ? CXA::RIGHT_CLIPPED : CXA::LEFT_CLIPPED;
    //            }
    //            //We are setting the host side as interval 1 by convention,
    //            //as a result, any molecules which do not have both intervals
    //            //is non chimeric, it will also be incomplete
    //            ChimericFragment_t::IV_IDX ivIdx =  (isViral) ?
    //                                                ChimericFragment_t::IV2 :
    //                                                ChimericFragment_t::IV1;
    //            //Attempt addition of the alignment assuming it is:
    //            //  left-clipped - unless there are no left clipped bases
    //            //  right-clipped - unless there are no right clipped bases
    //            //  unclipped - unless it is already known to be clipped
    //            for( CXA::CLIP_SIDE side : 
    //                    {CXA::UNCLIPPED, CXA::LEFT_CLIPPED, CXA::RIGHT_CLIPPED} )
    //            {
    //                //Skip attempt UNCLIPPED if the alignment is a clip!
    //                if((side == CXA::UNCLIPPED) && (flag & 0x1)) { continue; }
    //                //Skip left clip if the alignment isn't left clipped
    //                if((side == CXA::LEFT_CLIPPED) && !(clipSide & CXA::LEFT_CLIPPED)) { continue; }
    //                //Skip right clip if the alignment isn't right clipped
    //                if((side == CXA::RIGHT_CLIPPED) && !(clipSide & CXA::RIGHT_CLIPPED)) { continue; }
    //                //Attempt to add the alignment interpretting it with the 
    //                //  current clip status and direction
    //                if (bLog) std::cerr << "Attempt Add Entry " << i << " Part " << part << " Side " << side << "\n";
    //                ChimericFragment_t frag = parentFrag;
    //                if(bLog) std::cerr << fragmentVecTmp.size() << ": PRE\t" <<frag.to_bedpe(true) << "\n";
    //                if(frag.add_alignment(  cxa, (flag & 0x1),
    //                                        (side != CXA::UNCLIPPED) && !(flag & 0x1),
    //                                        side == CXA::RIGHT_CLIPPED,
    //                                        flag & 0x4, clipSide, ivIdx))
    //                {
    //                    fragmentVecTmp.push_back(frag);
    //                    //bMod = true;
    //                }
    //                else {
    //                    if(bLog) std::cerr << "Failure to add\n";
    //                }
    //                if(bLog) std::cerr << fragmentVecTmp.size() << ": POST\t" << frag.to_bedpe(true) << "\n";
    //            }
    //        }
    //        //Retain the parent if no alignments could be added
    //        //if(!bMod){
    //        //Retain the parent so that fragments which don't require all alignments can be considered
    //            fragmentVecTmp.push_back(parentFrag);
    //        //}
    //        fragmentVec.pop_back();
    //    }
    //    std::swap(fragmentVec,fragmentVecTmp);
    //}
    bool bLog = true;
    //if(bLog) std::cerr << fragmentVec.size() << "\n";
    bool bValid = true;
    //Perform check that
    //all fragments are chimeric
    for(const ChimericFragment_t & frag : fragmentVec ){
        //if(bLog) std::cerr << "TEST\n" << frag.to_bedpe(true) << "\n";
        if(frag.not_chimeric()){
            //if(bLog) std::cerr << "FAILURE\n";
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
            //std::cerr << "PARTIAL\t"<< frag.to_bedpe(true) << "\n";
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
        //std::unordered_set<std::string> knownFragments;
        for( const ChimericFragment_t & frag : domFragList ){
            std::string bedpeStr = frag.to_bedpe();
            //auto pair = knownFragments.insert(bedpeStr);
            //if(pair.second){
            //std::cerr << "FINAL\t"<< frag.to_bedpe(true) << "\n";
            mu.lock();
            outbed << bedpeStr << "\n";
            mu.unlock();
            //} else {
            //    std::cerr << "Duplication occured\n";
            //}
        }
    }
    DestroyAlnVector(alnVecPtr);
    //if(bLog) std::cerr << "Terminate Proc Aln\n";
}

//void ProcessPair(   bam1_t *r1, bam1_t *r2, std::string cname1,
//		    std::string cname2, std::ofstream & outbed){
//    std::string qname = bam_get_qname(r1);
//    //Clipped Reads get priority
//    if(GoodClipSet.count(qname)) return;
//    //Small optimization to quickly exclude non-chimeric reads
//    if(r1->core.tid == r2->core.tid) return;
//    bool r1IsVirus = (VirusNameSet.count(cname1));
//    bool r2IsVirus = (VirusNameSet.count(cname2));
//    //Again can quickly skip non-chimerics (Double host or double virus)
//    if(r1IsVirus == r2IsVirus) return;
//    //Exclude Homopolymers
//    if(is_poly_ACGT(r1) || is_poly_ACGT(r2)) return;
//
//
//
//    std::unordered_set<std::string> potentialEntries;
//
//    std::vector<CXA> r1Mappings;
//    std::vector<CXA> r2Mappings;
//    ParseReadXA(r1,cname1,r1Mappings);
//    ParseReadXA(r2,cname2,r2Mappings);
//
//    for(size_t i = 0; i < r1Mappings.size(); i++){
//	const CXA & r1Map = r1Mappings[i];
//	for(size_t j = 0; j < r2Mappings.size(); j++){
//	    const CXA & r2Map = r2Mappings[j];
//	    //All alignments for a segment mapping to host must be to host
//	    //All alignments for a segment mapping to virus must be to virus
//	    if(r1IsVirus != bool(VirusNameSet.count(r1Map.chr))) return;
//	    if(r2IsVirus != bool(VirusNameSet.count(r2Map.chr))) return;
//            //Determine junction orientation
//	    std::array<char,2> strands = DeterminePairedJunctionOrientation(
//		    r1IsVirus,r1Map.bRev,r2Map.bRev);
//            std::array<hts_pos_t,2> proxPos = DeterminePairedJuncRelPos(
//                    r1IsVirus,r1Map.bRev,r2Map.bRev,
//                    r1Map.pos,r1Map.endpos(),
//                    r2Map.pos,r2Map.endpos());
//            //Case: H + V + r1is Host
//	    //Construct Entries
//	    std::string hostChr = r1Map.chr;
//	    std::string virChr = r2Map.chr;
//            if(r1IsVirus){ // r2 is Host Side
//		std::swap(hostChr,virChr);
//            }
//	    hts_pos_t hostPos = proxPos[0];
//	    hts_pos_t virPos = proxPos[1];
//	    //Store the unique entries
//            potentialEntries.insert(
//                    ConstructCandidateString(   hostChr,hostPos,qname,
//                                                strands.front()));
//            potentialEntries.insert(
//                    ConstructCandidateString(   virChr,virPos,qname,
//                                                strands.back()));
//	}
//    }
//
//    //The pair passed all filters, can output the entries now
//    for(auto entry : potentialEntries){
//	outbed << entry << "\n";
//    }
//
//}
//
////Opens a given bam file and outputs all of the host/virus side junctions
////each properly mapped pair of reads suppports, only read pairs which
////only map in chimeric configurations are accepted
////Inputs - a path to a mapped clip file
////	 - an ofstream object to which to write
////	 - Also uses the Global Good Clips Set for filtering
////Output - None, writes to outbed
//void ProcessPairs(std::string fname, std::ofstream & outbed){
//    open_samFile_t* bam_file = open_samFile(fname.c_str(), false, false);
//    bam1_t* read1 = bam_init1();
//    bam1_t* read2 = bam_init1();
//
//
//    while (sam_read1(bam_file->file, bam_file->header, read1) >= 0 &&
//	   sam_read1(bam_file->file, bam_file->header, read2) >= 0 ) {
//	if(read1->core.flag & (BAM_FUNMAP | BAM_FMUNMAP)) continue;
//	std::string qname = bam_get_qname(read1);
//	if(qname != bam_get_qname(read2)){
//	    throw std::invalid_argument("Mates not adjacent in paired reads bam");
//	}
//	std::string cname1 = sam_hdr_tid2name(bam_file->header,read1->core.tid);
//	std::string cname2 = sam_hdr_tid2name(bam_file->header,read2->core.tid);
//	//Provide the pair in segment 1 segment 2 order
//	if(read1->core.flag & BAM_FREAD1){
//	    ProcessPair(read1,read2,cname1,cname2,outbed);
//	} else {
//	    ProcessPair(read2,read1,cname2,cname1,outbed);
//	}
//    }
//
//
//    close_samFile(bam_file);
//    bam_destroy1(read1);
//    bam_destroy1(read2);
//}
//
//void ProcessSplitRead(	bam1_t *anchor, bam1_t *clip, int jSide, 
//			std::string primaryContig, std::string clipCName,
//			int lrIdx, std::ofstream & outbed){
//
//    bool bViralAnchor = (jSide == JS_VIRUS);
//    bool isLeftClip = (lrIdx % 2 == 0);
//    bool isR1 = (lrIdx < 2);
//    //Load Clip Alts
//    std::vector<CXA> anchorMappings;
//    std::vector<CXA> clipMappings;
//    ParseReadXA(anchor,primaryContig,anchorMappings);
//    ParseReadXA(clip,clipCName,clipMappings);
//    
//    //Process Primary
//    
//    std::string qname = bam_get_qname(anchor);
//    qname += (isR1) ? "_1" : "_2";
//
//
//    //Iterate over all pairs of anchor and clips
//    //And note all uniq breakpoints this read supports
//    std::unordered_set<std::string> uniqBPStrSet;
//    for(size_t i = 0; i < anchorMappings.size(); i++){
//	const CXA & anchorMap = anchorMappings[i];
//        bool isViralAnchorXA = (VirusNameSet.count(anchorMap.chr));
//        //Skip alt anchors which don't have the same virus status as the primary
//        // anchor alignment
//        if(bViralAnchor != isViralAnchorXA){
//            continue;
//        }
//	for(size_t j = 0; j < clipMappings.size(); j++){
//	    const CXA & clipMap = clipMappings[j];
//	    hts_pos_t anchorPos = (isLeftClip) ? anchorMap.pos : anchorMap.endpos();
//	    hts_pos_t clipPos = (isLeftClip) ? clipMap.endpos() : clipMap.pos;
//	    std::array<char,2> strands = DetermineJunctionOrientation(bViralAnchor,
//				    isLeftClip,anchorMap.bRev,clipMap.bRev,isR1);
//            uniqBPStrSet.insert(ConstructCandidateString(anchorMap.chr,
//                                                         anchorPos,
//                                                         qname,
//                                                         strands.front()));
//            uniqBPStrSet.insert(ConstructCandidateString(clipMap.chr,
//                                                         clipPos,
//                                                         qname,
//                                                         strands.back()));
//	}
//    }
//    //Ouput each breakpoint
//    for(const std::string & bpStr : uniqBPStrSet){
//        outbed << bpStr << "\n";
//    }
//}
//
////Opens a given bam file and outputs all of the host/virus side junctions
////each clip or anchor supports
////Inputs - paths to both anchor and clip files the latter is only used for
////	    its header
////	 - a boolean indicating wether a clip or anchor file has been provided
////	 - an ofstream object to which to write
////Output - None, writes to outbed
//void ProcessSplitReads(	std::string anchor_fname, std::string clip_fname,
//			int jSide, std::ofstream & outbed){
//    open_samFile_t* anchors_file = open_samFile(anchor_fname.c_str(),false,false);
//    //Only needed for its header
//    open_samFile_t* clips_file = open_samFile(clip_fname.c_str(),false,false);
//    bam1_t* read = bam_init1();
//
//    while (sam_read1(anchors_file->file, anchors_file->header, read) >= 0) {
//        std::string qname = bam_get_qname(read);
//	std::string cname = sam_hdr_tid2name(	anchors_file->header,
//						read->core.tid);
//	//Make sure there are any good clips for this anchor
//	if(!GoodClipMap.count(qname)) continue;
//	bool isR1 = (read->core.flag & BAM_FREAD1);
//	int leftIdx = (isR1) ? GCT_R1L : GCT_R2L;
//	int rightIdx = (isR1) ? GCT_R1R : GCT_R2R;
//	//Iterate over both left and right clips
//	for(int idx = leftIdx; idx <= rightIdx; idx++){
//	    //Skip there isn't a matching read
//	    if(!GoodClipMap[qname][idx]) continue;
//	    //No good clip exists for this segment on this side;
//	    std::string clip_cname = sam_hdr_tid2name(	clips_file->header,
//						GoodClipMap[qname][idx]->core.tid);
//	    ProcessSplitRead(	read,GoodClipMap[qname][idx],jSide,cname,
//				clip_cname,idx,outbed);
//	}
//    }
//    close_samFile(anchors_file);
//    close_samFile(clips_file);
//    bam_destroy1(read);
//}



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
        //std::cerr << alnVector->size() << "\t"  << qName << "\t" << curQName << "\t" << int(flag) << "\n";
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
