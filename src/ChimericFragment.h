#ifndef CHIMERIC_FRAGMENT_H
#define CHIMERIC_FRAGMENT_H

#include <string>
#include <htslib/sam.h>
#include "sam_utils.h"
#include <stdexcept>
#include <bitset>

//Structure for building and tracking potentially chimeric fragments
//Considered complete if distal ends of the fragment correspond to read termini
//  are mapped
//Considered chimeric if the host status of the up and downstream intervals are
//  not equal
struct ChimericFragment_t {
    //Static Members
    public:
    enum IV_IDX {
        IV1 = 0, IV2 = 1
    };
    enum INFOFLAGBIT {
        OPENS_LEFT = 0x1,
        IS_SPLIT = 0x2,
        HAS_INTERVAL = 0x4,
        DISTAL_IS_TERMINAL = 0x8,
        PROXIMAL_IS_TERMINAL = 0x10,
        FROM_R1 = 0x20,
        FROM_R2 = 0x40,
    };
    static const size_t MATE_SHIFT = 5;
    static const size_t TERMINAL_SHIFT = 3;
    static const size_t FLAG_BITS = 7;
    //static const size_t IV1 = 0;
    //static const size_t IV2 = 1;
    //Members
    protected:
    std::string name;
    std::array<std::string,2> chr;
    std::array<size_t,2> off;
    std::array<size_t,2> end;
    std::array<uint16_t,2> flag;
    //std::array<char,2> bOpensLeft;
    //std::array<std::string,2> cigar;
    //std::array<bool,2> bDistalIsTerminal;
    //std::array<bool,2> bHasInterval;
    //bool bSplit;
    //Con-/Destruction
    public:
    ChimericFragment_t(const std::string name = "") :
        name(name),
        chr({"",""}), off({0,0}), end({0,0}), flag({0,0})
        //bOpensLeft({false,false}),
        //cigar({"",""}), bDistalIsTerminal({false,false}),
        //bHasInterval({false,false}),
        //bSplit(false)
    {}
    //Accessors
    protected:
    bool both(INFOFLAGBIT bit) const {
        return (flag[IV1] & bit) && (flag[IV2] & bit);
    }
    bool either(INFOFLAGBIT bit) const {
        return (flag[IV1] & bit) || (flag[IV2] & bit);
    }
    bool neither(INFOFLAGBIT bit) const {
        return !((flag[IV1] & bit) || (flag[IV2] & bit));
    }
    public:
    bool is_complete() const { return this->both(DISTAL_IS_TERMINAL); }
    bool is_chimeric() const { return this->both(HAS_INTERVAL); }
    bool not_chimeric() const;
    //bool is_chimeric() const {
    //    return this->bHost[IV1] != this->bHost[IV2];
    //}
    bool is_split(IV_IDX ivIdx) const { return this->flag[ivIdx] & IS_SPLIT; }
    
    //std::string get_cigar(bool bUp) const { return this->cigar[bUp]; }
    //Mutators
    protected:
    static void set_bit(uint16_t & flag, INFOFLAGBIT bit) {
        flag |= (~0 & bit);
    }
    void set_bit(IV_IDX ivIdx, INFOFLAGBIT bit) {
        set_bit(this->flag[ivIdx], bit);
    }
    //Methods
    public:
    bool add_alignment( const CXA & aln, bool isClip,
                        bool isLeft, bool isR1,
                        ChimericFragment_t::IV_IDX ivIdx);
    int dominant_comparison (const ChimericFragment_t & other);
    std::string to_bedpe(bool bitflag = false) const ;
};

//Output - true if the alignment was successfully added, false otherwise
//          an alignment fails to be added if a fragment cannot be generated
//          which is consistent with the existing fragment and the alignment
bool ChimericFragment_t::add_alignment( const CXA & aln, bool isClip,
                                        bool isLeft, bool isR1,
                                        ChimericFragment_t::IV_IDX ivIdx) 
{
    //std::cerr << "\tInit Add\n";
    uint16_t & flag = this->flag[ivIdx];
    //Step one: Assign up and down based on the host sequence
    //  this collapses H+V+ and V-H-, etc. together in the end
    //aln coords are zero indexed
    //bed uses half open intervals
    //off should be inclusive 0 indexed (or exclusive 1 indexed)
    //end should be exclusive 0 indexed (i.e inclusive 1 indexed)
    size_t off = (aln.pos > 0) ? (aln.pos) : 0;
    size_t end = (aln.endpos());
    //OpensRight to avoid confusion with member variable OpensLeft
    bool bOpensRight;
    if(isClip){
        bOpensRight = isLeft; 
    } else {
        bOpensRight = (aln.bRev) ? false : true;
    }
    //Get the information on whehter the terminal bases of the read are mapped
    uint8_t clipSide = aln.clipSide();
    if(isClip){
        clipSide |= (isLeft) ? CXA::RIGHT_CLIPPED : CXA::LEFT_CLIPPED;
    }
    //Named to avoid confusion with member variable
    //FIXME: terminalSide isn't working as intended
    CXA::CLIP_SIDE terminalSide = (aln.bRev) ? CXA::RIGHT_CLIPPED : CXA::LEFT_CLIPPED;
    CXA::CLIP_SIDE distalSide = bOpensRight ? CXA::LEFT_CLIPPED : CXA::RIGHT_CLIPPED;
    bool bDistalNotTerminal = (terminalSide != distalSide) || (clipSide & distalSide);
    //std::cerr << "\t" << this->to_bedpe() << "\n";
    if(name != "New"){
        ChimericFragment_t tmp("New");
        tmp.add_alignment(aln,isClip,isLeft,isR1,ivIdx);
        std::cerr << "\t" << tmp.to_bedpe(true) << "\n";
    }
    //std::cerr << "\t" << aln.chr << "\t" << "\t" << off << "-" << end << "\t" << aln.bRev << "\t" << bOpensRight << "\t" << isClip << "\t" << isLeft << "\t" << isR1 << "\t" << ivIdx <<"\t" << int(clipSide) << "\t"  << bDistalNotTerminal <<"\n";
    //std::cerr << "\t*" << this->chr[ivIdx] << "*\t*" << (flag & OPENS_LEFT) << "*\n";
    if(!(flag & HAS_INTERVAL)){
        if(name != "New") std::cerr << "\tNEW\n";
        //First alignment on this side, take it as is
        this->chr[ivIdx] = aln.chr;
        if(!bOpensRight) { set_bit(flag,OPENS_LEFT); }
        this->off[ivIdx] = off;
        this->end[ivIdx] = end;
        set_bit(flag,isR1 ? FROM_R1 : FROM_R2);
        //If this alignment is a clip, then the fragment has a split read
        //Otherwise no update
        if(isClip) {set_bit(flag, IS_SPLIT); }
        set_bit(flag,HAS_INTERVAL);
        if(!bDistalNotTerminal) {set_bit(flag,DISTAL_IS_TERMINAL); }
        return true;
    }
    std::cerr << "\tUPDATE";
    //There is something on this side already
    if(aln.chr != this->chr[ivIdx]) {
        std::cerr << "\tDIFF Contig\n";
        return false; //This alignment is completely inconsistent
    }
    bool bChange = false;
    bool bProperPair = false;
    //Check if the new alignment disagrees on the direction of the breakpoint
    if(bOpensRight == bool(flag & OPENS_LEFT)) {
        std::cerr << "\tConflicting BP side";
        //Check if the new alignment forms a discordant pair
        if( (bOpensRight && (end >= this->end[ivIdx])) || //New says open right and is right
            (!bOpensRight && (off <= this->off[ivIdx])) ) //New says open left and is left
        { 
            std::cerr << "\tDiscordant\n";
            //TODO: Account for clipping (the aligned ends of the old might be mapped in the new)
            return false;
        }
        std::cerr << "\tConcordant";
        //The new alignment is concordant
        bProperPair = true;
        //The distal end of this alignment is considered proximal by the existing fragment
        //If the distal end is terminal, then the fragment's proximal end is terminal
        //Indicating a non-chimeric fragment
        if(!bDistalNotTerminal && !(flag & PROXIMAL_IS_TERMINAL)){
            set_bit(flag,PROXIMAL_IS_TERMINAL);
            bChange = true;
        }
    }
    std::cerr << "\tAdding\n";
    //track whether incorporating the alignment changes the fragment
    //Update information with the added alignment
    if(!(flag & IS_SPLIT) & isClip){
        set_bit(flag, IS_SPLIT);
        bChange = true;
    }
    //If adding another read doesn't otherwise change the fragment
    // it doesn't need to be carried along
    set_bit(flag,isR1 ? FROM_R1 : FROM_R2);
    if(off <= this->off[ivIdx]){ //Extend to the left / check for terminality
        if(off < this->off[ivIdx]) {
            this->off[ivIdx] = off;
            bChange = true;
        }
        if(!bProperPair){//proper paired alignments can't inform on the distal end
            if(!(flag & OPENS_LEFT) && !bDistalNotTerminal) {
                set_bit(flag, DISTAL_IS_TERMINAL);
                bChange = true;
            }
        }
    }
    if(end >= this->end[ivIdx]){ //Extend to the right / check for terminality
        if(end >= this->end[ivIdx]){
            this->end[ivIdx] = end;
            bChange = true;
        }
        if(!bProperPair){ //proper paired alignments can't inform on the distal end
            if((flag & OPENS_LEFT) && !bDistalNotTerminal) {
                set_bit(flag, DISTAL_IS_TERMINAL);
                bChange = true;
            }
        }
    }
    //std::cerr << "\tTerm Add\n";
    return bChange;
}


//Determine if a fragment is strictly better than another
//  If two fragments have the same configuration, and the same distal positions,
//  the fragment with the more proximal positions (up to being split) is dominant
//  Given a tie, the fragment with more sources is better
//Inputs - A chimericFragment to which to compare
//Output -  -1 if this fragment is dominated by the other
//          0 if incomparable or neither dominates
//          1 if this fragment dominates the other
int ChimericFragment_t::dominant_comparison (const ChimericFragment_t & other) {
    //TODO: Implement
    return 0;
}

bool ChimericFragment_t::not_chimeric() const {
    for(IV_IDX ivIdx : {IV1 , IV2}) {
        std::cerr << bool(flag[ivIdx] & (HAS_INTERVAL)) << "\t"  << bool(flag[ivIdx] & (PROXIMAL_IS_TERMINAL)) << "\t" << (flag[ivIdx] & (HAS_INTERVAL | PROXIMAL_IS_TERMINAL)) << "\t" << HAS_INTERVAL << flag[ivIdx] << "\n";
        //Check if the interval is present and the proximal position is terminal
        if((flag[ivIdx] & (HAS_INTERVAL | PROXIMAL_IS_TERMINAL)) > HAS_INTERVAL){
            return true;
        }
    }
    return false;
}

std::string ChimericFragment_t::to_bedpe(bool bitflag) const {
    std::string line("");
    ////Empty bedpe entry for incomplete fragments
    //if(bEmptyIncomplete && !this->is_complete()){
    //    return line;
    //}
    for(size_t ivIdx : {IV1, IV2}) {
        if(!(flag[ivIdx] & HAS_INTERVAL)){
            line += ".\t.\t.\t";
            continue;
        }
        line += chr[ivIdx] + "\t" + std::to_string(off[ivIdx]) + "\t" +
                std::to_string(end[ivIdx]) + "\t";
    }
    line += name + "\t";
    uint32_t comboFlag = (flag[IV1] << (FLAG_BITS + 1)) | flag[IV2];
    line += std::to_string(comboFlag) + "\t";
    //line += std::to_string(comboFlag) + ":" +std::to_string(flag[IV1]) + "," + std::to_string(flag[IV2]) + "\t";
    for(size_t ivIdx : {IV1, IV2}) {
        std::string strand = ".";
        if((flag[ivIdx] & HAS_INTERVAL)){
            strand = bool(flag[ivIdx] & OPENS_LEFT) ? "-" : "+";
        }
        line += strand;
        if(ivIdx != IV2){
            line += "\t";
        }
    }
    if(bitflag){
        line += "\t" + std::bitset<FLAG_BITS>(flag[IV1]).to_string() + "." +
                std::bitset<FLAG_BITS>(flag[IV2]).to_string();
    }
    return line;
}

#endif //CHIMERIC_FRAGMENT_H
