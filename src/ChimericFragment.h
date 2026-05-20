#ifndef CHIMERIC_FRAGMENT_H
#define CHIMERIC_FRAGMENT_H

#include <string>
#include <htslib/sam.h>
#include "sam_utils.h"
#include <stdexcept>

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
        DISTAL_IS_TERMINAL = 0x2,
        HAS_INTERVAL = 0x4,
        FROM_R1 = 0x8,
        FROM_R2 = 0x10,
    };
    static const size_t MATE_SHIFT = 3;
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
    bool bSplit;
    //Con-/Destruction
    public:
    ChimericFragment_t(const std::string name = "") :
        name(name),
        chr({"",""}), off({0,0}), end({0,0}), flag({0,0}),
        //bOpensLeft({false,false}),
        //cigar({"",""}), bDistalIsTerminal({false,false}),
        //bHasInterval({false,false}),
        bSplit(false)
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
    } public:
    bool is_complete() const { return this->both(DISTAL_IS_TERMINAL); }
    bool is_chimeric() const { return this->both(HAS_INTERVAL); }
    //bool is_chimeric() const {
    //    return this->bHost[IV1] != this->bHost[IV2];
    //}
    bool is_split() const { return this->bSplit; }
    
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
    std::string to_bedpe(bool bEmptyIncomplete = false) const ;
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
    CXA::CLIP_SIDE distalSide = bOpensRight ? CXA::LEFT_CLIPPED : CXA::RIGHT_CLIPPED;
    bool bDistalNotTerminal = clipSide & distalSide;
    //std::cerr << "\t" << this->to_bedpe() << "\n";
    //std::cerr << "\t" << aln.chr << "\t" << "\t" << off << "-" << end << "\t" << aln.bRev << "\t" << bOpensRight << "\t" << isClip << "\t" << isLeft << "\t" << isR1 << "\t" << ivIdx <<"\t" << int(clipSide) << "\t"  << bDistalNotTerminal <<"\n";
    //std::cerr << "\t*" << this->chr[ivIdx] << "*\t*" << (flag & OPENS_LEFT) << "*\n";
    if(!(flag & HAS_INTERVAL)){
        //std::cerr << "\tnostream\n";
        //First alignment on this side, take it as is
        this->chr[ivIdx] = aln.chr;
        if(!bOpensRight) { set_bit(flag,OPENS_LEFT); }
        this->off[ivIdx] = off;
        this->end[ivIdx] = end;
        set_bit(flag,isR1 ? FROM_R1 : FROM_R2);
        //If this alignment is a clip, then the fragment has a slip read
        //Otherwise no update
        this->bSplit |= isClip;
        set_bit(flag,HAS_INTERVAL);
        if(!bDistalNotTerminal) {set_bit(flag,DISTAL_IS_TERMINAL); }
        return true;
    }
    //There is something on this side already
    if( (aln.chr != this->chr[ivIdx]) || 
        (bOpensRight == bool(flag & OPENS_LEFT)) )
    {
        //The alignment is inconsistent with the breakpoint implied by
        // the current fragment information
        return false;
    }
    //track whether incorporating the alignment changes the fragment
    bool bChange = false;
    //Update information with the added alignment
    if(!this->bSplit & isClip){
        this->bSplit |= isClip;
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
        if(!(flag & OPENS_LEFT) && !bDistalNotTerminal) {
            set_bit(flag, DISTAL_IS_TERMINAL);
            bChange = true;
        }
    }
    if(end >= this->end[ivIdx]){ //Extend to the right / check for terminality
        if(end >= this->end[ivIdx]){
            this->end[ivIdx] = end;
            bChange = true;
        }
        if((flag & OPENS_LEFT) && !bDistalNotTerminal) {
            set_bit(flag, DISTAL_IS_TERMINAL);
            bChange = true;
        }
    }
    //std::cerr << "\tTerm Add\n";
    return bChange;
}

std::string ChimericFragment_t::to_bedpe(bool bEmptyIncomplete) const {
    std::string line("");
    //Empty bedpe entry for incomplete fragments
    if(bEmptyIncomplete && !this->is_complete()){
        return line;
    }
    for(size_t ivIdx : {IV1, IV2}) {
        if(!(flag[ivIdx] & HAS_INTERVAL)){
            line += ".\t.\t.\t";
            continue;
        }
        line += chr[ivIdx] + "\t" + std::to_string(off[ivIdx]) + "\t" +
                std::to_string(end[ivIdx]) + "\t";
    }
    line += name + "\t" + ((bSplit) ? "1" : "0") + "\t";
    for(size_t ivIdx : {IV1, IV2}) {
        if(!(flag[ivIdx] & HAS_INTERVAL)){
            line += ".\t";
            continue;
        }
        line += bool(flag[ivIdx] & OPENS_LEFT) ? "-" : "+";
        line += "\t";
    }
    for(size_t ivIdx : {IV1, IV2}) {
        if(!(flag[ivIdx] & HAS_INTERVAL)){
            line += ".\t";
            continue;
        }
        line += std::to_string((flag[ivIdx] & (FROM_R1 | FROM_R2)) >> MATE_SHIFT);
        if(ivIdx != IV2) {
            line += "\t";
        }
    }
    return line;
}

#endif //CHIMERIC_FRAGMENT_H
