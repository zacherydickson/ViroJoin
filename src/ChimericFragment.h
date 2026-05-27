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
class ChimericFragment_t {
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
    //Members
    protected:
    std::string name;
    std::array<std::string,2> chr;
    std::array<size_t,2> off;
    std::array<size_t,2> end;
    std::array<uint16_t,2> flag;
    //Con-/Destruction
    public:
    ChimericFragment_t(const std::string name = "") :
        name(name),
        chr({"",""}), off({0,0}), end({0,0}), flag({0,0})
    {}
    ChimericFragment_t(const ChimericFragment_t & other) :
        name(other.name), chr(other.chr), off(other.off),
        end(other.end), flag(other.flag) 
    {
    }
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
    size_t n_reads(IV_IDX ivIdx) const {
        return bool(flag[ivIdx] & FROM_R1) + bool(flag[ivIdx] & FROM_R2);
    }
    size_t n_split() const;
    public:
    size_t distal_pos(IV_IDX ivIdx) const {
        return (flag[ivIdx] & OPENS_LEFT) ? end[ivIdx] : off[ivIdx] + 1;
    }
    std::array<size_t,2> distal_pos() const {
        return {this->distal_pos(IV1),this->distal_pos(IV2)};
    }
    bool is_complete() const { return this->both(DISTAL_IS_TERMINAL); }
    bool is_chimeric() const { return this->both(HAS_INTERVAL); }
    bool is_split(IV_IDX ivIdx) const { return this->flag[ivIdx] & IS_SPLIT; }
    bool not_chimeric() const;
    size_t proximal_pos(IV_IDX ivIdx) const { 
        return (flag[ivIdx] & OPENS_LEFT) ? off[ivIdx]+1 : end[ivIdx];
    }
    std::array<size_t,2> proximal_pos() const {
        return {this->proximal_pos(IV1),this->proximal_pos(IV2)};
    }
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
    bool add_alignment( const CXA & aln, bool isClip, bool isAnchor,
                        bool isLeft, bool isR1, uint8_t clipSide,
                        ChimericFragment_t::IV_IDX ivIdx);
    int dominant_comparison (const ChimericFragment_t & other) const;
    std::string to_bedpe(bool bitflag = false) const ;
    //Static Methods
    ChimericFragment_t from_bedpe(const std::string bedpe) {
        auto fields = strsplit(bedpe,'\t');
        if(fields.size() < 8){
            throw std::invalid_argument("Attempt to construct Chimeric Fragment from incomplete bedpe entry");
        }
        ChimericFragment_t frag(fields[6]);
        frag.chr = {fields[0],fields[3]};
        size_t vals[5];
        int idx = 0;
        for(int i : {1,2,4,5,7}) {
            vals[idx++] = (fields[i] == ".") ? 0 : std::stoul(fields[i]);
        }
        frag.off = {vals[0],vals[1]};
        frag.end = {vals[2],vals[3]};
        frag.flag = {   uint16_t(vals[4] >> (FLAG_BITS + 1)),
                        uint16_t(vals[4] & ((1 << FLAG_BITS) - 1))
        };
        return frag;
    }
};

//Output - true if the alignment was successfully added, false otherwise
//          an alignment fails to be added if a fragment cannot be generated
//          which is consistent with the existing fragment and the alignment
bool ChimericFragment_t::add_alignment( const CXA & aln, bool isClip, bool isAnchor,
                                        bool isLeft, bool isR1, uint8_t clipSide,
                                        ChimericFragment_t::IV_IDX ivIdx) 
{
    uint16_t & flag = this->flag[ivIdx];
    //Step one: Assign up and down based on the host sequence
    //  this collapses H+V+ and V-H-, etc. together in the end
    //aln coords are zero indexed
    //bed uses half open intervals
    //off should be inclusive 0 indexed (or exclusive 1 indexed)
    //end should be exclusive 0 indexed (i.e inclusive 1 indexed)
    size_t off = (aln.pos > 0) ? (aln.pos) : 0;
    size_t end = (aln.endpos());
    //Gt information on where the breakpoint is relative to the alignment
    bool bOpensRight;
    if(isClip || isAnchor){
        bOpensRight = isLeft; 
        if(isClip && aln.bRev) {
            bOpensRight = !bOpensRight;
        }
    } else {
        bOpensRight = (aln.bRev) ? false : true;
    }
    //Get the information on whehter the terminal bases of the read are mapped
    CXA::CLIP_SIDE terminalSide;
    if(isClip){ 
        terminalSide = (bOpensRight) ? CXA::LEFT_CLIPPED : CXA::RIGHT_CLIPPED;
    } else {
        terminalSide = (aln.bRev) ? CXA::RIGHT_CLIPPED : CXA::LEFT_CLIPPED;
    }
    CXA::CLIP_SIDE distalSide = bOpensRight ? CXA::LEFT_CLIPPED : CXA::RIGHT_CLIPPED;
    //Named in a prior state to minimize confusion :P
    bool bDistalNotTerminal = (terminalSide != distalSide) || (clipSide & distalSide);
    if(!(flag & HAS_INTERVAL)){
        //First alignment on this side, take it as is
        this->chr[ivIdx] = aln.chr;
        if(!bOpensRight) { set_bit(flag,OPENS_LEFT); }
        this->off[ivIdx] = off;
        this->end[ivIdx] = end;
        set_bit(flag,isR1 ? FROM_R1 : FROM_R2);
        //If this alignment is a clip, then the fragment has a split read
        //Otherwise no update
        if(isClip || isAnchor) {set_bit(flag, IS_SPLIT); }
        set_bit(flag,HAS_INTERVAL);
        if(!bDistalNotTerminal) {set_bit(flag,DISTAL_IS_TERMINAL); }
        return true;
    }
    //There is something on this side already
    if(aln.chr != this->chr[ivIdx]) {
        return false; //This alignment is completely inconsistent
    }
    bool bChange = false;
    bool bProperPair = false;
    //Check if the new alignment disagrees on the direction of the breakpoint
    if(bOpensRight == bool(flag & OPENS_LEFT)) {
        //Check if the new alignment forms a discordant pair
        if( (bOpensRight && (end > this->end[ivIdx])) || //New says open right and is right
            (!bOpensRight && (off < this->off[ivIdx])) ) //New says open left and is left
        { 
            //TODO: Account for clipping (the aligned ends of the old might be mapped in the new)
            return false;
        }
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
    //track whether incorporating the alignment changes the fragment
    //Update information with the added alignment
    if(!(flag & IS_SPLIT) & (isClip || isAnchor)){
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
    return bChange;
}


//Determine if a fragment is strictly better than another
//  If two fragments have the same configuration, and the same distal positions,
//  the fragment with the more proximal positions (up to being split) is dominant
//  Given a tie, the fragment with more sources is better
//Inputs - A chimericFragment to which to compare
//Output -  -1 if this fragment is dominated by the other
//          0 if incomparable 
//          1 if this fragment dominates the other, or the fragments are equal
//NOTE: If either fragment is not complete (!is_complete) or
//  not chimeric (not_chimeric), behaviour is undefined
int ChimericFragment_t::dominant_comparison (const ChimericFragment_t & other) const 
{
    for(IV_IDX ivIdx : {IV1, IV2}){
        //Chr must match to be compared
        if(this->chr[ivIdx] != other.chr[ivIdx]){ return 0; }
        //Distal positions must match to be compared
        if(this->distal_pos(ivIdx) != other.distal_pos(ivIdx)){ return 0; }
        //Configurations must be the same to be be compared
        if( (this->flag[ivIdx] & OPENS_LEFT) !=
            (other.flag[ivIdx] & OPENS_LEFT) )
        {
            return 0;
        }
    }
    size_t nSplits[2] = { this->n_split(), other.n_split() };
    //The fragments are comparable
    //If one fragment has more split intervals than the other
    if(nSplits[0] != nSplits[1]) {
        //The one with more splits is domininant
        return (nSplits[0] > nSplits[1]) ? 1 : -1;
    }
    //The fragments have the same number of split intervals
    //If the fragments don't have the same split status on IV1
    if((this->flag[IV1] & IS_SPLIT) != (other.flag[IV1] & IS_SPLIT)){
        //then each fragment is split on opposite intervals
        //These are incomparable
        return 0;
        //Note: If a set of alignments can produce both such fragments,
        //Then it should be able to make the combined fragment which
        //would dominate both of these
    }
    //The two fragments have exactly the same split statuses
    //Assess whether the proximal positions for split intervals match
    for(IV_IDX ivIdx : {IV1, IV2}){
        //Skip unsplit intervals
        if(!(flag[ivIdx] & IS_SPLIT)) { continue; }
        //If the proximal positions do not match
        if(this->proximal_pos(ivIdx) != other.proximal_pos(ivIdx)){
            //They imply different breakpoints and cannot be compared
            return 0;
        }
    }
    //Proximal positions on any split intervals match
    // Assess which fragment has more R1,R2 support on each interval
    bool OtherHasMoreSupport = !(   (this->n_reads(IV1) + this->n_reads(IV2)) >= 
                                    (other.n_reads(IV1) + other.n_reads(IV2)) );
    //If the fragments are both double split
    if(nSplits[0] == 2){
        //If this fragment has less support, it is dominated, otherwise it wins
        return (OtherHasMoreSupport) ?  -1 : 1;
    }
    //The two fragments each have at least one unsplit interval
    //Assess both intervals to see which has proximal positions closer to the junction
    int betterIV[2] = {0, 0}; //Assuming the two are equal
    for(IV_IDX ivIdx : {IV1, IV2}){
        //Skip any split intervals
        if(flag[ivIdx] & IS_SPLIT) { continue; }
        //Get Distal-proximal distances
        int dpDist[2] = { 
            std::abs(int(this->distal_pos(ivIdx)) - int(this->proximal_pos(ivIdx))), 
            std::abs(int(other.distal_pos(ivIdx)) - int(other.proximal_pos(ivIdx)))
        };
        //If the distances are different
        if(dpDist[0] != dpDist[1]){
            //The greater distance is more coverage
            betterIV[ivIdx] = (dpDist[0] > dpDist[1]) ? 1 : -1;
        }
    }
    // If The better interval is the same for all unsplit intervals
    if(betterIV[IV1] == betterIV[IV2]){
        //If a better interval was called
        if(betterIV[IV1] != 0){
            //Then the better fragment is found
            return betterIV[IV1];
        }
        //No calls were made
        return (OtherHasMoreSupport) ?  -1 : 1;
    }
    //The two intervals disagree
    //If at least one interval didn't make a call
    if(betterIV[IV1] * betterIV[IV2] == 0) {
        //The dominant IV is the one which was called
        return (betterIV[IV1] == 0) ? betterIV[IV2] : betterIV[IV1];
    }
    //They are incomparable
    return 0;
}

size_t ChimericFragment_t::n_split() const {
    size_t n = 0;
    for(const uint16_t & f : flag){
        if(f & IS_SPLIT) { n++; }
    }
    return n;
}

bool ChimericFragment_t::not_chimeric() const {
    for(IV_IDX ivIdx : {IV1 , IV2}) {
        //Check if the interval is present and botht he proximal and distal positions are terminal
        //I.e. the implied fragment doesn't span a breakpoint
        uint16_t mask = (HAS_INTERVAL | PROXIMAL_IS_TERMINAL | DISTAL_IS_TERMINAL);
        if((flag[ivIdx] & mask)  == mask){
            return true;
        }
    }
    return false;
}

std::string ChimericFragment_t::to_bedpe(bool bitflag) const {
    std::string line("");
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
    for(size_t ivIdx : {IV1, IV2}) {
        std::string strand = ".";
        if((flag[ivIdx] & HAS_INTERVAL)){
            strand = (bool(flag[ivIdx] & OPENS_LEFT) == (ivIdx == IV1)) ? "-" : "+";
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
