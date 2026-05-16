#ifndef CHIMERIC_FRAGMENT_H
#define CHIMERIC_FRAGMENT_H

#include <string>
#include <htslib/sam.h>
#include "sam_utils.h"
#include <stdexcept>

//Structure for building and tracking potentially chimeric fragments
//Considered complete if the 5` end of both the forward and reverse reads
//  are mapped
//Considered chimeric if the host status of the up and downstream intervals are
//  not equal
struct ChimericFragment_t {
    //Static Members
    public:
    static const size_t UPSTREAM = 0;
    static const size_t DOWNSTREAM = 1;
    //Members
    public:
    std::string name;
    std::array<std::string,2> chr;
    std::array<size_t,2> off;
    std::array<size_t,2> end;
    std::array<char,2> strand;
    protected:
        std::array<std::string,2> cigar;
        std::array<bool,2> bHost;
        std::array<bool,2> bComplete;
        bool bSplit;
    //Con-/Destruction
    public:
    ChimericFragment_t() :
        chr({"",""}), off({0,0}), end({0,0}), strand({'\0','\0'}),
        cigar({"",""}), bHost({false,false}), bComplete({false,false}),
        bSplit(false)
    {}
    //Accessors
    public:
    bool is_complete() const {
        return (this->bComplete[UPSTREAM] && this->bComplete[DOWNSTREAM]);
    }
    bool is_chimeric() const {
        return this->bHost[UPSTREAM] != this->bHost[DOWNSTREAM];
    }
    bool is_split() const { return this->bSplit; }
    bool has_stream(bool bUp) const {
        return (this->chr[bUp] != "" && this->strand[bUp] != '\0');
    }
    std::string get_cigar(bool bUp) const { return this->cigar[bUp]; }
    //Mutators
    //Methods
    public:
    bool add_alignment( const CXA & aln, bool isClip,
                        bool isLeft, bool isR1, bool isHost);
    std::string to_bedpe() const ;
};

//Output - true if the alignment was successfully added, false otherwise
//          an alignment fails to be added if a fragment cannot be generated
//          which is consistent with the existing fragment and the alignment
bool ChimericFragment_t::add_alignment( const CXA & aln, bool isClip,
                                        bool isLeft, bool isR1, bool isHost) 
{
    //Step one: Assign up and down based on the host sequence
    //  this collapses H+V+ and V-H-, etc. together in the end
    bool bUpstream = (isHost);
    size_t off = (aln.pos > 0) ? (aln.pos - 1) : 0;
    size_t end = (aln.endpos());
    uint8_t clipFlag = aln.clipSide();
    if(isClip){
        clipFlag = (isLeft) ? 0b10 : 0b01;
    }
    char strand = (aln.bRev) ? '+' : '-';
    if(this->has_stream(bUpstream)){
        //First alignment on this side, take it as is
        this->chr[bUpstream] = aln.chr;
        this->strand[bUpstream] = strand;
        this->off[bUpstream] = off;
        this->end[bUpstream] = end;
        this->cigar[bUpstream] = get_cigar_code(aln.cigar,aln.nCigar);
        this->bHost[bUpstream] = isHost;
        this->bComplete[bUpstream] = (bUpstream) ?  !(clipFlag && 0b01) :
                                                    !(clipFlag && 0b10);
        this->bSplit = clipFlag;
        return true;
    }
    //There is something on this side already
    if(aln.chr != this->chr[bUpstream] || strand != this->strand[bUpstream]) {
        //The alignment on the side of the junction is not 
        // on the same contig and strand as previous data
        return false;
    }
    //The previous data is on the same contig+strand 
    if(off < this->off[bUpstream]){
        //TODO: Update CIGAR and offset
    }
    //TODO: Implement
    return true;
}

std::string ChimericFragment_t::to_bedpe() const {
    std::string line("");
    //Empty bedpe entry for incomplete fragments
    if(!this->is_complete()){
        return line;
    }
    for(size_t stream : {UPSTREAM, DOWNSTREAM}) {
        line = chr[stream] + "\t" + std::to_string(off[stream]) + "\t" +
                std::to_string(end[stream]) + "\t";
    }
    line += name + "\t" + ((bSplit) ? "1" : "0") + "\t";
    for(size_t stream : {UPSTREAM, DOWNSTREAM}) {
        line += strand[stream];
    }
    for(size_t stream : {UPSTREAM, DOWNSTREAM}) {
        line += cigar[stream];
    }
    line += ((bHost[UPSTREAM]) ? "1" : "0");
    return line;
}

#endif //CHIMERIC_FRAGMENT_H
