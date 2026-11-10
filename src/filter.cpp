#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <htslib/sam.h>

#include "config.h"
#include "sam_utils.h"
#include <ssw.h>
#include <ssw_cpp.h>
#include "utils.h"

// DECLARATIONS

using CallVector_t = std::vector<call_t>;

// GLOBAL VARIABLES

static double MaxAltPropDistance = 0.001;
static int MaxAltRawDistance = -1;
std::string ContigLengthsFileName;

// FUNCTIONS

void print_calls(CallVector_t & calls) {
    for (size_t i = 0; i < calls.size(); i++) {
        call_t call = calls[i];
        if (!call.removed) {
            std::cout << call.to_human_string() << std::endl;
        }
    }
}

StripedSmithWaterman::Alignment align(StripedSmithWaterman::Aligner& aligner, std::string& ref, std::string& query, bool align_rc) {
    std::string rc_query = query;
    get_rc(rc_query);

    StripedSmithWaterman::Alignment aln;
    StripedSmithWaterman::Filter filter;
    if (align_rc) {
        aligner.Align(rc_query.c_str(), ref.c_str(), ref.length(), filter, &aln, 0);
    } else {
        aligner.Align(query.c_str(), ref.c_str(), ref.length(), filter, &aln, 0);
    }
    return aln;
}


//Function to identify 'accepted' calls which are likely to be alternate mapping of other calls
//  This is done with the assumption that a particular genomic location is likely to to have only
//  one junction with a given configuration, so additional junctions in the same vicinity may be
//  alternate mappings with other reads
//Inputs - a vector of calls to filter, and a vector of calls wihin which to store identified junctions
//       - Also uses global MaxAltPropDistance and MaxAltDistance
//Output - None, modifies the inputs
void IdentifyAlternates(CallVector_t & accepted_calls,
                        CallVector_t & alternate_calls)
{
    bool rawFilter = (MaxAltRawDistance >= 0);
    bool propFilter = (MaxAltPropDistance >= 0);
    //The Filter is disabled if maxAlts are both negative
    if(!rawFilter && !propFilter) { return; }
    //Need the contig lengths in order to calculate proportions
    std::unordered_map<std::string,size_t> contigLengthMap;
    if(propFilter) {
        std::ifstream in(ContigLengthsFileName);
        std::string contigName;
        size_t length;
        while(in >> contigName >> length) {
            contigLengthMap[contigName] = length;
        }
    }
    //The First call is always accepted, so we start at the second call
    for(size_t i = 1; i < accepted_calls.size(); i++) {
        call_t & call_i = accepted_calls[i];
        if(call_i.removed) { continue; }
        //Iterate over all prior (non-alt) calls
        for(size_t j = 0; j < i; j++) {
            call_t & call_j = accepted_calls[j];
            //Skip removed prior calls
            // NOTE: As we use this to remove alts, this can lead to
            // A is alt of B, and C is alt of A, but C is not an alt of B, so C is allowed through
            if(call_j.removed) { continue; }
            //Determine the raw and proportional distance(s) between
            // bp positions on the same contig strand
            std::vector<double> rawDists, propDists;
            rawDists.reserve(2);
            if (propFilter) {propDists.reserve(2); }
            using BPPair_t = std::pair<breakpoint_t,breakpoint_t>;
            for(const BPPair_t & bpPair :
                    {   BPPair_t{call_i.host_bp,call_j.host_bp},
                        BPPair_t{call_i.virus_bp,call_j.virus_bp} } )
            {
                if( (bpPair.first.chr != bpPair.second.chr) ||
                    (bpPair.first.rev != bpPair.second.rev) )
                { continue; }
                double dist = std::abs(bpPair.first.start - bpPair.second.start);
                // no need to save the distance if we aren't raw filtering
                if( rawFilter) { rawDists.push_back(dist); }
                if ( !propFilter) { continue; }
                double propDist = dist / contigLengthMap.at(bpPair.first.chr);
                propDists.push_back(propDist);
            }
            //Process both the raw and proportional distances
            // Note: An empty vector is the result of either the filter being disabled
            //  or the calls being incomparable
            //Both filters must say filter for the junction to be filtered
            bool bFilter = false;
            for(std::vector<double> * dists_ptr : {&rawDists,&propDists}) {
                //If the vector is empty, skip that filter
                if(!dists_ptr->size()) { continue; }
                bFilter = true;
                //Get the minimum distance
                std::sort(dists_ptr->begin(),dists_ptr->end());
                double dist = dists_ptr->front();
                //Determine which threshold is being used
                double thresh = (dists_ptr == &rawDists) ?  MaxAltRawDistance :
                                                            MaxAltPropDistance;
                if(dist > thresh) {
                    bFilter = false;
                    break;
                }
            }
            //Move along if filtering isn't requested
            if( !bFilter ) { continue; }
            //Filter Call_J
            //Copy Call_I as an alt of call_J, then mark the original removed
            alternate_calls.push_back(call_i);
            alternate_calls.back().paired_with = call_j.id;
            call_i.removed = true;
            //An alt was found for i, so we don't need to consider any additional j's
            break;
        }
    }
}




int main(int argc, char* argv[]) {

    std::string workdir = argv[1];

    ContigLengthsFileName = workdir + "/contig_lengths.tab";
    config_t config = parse_config(workdir + "/config.txt");
    double min_pairs = config.min_frag;
    double split_bonus = config.split_bonus;
    double min_seq_len = double(config.read_len)/2.0;

    std::vector<std::string> args;
    bool    print_rejected = false, accept_clipped = false,
            remapped = false, print_alternates = false;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--print-rejected") == 0) {
            print_rejected = true;
        } else if (strcmp(argv[i], "--print-alternates") == 0) {
            print_alternates = true;
        } else if (strcmp(argv[i], "--accept-all-clipped") == 0) {
            accept_clipped = true;
        } else if (strcmp(argv[i], "--remapped") == 0) {
            remapped = true;
        } else {
            args.push_back(argv[i]);
        }
    }

    std::ifstream fin(workdir + (remapped ? "/results.remapped.txt" : "/results.txt"));


    double min_host_pbs = 0.8;
    if (args.size() >= 2) min_host_pbs = std::stod(args[2]);

    double min_cov_for_accept = 0.5;
    if (args.size() >= 3) min_cov_for_accept = std::stod(args[3]);


    // read host and virus sequences
    std::string host_bp_seqs_fname = workdir + "/host_bp_seqs.fa", virus_bp_seqs_fname = workdir + "/virus_bp_seqs.fa";
    FILE* host_bp_fasta = fopen(host_bp_seqs_fname.c_str(), "r"), *virus_bp_fasta = fopen(virus_bp_seqs_fname.c_str(), "r");
    kseq_t* hseq = kseq_init(fileno(host_bp_fasta)), *vseq = kseq_init(fileno(virus_bp_fasta));
    std::unordered_map<size_t,std::pair<std::string, std::string> > host_virus_bp_seqs;
    std::unordered_map<size_t,std::pair<int, int> > host_virus_bp_Ns;
    while (kseq_read(hseq) >= 0 && kseq_read(vseq) >= 0) {
        std::string name(hseq->name.s);
        size_t callid = std::stoul(name.substr(0,name.length()-2));
        host_virus_bp_seqs[callid] = {hseq->seq.s, vseq->seq.s};
        host_virus_bp_Ns[callid] = {0, 0};
    }

    std::string host_bp_seqs_masks_fname = workdir + "/host_bp_seqs.masked.bed",
                virus_bp_seqs_masks_fname = workdir + "/virus_bp_seqs.masked.bed";
    std::ifstream host_bp_masks_fin(host_bp_seqs_masks_fname), virus_bp_masks_fin(virus_bp_seqs_masks_fname);
    std::vector<std::pair<int, int> > host_virus_bp_seqs_masked;
    std::string id;
    int start, end;
    while (host_bp_masks_fin >> id >> start >> end) {
        size_t callid = std::stoul(id.substr(0,id.length()-2));
        host_virus_bp_Ns[callid].first = end-start+1;
    }
    while (virus_bp_masks_fin >> id >> start >> end) {
        size_t callid = std::stoul(id.substr(0,id.length()-2));
        host_virus_bp_Ns[callid].second = end-start+1;
    }

    CallVector_t accepted_calls, rejected_calls, alternate_calls;
    std::string line;
    while (std::getline(fin, line)) {
        call_t call(line);

        bool accept = true;
        if (call.host_pbs < min_host_pbs) accept = false;
        if (    (call.good_pairs < min_pairs-split_bonus) ||
                (call.good_pairs < min_pairs && call.split_reads == 0) )
        { accept = false; }

        std::string host_seq = host_virus_bp_seqs[call.id].first;
        std::string virus_seq = host_virus_bp_seqs[call.id].second;
        if (host_seq.length() < min_seq_len || virus_seq.length() < min_seq_len) {
            accept = false;
        }

        if (call.coverage() >= min_cov_for_accept && call.good_pairs > 2) accept = true;
        if (accept_clipped && call.split_reads > 0) accept = true;

        int host_Ns = host_virus_bp_Ns[call.id].first;
        int virus_Ns = host_virus_bp_Ns[call.id].second;
        if (double(host_Ns)/host_seq.length() > 0.80 || double(virus_Ns)/virus_seq.length() > 0.80) {
            accept = false;
        }

        if (accept) {
            accepted_calls.push_back(call);
        } else {
            rejected_calls.push_back(call);
        }
    }


    // find pairs
    for (size_t i = 0; i < accepted_calls.size(); i++) {
        call_t a_call1 = accepted_calls[i];

        int pair_j = -1, min_dist = INT32_MAX;
        for (size_t j = i+1; j < accepted_calls.size(); j++) {
            call_t a_call2 = accepted_calls[j];
            if (a_call2.is_paired()) continue;

            int dist = pair_dist(a_call1, a_call2);
            if (abs(dist) < min_dist) {
                pair_j = j;
                min_dist = abs(dist);
            }
        }

        if (pair_j != -1) {
            accepted_calls[i].paired_with = accepted_calls[pair_j].id;
            accepted_calls[pair_j].paired_with = accepted_calls[i].id;
        }
    }

    for (size_t i = 0; i < accepted_calls.size(); i++) {
        call_t a_call = accepted_calls[i];
        if (a_call.is_paired()) continue;

        int pair_j = -1, min_dist = INT32_MAX;
        for (size_t j = 0; j < rejected_calls.size(); j++) {
            if (rejected_calls[j].is_paired()) continue;

            call_t r_call = rejected_calls[j];
            int dist = pair_dist(a_call, r_call);
            if (abs(dist) < min_dist) {
                pair_j = j;
                min_dist = abs(dist);
            }
        }

        if (pair_j != -1) {
            a_call.paired_with = rejected_calls[pair_j].id;
            rejected_calls[pair_j].paired_with = a_call.id;
            accepted_calls.push_back(rejected_calls[pair_j]);
        }
    }


    // find duplicates
    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
    for (size_t j = 0; j < accepted_calls.size(); j++) {
        std::string hseq_j = host_virus_bp_seqs[accepted_calls[j].id].first;
        std::string vseq_j = host_virus_bp_seqs[accepted_calls[j].id].second;

        for (size_t i = 0; i < j; i++) {
            std::string hseq_i = host_virus_bp_seqs[accepted_calls[i].id].first;
            std::string vseq_i = host_virus_bp_seqs[accepted_calls[i].id].second;

            std::string h_ref = hseq_i.length() >= hseq_j.length() ? hseq_i : hseq_j;
            std::string h_query = hseq_i.length() >= hseq_j.length() ? hseq_j : hseq_i;
            std::string v_ref = vseq_i.length() >= vseq_j.length() ? vseq_i : vseq_j;
            std::string v_query = vseq_i.length() >= vseq_j.length() ? vseq_j : vseq_i;

            // TODO: if the two calls are on opposite strand, align rc, otherwise, no
            bool h_rc = accepted_calls[i].host_bp.rev != accepted_calls[j].host_bp.rev;
            bool v_rc = accepted_calls[i].virus_bp.rev != accepted_calls[j].virus_bp.rev;
            StripedSmithWaterman::Alignment h_aln = align(aligner, h_ref, h_query, h_rc),
                                            v_aln = align(aligner, v_ref, v_query, v_rc);
            bool same_hseq = h_aln.query_end-h_aln.query_begin >= h_query.length()*0.8;
            bool same_vseq = v_aln.query_end-v_aln.query_begin >= v_query.length()*0.8;

            if (same_hseq && same_vseq) {
                accepted_calls[j].removed = true;
                break;
            }
        }
    }

    // No additional work if printing rejected

    if (print_rejected) {
        print_calls(rejected_calls);
        return 0;
    }

    // find alternates

    IdentifyAlternates(accepted_calls,alternate_calls);

    if (print_alternates) {
        print_calls(alternate_calls);
    } else {
        print_calls(accepted_calls);
    }
}
