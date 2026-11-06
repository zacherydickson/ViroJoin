#ifndef SURVEYOR_CONFIG_H
#define SURVEYOR_CONFIG_H

#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

const int MIN_CLIP_LEN = 5;

struct config_t {
    int threads;
    int min_sc_size, max_sc_dist, max_alt_aln;
    int min_frag;
    int split_bonus;
    int read_len;
    std::string cram_reference;
    bool explore;
    bool clique;
    double max_used_frac;
};

config_t parse_config(std::string file) {
    std::unordered_map<std::string, std::string> config_params;
    std::ifstream fin(file);
    std::string name, value;
    while (fin >> name >> value) {
        config_params[name] = value;
    }
    fin.close();

    config_t config;
    //Integer configs
    std::unordered_map<std::string,int*> integerParams = {
        {"threads", &config.threads},
        {"min_sc_size", &config.min_sc_size},
        {"max_sc_dist", &config.max_sc_dist},
        {"max_alt_aln", &config.max_alt_aln},
        {"min_frag", &config.min_frag},
        {"split_bonus", &config.split_bonus},
        {"read_len", &config.read_len},
    };
    for(const auto & pair : integerParams) {
        auto it = config_params.find(pair.first);
        if(it == config_params.end()) {
            throw std::invalid_argument("Provided config file is missing " + pair.first);
        }
        *(pair.second) = stoi(it->second);
    }
    if (config_params.count("cram_reference")) {
        config.cram_reference = config_params["cram_reference"];
    }
    config.explore = config_params["explore"] == "True";
    config.clique = config_params["clique"] == "True";
    config.max_used_frac = stod(config_params["max_used_frac"]);
    return config;
};

int parse_config_threads(const std::string & file){
    return parse_config(file).threads;
}


struct stats_t {
    int max_is;
};

stats_t parse_stats(std::string file) {
    std::unordered_map<std::string, std::string> config_params;
    std::ifstream fin(file);
    std::string name, value;
    while (fin >> name >> value) {
        config_params[name] = value;
    }
    fin.close();

    stats_t stats;
    stats.max_is = stoi(config_params["max_is"]);
    return stats;
}

#endif //SURVEYOR_CONFIG_H
