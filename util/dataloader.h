#ifndef DATALOADER
#define DATALOADER

#include <fstream>
#include <iostream>
#include <random>
#include <chrono>
#include <assert.h>
#include <filesystem>
#include <algorithm>
#include "key.h"
#include <arpa/inet.h>
#include <malloc.h>
#include "genzipf.h"
#include "murmurhash/MurmurHash3.h"
// #include <set>
#include <unordered_map>
#include <set>
#include <unordered_set>
// Dataset directory
//#define SHALLA_PATH  "../data/dataset_hosefire_long/hosefile0_25"

//#define SHALLA_PATH  "../data/dataset_hosefire/dataset_hosefire_0_5/hosefile"
#define SHALLA_PATH  "../data/dataset_hosefire/dataset_hosefire_1_5/hosefire1_5_data"
#define HOSEFIRE_PATH  "../data/dataset_hosefire/dataset_hosefire_1_0/hosefire1_0_data"

#define YCSB_PATH  "../data/dataset_ycsb/"


#define SHALLA_PATH_SK  "../data/dataset_shalla/shalla_1_5/shalla_cost_1.5-"
// #define SHALLA_PATH_SK  "../data/dataset_skewness/skewness1_0/shalla_1_0-"
//#define SHALLA_PATH_SK  "../data/dataset_skewness/tenDataset/shalla_cost_1.5-"
// #define SHALLA_PATH_SK "../data/dataset_skewness/temp/shalla_cost_2.0-"
#define CAIDA_PATH  "../data/dataset_caida/"

#define YCSB_PATH_SK  "../data/dataset_shalla/shalla_2_5/shalla_cost2.0-"

#define HOSEFIRE_PATH_1  "../data/dataset_hosefire_0_25/hosefile0_25_data1"
#define HOSEFIRE_PATH_2  "../data/dataset_hosefire_0_25/hosefile0_25_data2"
#define HOSEFIRE_PATH_3  "../data/dataset_hosefire_0_25/hosefile0_25_data3"
#define HOSEFIRE_PATH_4  "../data/dataset_hosefire_0_25/hosefile0_25_data4"
#define HOSEFIRE_PATH_5  "../data/dataset_hosefire_0_25/hosefile0_25_data5"




#define RANDOM_KEYSTR_PATH "../util/randomKeyStr.txt"
#define RANDOM_COST_TYPE zipf

enum {uniform, hotcost, normal, zipf};

// data loader for datasets and random data
class dataloader{
    public:
        std::vector<Slice *> pos_keys_; // positive keys
        std::vector<Slice *> neg_keys_; // negative keys
        std::vector<Slice *> query_keys; // query keys
        std::vector<size_t> pos_keys_ip_; // positive keys
        std::vector<size_t> query_keys_ip_; // query keys
        std::set<uint32_t> pos_keys_ipaddr_set;
        ~dataloader();
        bool load(std::string data_name_, bool using_cost_);
        bool load(std::string data_name_, bool using_cost_, std::string epcho);
        bool loadRandomKey(int positives_number, int negatives_number, bool using_cost_);
        bool loadYCSB(bool using_cost_, std::string epcho);
        bool loadShalla(bool using_cost_, std::string epcho);
        bool loadHosefire(bool using_cost_, std::string epcho);
        bool loadShalla_sk(bool using_cost_, std::string epcho);
        bool loadNegtive();
        bool loadShalla_vulNeg(std::string name, bool using_cost_, std::string epcho);
        bool loadCAIDA(bool using_cost_, std::string epcho);
        bool loadCAIDA_adp(bool using_cost_, std::string epcho, double load_factor);
        bool loadZipf_adp(bool using_cost_, std::string epcho, double load_factor, double sknewness);
        bool loadZipf_semi(std::string epcho, double sknewness);
        bool handleCAIDA_adp(bool using_cost_, std::string epcho, double load_factor);
        bool genarateZipf(std::string epcho, double load_factor, double sknewness);
        bool handleCAIDA_adp_A_S(bool using_cost_, std::string epcho, double load_factor, size_t a);
        bool loadCAIDA_adp_A_S(bool using_cost_, std::string epcho, double load_factor, size_t a);
        bool testYCSB_adp(std::string epcho, double load_factor, double sknewness);
        bool loadYCSB_adp(std::string epcho, double load_factor, double sknewness);
        bool loadYCSB(std::string epcho, double sknewness);
        bool loadYCSB_semi(std::string epcho, double sknewness);
        bool loadShalla_semi(std::string epcho, double sknewness);

        private:
        
        
        bool loadHosefire_1(bool using_cost_, std::string epcho);
        bool loadHosefire_2(bool using_cost_, std::string epcho);
        bool loadHosefire_3(bool using_cost_, std::string epcho);
        bool loadHosefire_4(bool using_cost_, std::string epcho);
        bool loadHosefire_5(bool using_cost_, std::string epcho);
        

};

// generate random keys and costs
class KeyBuilder{
    public:
        KeyBuilder();
        std::string GetKeyStr();
        bool ReadKeys(std::vector<Slice *> &v, int start_position_);
        void GenKeyStrAndToFile();
        void GenKeysUniformCosts(std::vector<Slice *> &keys, int interval);
        void GenKeysHotCosts(std::vector<Slice *> &keys, double hotNumberpro, int hotcost, int coldcost);
        void GenKeysNormalCosts(std::vector<Slice *> &keys, int u, int d);
        void GenKeysZipfCosts(std::vector<Slice *> &keys, double a, double c);
    private:
        std::vector<std::string> key_strs;
};

uint32_t ipv4ToUint32(const std::string& ip) {
    struct in_addr addr;
    if (inet_pton(AF_INET, ip.c_str(), &addr) != 1) {
        throw std::invalid_argument("Invalid IPv4 address");
    }
    return ntohl(addr.s_addr);
}

bool dataloader::loadCAIDA_adp(bool using_cost_, std::string epcho, double load_factor){
    std::cout << "caida reading..."  << std::endl;
    int lf = load_factor * 100;
    std::ifstream is(std::string(CAIDA_PATH) + std::__cxx11::to_string(lf) + "/" + "0" + epcho + "query.txt", std::ios::binary);
    size_t pos_num = (1<<18) * load_factor;
    if(is){
        size_t ip;
        // size_t totalIpPair = 0;
        // size_t posIpPair = 0;
        // std::unordered_map<uint32_t, size_t> ipset = std::unordered_map<uint32_t, size_t>();
        while (is >> ip) {  
            try {
                // if (ipset.find(ip) == ipset.end()) {
                //     ipset.insert(std::make_pair(ip, 1));
                // } else {
                //     ipset[ip]++;                
                // }
                if (pos_keys_ip_.size() < pos_num) {
                    pos_keys_ip_.push_back(ip);
                } else {
                    query_keys_ip_.push_back(ip);
                }
            } catch (const std::exception& e) {
            }
        }
        pos_keys_ipaddr_set = std::set<uint32_t>(pos_keys_ip_.begin(), pos_keys_ip_.end());
        // auto it = ipset.begin();
        // // std::advance(it, ipset.size() / 2);
        // std::advance(it, pos_num);
        // for (auto it2 = ipset.begin(); it2 != it; ++it2) {
        //     Slice * key = new Slice();
        //     key->ipaddr = it2->first;
        //     key->str = std::to_string(it2->first);
        //     pos_keys_.push_back(key);
        // }
        // for (auto it2 = it; it2 != ipset.end(); ++it2) {
        //     Slice * key = new Slice();
        //     key->ipaddr = it2->first;
        //     key->str = std::to_string(it2->first);
        //     key->cost = using_cost_ ? it2->second : 1;
        //     neg_keys_.push_back(key);
        // }
        // std::cout << "ipv4: " << ipset.size() << std::endl;
        std::cout << "pos_keys_: " << pos_keys_ip_.size() << std::endl;
        // std::cout << "totalIpPair: " << totalIpPair << std::endl;
        std::cout << "query_keys: " << query_keys_ip_.size() << std::endl;
        // std::cout << "posIpPair: " << posIpPair << std::endl;
        // std::cout << "neg_keys_: " << neg_keys_.size() << std::endl;
        // for (const auto& slice : pos_keys_) {
        //     pos_keys_ipaddr_set.insert(slice->ipaddr);
        // }
        // query_keys.erase(
        //     std::remove_if(query_keys.begin(), query_keys.end(),
        //                 [&](Slice* key) {
        //                     return pos_keys_ipaddr_set.find(key->ipaddr) != pos_keys_ipaddr_set.end();
        //                 }),
        // query_keys.end());
        // ipset.clear();    
        // std::unordered_map<uint32_t, size_t>().swap(ipset);

        // std::string outputFilename = std::string(CAIDA_PATH) + "0" + epcho + "query.txt";
        // std::ofstream os(outputFilename);
        // if (!os.is_open()) {
        //     std::cerr << "Failed to create file: " << outputFilename << std::endl;
        //     return false;
        // }

        // for (const auto& key : query_keys) {
        //     os << key->ipaddr << std::endl;
        // }

        // os.close();
        // std::cout << "Query keys written to file: " << outputFilename << std::endl;
        is.close();
        return true;
    } else {
        std::cout << "Failed to open file: " << std::string(CAIDA_PATH) + std::__cxx11::to_string(lf) + "/" + "0" + epcho + "query.txt" << std::endl;
    }
    is.close();
    return false;
}

bool dataloader::testYCSB_adp(std::string epcho, double load_factor, double sknewness){
    std::cout << "YCSB reading..."  << std::endl;
    int lf = load_factor * 100;
    int s = sknewness * 100;
    std::cout << "skewness_" + std::__cxx11::to_string(s) + "/load_factor_"
    + std::__cxx11::to_string(lf) + "/" + "query" + epcho + ".txt" << endl;
    std::ifstream is(std::string(YCSB_PATH) + "skewness_" + std::__cxx11::to_string(s) + "/load_factor_"
    + std::__cxx11::to_string(lf) + "/" + "query" + epcho + ".txt", std::ios::binary);
    size_t pos_num = (1<<20) * load_factor;
    if(is){
        size_t ip;
        while (is >> ip) {  
            try {

                if (pos_keys_ipaddr_set.size() < pos_num) {
                    pos_keys_ipaddr_set.insert(ip);
                } else {
                    query_keys_ip_.push_back(ip);
                }
            } catch (const std::exception& e) {
            }
        }
        // for (auto & ip : query_keys_ip_) {
        //     if (pos_keys_ipaddr_set.find(ip) != pos_keys_ipaddr_set.end()) {
        //         std::cout << "ip: " << ip << std::endl;
        //         break;
        //     }
        // }
        std::cout << "pos_keys_: " << pos_keys_ipaddr_set.size() << std::endl;
        std::cout << "query_keys: " << query_keys_ip_.size() << std::endl;

        // std::cout << "posIpPair: " << posIpPair << std::endl;
        // std::cout << "neg_keys_: " << neg_keys_.size() << std::endl;
        is.close();
        return true;
    } else {
        std::cout << "Failed to open file: " << std::endl;
    }
    is.close();
    return false;
}

bool dataloader::loadYCSB(std::string epcho, double sknewness){
    std::cout << "YCSB reading..."  << std::endl;
    int lf = 95;
    int s = sknewness * 100;
    // std::cout << "skewness_" + std::__cxx11::to_string(s) + "/load_factor_"
    // + std::__cxx11::to_string(lf) + "/" + "query" + epcho + ".txt" << endl;
    std::ifstream is(std::string(YCSB_PATH) + "skewness_" + std::__cxx11::to_string(s) + "/load_factor_"
    + std::__cxx11::to_string(lf) + "/" + "query" + epcho + ".txt", std::ios::binary);
    size_t pos_num = (1<<20) * 0.95;
    std::unordered_map<uint32_t, size_t> ipset = std::unordered_map<uint32_t, size_t>();
    if(is){
        size_t ip;
        while (is >> ip) {  
            try {
                if (pos_keys_.size() < pos_num) {
                    Slice * key = new Slice();
                    key->ipaddr = ip;
                    key->str = std::to_string(ip);
                    pos_keys_.push_back(key);
                } else {
                    if (ipset.find(ip) == ipset.end()) {
                        ipset.insert(std::make_pair(ip, 1));
                    } else {
                        ipset[ip]++;                
                    }
                }
            } catch (const std::exception& e) {
            }
        }
        for (auto it = ipset.begin(); it != ipset.end(); ++it) {
            Slice * key = new Slice();
            key->ipaddr = it->first;
            key->str = std::to_string(it->first);
            key->cost = it->second;
            neg_keys_.push_back(key);
        }
        // pos_keys_ipaddr_set = std::set<uint32_t>(pos_keys_ip_.begin(), pos_keys_ip_.end());
        // for (auto & ip : query_keys_ip_) {
        //     if (pos_keys_ipaddr_set.find(ip) != pos_keys_ipaddr_set.end()) {
        //         std::cout << "ip: " << ip << std::endl;
        //         break;
        //     }
        // }
        // std::cout << "pos_keys_: " << pos_keys_ipaddr_set.size() << std::endl;
        // std::cout << "query_keys: " << query_keys_ip_.size() << std::endl;

        // std::cout << "posIpPair: " << posIpPair << std::endl;
        // std::cout << "neg_keys_: " << neg_keys_.size() << std::endl;
        is.close();
        return true;
    } else {
        std::cout << "Failed to open file: " << std::endl;
    }
    is.close();
    return false;
}

bool dataloader::loadYCSB_semi(std::string epcho, double sknewness){
    std::cout << "YCSB reading..."  << std::endl;
    int s = sknewness * 100;
    // std::cout << "skewness_" + std::__cxx11::to_string(s) + "/load_factor_"
    // + std::__cxx11::to_string(lf) + "/" + "query" + epcho + ".txt" << endl;
    std::ifstream is(std::string(YCSB_PATH) + "/semi/skewness_" + std::__cxx11::to_string(s) + "/query" + epcho + ".txt", std::ios::binary);
    std::unordered_map<uint32_t, size_t> ipset = std::unordered_map<uint32_t, size_t>();
    if(is){
        size_t ip;
        size_t pos_num;
        is>> pos_num;
        while (is >> ip) {  
            try {
                if (pos_keys_.size() < pos_num) {
                    Slice * key = new Slice();
                    key->ipaddr = ip;
                    key->str = std::to_string(ip);
                    pos_keys_.push_back(key);
                } else {
                    if (ipset.find(ip) == ipset.end()) {
                        ipset.insert(std::make_pair(ip, 1));
                    } else {
                        ipset[ip]++;                
                    }
                }
            } catch (const std::exception& e) {
            }
        }
        for (auto it = ipset.begin(); it != ipset.end(); ++it) {
            Slice * key = new Slice();
            key->ipaddr = it->first;
            key->str = std::to_string(it->first);
            key->cost = it->second;
            neg_keys_.push_back(key);
        }
        // pos_keys_ipaddr_set = std::set<uint32_t>(pos_keys_ip_.begin(), pos_keys_ip_.end());
        // for (auto & ip : query_keys_ip_) {
        //     if (pos_keys_ipaddr_set.find(ip) != pos_keys_ipaddr_set.end()) {
        //         std::cout << "ip: " << ip << std::endl;
        //         break;
        //     }
        // }
        // std::cout << "pos_keys_: " << pos_keys_ipaddr_set.size() << std::endl;
        // std::cout << "query_keys: " << query_keys_ip_.size() << std::endl;

        // std::cout << "posIpPair: " << posIpPair << std::endl;
        // std::cout << "neg_keys_: " << neg_keys_.size() << std::endl;
        is.close();
        return true;
    } else {
        std::cout << "Failed to open file: " << std::endl;
    }
    is.close();
    return false;
}

bool dataloader::loadYCSB_adp(std::string epcho, double load_factor, double sknewness){
    std::cout << "YCSB reading..."  << std::endl;
    int lf = load_factor * 100;
    int s = sknewness * 100;
    // std::cout << "skewness_" + std::__cxx11::to_string(s) + "/load_factor_"
    // + std::__cxx11::to_string(lf) + "/" + "query" + epcho + ".txt" << endl;
    std::ifstream is(std::string(YCSB_PATH) + "skewness_" + std::__cxx11::to_string(s) + "/load_factor_"
    + std::__cxx11::to_string(lf) + "/" + "query" + epcho + ".txt", std::ios::binary);
    size_t pos_num = (1<<20) * load_factor;
    if(is){
        size_t ip;
        while (is >> ip) {  
            try {

                if (pos_keys_ip_.size() < pos_num) {
                    pos_keys_ip_.push_back(ip);
                } else {
                    query_keys_ip_.push_back(ip);
                }
            } catch (const std::exception& e) {
            }
        }
        pos_keys_ipaddr_set = std::set<uint32_t>(pos_keys_ip_.begin(), pos_keys_ip_.end());
        // for (auto & ip : query_keys_ip_) {
        //     if (pos_keys_ipaddr_set.find(ip) != pos_keys_ipaddr_set.end()) {
        //         std::cout << "ip: " << ip << std::endl;
        //         break;
        //     }
        // }
        std::cout << "pos_keys_: " << pos_keys_ipaddr_set.size() << std::endl;
        std::cout << "query_keys: " << query_keys_ip_.size() << std::endl;

        // std::cout << "posIpPair: " << posIpPair << std::endl;
        // std::cout << "neg_keys_: " << neg_keys_.size() << std::endl;
        is.close();
        return true;
    } else {
        std::cout << "Failed to open file: " << std::endl;
    }
    is.close();
    return false;
}

bool dataloader::loadZipf_semi(std::string epcho, double sknewness){
    double load_factor = 0.95;
    std::cout << "zipf reading..."  << std::endl;
    int lf = load_factor * 100;
    int s = sknewness * 100;
    std::ifstream is(std::string(CAIDA_PATH) + "skewness_" + std::__cxx11::to_string(s) + "/lf_"
    + std::__cxx11::to_string(lf) + "/" + "0" + epcho + "query.txt", std::ios::binary);
    size_t pos_num = (1<<18) * load_factor;
    std::unordered_map<uint32_t, size_t> ipset = std::unordered_map<uint32_t, size_t>();
    if(is){
        size_t ip;
        // size_t totalIpPair = 0;
        // size_t posIpPair = 0;
        // std::unordered_map<uint32_t, size_t> ipset = std::unordered_map<uint32_t, size_t>();

        while (is >> ip) {  
            try {
                if (pos_keys_.size() < pos_num) {
                    Slice * key = new Slice();
                    key->ipaddr = ip;
                    key->str = std::to_string(ip);
                    pos_keys_.push_back(key);
                } else {
                    if (ipset.find(ip) == ipset.end()) {
                        ipset.insert(std::make_pair(ip, 1));
                    } else {
                        ipset[ip]++;                
                    }
                }
            } catch (const std::exception& e) {
            }
        }
        for (auto it = ipset.begin(); it != ipset.end(); ++it) {
            Slice * key = new Slice();
            key->ipaddr = it->first;
            key->str = std::to_string(it->first);
            key->cost = it->second;
            neg_keys_.push_back(key);
        } // auto it = ipset.begin();
        // // std::advance(it, ipset.size() / 2);
        // std::advance(it, pos_num);
        // for (auto it2 = ipset.begin(); it2 != it; ++it2) {
        //     Slice * key = new Slice();
        //     key->ipaddr = it2->first;
        //     key->str = std::to_string(it2->first);
        //     pos_keys_.push_back(key);
        // }
        // for (auto it2 = it; it2 != ipset.end(); ++it2) {
        //     Slice * key = new Slice();
        //     key->ipaddr = it2->first;
        //     key->str = std::to_string(it2->first);
        //     key->cost = using_cost_ ? it2->second : 1;
        //     neg_keys_.push_back(key);
        // }
        // std::cout << "ipv4: " << ipset.size() << std::endl;
        // std::cout << "pos_keys_: " << pos_keys_ip_.size() << std::endl;
        // // std::cout << "totalIpPair: " << totalIpPair << std::endl;
        // std::cout << "query_keys: " << query_keys_ip_.size() << std::endl;
        // std::cout << "posIpPair: " << posIpPair << std::endl;
        // std::cout << "neg_keys_: " << neg_keys_.size() << std::endl;
        // for (const auto& slice : pos_keys_) {
        //     pos_keys_ipaddr_set.insert(slice->ipaddr);
        // }
        // query_keys.erase(
        //     std::remove_if(query_keys.begin(), query_keys.end(),
        //                 [&](Slice* key) {
        //                     return pos_keys_ipaddr_set.find(key->ipaddr) != pos_keys_ipaddr_set.end();
        //                 }),
        // query_keys.end());
        // ipset.clear();    
        // std::unordered_map<uint32_t, size_t>().swap(ipset);

        // std::string outputFilename = std::string(CAIDA_PATH) + "0" + epcho + "query.txt";
        // std::ofstream os(outputFilename);
        // if (!os.is_open()) {
        //     std::cerr << "Failed to create file: " << outputFilename << std::endl;
        //     return false;
        // }

        // for (const auto& key : query_keys) {
        //     os << key->ipaddr << std::endl;
        // }

        // os.close();
        // std::cout << "Query keys written to file: " << outputFilename << std::endl;
        is.close();
        return true;
    } else {
        std::cout << "Failed to open file: " << std::endl;
    }
    is.close();
    return false;
}

bool dataloader::loadZipf_adp(bool using_cost_, std::string epcho, double load_factor, double sknewness){
    std::cout << "zipf reading..."  << std::endl;
    int lf = load_factor * 100;
    int s = sknewness * 100;
    std::ifstream is(std::string(CAIDA_PATH) + "skewness_" + std::__cxx11::to_string(s) + "/lf_"
    + std::__cxx11::to_string(lf) + "/" + "0" + epcho + "query.txt", std::ios::binary);
    size_t pos_num = (1<<18) * load_factor;
    if(is){
        size_t ip;
        // size_t totalIpPair = 0;
        // size_t posIpPair = 0;
        // std::unordered_map<uint32_t, size_t> ipset = std::unordered_map<uint32_t, size_t>();
        while (is >> ip) {  
            try {
                // if (ipset.find(ip) == ipset.end()) {
                //     ipset.insert(std::make_pair(ip, 1));
                // } else {
                //     ipset[ip]++;                
                // }
                if (pos_keys_ip_.size() < pos_num) {
                    pos_keys_ip_.push_back(ip);
                } else {
                    query_keys_ip_.push_back(ip);
                }
            } catch (const std::exception& e) {
            }
        }
        pos_keys_ipaddr_set = std::set<uint32_t>(pos_keys_ip_.begin(), pos_keys_ip_.end());
        // auto it = ipset.begin();
        // // std::advance(it, ipset.size() / 2);
        // std::advance(it, pos_num);
        // for (auto it2 = ipset.begin(); it2 != it; ++it2) {
        //     Slice * key = new Slice();
        //     key->ipaddr = it2->first;
        //     key->str = std::to_string(it2->first);
        //     pos_keys_.push_back(key);
        // }
        // for (auto it2 = it; it2 != ipset.end(); ++it2) {
        //     Slice * key = new Slice();
        //     key->ipaddr = it2->first;
        //     key->str = std::to_string(it2->first);
        //     key->cost = using_cost_ ? it2->second : 1;
        //     neg_keys_.push_back(key);
        // }
        // std::cout << "ipv4: " << ipset.size() << std::endl;
        std::cout << "pos_keys_: " << pos_keys_ip_.size() << std::endl;
        // std::cout << "totalIpPair: " << totalIpPair << std::endl;
        std::cout << "query_keys: " << query_keys_ip_.size() << std::endl;
        // std::cout << "posIpPair: " << posIpPair << std::endl;
        // std::cout << "neg_keys_: " << neg_keys_.size() << std::endl;
        // for (const auto& slice : pos_keys_) {
        //     pos_keys_ipaddr_set.insert(slice->ipaddr);
        // }
        // query_keys.erase(
        //     std::remove_if(query_keys.begin(), query_keys.end(),
        //                 [&](Slice* key) {
        //                     return pos_keys_ipaddr_set.find(key->ipaddr) != pos_keys_ipaddr_set.end();
        //                 }),
        // query_keys.end());
        // ipset.clear();    
        // std::unordered_map<uint32_t, size_t>().swap(ipset);

        // std::string outputFilename = std::string(CAIDA_PATH) + "0" + epcho + "query.txt";
        // std::ofstream os(outputFilename);
        // if (!os.is_open()) {
        //     std::cerr << "Failed to create file: " << outputFilename << std::endl;
        //     return false;
        // }

        // for (const auto& key : query_keys) {
        //     os << key->ipaddr << std::endl;
        // }

        // os.close();
        // std::cout << "Query keys written to file: " << outputFilename << std::endl;
        is.close();
        return true;
    } else {
        std::cout << "Failed to open file: " << std::endl;
    }
    is.close();
    return false;
}

bool dataloader::genarateZipf(std::string epcho, double load_factor, double sknewness){
    size_t pos_num = (1<<18) * load_factor;
    uint32_t key_len = 4;
    uint32_t flow_num = 36E4;
    uint32_t packet_num = 32E6;
    // uint32_t hash_seed = epcho[0] - '0';
    uint32_t hash_seed = std::chrono::steady_clock::now().time_since_epoch().count();
    int lf = load_factor * 100;
    int s = sknewness * 100;
    while (pos_keys_ipaddr_set.size() < pos_num) {
        size_t rank = rand() % flow_num + 1;
        uint32_t key;
        MurmurHash3_x86_32(&rank, key_len, hash_seed, &key);
        pos_keys_ipaddr_set.insert(key);
    }
    for (int i = 0; i < packet_num; ++i)
    {
        uint32_t rand_num = zipf_(sknewness, flow_num, i == 0);
        uint32_t key;
        MurmurHash3_x86_32(&rand_num, key_len, hash_seed, &key);
        if (pos_keys_ipaddr_set.find(key) == pos_keys_ipaddr_set.end()) query_keys_ip_.push_back(key);

    }
        std::cout << "pos_keys_: " << pos_keys_ipaddr_set.size() << std::endl;
        std::cout << "query_keys: " << query_keys_ip_.size() << std::endl;
    std::string outputFilename = std::string(CAIDA_PATH) + "skewness_" + std::__cxx11::to_string(s) + "/lf_"
    + std::__cxx11::to_string(lf) + "/" + "0" + epcho + "query.txt";
    std::filesystem::create_directories(outputFilename.substr(0, outputFilename.find_last_of('/')));
    std::ofstream os(outputFilename);
    if (!os.is_open()) {
        std::cerr << "Failed to create file: " << outputFilename << std::endl;
        return false;
    }
    for (const auto& ip : pos_keys_ipaddr_set) {
        os << ip << std::endl;
    }
    for (const auto& ip : query_keys_ip_) {
        os << ip << std::endl;
    }
    os.close();
    std::cout << "Query keys written to file: " << outputFilename << std::endl;
    return true;

}

bool dataloader::handleCAIDA_adp(bool using_cost_, std::string epcho, double load_factor){
    std::cout << "caida reading..."  << std::endl;
    std::ifstream is(std::string(CAIDA_PATH) + "0" + epcho + ".txt", std::ios::binary);
    size_t pos_num = (1<<18) * load_factor;
    if(is){
        std::string sor, des;
        size_t totalIpPair = 0;
        size_t posIpPair = 0;
        // std::unordered_map<uint32_t, size_t> ipset = std::unordered_map<uint32_t, size_t>();
        while (is >> sor >> des) {  
            try {
                uint32_t ip = ipv4ToUint32(sor);
                // if (ipset.find(ip) == ipset.end()) {
                //     ipset.insert(std::make_pair(ip, 1));
                // } else {
                //     ipset[ip]++;                
                // }
                Slice * key = new Slice();
                key->ipaddr = ip;
                key->str = std::to_string(ip);
                if (pos_keys_ipaddr_set.size() < pos_num) {
                    // pos_keys_.push_back(key);
                    pos_keys_ipaddr_set.insert(ip);
                    posIpPair++;
                } else {
                    if (pos_keys_ipaddr_set.find(ip) == pos_keys_ipaddr_set.end()) {
                        query_keys.push_back(key);
                    } else posIpPair++;
                }
                totalIpPair++;
            } catch (const std::exception& e) {
            }
        }
        std::cout << "pos_keys_: " << pos_keys_ipaddr_set.size() << std::endl;
        std::cout << "totalIpPair: " << totalIpPair << std::endl;
        std::cout << "query_keys: " << query_keys.size() << std::endl;
        std::cout << "posIpPair: " << posIpPair << std::endl;
        // std::cout << "neg_keys_: " << neg_keys_.size() << std::endl;
        // for (const auto& slice : pos_keys_) {
        //     pos_keys_ipaddr_set.insert(slice->ipaddr);
        // }
        // query_keys.erase(
        //     std::remove_if(query_keys.begin(), query_keys.end(),
        //                 [&](Slice* key) {
        //                     return pos_keys_ipaddr_set.find(key->ipaddr) != pos_keys_ipaddr_set.end();
        //                 }),
        // query_keys.end());
        // ipset.clear();    
        // std::unordered_map<uint32_t, size_t>().swap(ipset);
        int lf = load_factor * 100;
        std::string outputFilename = std::string(CAIDA_PATH) + std::__cxx11::to_string(lf) + "/" + "0" + epcho + "query.txt";
        std::ofstream os(outputFilename);
        if (!os.is_open()) {
            std::cerr << "Failed to create file: " << outputFilename << std::endl;
            return false;
        }
        for (const auto& ip : pos_keys_ipaddr_set) {
            os << ip << std::endl;
        }
        for (const auto& key : query_keys) {
            os << key->ipaddr << std::endl;
        }

        os.close();
        std::cout << "Query keys written to file: " << outputFilename << std::endl;
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::handleCAIDA_adp_A_S(bool using_cost_, std::string epcho, double load_factor, size_t a) {
    std::cout << "caida reading..."  << std::endl;
    std::ifstream is(std::string(CAIDA_PATH) + "0" + epcho + ".txt", std::ios::binary);
    size_t pos_num = (1<<12) * load_factor;
    size_t unique_query_num = pos_num * a;
    if(is){
        std::string sor, des;
        size_t totalIpPair = 0;
        size_t posIpPair = 0;
        std::unordered_set< size_t> querySet = std::unordered_set< size_t>();
        std::vector<size_t> query1;
        while (is >> sor >> des) {  
            try {
                uint32_t ip = ipv4ToUint32(sor);
                // if (ipset.find(ip) == ipset.end()) {
                //     ipset.insert(std::make_pair(ip, 1));
                // } else {
                //     ipset[ip]++;                
                // }
                // Slice * key = new Slice();
                // key->ipaddr = ip;
                // key->str = std::to_string(ip);
                if (pos_keys_ipaddr_set.size() < pos_num) {
                    // pos_keys_.push_back(key);
                    pos_keys_ipaddr_set.insert(ip);
                    posIpPair++;
                } else {
                    if (pos_keys_ipaddr_set.find(ip) == pos_keys_ipaddr_set.end()) {
                        if (querySet.size() < unique_query_num) {
                            querySet.insert(ip);
                            // query_keys.push_back(key);
                        }
                        // query_keys.push_back(key);
                        query1.push_back(ip);
                    } else posIpPair++;
                }
                totalIpPair++;
            } catch (const std::exception& e) {
            }
        }
        // std::cout << "keys_: " << ipset.size() << std::endl;    
        std::cout << "pos_keys_: " << pos_keys_ipaddr_set.size() << std::endl;
        std::cout << "querySe:  " << querySet.size() << std::endl;
        // std::cout << "totalIpPair: " << totalIpPair << std::endl;

        // std::cout << "posIpPair: " << posIpPair << std::endl;
        // std::cout << "neg_keys_: " << neg_keys_.size() << std::endl;
        // for (const auto& slice : pos_keys_) {
        //     pos_keys_ipaddr_set.insert(slice->ipaddr);
        // }
        std::vector<size_t> query2;
        for (size_t i = 0; i < query1.size(); i++) {
            if (querySet.find(query1[i]) != querySet.end()) {
                query2.push_back(query1[i]);
            }
        }
        // query_keys.erase(
        //     std::remove_if(query_keys.begin(), query_keys.end(),
        //                 [&](Slice* key) {
        //                     return querySet.find(key->ipaddr) == querySet.end();
        //                     // pos_keys_ipaddr_set.find(key->ipaddr) != pos_keys_ipaddr_set.end();
        //                 }),
        // query_keys.end());
        std::cout << "query_keys: " << query2.size() << std::endl;
        // ipset.clear();    
        // std::unordered_map<uint32_t, size_t>().swap(ipset);

        int lf = load_factor * 100;
        std::string outputFilename = std::string(CAIDA_PATH) + std::__cxx11::to_string(lf) + "/" 
        + "A_" + std::__cxx11::to_string(a) + "/" + "0" + epcho + "query.txt";
        std::filesystem::create_directories(outputFilename.substr(0, outputFilename.find_last_of('/')));
        // std::filesystem::create_directories(outputDir);
        std::ofstream os(outputFilename);
        if (!os.is_open()) {
            std::cerr << "Failed to create file: " << outputFilename << std::endl;
            return false;
        }
        for (const auto& ip : pos_keys_ipaddr_set) {
            os << ip << std::endl;
        }
        for (const auto& ip : query2) {
            os << ip << std::endl;
        }

        os.close();
        std::cout << "Query keys written to file: " << outputFilename << std::endl;
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadCAIDA_adp_A_S(bool using_cost_, std::string epcho, double load_factor, size_t a) {
    std::cout << "caida reading..."  << std::endl;
    int lf = load_factor * 100;
    std::ifstream is(std::string(CAIDA_PATH) + std::__cxx11::to_string(lf) + "/" 
        + "A_" + std::__cxx11::to_string(a) + "/" + "0" + epcho + "query.txt", std::ios::binary);
    size_t pos_num = (1<<12) * load_factor;
    if(is){
        size_t ip;
        // size_t totalIpPair = 0;
        // size_t posIpPair = 0;
        // std::unordered_map<uint32_t, size_t> ipset = std::unordered_map<uint32_t, size_t>();
        while (is >> ip) {  
            try {
                // if (ipset.find(ip) == ipset.end()) {
                //     ipset.insert(std::make_pair(ip, 1));
                // } else {
                //     ipset[ip]++;                
                // }
                if (pos_keys_ip_.size() < pos_num) {
                    pos_keys_ip_.push_back(ip);
                } else {
                    query_keys_ip_.push_back(ip);
                }
            } catch (const std::exception& e) {
            }
        }
        pos_keys_ipaddr_set = std::set<uint32_t>(pos_keys_ip_.begin(), pos_keys_ip_.end());
        // auto it = ipset.begin();
        // // std::advance(it, ipset.size() / 2);
        // std::advance(it, pos_num);
        // for (auto it2 = ipset.begin(); it2 != it; ++it2) {
        //     Slice * key = new Slice();
        //     key->ipaddr = it2->first;
        //     key->str = std::to_string(it2->first);
        //     pos_keys_.push_back(key);
        // }
        // for (auto it2 = it; it2 != ipset.end(); ++it2) {
        //     Slice * key = new Slice();
        //     key->ipaddr = it2->first;
        //     key->str = std::to_string(it2->first);
        //     key->cost = using_cost_ ? it2->second : 1;
        //     neg_keys_.push_back(key);
        // }
        // std::cout << "ipv4: " << ipset.size() << std::endl;
        std::cout << "pos_keys_: " << pos_keys_ip_.size() << std::endl;
        // std::cout << "totalIpPair: " << totalIpPair << std::endl;
        std::cout << "query_keys: " << query_keys_ip_.size() << std::endl;
        // std::cout << "posIpPair: " << posIpPair << std::endl;
        // std::cout << "neg_keys_: " << neg_keys_.size() << std::endl;
        // for (const auto& slice : pos_keys_) {
        //     pos_keys_ipaddr_set.insert(slice->ipaddr);
        // }
        // query_keys.erase(
        //     std::remove_if(query_keys.begin(), query_keys.end(),
        //                 [&](Slice* key) {
        //                     return pos_keys_ipaddr_set.find(key->ipaddr) != pos_keys_ipaddr_set.end();
        //                 }),
        // query_keys.end());
        // ipset.clear();    
        // std::unordered_map<uint32_t, size_t>().swap(ipset);

        // std::string outputFilename = std::string(CAIDA_PATH) + "0" + epcho + "query.txt";
        // std::ofstream os(outputFilename);
        // if (!os.is_open()) {
        //     std::cerr << "Failed to create file: " << outputFilename << std::endl;
        //     return false;
        // }

        // for (const auto& key : query_keys) {
        //     os << key->ipaddr << std::endl;
        // }

        // os.close();
        // std::cout << "Query keys written to file: " << outputFilename << std::endl;
        is.close();
        return true;
    } else {
        std::cout << "Failed to open file: " << std::endl;
    }
    is.close();
    return false;
}


bool dataloader::loadCAIDA(bool using_cost_, std::string epcho){
    std::cout << "caida reading..."  << std::endl;
    std::ifstream is(std::string(CAIDA_PATH) + "0" + epcho + ".txt", std::ios::binary);
    if(is){
        std::string sor, des;
        std::unordered_map<uint32_t, size_t> ipset = std::unordered_map<uint32_t, size_t>();
        while (is >> sor >> des) {  
            try {
                uint32_t ip = ipv4ToUint32(sor);
                if (ipset.find(ip) == ipset.end()) {
                    ipset.insert(std::make_pair(ip, 1));
                } else {
                    ipset[ip]++;                
                }
            } catch (const std::exception& e) {
            }
        }
        auto it = ipset.begin();
        std::advance(it, ipset.size() / 2);
        for (auto it2 = ipset.begin(); it2 != it; ++it2) {
            Slice * key = new Slice();
            key->ipaddr = it2->first;
            key->str = std::to_string(it2->first);
            pos_keys_.push_back(key);
        }
        for (auto it2 = it; it2 != ipset.end(); ++it2) {
            Slice * key = new Slice();
            key->ipaddr = it2->first;
            key->str = std::to_string(it2->first);
            key->cost = using_cost_ ? it2->second : 1;
            neg_keys_.push_back(key);
        }
        std::cout << "ipv4: " << ipset.size() << std::endl;
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::load(std::string data_name_, bool using_cost_){
    if(data_name_ == "shalla_0") return loadShalla_sk(using_cost_, "0");
    else if(data_name_ == "shalla_1") return loadShalla_sk(using_cost_, "1");
    else if(data_name_ == "shalla_2") return loadShalla_sk(using_cost_, "2");
    else if(data_name_ == "shalla_3") return loadShalla_sk(using_cost_, "3");
    else if(data_name_ == "shalla_4") return loadShalla_sk(using_cost_, "4");
    else if(data_name_ == "ycsb") return loadYCSB(using_cost_, "");
    else if(data_name_ == "hosefire_1") return loadHosefire_1(using_cost_, "");
    else if(data_name_ == "hosefire_2") return loadHosefire_2(using_cost_, "");
    else if(data_name_ == "hosefire_3") return loadHosefire_3(using_cost_, "");
    else if(data_name_ == "hosefire_4") return loadHosefire_4(using_cost_, "");
    else if(data_name_ == "hosefire_5") return loadHosefire_5(using_cost_, "");
    

    else return false;
}

bool dataloader::load(std::string data_name_, bool using_cost_, std::string epcho){
    if(data_name_ == "shalla") return loadShalla(using_cost_, epcho);
    else if(data_name_ == "ycsb") return loadYCSB(using_cost_, epcho);
    else return false;
}
bool dataloader::loadShalla_sk(bool using_cost_, std::string epcho){
    std::cout << "shalla reading..."  << std::endl;
    std::ifstream is(SHALLA_PATH_SK+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        double max_cost = 0;
        size_t total_cost = 0;
        // size_t neg_cnt = 0;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                max_cost = max_cost > key->cost ? max_cost : key->cost;
                neg_keys_.push_back(key); 
                total_cost += key->cost;
            }
        }
        std::cout << "max_cost: " << max_cost << std::endl;
        std::cout << "total_cost: " << total_cost << std::endl;
        std::cout << "neg_keys_: " << neg_keys_.size() << std::endl;
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadShalla_semi(std::string epcho, double sknewness){
    std::cout << "shalla reading..."  << std::endl;
    std::string path;
     if (sknewness == 1.0) {
        path = "../data/dataset_shalla/shalla_1_0/shalla_cost_1.0-";
    } else if (sknewness == 1.5) {
        path = "../data/dataset_shalla/shalla_1_5/shalla_cost_1.5-";
    } else if (sknewness == 2.0) {
        path = "../data/dataset_shalla/shalla_2_0/shalla_cost_2.0-";
    } else if (sknewness == 2.5) {
        path = "../data/dataset_shalla/shalla_2_5/shalla_cost_2.5-";
    }
    // SHALLA_PATH_SK  "../data/dataset_shalla/shalla_1_5/shalla_cost_1.5-"
    std::ifstream is(path+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        double max_cost = 0;
        size_t total_cost = 0;
        // size_t neg_cnt = 0;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = cost;
                max_cost = max_cost > key->cost ? max_cost : key->cost;
                neg_keys_.push_back(key); 
                total_cost += key->cost;
            }
        }
        std::cout << "max_cost: " << max_cost << std::endl;
        std::cout << "total_cost: " << total_cost << std::endl;
        std::cout << "neg_keys_: " << neg_keys_.size() << std::endl;
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadShalla_vulNeg(std::string name, bool using_cost_, std::string epcho){
    std::cout << "shalla reading..."  << std::endl;
    std::ifstream is(name+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadNegtive(){
    std::cout << "negative item reading..."  << std::endl;
    std::ifstream is("../data/addNegtive.txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = cost;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadShalla(bool using_cost_, std::string epcho){
    std::cout << "shalla reading..."  << std::endl;
    std::ifstream is(SHALLA_PATH+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}
bool dataloader::loadHosefire(bool using_cost_, std::string epcho){
    std::cout << "hosefire reading..."  << std::endl;
    std::ifstream is(HOSEFIRE_PATH+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadHosefire_1(bool using_cost_, std::string epcho){
    std::cout << "hosefire1 reading..."  << std::endl;
    std::ifstream is(HOSEFIRE_PATH_1+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadHosefire_2(bool using_cost_, std::string epcho){
    std::cout << "hosefire2 reading..."  << std::endl;
    std::ifstream is(HOSEFIRE_PATH_2+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadHosefire_3(bool using_cost_, std::string epcho){
    std::cout << "hosefire3 reading..."  << std::endl;
    std::ifstream is(HOSEFIRE_PATH_3+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadHosefire_4(bool using_cost_, std::string epcho){
    std::cout << "hosefire4 reading..."  << std::endl;
    std::ifstream is(HOSEFIRE_PATH_4+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadHosefire_5(bool using_cost_, std::string epcho){
    std::cout << "hosefire5 reading..."  << std::endl;
    std::ifstream is(HOSEFIRE_PATH_5+epcho+".txt", std::ios::binary);
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "1") pos_keys_.push_back(key);
            else if(optype == "0"){    
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}



bool dataloader::loadYCSB(bool using_cost_, std::string epcho){
    std::cout << "ycsb reading..."  << std::endl;
    std::ifstream is(YCSB_PATH+epcho+".txt");
    if(is){
        std::string optype, keystr;
        double cost;
        while (is >> optype >> keystr >> cost) 
        { 
            Slice * key = new Slice();    
            key->str = keystr;   
            if(optype == "FILTERKEY" || optype == "1") pos_keys_.push_back(key);
            else if(optype == "OTHERKEY" || optype == "0"){          
                key->cost = using_cost_ ? cost : 1;
                neg_keys_.push_back(key); 
            }
        }
        is.close();
        return true;
    }
    is.close();
    return false;
}

bool dataloader::loadRandomKey(int positives_number, int negatives_number, bool using_cost_){
    KeyBuilder kb;
    for(int i=0; i<positives_number; i++){
        Slice *key = new Slice();
        pos_keys_.push_back(key);
    }
    for(int i=0; i<negatives_number; i++){
        Slice *key = new Slice();
        neg_keys_.push_back(key);
    }
    if(!kb.ReadKeys(pos_keys_, 0)) return false; 
    if(!kb.ReadKeys(pos_keys_, 0)) return false; 
    if(using_cost_)
        switch (RANDOM_COST_TYPE)
        {
        case uniform:
            kb.GenKeysUniformCosts(neg_keys_, 5);
            break;
        case hotcost:
            kb.GenKeysHotCosts(neg_keys_, 0.01, 100, 1);
            break;
        case normal:
            kb.GenKeysNormalCosts(neg_keys_, 20, 10);
            break;
        case zipf:
            kb.GenKeysZipfCosts(neg_keys_, 1.25, 1.0);
            break;
        default:
            break;
        }
    return true;
}
        // std::vector<Slice *> pos_keys_; // positive keys
        // std::vector<Slice *> neg_keys_; // negative keys
        // std::vector<Slice *> query_keys; // query keys
        // std::set<uint32_t> pos_keys_ipaddr_set;
dataloader::~dataloader() {
    for (Slice *key : pos_keys_) {
        delete key;
    }
    for (Slice *key : neg_keys_) {
        delete key;
    }
    for (Slice *key : query_keys) {
        delete key;
    }
    std::cout << "dataloader delete..." << std::endl;
    pos_keys_ipaddr_set.clear();
    std::set<uint32_t>().swap(pos_keys_ipaddr_set); // 释放内存

    pos_keys_.clear();
    std::vector<Slice*>().swap(pos_keys_); // 释放内存

    neg_keys_.clear();
    std::vector<Slice*>().swap(neg_keys_); // 释放内存

    query_keys.clear();
    std::vector<Slice*>().swap(query_keys); // 释放内存
    malloc_trim(0);
}

KeyBuilder::KeyBuilder(){
    std::fstream ifs(RANDOM_KEYSTR_PATH);
    if(ifs.is_open()){
        std::string s;
        while(std::getline(ifs,s))
            key_strs.push_back(s);
    }else{
        std::cout << "Keystr file not exists, generate again..." << std::endl;
        GenKeyStrAndToFile();
    }
    ifs.close();
}
std::string KeyBuilder::GetKeyStr(){
    int k=rand()%10+1;
    char arr[10];
    for(int i=1;i<=k;i++){
        int x,s;                         
        s=rand()%2;                     
        if(s==1) x=rand()%('Z'-'A'+1)+'A';        
        else x=rand()%('z'-'a'+1)+'a';      
        arr[i-1] = x;                  
    }
    return std::string(arr,k);
}

bool KeyBuilder::ReadKeys(std::vector<Slice *> &v, int start_position_){
    int size_ = v.size();
    if(start_position_ + size_ >= key_strs.size()) return false;
    for(int j=0; j<size_; j++)
        v[j]->str = key_strs[start_position_+j];
    return true;
}

void KeyBuilder::GenKeyStrAndToFile(){
    int gen_key_size_ = 200000;
    int i = key_strs.size();
    std::ofstream ofs(RANDOM_KEYSTR_PATH);
    while(i < gen_key_size_){
        if(0 == i%10000) std::cout << i << "keys have been created..." << std::endl;
        std::string str = GetKeyStr();
        if(std::find(key_strs.begin(), key_strs.end(), str) == key_strs.end()){
            key_strs.push_back(str);
            ofs << str << std::endl;
            i++;
        }
    }
    ofs.close();
}

void KeyBuilder::GenKeysUniformCosts(std::vector<Slice *> &keys, int interval){
    for(int i=0; i<keys.size(); i++)
        keys[i]->cost = 1+(i+1)*interval;
}
void KeyBuilder::GenKeysHotCosts(std::vector<Slice *> &keys, double hotNumberpro, int hotcost, int coldcost){
    int hotNumSize = hotNumberpro * keys.size();
    for(int i=0; i<keys.size(); i++){
        if(i <= hotNumSize) 
            keys[i]->cost = hotcost;
        else 
            keys[i]->cost = coldcost;
    }
}
void KeyBuilder::GenKeysNormalCosts(std::vector<Slice *> &keys, int u, int d){
    for(int i=0; i<keys.size(); i++){
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);
        std::normal_distribution<double> distribution(u, d);
        keys[i]->cost = distribution(generator);
    }
}

void KeyBuilder::GenKeysZipfCosts(std::vector<Slice *> &keys, double a, double c){
    int r = 10000;
    double pf[10000];
    double sum = 0.0;
    for (int i = 0; i < r; i++)        
        sum += c/pow((double)(i+2), a);  
    for (int i = 0; i < r; i++){ 
        if (i == 0)
            pf[i] = c/pow((double)(i+2), a)/sum;
        else
            pf[i] = pf[i-1] + c/pow((double)(i+2), a)/sum;
    }
     for (int i = 0; i < keys.size(); i++){
        int index = 0;
        double data = (double)rand()/RAND_MAX;  
        while (data > pf[index])  
            index++;
        keys[i]->cost = index;
    }
}

#endif