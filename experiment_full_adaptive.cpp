#include <cassert>
#include <cstdlib>
#include <iostream>
#include <time.h>
#include <vector>
#include <queue>
#include <utility> // for std::pair
#include <functional> // for std::greater
#include <assert.h>
#include <math.h>
#include<string>

#include "util/dataloader.h"
#include "nonlearnedfilter/countingBloom.h"
#include "nonlearnedfilter/cuckoofilter/src/perfectCF.h"
#include "nonlearnedfilter/cuckoofilter/src/cuckoofilter.h"
#include "nonlearnedfilter/cuckoofilter/src/RCF.h"
#include "nonlearnedfilter/cuckoofilter/src/adpPerfectCF.h"
// #include "../telescoping-filter/src/taf.h"
#include "nonlearnedfilter/telescoping-filter/src/taf.h"
// #include "../cuckoofilter/src/perfectCF.h"
// #include "../cuckoofilter/src/cuckoofilter.h"
// #include "../cuckoofilter/src/RCF.h"
// #include "../cuckoofilter/src/adpPerfectCF.h"
// #include "RF/adaPerfectCF.h"
#define SHALLA_NAME "shalla_"
#define YCSBT_NAME "ycsbt_"
#define HOSEFIRE_NAME "hosefire_"
const size_t datasetNum = 5; 

using namespace std;
using namespace std::chrono;
using cuckoofilter::PerfectCF;
using namespace std;

void TestPerfectCFCAIDA(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio);
void TestPerfectCFCAIDA_adp(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio);
void TestAdpPerfectCFCAIDA_adp(dataloader * dlSet, double * bits_per_key_sets);
void TestPerfectCFCAIDA_adp_potential(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio);

void TestPerfectCFCAIDA_adp_cf(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio);
void TestAdpPerfectCFCAIDA_adp_cf(dataloader * dlSet, double * bits_per_key_sets);
void TestPerfectCFCAIDA_adp_cf_AND_RCF_CAIDA(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio);
void TestPerfectCFCAIDA_adp_cf_potential(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio);
void TestTelescopingFilter(dataloader * dlSet, double * bits_per_key_sets);
void TestRCF_CAIDA(dataloader * dlSet, double * bits_per_key_sets);
void TestRCF(dataloader * dlSet);

bool compare_func(Slice * key1, Slice * key2);
// #define  dataSet 2

vector<vector<double>> wFPR = vector<vector<double>>(7);
vector<vector<double>> insert_latency = vector<vector<double>>(7);
vector<vector<double>> query_latency = vector<vector<double>>(7);
vector<vector<double>> delete_latency = vector<vector<double>>(7);
// vector<size_t> A_S = {1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50};
// vector<double> skewness = {0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5};
// vector<vector<double>> skewness = {{0.25, 1.0}, {}, {0.25, 0.99}};
vector<vector<double>> skewness = {{1.0}, {}, {}};
// 0->CAIDA; 2->YCSB 
vector<int> dataSets = {0};

void display() {
    vector<string> filters = {
        "RCF",
        "TAF",
        "RF(CF)",
        "RF(BF)"
    };

    std::cout << "REM_SIZE (FINGERPRINT SIZE IN TAF): " << REM_SIZE << endl;
    // display wFPR
    std::cout << "--------------wFPR testing--------------" << endl;
    std::cout << "Load Factor: 0.55,  0.65,  0.75,  0.85,  0.95" << endl;
    for (size_t i = 0; i < wFPR.size(); i++)
    {
        if (wFPR[i].empty()) continue;
        std::cout << filters[i] << ": ";
        for (size_t j = 0; j < wFPR[i].size(); j++)
        {
            std::cout << wFPR[i][j] << " ";
        }
        if (wFPR[i].size()) std::cout << endl;
        wFPR[i].clear();
    }
    //display load_factor
    // for (size_t i = 0; i < load_factor.size(); i++)
    // {
    //     std::cout << load_factor[i] << " ";
    // }
    //display insert_latency
    std::cout << "--------------insert_latency testing--------------" << endl;
    for (size_t i = 0; i < insert_latency.size(); i++)
    {
        if (insert_latency[i].empty()) continue;
        std::cout << filters[i] << ": ";
        for (size_t j = 0; j < insert_latency[i].size(); j++)
        {
            std::cout << insert_latency[i][j] << " ";
        }
        if (insert_latency[i].size()) std::cout << endl;
        insert_latency[i].clear();
    }
    //display query_latency
    std::cout << "--------------query_latency testing--------------" << endl;
    for (size_t i = 0; i < query_latency.size(); i++)
    {
        if (query_latency[i].empty()) continue;
        std::cout << filters[i] << ": ";
        for (size_t j = 0; j < query_latency[i].size(); j++)
        {
            std::cout << query_latency[i][j] << " ";
        }
        if (query_latency[i].size()) std::cout << endl;
        query_latency[i].clear();
    }
    //display delete_latency
    // std::cout << "delete_latency" << endl;
    // for (size_t i = 0; i < 5; i++)
    // {
    //     for (size_t j = 0; j < delete_latency[i].size(); j++)
    //     {
    //         std::cout << delete_latency[i][j] << " ";
    //     }
    //     std::cout << endl;
    // }
}

int main(){
    for (int dataSet : dataSets) {
        std::cout << "----------------------------" << endl;
        std::cout << "dataSet: " << dataSet << endl;
        for (double s : skewness[dataSet]) {
            std::cout.setstate(std::ios_base::failbit);
            for (double lf = 0.55; lf < 0.96; lf+=0.1) {
                dataloader * dlSet = new dataloader[datasetNum];
                double *bits_per_key_sets = new double[datasetNum];
                for (size_t i = 0; i < datasetNum; i++)
                {
                    if (dataSet == 0) {
                        if (!dlSet[i].loadZipf_adp(true, to_string(i), lf, s))
                        {
                            std::cout << "error!\n";
                        }

                    }
                    else if (dataSet == 2) {
                        // if (!dlSet[i].testYCSB_adp(to_string(i), lf, s))
                        // {
                        //     std::cout << "error!\n";
                        // }
                        if (!dlSet[i].loadYCSB_adp(to_string(i), lf, s))
                        {
                            std::cout << "error!\n";
                        }
                    }
                    // if (!dlSet[i].loadShalla_sk(true, to_string(i)))
                    // {
                    //     std::cout << "error!\n";
                    // }
                    bits_per_key_sets[i] = REM_SIZE + 4;
                }    

                for (int round = 0; round < 1 ; round++){
                std::cout << "--------------TestReinforcementCuckooFilter testing--------------" << endl;
                TestRCF_CAIDA(dlSet, bits_per_key_sets);
                std::cout << "--------------TestTelescopingFilter testing--------------" << endl;
                TestTelescopingFilter(dlSet,bits_per_key_sets);
                // std::cout << "--------------TestPerfectCFCAIDA_adp_cf testing----------------" << endl;
                // TestPerfectCFCAIDA_adp_cf(dlSet, bits_per_key_sets, 0.05);
                // std::cout << "--------------TestPerfectCFCAIDA_adp testing----------------" << endl;
                // TestPerfectCFCAIDA_adp(dlSet, bits_per_key_sets, 0.05);
                // std::cout << "--------------TestPerfectCFCAIDA_adp_cf_potential testing----------------" << endl;
                // TestPerfectCFCAIDA_adp_cf_potential(dlSet, bits_per_key_sets, 0.05);
                std::cout << "--------------TestAdaPerfectCF_CuckooFilter testing--------------" << endl;
                TestAdpPerfectCFCAIDA_adp_cf(dlSet, bits_per_key_sets);
                std::cout << "--------------TestAdaPerfectCF_BloomFilter testing--------------" << endl;
                TestAdpPerfectCFCAIDA_adp(dlSet, bits_per_key_sets);
                }
                delete[] dlSet;
                delete[] bits_per_key_sets;                
        }
            std::cout.clear();
            std::cout << "--------------display--------------" << endl;
            std::cout << "skewness: "  << s <<endl;
            display();
        }
    }
}

void TestRCF(dataloader * dlSet) {
    for (size_t i = 0; i < 1; i++)
    {
        // bits_per_key_sets[i] *= 0.9;
        size_t cnt = 0;
        std::set<size_t>* s = new std::set<size_t>();
        // cuckoofilter::RCF<size_t,64>* RCF = new cuckoofilter::RCF<size_t,64>(dlSet[i].pos_keys_.size(), 1, 8);
        cuckoofilter::RCF<size_t,32>* RCF = new cuckoofilter::RCF<size_t,32>(dlSet[i].pos_keys_ipaddr_set.size(), 6, 8);
        for (Slice *v : dlSet[i].pos_keys_ )
        {    
            // size_t temp=XXH3_128bits(&v->str[0], v->str.size()).low64;
            // RCF->Add(temp);
            // s->insert(temp);
            RCF->Add(v->ipaddr);
            s->insert(v->ipaddr);
        }

        for (Slice *v : dlSet[i].query_keys )
        {    
            size_t temp=XXH3_128bits(&v->str[0], v->str.size()).low64;
            // if (RCF->Contain(temp) == cuckoofilter::Ok)
            if (RCF->Contain(v->ipaddr) == cuckoofilter::Ok)
            {
                // if (s->find(temp) == s->end())
                if (s->find(v->ipaddr) == s->end())
                {
                    RCF->adp(v->ipaddr);
                    // cnt++;//7895
                    // if (RCF->adp(v->ipaddr)){
                        // cnt++;
                        for (Slice *v2 : dlSet[i].pos_keys_ )
                        {
                            if (RCF->Contain(v2->ipaddr) != cuckoofilter::Ok)
                            {
                                RCF->testGenerateIndexTagHash(v->ipaddr);
                                RCF->testshow(v->ipaddr);
                                std::cout<<endl;
                                RCF->testGenerateIndexTagHash(v2->ipaddr);
                                RCF->testshow(v2->ipaddr);
                                 std::cout<<endl;
                                 RCF->testContain(v2->ipaddr);
                                // std::cout<<"\n";
                                cnt++;
                                break;
                            }
                        }  
                        if (cnt >= 3) break;
                    // }
                }
            }
        }
        std::cout << "adp cnt: " << cnt << "\n";
        cnt = 0;
        for (Slice *v : dlSet[i].pos_keys_ )
        {    
            // size_t temp=XXH3_128bits(&v->str[0], v->str.size()).low64;
            // RCF->Add(temp);
            // s->insert(temp);
            // RCF->Add(v->ipaddr);
            // s->insert(v->ipaddr);
            
            if (RCF->Contain(v->ipaddr) != cuckoofilter::Ok)
            {
                // if (s->find(temp) == s->end())
                // std::cout<<"query hou error!\n";
                cnt++;
            }
            
        }
        std::cout<<"query hou error    " <<cnt<<  "\n";
        for (Slice *v : dlSet[i].pos_keys_ )
        {
            // size_t temp=XXH3_128bits(&v->str[0], v->str.size()).low64;
            // RCF->Delete(temp);
            if (!RCF->Delete(v->ipaddr)) {
                // std::cout<<"Delete error!\n";
                // RCF->testContain(v->ipaddr);
            }
        }
        for (Slice *v : dlSet[i].neg_keys_ )
        {
            // size_t temp=XXH3_128bits(&v->str[0], v->str.size()).low64;
            // if (RCF->Contain(temp) == cuckoofilter::Ok)
            if (RCF->Contain(v->ipaddr) == cuckoofilter::Ok)
            {
                std::cout<<"error!\n";
                // RCF->testContain(v->ipaddr);
                // RCF->Delete(v->ipaddr);
                // RCF->testContain(v->ipaddr);
            }
        }
    //     perfectCF[i] = new cuckoofilter::PerfectCF<size_t>(dlSet[i].pos_keys_.size()*vulNegKeyRatio);
    //     double alpha = 1 - perfectCF[i]->SizeInBytes() * 8.0 / (bits_per_key_sets[0] * dlSet[0].pos_keys_.size());
    //     std::cout<<"alpha: "<<alpha<<endl;
    //     nfuncs = bits_per_key_sets[i] * alpha/ 4 * log(2);
    //     cbfArr[i] = new countingBloom::CountingBloom(bits_per_key_sets[0] *alpha* dlSet[0].pos_keys_.size() /bitPerCounter, bitPerCounter, nfuncs);
    }
}
void TestRCF_CAIDA(dataloader * dlSet, double * bits_per_key_sets) 
{

    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;
    cuckoofilter::RCF<size_t,32> * RCF[datasetNum];
    // set<size_t> * pos_keys_ipaddr_set[5];
    for (size_t i = 0; i < datasetNum; i++)
    {
        // bits_per_key_sets[i] *= 0.9;
        // if (REM_SIZE == 4) {
        //     RCF[i] = new cuckoofilter::RCF<size_t, 32>(dlSet[i].pos_keys_.size(), bits_per_key_sets[i]-2, 11);
        //     // perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(8192);//REM_SIZE 4
        // } else if (REM_SIZE == 8) {
        //     RCF[i] = new cuckoofilter::RCF<size_t, 32>(dlSet[i].pos_keys_.size(), bits_per_key_sets[i]-2, 8);
        //     // perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(1024);//REM_SIZE 8
        // } else if (REM_SIZE == 16) {
        //     RCF[i] = new cuckoofilter::RCF<size_t, 32>(dlSet[i].pos_keys_.size(), bits_per_key_sets[i]-2, 5);
        //     // RCF[i] = new cuckoofilter::RCF<size_t, 32>(dlSet[i].pos_keys_.size(), bits_per_key_sets[i]-2, 5);
        //     // perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(128);//REM_SIZE 16
        // } 

        size_t totalBits = dlSet[i].pos_keys_ipaddr_set.size() * bits_per_key_sets[i];
        size_t h1 = std::ceil(std::log2(dlSet[i].pos_keys_ipaddr_set.size())) - 2;
        size_t h2 = bits_per_key_sets[i]-2 <= 32 - 1 - h1 ? bits_per_key_sets[i]-2 : 32 - 1 - h1;
        size_t restBits = totalBits > (1ULL << h1) * 4 * (h2+1) ? totalBits - (1ULL << h1) * 4 * (h2+1) : (1ULL << h1) * 4 * (h2+1) * 0.01;
        restBits = max(restBits, static_cast<size_t>((1ULL << h1) * 4 * (h2+1) * 0.01));
        int is = 1;
        while ((1ULL << (is+2)) * (32 - h2 + 1 - is) <= restBits) is++;
        is--;
        size_t size = 1ULL << (is + 2);
        std::cout<<"size: "<<size<<endl;
        RCF[i] = new cuckoofilter::RCF<size_t, 32>(dlSet[i].pos_keys_ipaddr_set.size(), bits_per_key_sets[i]-2, is);

    }
    for (size_t i = 0; i < datasetNum; i++)
    {
        dlSet[i].pos_keys_ipaddr_set.clear();
        // pos_keys_ipaddr_set[i] = new set<size_t>();
        auto t1 = steady_clock::now(); 
        for (auto& ip : dlSet[i].pos_keys_ip_)
        {    
            RCF[i]->Add(ip);
            dlSet[i].pos_keys_ipaddr_set.insert(ip);
            // if (RCF[i]->Add(ip) != cuckoofilter::Ok)
            // {
            //     std::cout<<"add error!\n";
            // }
            // pos_keys_ipaddr_set[i]->insert(v->ipaddr);
        }
        auto t2 = steady_clock::now();
        std::cout<<RCF[i]->Info()<<endl;
        // std::cout<<"set size: "<<s[i]->size()<<endl;
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_ip_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;
    }

    std::cout<<"Average Insertion Latency (Time(ns)/Key): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";

    double sumFNR = 0;
    double sumQueryTime = 0;
    double minQueryTime = MAXFLOAT;
    double maxQueryTime = 0;
    double sumWeightedFPR = 0;
    double sumFPR = 0;
    double max_weightedFpr = 0;
    double min_weightedFpr = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0;
        double wfpr_c = 0;
        size_t fpr_Num = 0;
        auto t1 = steady_clock::now();
        for (auto& ip : dlSet[i].pos_keys_ip_)
        {
            // if (RCF[i]->Contain(v->ipaddr) != cuckoofilter::Ok) {
            //     fnr_c++;
            // }
            if (RCF[i]->Contain(ip) == cuckoofilter::Ok) {
                assert(dlSet[i].pos_keys_ipaddr_set.find(ip) != dlSet[i].pos_keys_ipaddr_set.end());
            } else fnr_c++;
        }

        for (auto& ip: dlSet[i].query_keys_ip_)
        {
            if(RCF[i]->Contain(ip) == cuckoofilter::Ok)
            {
                wfpr_c ++;
                fpr_Num ++;
                // if (pos_keys_ipaddr_set[i]->find(v->ipaddr) == pos_keys_ipaddr_set[i]->end())
                if (dlSet[i].pos_keys_ipaddr_set.find(ip) == dlSet[i].pos_keys_ipaddr_set.end())
                {
                    RCF[i]->adp(ip);
                }
            }

        }
        auto t2 = steady_clock::now(); 
        std::cout<<"fpr_Num: "<<fpr_Num<<endl;
        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_ip_.size()+dlSet[i].query_keys_ip_.size());
        sumQueryTime += tempV;
        minQueryTime = tempV < minQueryTime ? tempV : minQueryTime;
        maxQueryTime = tempV > maxQueryTime ? tempV : maxQueryTime;
    
        sumFPR += (double)(fpr_Num*100) / (double)dlSet[i].query_keys_ip_.size();

        sumWeightedFPR += (double)(wfpr_c*100) / dlSet[i].query_keys_ip_.size();
        if((double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size() > max_weightedFpr){
            max_weightedFpr = (double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size();
        }
        if((double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size() < min_weightedFpr){
            min_weightedFpr = (double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size();
        }
        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_ip_.size();
    }
    
    std::cout<<"Average Query Latency (Time(ns)/Key): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"sumFPR: "<< sumFPR <<endl;
    std::cout<<"datasetNum: "   << datasetNum <<endl;
    std::cout<<"sumFNR: "<< sumFNR <<endl;
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;
    wFPR[0].push_back((double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2));
    insert_latency[0].push_back((sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2));
    query_latency[0].push_back((sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2));
    // load_factor.push_back(perfectCF[0]->LoadFactor());
    // double sumDeleteTime = 0;
    // double minDeleteTime = MAXFLOAT;
    // double maxDeleteTime = 0;
    // for (size_t i = 0; i < datasetNum; i++)
    // {

        
    //     auto t1 = steady_clock::now();
    //     for (auto& ip: dlSet[i].pos_keys_ip_ )
    //     {
    //         RCF[i]->Delete(ip);
    //         // if(RCF[i]->Contain(v->ipaddr) == cuckoofilter::Ok) std::cout<<"--------error----------\n";
    //     }
    //     auto t2 = steady_clock::now();       
    //     double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_ip_.size());  
    //     sumDeleteTime += tempV;
    //     minDeleteTime = tempV < minDeleteTime ? tempV : minDeleteTime;
    //     maxDeleteTime = tempV > maxDeleteTime ? tempV : maxDeleteTime;
    // }  
    // std::cout<<"Average Deletion Latency (Time(ns)/Key) : "<<(sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2)<<endl;
    // // delete_latency[0].push_back((sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2));
    // double sumFNR_afterDelete=0;
    // for (size_t i = 0; i < datasetNum; i++)
    // {
    //     size_t fnr_c = 0;
    //     auto t1 = steady_clock::now();
    //     for (auto& ip : dlSet[i].pos_keys_ip_ )
    //     {
    //         // RCF[i]->Delete(v->ipaddr);
    //         if(RCF[i] ->Contain(ip) == cuckoofilter::Ok){
    //             fnr_c++;
    //             // RCF[i]->Delete(ip);
    //             // // size_t i;
    //             // // size_t tag;
    //             // // size_t hs;
    //             // // RCF[i] ->TESTHASH(v->ipaddr, &i, &tag, &hs);
    //             // // std::cout<<i<<"   "<<tag<<"   "<<hs<<endl;
    //             // if(RCF[i]->Contain(v->ipaddr) == cuckoofilter::Ok) std::cout<<"error\n";
    //             // else std::cout<<"correct\n";
                
    //         } 
    //     }
    
    //     sumFNR_afterDelete += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_ip_.size();
    // }
    // std::cout<<"FNR(After delete): "<< (double)sumFNR_afterDelete / datasetNum << "%" <<endl; 
    for (size_t i = 0; i < datasetNum; i++)
    {
        delete RCF[i];
    }
}

void TestTelescopingFilter(dataloader * dlSet, double * bits_per_key_sets)
{
    TAF * TAFArr[datasetNum];
    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;
    for (size_t i = 0; i < datasetNum; i++)
    {
        // 0.0796877%
        TAFArr[i] = (TAF *)malloc(sizeof(TAF));
        // taf_init(TAFArr[i], (size_t)(dlSet[i].pos_keys_.size() * bits_per_key_sets[i]/(bits_per_key_sets[i]-1)), 32776517);
        taf_init(TAFArr[i], (size_t)(dlSet[i].pos_keys_ip_.size() ), 32776517);
    }

    for (size_t i = 0; i < datasetNum; i++)
    {
        auto t1 = steady_clock::now(); 
        for (auto& ip : dlSet[i].pos_keys_ip_ )
        {    
            taf_insert(TAFArr[i], ip);
        }
        

        auto t2 = steady_clock::now();
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_ip_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;
    }

    std::cout<<"Average Insertion Latency (Time(ns)/Key): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";

    double sumFNR = 0;
    double sumQueryTime = 0;
    double minQueryTime = MAXFLOAT;
    double maxQueryTime = 0;
    double sumWeightedFPR = 0;
    double max_weightedFpr = 0;
    double min_weightedFpr = MAXFLOAT;
    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0;
        size_t wfpr_c = 0;
        auto t1 = steady_clock::now();
        for (auto& ip : dlSet[i].pos_keys_ip_ )
        {
            if (!taf_lookup(TAFArr[i], ip))
            {
                fnr_c++;
            }

        }
        
        for (auto& ip : dlSet[i].query_keys_ip_)
        {
            if (taf_lookup(TAFArr[i], ip))
            {
                wfpr_c ++;
            }

        }
        auto t2 = steady_clock::now();   
        // print_taf_stats(TAFArr[i]);
        // print_taf_metadata(TAFArr[i]);
        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_ip_.size()+dlSet[i].query_keys_ip_.size());
        sumQueryTime += tempV;
        minQueryTime = tempV < minQueryTime ? tempV : minQueryTime;
        maxQueryTime = tempV > maxQueryTime ? tempV : maxQueryTime;
    


        sumWeightedFPR += (double)(wfpr_c*100) / dlSet[i].query_keys_ip_.size();
        if((double)(wfpr_c*100) / dlSet[i].query_keys_ip_.size() > max_weightedFpr){
            max_weightedFpr = (double)(wfpr_c*100) / dlSet[i].query_keys_ip_.size();
        }
        if((double)(wfpr_c*100) / dlSet[i].query_keys_ip_.size() < min_weightedFpr){
            min_weightedFpr = (double)(wfpr_c*100) / dlSet[i].query_keys_ip_.size();
        }
        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_ip_.size();
    }
    
    std::cout<<"Average Query Latency (Time(ns)/Key): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"datasetNum: "   << datasetNum <<endl;
    // std::cout<<"sumFNR: "<< sumFNR <<endl;
    // std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;
    wFPR[1].push_back((double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2));
    insert_latency[1].push_back((sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2));
    query_latency[1].push_back((sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2));   
    for (size_t i = 0; i < datasetNum; i++)
    {
        taf_destroy(TAFArr[i]);
        // free(TAFArr[i]);
    }
}
// 定义一个比较函数，用于根据 cnt 的大小进行排序
struct CompareCnt {
    bool operator()(const std::pair<size_t, size_t>& a, const std::pair<size_t, size_t>& b) {
        return a.second < b.second; // 最大堆：cnt 较大的元素优先
    }
};
void TestPerfectCFCAIDA_adp_cf(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio) 
{
    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;
    cuckoofilter::PerfectCF<size_t, 32> * perfectCF[datasetNum];
    double sumFNR = 0;
    double sumQueryTime = 0;
    double minQueryTime = MAXFLOAT;
    double maxQueryTime = 0;
    double sumWeightedFPR = 0;
    double sumFPR = 0;
    double max_weightedFpr = 0;
    double min_weightedFpr = MAXFLOAT;
    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t totalBits = dlSet[i].pos_keys_ip_.size() * bits_per_key_sets[i];
        cuckoofilter::CuckooFilter<size_t, REM_SIZE+4-1>* cf = new cuckoofilter::CuckooFilter<size_t, REM_SIZE+4-1>(dlSet[i].pos_keys_ip_.size());
        size_t restBits = totalBits > cf->SizeInBytes() * 8 ? totalBits - cf->SizeInBytes() * 8 : cf->SizeInBytes() * 8 * 0.01;
        restBits = max(restBits, static_cast<size_t>(cf->SizeInBytes() * 8 * 0.01));
        std::cout<<"restBits: "<<restBits<<endl;
        int is = 1;
        while ((1ULL << (is+2)) * (33-is) <= restBits) is++;
        is--;
        size_t size = 1ULL << (is + 2);
        std::cout<<"size: "<<size<<endl;
        perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(size);
        // if (REM_SIZE == 4) {
        //     perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(8192);//REM_SIZE 4
        // } else if (REM_SIZE == 8) {
        //     perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(1024);//REM_SIZE 8
        // } else if (REM_SIZE == 16) {
        //     perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(128);//REM_SIZE 16
        // } 
        // double alpha = 1 - perfectCF[i]->SizeInBytes() * 8.0 / (bits_per_key_sets[0] * dlSet[0].pos_keys_ip_.size());
        // std::cout<<"alpha: "<<alpha<<endl;
        // nfuncs = bits_per_key_sets[i] * alpha * log(2);
        
        
        // cbfArr[i] = new countingBloom::CountingBloom(bits_per_key_sets[0] *alpha* dlSet[0].pos_keys_.size(), bitPerCounter, nfuncs);
        dlSet[i].pos_keys_ipaddr_set.clear();
        auto t1 = steady_clock::now(); 
        for (auto& ip: dlSet[i].pos_keys_ip_ )
        {    
            cf->Add(ip);
            dlSet[i].pos_keys_ipaddr_set.insert(ip);
            // cfArr[i]->counting_bloom_add( *v);
        }
        std::cout<< "pos_keys_ipaddr_set size: " << dlSet[i].pos_keys_ipaddr_set.size() << endl;
        // 创建一个最大堆，元素是 <size_t ipaddr, size_t cnt>
        // std::priority_queue<std::pair<size_t, size_t>, std::vector<std::pair<size_t, size_t>>, CompareCnt> maxHeap;
        // size_t isOneNum = dlSet[i].pos_keys_.size() * vulNegKeyRatio;
        // for (int j=0; j<isOneNum; j++)
        // {
        //     Slice *v = dlSet[i].neg_keys_[j];
        //     // s[i]->insert(XXH3_128bits(&v->str[0], v->str.size()).low64);
        //     perfectCF[i]->Add32(v->ipaddr);
        // }
        // std::cout<<perfectCF[i]->Info()<<endl;
        // std::cout<<"set size: "<<s[i]->size()<<endl;

        auto t2 = steady_clock::now();
        std::cout<<cf->Info()<<endl;
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_ip_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;

        size_t fnr_c = 0;
        double wfpr_c = 0;
        size_t fpr_Num = 0;
        t1 = steady_clock::now();
        for (auto& ip: dlSet[i].pos_keys_ip_ )
        {
            // if (cf->Contain(ip) != cuckoofilter::Ok || perfectCF[i]->Contain32(ip)== cuckoofilter::Ok)
            // {
            //     fnr_c++;
            // }
            if(cf->Contain(ip) == cuckoofilter::Ok && perfectCF[i]->Contain32(ip) != cuckoofilter::Ok) {
                assert(dlSet[i].pos_keys_ipaddr_set.find(ip) != dlSet[i].pos_keys_ipaddr_set.end());
            } else fnr_c++;

        }
        size_t cnt = 0;
        size_t cnt2 = 0;
        for (auto& ip: dlSet[i].query_keys_ip_)
        {
            if(cf->Contain(ip) == cuckoofilter::Ok && perfectCF[i]->Contain32(ip) != cuckoofilter::Ok)
            {
                wfpr_c ++;
                fpr_Num ++;
                if (dlSet[i].pos_keys_ipaddr_set.find(ip) == dlSet[i].pos_keys_ipaddr_set.end())
                {
                    perfectCF[i]->Add32(ip);
                    // if (perfectCF[i]->LoadFactor() < 0.95) {
                    //     perfectCF[i]->Add32(ip);
                    //     cnt2++;
                    // } else cnt++;
                }
            }

        }
        t2 = steady_clock::now();   
        std::cout<<"cnt: "<<cnt<<endl;
        std::cout<<"cnt2: "<<cnt2<<endl;
        std::cout<<"fpr_Num: "<<fpr_Num<<endl;
        std::cout<<"load factor: "<<perfectCF[i]->LoadFactor()<<endl;
        
        tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_ip_.size()+dlSet[i].query_keys_ip_.size());
        sumQueryTime += tempV;
        minQueryTime = tempV < minQueryTime ? tempV : minQueryTime;
        maxQueryTime = tempV > maxQueryTime ? tempV : maxQueryTime;
    
        sumFPR += (double)(fpr_Num*100) / (double)dlSet[i].query_keys_ip_.size();

        sumWeightedFPR += (double)(wfpr_c*100) / dlSet[i].query_keys_ip_.size();
        if((double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size() > max_weightedFpr){
            max_weightedFpr = (double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size();
        }
        if((double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size() < min_weightedFpr){
            min_weightedFpr = (double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size();
        }
        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_ip_.size();
        delete cf;
    }
    std::cout<<"Average Insertion Latency (Time(ns)/Key): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";
    std::cout<<"Average Query Latency (Time(ns)/Key): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"sumFPR: "<< sumFPR <<endl;
    std::cout<<"datasetNum: "   << datasetNum <<endl;
    std::cout<<"sumFNR: "<< sumFNR <<endl;
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;
    wFPR[2].push_back((double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2));
    insert_latency[2].push_back((sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2));
    query_latency[2].push_back((sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2));
    for (size_t i = 0; i < datasetNum; i++)
    {
        delete perfectCF[i];
    }
}

void TestAdpPerfectCFCAIDA_adp_cf(dataloader * dlSet, double * bits_per_key_sets) 
{
    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;
    // cuckoofilter::PerfectCF<size_t, 32> * perfectCF[datasetNum];
    const size_t cntSize = 2;
    cuckoofilter::AdpPerfectCF<size_t, 32, cntSize> * adpPerfectCF[datasetNum];
    // cuckoofilter::AdaPerfectCF<size_t, 32, cntSize> * adaPerfectCF[datasetNum];
    double sumFNR = 0;
    double sumQueryTime = 0;
    double minQueryTime = MAXFLOAT;
    double maxQueryTime = 0;
    double sumWeightedFPR = 0;
    double sumFPR = 0;
    double max_weightedFpr = 0;
    double min_weightedFpr = MAXFLOAT;
    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t totalBits = dlSet[i].pos_keys_ip_.size() * bits_per_key_sets[i];
        cuckoofilter::CuckooFilter<size_t, REM_SIZE+4-1>* cf = new cuckoofilter::CuckooFilter<size_t, REM_SIZE+4-1>(dlSet[i].pos_keys_ip_.size());
        size_t restBits = totalBits > cf->SizeInBytes() * 8 ? totalBits - cf->SizeInBytes() * 8 : cf->SizeInBytes() * 8 * 0.01;
        restBits = max(restBits, static_cast<size_t>(cf->SizeInBytes() * 8 * 0.01));
        std::cout<<"restBits: "<<restBits<<endl;
        int is = 1;
        while ((1ULL << (is+2)) * (33 + cntSize -is) <= restBits) is++;
        is--;
        size_t size = 1ULL << (is + 2);
        std::cout<<"size: "<<size<<endl;
        adpPerfectCF[i] = new cuckoofilter::AdpPerfectCF<size_t, 32, cntSize>(size);
        // if (REM_SIZE == 4) {
        //     perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(8192);//REM_SIZE 4
        // } else if (REM_SIZE == 8) {
        //     perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(1024);//REM_SIZE 8
        // } else if (REM_SIZE == 16) {
        //     perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(128);//REM_SIZE 16
        // } 
        // double alpha = 1 - perfectCF[i]->SizeInBytes() * 8.0 / (bits_per_key_sets[0] * dlSet[0].pos_keys_ip_.size());
        // std::cout<<"alpha: "<<alpha<<endl;
        // nfuncs = bits_per_key_sets[i] * alpha * log(2);
        
        
        // cbfArr[i] = new countingBloom::CountingBloom(bits_per_key_sets[0] *alpha* dlSet[0].pos_keys_.size(), bitPerCounter, nfuncs);
        dlSet[i].pos_keys_ipaddr_set.clear();
        auto t1 = steady_clock::now(); 
        for (auto& ip: dlSet[i].pos_keys_ip_ )
        {    
            cf->Add(ip);
            dlSet[i].pos_keys_ipaddr_set.insert(ip);
            // cfArr[i]->counting_bloom_add( *v);
        }
        std::cout<< "pos_keys_ipaddr_set size: " << dlSet[i].pos_keys_ipaddr_set.size() << endl;
        // 创建一个最大堆，元素是 <size_t ipaddr, size_t cnt>
        // std::priority_queue<std::pair<size_t, size_t>, std::vector<std::pair<size_t, size_t>>, CompareCnt> maxHeap;
        // size_t isOneNum = dlSet[i].pos_keys_.size() * vulNegKeyRatio;
        // for (int j=0; j<isOneNum; j++)
        // {
        //     Slice *v = dlSet[i].neg_keys_[j];
        //     // s[i]->insert(XXH3_128bits(&v->str[0], v->str.size()).low64);
        //     perfectCF[i]->Add32(v->ipaddr);
        // }
        // std::cout<<perfectCF[i]->Info()<<endl;
        // std::cout<<"set size: "<<s[i]->size()<<endl;

        auto t2 = steady_clock::now();
        std::cout<<cf->Info()<<endl;
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_ip_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;

        size_t fnr_c = 0;
        double wfpr_c = 0;
        size_t fpr_Num = 0;
        t1 = steady_clock::now();
        for (auto& ip: dlSet[i].pos_keys_ip_ )
        {
            // if (cf->Contain(ip) != cuckoofilter::Ok || perfectCF[i]->Contain32(ip)== cuckoofilter::Ok)
            // {
            //     fnr_c++;
            // }
            if(cf->Contain(ip) == cuckoofilter::Ok && adpPerfectCF[i]->Contain32(ip) != cuckoofilter::Ok) {
                assert(dlSet[i].pos_keys_ipaddr_set.find(ip) != dlSet[i].pos_keys_ipaddr_set.end());
            } else fnr_c++;

        }
        size_t cnt = 0;
        size_t cnt2 = 0;
        std::cout<<"123 pos_keys_ipaddr_set size: "<<dlSet[i].pos_keys_ipaddr_set.size()<<endl;
        for (auto& ip: dlSet[i].query_keys_ip_)
        {
            if(cf->Contain(ip) == cuckoofilter::Ok && adpPerfectCF[i]->Contain32(ip) != cuckoofilter::Ok)
            {
                wfpr_c ++;
                fpr_Num ++;
                if (dlSet[i].pos_keys_ipaddr_set.find(ip) == dlSet[i].pos_keys_ipaddr_set.end())
                {
                    adpPerfectCF[i]->Add32(ip);
                    // std::cout<<"123\n";
                    // if (perfectCF[i]->LoadFactor() < 0.95) {
                    //     perfectCF[i]->Add32(ip);
                    //     cnt2++;
                    // } else cnt++;
                }
            }

        }
        t2 = steady_clock::now();   
        std::cout<<"cnt: "<<cnt<<endl;
        std::cout<<"cnt2: "<<cnt2<<endl;
        std::cout<<"fpr_Num: "<<fpr_Num<<endl;
        std::cout<<"load factor: "<<adpPerfectCF[i]->LoadFactor()<<endl;
        
        tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_ip_.size()+dlSet[i].query_keys_ip_.size());
        sumQueryTime += tempV;
        minQueryTime = tempV < minQueryTime ? tempV : minQueryTime;
        maxQueryTime = tempV > maxQueryTime ? tempV : maxQueryTime;
    
        sumFPR += (double)(fpr_Num*100) / (double)dlSet[i].query_keys_ip_.size();

        sumWeightedFPR += (double)(wfpr_c*100) / dlSet[i].query_keys_ip_.size();
        if((double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size() > max_weightedFpr){
            max_weightedFpr = (double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size();
        }
        if((double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size() < min_weightedFpr){
            min_weightedFpr = (double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size();
        }
        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_ip_.size();
        delete cf;
    }
    std::cout<<"Average Insertion Latency (Time(ns)/Key): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";
    std::cout<<"Average Query Latency (Time(ns)/Key): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"sumFPR: "<< sumFPR <<endl;
    std::cout<<"datasetNum: "   << datasetNum <<endl;
    std::cout<<"sumFNR: "<< sumFNR <<endl;
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;
    wFPR[2].push_back((double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2));
    insert_latency[2].push_back((sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2));
    query_latency[2].push_back((sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2));
    for (size_t i = 0; i < datasetNum; i++)
    {
        delete adpPerfectCF[i];
    }
}

void TestPerfectCFCAIDA_adp_cf_AND_RCF_CAIDA(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio) 
{
    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;
    cuckoofilter::PerfectCF<size_t, 32> * perfectCF[datasetNum];
    cuckoofilter::RCF<size_t,32> * RCF[datasetNum];
    double sumFNR = 0;
    double sumQueryTime = 0;
    double minQueryTime = MAXFLOAT;
    double maxQueryTime = 0;
    double sumWeightedFPR = 0;
    double sumFPR = 0;
    double max_weightedFpr = 0;
    double min_weightedFpr = MAXFLOAT;
    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t totalBits = dlSet[i].pos_keys_ip_.size() * bits_per_key_sets[i];
        cuckoofilter::CuckooFilter<size_t, REM_SIZE+4-1>* cf = new cuckoofilter::CuckooFilter<size_t, REM_SIZE+4-1>(dlSet[i].pos_keys_ip_.size());
        size_t restBits = totalBits > cf->SizeInBytes() * 8 ? totalBits - cf->SizeInBytes() * 8 : cf->SizeInBytes() * 8 * 0.01;
        restBits = max(restBits, static_cast<size_t>(cf->SizeInBytes() * 8 * 0.01));
        std::cout<<"restBits: "<<restBits<<endl;
        int is = 1;
        while ((1ULL << (is+2)) * (33-is) <= restBits) is++;
        is--;
        size_t size = 1ULL << (is + 2);
        std::cout<<"size: "<<size<<endl;
        perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(size);
        size_t h1 = std::ceil(std::log2(dlSet[i].pos_keys_ipaddr_set.size())) - 2;
        size_t h2 = bits_per_key_sets[i]-2 <= 32 - 1 - h1 ? bits_per_key_sets[i]-2 : 32 - 1 - h1;
        restBits = totalBits > (1ULL << h1) * 4 * (h2+1) ? totalBits - (1ULL << h1) * 4 * (h2+1) : (1ULL << h1) * 4 * (h2+1) * 0.01;
        restBits = max(restBits, static_cast<size_t>((1ULL << h1) * 4 * (h2+1) * 0.01));
        is = 1;
        while ((1ULL << (is+2)) * (32 - h2 + 1 - is) <= restBits) is++;
        is--;
        size = 1ULL << (is + 2);
        std::cout<<"size: "<<size<<endl;
        RCF[i] = new cuckoofilter::RCF<size_t, 32>(dlSet[i].pos_keys_ipaddr_set.size(), bits_per_key_sets[i]-2, is);
    

        auto t1 = steady_clock::now(); 
        size_t cntt = 0;
        for (auto& ip: dlSet[i].pos_keys_ip_ )
        {    
            // if (++cntt > 10) break;
            cf->Add(ip);
            RCF[i]->Add(ip);
            
            // cfArr[i]->counting_bloom_add( *v);
        }
        // 创建一个最大堆，元素是 <size_t ipaddr, size_t cnt>
        // std::priority_queue<std::pair<size_t, size_t>, std::vector<std::pair<size_t, size_t>>, CompareCnt> maxHeap;
        // size_t isOneNum = dlSet[i].pos_keys_.size() * vulNegKeyRatio;
        // for (int j=0; j<isOneNum; j++)
        // {
        //     Slice *v = dlSet[i].neg_keys_[j];
        //     // s[i]->insert(XXH3_128bits(&v->str[0], v->str.size()).low64);
        //     perfectCF[i]->Add32(v->ipaddr);
        // }
        // std::cout<<perfectCF[i]->Info()<<endl;
        // std::cout<<"set size: "<<s[i]->size()<<endl;

        auto t2 = steady_clock::now();
        std::cout<<cf->Info()<<endl;
        std::cout<<RCF[i]->Info()<<endl;
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_ip_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;

        size_t fnr_c = 0;
        double wfpr_c = 0;
        size_t fpr_Num_CF = 0;
        size_t fpr_Num_RCF = 0;
        t1 = steady_clock::now();

        size_t cnt = 0;
        size_t cnt2 = 0;
        // std::set<size_t> tempSet(dlSet[i].query_keys_ip_.begin(), dlSet[i].query_keys_ip_.end());
        // cntt = 0;
        // for (auto& ip: tempSet)
        // {
        //     // if (++cntt > 10) break; 
        //     if(cf->Contain(ip) == cuckoofilter::Ok)
        //     {
        //         fpr_Num_CF ++;
        //     }
        //     if(RCF[i]->Contain(ip) == cuckoofilter::Ok)
        //     {
        //         fpr_Num_RCF ++;
        //     }
        // }
        // t2 = steady_clock::now();   
        // std::cout<<"fpr_Num_CF: "<<fpr_Num_CF<<endl;
        // std::cout<<"fpr_Num_RCF: "<<fpr_Num_RCF<<endl;
        
    } 
        
    for (size_t i = 0; i < datasetNum; i++)
    {
        delete perfectCF[i];
        delete RCF[i];
    }
}

void TestPerfectCFCAIDA_adp_cf_potential(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio) 
{
    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;
    cuckoofilter::PerfectCF<size_t, 32> * perfectCF[datasetNum];
    double sumFNR = 0;
    double sumQueryTime = 0;
    double minQueryTime = MAXFLOAT;
    double maxQueryTime = 0;
    double sumWeightedFPR = 0;
    double sumFPR = 0;
    double max_weightedFpr = 0;
    double min_weightedFpr = MAXFLOAT;
    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t totalBits = dlSet[i].pos_keys_ip_.size() * bits_per_key_sets[i];
        cuckoofilter::CuckooFilter<size_t, REM_SIZE+4-1>* cf = new cuckoofilter::CuckooFilter<size_t, REM_SIZE+4-1>(dlSet[i].pos_keys_ip_.size());
        size_t restBits = totalBits > cf->SizeInBytes() * 8 ? totalBits - cf->SizeInBytes() * 8 : cf->SizeInBytes() * 8 * 0.01;
        restBits = max(restBits, static_cast<size_t>(cf->SizeInBytes() * 8 * 0.01));
        std::cout<<"restBits: "<<restBits<<endl;
        int is = 1;
        while ((1ULL << (is+2)) * (33-is) <= restBits) is++;
        is--;
        size_t size = 1ULL << (is + 2);
        std::cout<<"size: "<<size<<endl;
        perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(size);
        vector<size_t>* potential[5];
        size_t thr = 12;
        double theta = 0.8;
        size_t cntBits = pow(2, std::ceil(std::log2(thr)));
        size_t potentialNum = (restBits - (1ULL << (is+2)) * (33-is)) / (32 + cntBits);
        std::cout<<"potentialNum: "<<potentialNum<<endl;
        potential[i] = new vector<size_t>;
        for (size_t j = 0; j < potentialNum; j++)
        {
            // size_t temp = 1ULL << 32;
            potential[i]->push_back(0);
        }
         dlSet[i].pos_keys_ipaddr_set.clear();
        auto t1 = steady_clock::now(); 
        for (auto& ip: dlSet[i].pos_keys_ip_ )
        {    
            cf->Add(ip);
            dlSet[i].pos_keys_ipaddr_set.insert(ip);
            // cfArr[i]->counting_bloom_add( *v);
        }
        std::cout<< "pos_keys_ipaddr_set size: " << dlSet[i].pos_keys_ipaddr_set.size() << endl;
        // 创建一个最大堆，元素是 <size_t ipaddr, size_t cnt>
        // std::priority_queue<std::pair<size_t, size_t>, std::vector<std::pair<size_t, size_t>>, CompareCnt> maxHeap;
        // size_t isOneNum = dlSet[i].pos_keys_.size() * vulNegKeyRatio;
        // for (int j=0; j<isOneNum; j++)
        // {
        //     Slice *v = dlSet[i].neg_keys_[j];
        //     // s[i]->insert(XXH3_128bits(&v->str[0], v->str.size()).low64);
        //     perfectCF[i]->Add32(v->ipaddr);
        // }
        // std::cout<<perfectCF[i]->Info()<<endl;
        // std::cout<<"set size: "<<s[i]->size()<<endl;

        auto t2 = steady_clock::now();
        std::cout<<cf->Info()<<endl;
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_ip_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;

        size_t fnr_c = 0;
        double wfpr_c = 0;
        size_t fpr_Num = 0;
        t1 = steady_clock::now();
        for (auto& ip: dlSet[i].pos_keys_ip_ )
        {
            // if (cf->Contain(ip) != cuckoofilter::Ok || perfectCF[i]->Contain32(ip)== cuckoofilter::Ok)
            // {
            //     fnr_c++;
            // }

            // if(cf->Contain(ip) == cuckoofilter::Ok && perfectCF[i]->Contain32(ip) != cuckoofilter::Ok) {
            //     assert(dlSet[i].pos_keys_ipaddr_set.find(ip) != dlSet[i].pos_keys_ipaddr_set.end());
            // } else fnr_c++;

            if(cf->Contain(ip) == cuckoofilter::Ok && perfectCF[i]->Contain32(ip) != cuckoofilter::Ok)
            {
                bool isPotential = false;
                bool isEmpty = false;
                size_t emptyIdx = -1;
                size_t smallestNum = 1 << cntBits;
                size_t smallestIdx = -1;
                for (size_t j = 0; j < potentialNum; j++)
                {
                    if (ip == (potential[i]->at(j) & 0xFFFFFFFF))
                    {
                        if (potential[i]->at(j) >> 32 == thr - 1) {
                            potential[i]->at(j) = 0;
                            perfectCF[i]->Add32(ip);
                        } else {
                            potential[i]->at(j) = ((((potential[i]->at(j) >> 32) + 1)) << 32) + ip;
                        }
                        isPotential = true;
                        break;
                    }
                    if (potential[i]->at(j) == 0)
                    {
                        isEmpty = true;
                        emptyIdx = j;
                    }
                    if ((potential[i]->at(j) >> 32) < smallestNum)
                    {
                        smallestNum = potential[i]->at(j) >> 32;
                        smallestIdx = j;
                    }
                }
                // std::cout<<"smallestNum: "<<smallestNum<<endl;
                // std::cout<<"smallestIdx: "<<smallestIdx<<endl;
                if (isPotential) continue;
                assert(dlSet[i].pos_keys_ipaddr_set.find(ip) != dlSet[i].pos_keys_ipaddr_set.end());
                if (dlSet[i].pos_keys_ipaddr_set.find(ip) == dlSet[i].pos_keys_ipaddr_set.end())
                {
                    if (isEmpty) {
                        // std::cout<<"emptyIdx: "<<emptyIdx<<endl;
                        potential[i]->at(emptyIdx) += ip + (1ULL << 32);
                        if ((potential[i]->at(emptyIdx) >> 32) == thr) {
                            perfectCF[i]->Add32(potential[i]->at(emptyIdx) & 0xFFFFFFFF);
                            potential[i]->at(emptyIdx) = 0;
                        } 
                    } else {
                        // std::cout<<"smallestIdx: "<<smallestIdx<<endl;
                        size_t randNum = rand() % (smallestNum + 1);
                        if (randNum < 1) {
                            // std::cout<<"cnt2: "<<cnt2<<endl;
                            potential[i]->at(smallestIdx) = ((smallestNum + 1) << 32) + ip;
                        } else potential[i]->at(smallestIdx) = ((smallestNum + 1) << 32) + (potential[i]->at(smallestIdx) & 0xFFFFFFFF);
                        if ((potential[i]->at(smallestIdx) >> 32) == thr) {
                            perfectCF[i]->Add32(potential[i]->at(smallestIdx) & 0xFFFFFFFF);
                            potential[i]->at(smallestIdx) = 0;
                        }

                    }
                    // if (perfectCF[i]->LoadFactor() < 0.95) {
                    //     perfectCF[i]->Add32(ip);
                    //     cnt2++;
                    // } else cnt++;
                }
            } else {
                fnr_c++;
            }

        }
        size_t cnt = 0;
        size_t cnt2 = 0;
        for (auto& ip: dlSet[i].query_keys_ip_)
        {
            if(cf->Contain(ip) == cuckoofilter::Ok && perfectCF[i]->Contain32(ip) != cuckoofilter::Ok)
            {
                // bool isPotential = false;
                // bool isEmpty = true;
                // size_t emptyIdx = 0;
                // size_t smallestNum = 1 << cntBits;
                // size_t smallestIdx = -1;
                bool isPotential = false;
                bool isEmpty = false;
                size_t emptyIdx = -1;
                size_t smallestNum = 1 << cntBits;
                size_t smallestIdx = -1;
                for (size_t j = 0; j < potentialNum; j++)
                {
                    if (ip == (potential[i]->at(j) & 0xFFFFFFFF))
                    {
                        if (potential[i]->at(j) >> 32 == thr - 1) {
                            potential[i]->at(j) = 0;
                            perfectCF[i]->Add32(ip);
                        } else {
                            potential[i]->at(j) = ((((potential[i]->at(j) >> 32) + 1)) << 32) + ip;
                        }
                        isPotential = true;
                        break;
                    }
                    if (potential[i]->at(j) == 0)
                    {
                        isEmpty = true;
                        emptyIdx = j;
                    }
                    if ((potential[i]->at(j) >> 32) < smallestNum)
                    {
                        smallestNum = potential[i]->at(j) >> 32;
                        smallestIdx = j;
                    }
                }
                // std::cout<<"smallestNum: "<<smallestNum<<endl;
                // std::cout<<"smallestIdx: "<<smallestIdx<<endl;
                if (isPotential) continue;
                wfpr_c ++;
                fpr_Num ++;
                assert(dlSet[i].pos_keys_ipaddr_set.find(ip) == dlSet[i].pos_keys_ipaddr_set.end());
                if (dlSet[i].pos_keys_ipaddr_set.find(ip) == dlSet[i].pos_keys_ipaddr_set.end())
                {
                    cnt++;
                    if (perfectCF[i]->LoadFactor() < theta) {
                        perfectCF[i]->Add32(ip);
                    } else {
                        if (isEmpty) {
                            // std::cout<<"emptyIdx: "<<emptyIdx<<endl;
                            potential[i]->at(emptyIdx) = ip + (1ULL << 32);
                            if ((potential[i]->at(emptyIdx) >> 32) == thr) {
                                cnt2++;
                                perfectCF[i]->Add32(potential[i]->at(emptyIdx) & 0xFFFFFFFF);
                                potential[i]->at(emptyIdx) = 0;
                            } 
                        } else {
                            
                            // std::cout<<"smallestIdx: "<<smallestIdx<<endl;
                            size_t randNum = rand() % (smallestNum + 1);
                            if (randNum < 1) {
                                
                                // std::cout<<"cnt2: "<<cnt2<<endl;
                                potential[i]->at(smallestIdx) = ((smallestNum + 1) << 32) + ip;
                            } else potential[i]->at(smallestIdx) = ((smallestNum + 1) << 32) + (potential[i]->at(smallestIdx) & 0xFFFFFFFF);
                            if ((potential[i]->at(smallestIdx) >> 32) == thr) {
                                perfectCF[i]->Add32(potential[i]->at(smallestIdx) & 0xFFFFFFFF);
                                potential[i]->at(smallestIdx) = 0;
                            }
                        }
                    }
                    // if (perfectCF[i]->LoadFactor() < 0.95) {
                    //     perfectCF[i]->Add32(ip);
                    //     cnt2++;
                    // } else cnt++;
                }
            }

        }
        t2 = steady_clock::now();   
        std::cout<<"cnt: "<<cnt<<endl;
        std::cout<<"cnt2: "<<cnt2<<endl;
        std::cout<<"fpr_Num: "<<fpr_Num<<endl;
        std::cout<<"load factor: "<<perfectCF[i]->LoadFactor()<<endl;
        
        tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_ip_.size()+dlSet[i].query_keys_ip_.size());
        sumQueryTime += tempV;
        minQueryTime = tempV < minQueryTime ? tempV : minQueryTime;
        maxQueryTime = tempV > maxQueryTime ? tempV : maxQueryTime;
    
        sumFPR += (double)(fpr_Num*100) / (double)dlSet[i].query_keys_ip_.size();

        sumWeightedFPR += (double)(wfpr_c*100) / dlSet[i].query_keys_ip_.size();
        if((double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size() > max_weightedFpr){
            max_weightedFpr = (double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size();
        }
        if((double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size() < min_weightedFpr){
            min_weightedFpr = (double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size();
        }
        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_ip_.size();
        delete cf;
    }
    std::cout<<"Average Insertion Latency (Time(ns)/Key): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";
    std::cout<<"Average Query Latency (Time(ns)/Key): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"sumFPR: "<< sumFPR <<endl;
    std::cout<<"datasetNum: "   << datasetNum <<endl;
    std::cout<<"sumFNR: "<< sumFNR <<endl;
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;
    wFPR[4].push_back((double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2));
    insert_latency[4].push_back((sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2));
    query_latency[4].push_back((sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2));
    for (size_t i = 0; i < datasetNum; i++)
    {
        delete perfectCF[i];
    }
}

void TestPerfectCFCAIDA_adp(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio) 
{
    size_t nfuncs;
    size_t bitPerCounter = 4;
    countingBloom::CountingBloom * cbfArr[datasetNum];
    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;
    cuckoofilter::PerfectCF<size_t, 32> * perfectCF[datasetNum];
    for (size_t i = 0; i < datasetNum; i++)
    {

        // perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(dlSet[i].pos_keys_.size()*vulNegKeyRatio*2);
        if (REM_SIZE == 4) {
            perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(4096);//REM_SIZE 4
        } else if (REM_SIZE == 8) {
            perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(512);//REM_SIZE 8
        } else if (REM_SIZE == 16) {
            perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(16);//REM_SIZE 16
        } else if (REM_SIZE == 2) {
            perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(8192);//REM_SIZE 32
        } 
        
    
        // size_t totalBits = pow(2, std::ceil(std::log2(dlSet[i].pos_keys_.size()))) * (bits_per_key_sets[i] - 1);
        size_t totalBits = max(dlSet[i].pos_keys_ip_.size() * bits_per_key_sets[i],
                               pow(2, std::ceil(std::log2(dlSet[i].pos_keys_ip_.size()))) * (bits_per_key_sets[i] - 1));   
        std::cout<<"totalBits: "<<totalBits<<endl;
        double new_bits_per_key_sets = (totalBits - perfectCF[i]->SizeInBytes() * 8.0) / dlSet[i].pos_keys_ip_.size();
        std::cout<<"new_bits_per_key_sets: "<<new_bits_per_key_sets<<endl;
        nfuncs = new_bits_per_key_sets * log(2);
        cbfArr[i] = new countingBloom::CountingBloom(totalBits - perfectCF[i]->SizeInBytes() * 8, bitPerCounter, nfuncs);

    }

    for (size_t i = 0; i < datasetNum; i++)
    {
        dlSet[i].pos_keys_ipaddr_set.clear();
        auto t1 = steady_clock::now(); 
        for (auto& ip : dlSet[i].pos_keys_ip_ )
        {    
            cbfArr[i]->counting_bloom_add_IP32(ip);
            dlSet[i].pos_keys_ipaddr_set.insert(ip);
            // Slice * key = new Slice();
            // key->ipaddr = ip;
            // key->str = std::to_string(ip);
            // cbfArr[i]->counting_bloom_add(*key);
        }
        // 创建一个最大堆，元素是 <size_t ipaddr, size_t cnt>
        // std::priority_queue<std::pair<size_t, size_t>, std::vector<std::pair<size_t, size_t>>, CompareCnt> maxHeap;
        // size_t isOneNum = dlSet[i].pos_keys_.size() * vulNegKeyRatio;
        // for (int j=0; j<isOneNum; j++)
        // {
        //     Slice *v = dlSet[i].neg_keys_[j];
        //     // s[i]->insert(XXH3_128bits(&v->str[0], v->str.size()).low64);
        //     perfectCF[i]->Add32(v->ipaddr);
        // }
        // std::cout<<perfectCF[i]->Info()<<endl;
        // std::cout<<"set size: "<<s[i]->size()<<endl;

        auto t2 = steady_clock::now();
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_ip_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;
    }

    std::cout<<"Average Insertion Latency (Time(ns)/Key): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";

    double sumFNR = 0;
    double sumQueryTime = 0;
    double minQueryTime = MAXFLOAT;
    double maxQueryTime = 0;
    double sumWeightedFPR = 0;
    double sumFPR = 0;
    double max_weightedFpr = 0;
    double min_weightedFpr = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0;
        double wfpr_c = 0;
        size_t fpr_Num = 0;

        auto t1 = steady_clock::now();
        for (auto& ip : dlSet[i].pos_keys_ip_)
        {
            // Slice * key = new Slice();
            // key->ipaddr = ip;
            // key->str = std::to_string(ip);
            // if (!cbfArr[i]->counting_bloom_check_IP32(ip) || perfectCF[i]->Contain32(ip)== cuckoofilter::Ok)
            // // if (!cbfArr[i]->counting_bloom_check(*key) || perfectCF[i]->Contain32(ip)== cuckoofilter::Ok)
            // {
            //     fnr_c++;
            // }
            if(cbfArr[i]->counting_bloom_check_IP32(ip) && perfectCF[i]->Contain32(ip) != cuckoofilter::Ok)
            // if(cbfArr[i]->counting_bloom_check(*key) && perfectCF[i]->Contain32(ip) != cuckoofilter::Ok)
            {
                assert(dlSet[i].pos_keys_ipaddr_set.find(ip) != dlSet[i].pos_keys_ipaddr_set.end());
            } else fnr_c++;

        }
        size_t cnt = 0;
        size_t cnt2 = 0;
        for (auto& ip : dlSet[i].query_keys_ip_)
        {
            // Slice * key = new Slice();
            // key->ipaddr = ip;
            // key->str = std::to_string(ip);
            if(cbfArr[i]->counting_bloom_check_IP32(ip) && perfectCF[i]->Contain32(ip) != cuckoofilter::Ok)
            // if(cbfArr[i]->counting_bloom_check(*key) && perfectCF[i]->Contain32(ip) != cuckoofilter::Ok)
            {
                wfpr_c ++;
                fpr_Num ++;
                if (dlSet[i].pos_keys_ipaddr_set.find(ip) == dlSet[i].pos_keys_ipaddr_set.end())
                {
                    if (perfectCF[i]->LoadFactor() < 0.95) {
                        perfectCF[i]->Add32(ip);
                        cnt2++;
                    } else cnt++;
                }
            }

        }
        auto t2 = steady_clock::now();   
        std::cout<<"cnt: "<<cnt<<endl;
        std::cout<<"cnt2: "<<cnt2<<endl;
        std::cout<<"fpr_Num: "<<fpr_Num<<endl;
        std::cout<<"load factor: "<<perfectCF[i]->LoadFactor()<<endl;
        
        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_ip_.size()+dlSet[i].query_keys_ip_.size());
        sumQueryTime += tempV;
        minQueryTime = tempV < minQueryTime ? tempV : minQueryTime;
        maxQueryTime = tempV > maxQueryTime ? tempV : maxQueryTime;
    
        sumFPR += (double)(fpr_Num*100) / (double)dlSet[i].query_keys_ip_.size();

        sumWeightedFPR += (double)(wfpr_c*100) / dlSet[i].query_keys_ip_.size();
        if((double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size() > max_weightedFpr){
            max_weightedFpr = (double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size();
        }
        if((double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size() < min_weightedFpr){
            min_weightedFpr = (double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size();
        }
        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_ip_.size();
    }
    
    std::cout<<"Average Query Latency (Time(ns)/Key): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"sumFPR: "<< sumFPR <<endl;
    std::cout<<"datasetNum: "   << datasetNum <<endl;
    std::cout<<"sumFNR: "<< sumFNR <<endl;
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;
    wFPR[3].push_back((double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2));
    insert_latency[3].push_back((sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2));
    query_latency[3].push_back((sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2));   
    // double sumDeleteTime = 0;
    // double minDeleteTime = MAXFLOAT;
    // double maxDeleteTime = 0;
    // for (size_t i = 0; i < datasetNum; i++)
    // {
    //     //random shuffle
    //     vector<Slice *> deleteItem(dlSet[i].pos_keys_.size());
    //     for(size_t j=0; j<dlSet[i].pos_keys_.size(); j++)
    //     {
    //         deleteItem[j] = dlSet[i].pos_keys_[j]; 
    //     }
    //     random_shuffle(deleteItem.begin(),deleteItem.end());
        
    //     auto t1 = steady_clock::now();
    //     for (Slice *v : deleteItem )
    //     {
    //         cbfArr[i]->counting_bloom_remove(*v);
    //     }
    //     auto t2 = steady_clock::now();       
    //     double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size());  
    //     sumDeleteTime += tempV;
    //     minDeleteTime = tempV < minDeleteTime ? tempV : minDeleteTime;
    //     maxDeleteTime = tempV > maxDeleteTime ? tempV : maxDeleteTime;
    //     cbfArr[i]->printHaveInsertNum();
    //     deleteItem.clear();
    // }  
    // std::cout<<"Average Deletion Latency (Time(ns)/Key) : "<<(sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2)<<endl;
    
    // double sumFNR_afterDelete=0;
    // for (size_t i = 0; i < datasetNum; i++)
    // {
    //     size_t fnr_c = 0;
    //     auto t1 = steady_clock::now();
    //     for (Slice *v : dlSet[i].pos_keys_ )
    //     {
    //         if(!cbfArr[i]->counting_bloom_check(*v)){
    //             fnr_c++;
    //         } 
    //     }
    
    //     sumFNR_afterDelete += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    // }
    // std::cout<<"FNR(After delete): "<< (double)sumFNR_afterDelete / datasetNum << "%" <<endl; 
    for (size_t i = 0; i < datasetNum; i++)
    {
        delete cbfArr[i];
        delete perfectCF[i];
    }
}

void TestAdpPerfectCFCAIDA_adp(dataloader * dlSet, double * bits_per_key_sets) 
{
    size_t nfuncs;
    size_t bitPerCounter = 4;
    countingBloom::CountingBloom * cbfArr[datasetNum];
    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;
    // const size_t cntSize = REM_SIZE >= 8 ? 0 : 1;
    const size_t cntSize = 1;
    // cuckoofilter::PerfectCF<size_t, 32> * perfectCF[datasetNum];
    cuckoofilter::AdpPerfectCF<size_t, 32, cntSize> * adpPerfectCF[datasetNum];
    for (size_t i = 0; i < datasetNum; i++)
    {

        // perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(dlSet[i].pos_keys_.size()*vulNegKeyRatio*2);
        // if (REM_SIZE == 4) {
        //     adpPerfectCF[i] = new cuckoofilter::AdpPerfectCF<size_t, 32, cntSize>(4096);//REM_SIZE 4
        // } else if (REM_SIZE == 8) {
        //     adpPerfectCF[i] = new cuckoofilter::AdpPerfectCF<size_t, 32, cntSize>(512);//REM_SIZE 8
        // } else if (REM_SIZE == 16) {
        //     adpPerfectCF[i] = new cuckoofilter::AdpPerfectCF<size_t, 32, cntSize>(16);//REM_SIZE 16
        // } else if (REM_SIZE == 2) {
        //     adpPerfectCF[i] = new cuckoofilter::AdpPerfectCF<size_t, 32, cntSize>(8192);//REM_SIZE 32
        // } 
        
    
        // size_t totalBits = pow(2, std::ceil(std::log2(dlSet[i].pos_keys_.size()))) * (bits_per_key_sets[i] - 1);
        size_t totalBits = max(dlSet[i].pos_keys_ip_.size() * bits_per_key_sets[i],
                               pow(2, std::ceil(std::log2(dlSet[i].pos_keys_ip_.size()))) * (bits_per_key_sets[i] - 1));   
        std::cout<<"totalBits: "<<totalBits<<endl;
        size_t adpPerfectCFSizeMax = totalBits * 0.3 / 32;
        if (REM_SIZE == 4) {
            size_t adpPerfectCFSizeMin = totalBits * 0.1 / 32;
            size_t adpPerfectCFSize = adpPerfectCFSizeMax > 4096 ? 4096 : adpPerfectCFSizeMax;
            adpPerfectCFSize = adpPerfectCFSize < adpPerfectCFSizeMin ? adpPerfectCFSizeMin : adpPerfectCFSize;
            adpPerfectCF[i] = new cuckoofilter::AdpPerfectCF<size_t, 32, cntSize>(adpPerfectCFSize);//REM_SIZE 4
        } else if (REM_SIZE == 8) {
            size_t adpPerfectCFSizeMin = totalBits * 0.01 / 32;
            size_t adpPerfectCFSize = adpPerfectCFSizeMax > 512 ? 512 : adpPerfectCFSizeMax;
            adpPerfectCFSize = adpPerfectCFSize < adpPerfectCFSizeMin ? adpPerfectCFSizeMin : adpPerfectCFSize;
            adpPerfectCF[i] = new cuckoofilter::AdpPerfectCF<size_t, 32, cntSize>(adpPerfectCFSize);//REM_SIZE 8
        } else if (REM_SIZE == 16) {
            size_t adpPerfectCFSizeMin = totalBits * 0.005 / 32;
            size_t adpPerfectCFSize = adpPerfectCFSizeMax > 16 ? 16 : adpPerfectCFSizeMax;
            adpPerfectCFSize = adpPerfectCFSize < adpPerfectCFSizeMin ? adpPerfectCFSizeMin : adpPerfectCFSize;
            adpPerfectCF[i] = new cuckoofilter::AdpPerfectCF<size_t, 32, cntSize>(adpPerfectCFSize);//REM_SIZE 16
        } else if (REM_SIZE == 2) {
            size_t adpPerfectCFSizeMin = totalBits * 0.2 / 32;
            size_t adpPerfectCFSize = adpPerfectCFSizeMax > 8192 ? 8192 : adpPerfectCFSizeMax;
            adpPerfectCFSize = adpPerfectCFSize < adpPerfectCFSizeMin ? adpPerfectCFSizeMin : adpPerfectCFSize;
            adpPerfectCF[i] = new cuckoofilter::AdpPerfectCF<size_t, 32, cntSize>(adpPerfectCFSize);//REM_SIZE 32
        } 
        double new_bits_per_key_sets = (totalBits - adpPerfectCF[i]->SizeInBytes() * 8.0) / dlSet[i].pos_keys_ip_.size();
        std::cout<<"new_bits_per_key_sets: "<<new_bits_per_key_sets<<endl;
        nfuncs = new_bits_per_key_sets * log(2);
        cbfArr[i] = new countingBloom::CountingBloom(totalBits - adpPerfectCF[i]->SizeInBytes() * 8, bitPerCounter, nfuncs);

        
    }

    for (size_t i = 0; i < datasetNum; i++)
    {
        dlSet[i].pos_keys_ipaddr_set.clear();
        auto t1 = steady_clock::now(); 
        for (auto& ip : dlSet[i].pos_keys_ip_ )
        {    
            cbfArr[i]->counting_bloom_add_IP32(ip);
            dlSet[i].pos_keys_ipaddr_set.insert(ip);
            // Slice * key = new Slice();
            // key->ipaddr = ip;
            // key->str = std::to_string(ip);
            // cbfArr[i]->counting_bloom_add(*key);
        }
        // 创建一个最大堆，元素是 <size_t ipaddr, size_t cnt>
        // std::priority_queue<std::pair<size_t, size_t>, std::vector<std::pair<size_t, size_t>>, CompareCnt> maxHeap;
        // size_t isOneNum = dlSet[i].pos_keys_.size() * vulNegKeyRatio;
        // for (int j=0; j<isOneNum; j++)
        // {
        //     Slice *v = dlSet[i].neg_keys_[j];
        //     // s[i]->insert(XXH3_128bits(&v->str[0], v->str.size()).low64);
        //     perfectCF[i]->Add32(v->ipaddr);
        // }
        // std::cout<<perfectCF[i]->Info()<<endl;
        // std::cout<<"set size: "<<s[i]->size()<<endl;

        auto t2 = steady_clock::now();
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_ip_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;
    }

    std::cout<<"Average Insertion Latency (Time(ns)/Key): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";

    double sumFNR = 0;
    double sumQueryTime = 0;
    double minQueryTime = MAXFLOAT;
    double maxQueryTime = 0;
    double sumWeightedFPR = 0;
    double sumFPR = 0;
    double max_weightedFpr = 0;
    double min_weightedFpr = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0;
        double wfpr_c = 0;
        size_t fpr_Num = 0;

        auto t1 = steady_clock::now();
        for (auto& ip : dlSet[i].pos_keys_ip_)
        {
            // Slice * key = new Slice();
            // key->ipaddr = ip;
            // key->str = std::to_string(ip);
            // if (!cbfArr[i]->counting_bloom_check_IP32(ip) || perfectCF[i]->Contain32(ip)== cuckoofilter::Ok)
            // // if (!cbfArr[i]->counting_bloom_check(*key) || perfectCF[i]->Contain32(ip)== cuckoofilter::Ok)
            // {
            //     fnr_c++;
            // }
            if(cbfArr[i]->counting_bloom_check_IP32(ip) && adpPerfectCF[i]->Contain32(ip) != cuckoofilter::Ok)
            // if(cbfArr[i]->counting_bloom_check(*key) && perfectCF[i]->Contain32(ip) != cuckoofilter::Ok)
            {
                assert(dlSet[i].pos_keys_ipaddr_set.find(ip) != dlSet[i].pos_keys_ipaddr_set.end());
            } else fnr_c++;

        }
        size_t cnt = 0;
        size_t cnt2 = 0;
        for (auto& ip : dlSet[i].query_keys_ip_)
        {
            // Slice * key = new Slice();
            // key->ipaddr = ip;
            // key->str = std::to_string(ip);
            if(cbfArr[i]->counting_bloom_check_IP32(ip) && adpPerfectCF[i]->Contain32(ip) != cuckoofilter::Ok)
            // if(cbfArr[i]->counting_bloom_check(*key) && perfectCF[i]->Contain32(ip) != cuckoofilter::Ok)
            {
                wfpr_c ++;
                fpr_Num ++;
                if (dlSet[i].pos_keys_ipaddr_set.find(ip) == dlSet[i].pos_keys_ipaddr_set.end())
                {
                    adpPerfectCF[i]->Add32(ip);
                    // if (adpPerfectCF[i]->LoadFactor() < 0.95) {
                        
                    //     cnt2++;
                    // } else cnt++;
                }
            }

        }
        auto t2 = steady_clock::now();   
        std::cout<<"cnt: "<<cnt<<endl;
        std::cout<<"cnt2: "<<cnt2<<endl;
        std::cout<<"fpr_Num: "<<fpr_Num<<endl;
        std::cout<<"load factor: "<<adpPerfectCF[i]->LoadFactor()<<endl;
        
        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_ip_.size()+dlSet[i].query_keys_ip_.size());
        sumQueryTime += tempV;
        minQueryTime = tempV < minQueryTime ? tempV : minQueryTime;
        maxQueryTime = tempV > maxQueryTime ? tempV : maxQueryTime;
    
        sumFPR += (double)(fpr_Num*100) / (double)dlSet[i].query_keys_ip_.size();

        sumWeightedFPR += (double)(wfpr_c*100) / dlSet[i].query_keys_ip_.size();
        if((double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size() > max_weightedFpr){
            max_weightedFpr = (double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size();
        }
        if((double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size() < min_weightedFpr){
            min_weightedFpr = (double)(wfpr_c*100) /  dlSet[i].query_keys_ip_.size();
        }
        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_ip_.size();
    }
    
    std::cout<<"Average Query Latency (Time(ns)/Key): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"sumFPR: "<< sumFPR <<endl;
    std::cout<<"datasetNum: "   << datasetNum <<endl;
    std::cout<<"sumFNR: "<< sumFNR <<endl;
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;
    wFPR[3].push_back((double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2));
    insert_latency[3].push_back((sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2));
    query_latency[3].push_back((sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2));   
    // double sumDeleteTime = 0;
    // double minDeleteTime = MAXFLOAT;
    // double maxDeleteTime = 0;
    // for (size_t i = 0; i < datasetNum; i++)
    // {
    //     //random shuffle
    //     vector<Slice *> deleteItem(dlSet[i].pos_keys_.size());
    //     for(size_t j=0; j<dlSet[i].pos_keys_.size(); j++)
    //     {
    //         deleteItem[j] = dlSet[i].pos_keys_[j]; 
    //     }
    //     random_shuffle(deleteItem.begin(),deleteItem.end());
        
    //     auto t1 = steady_clock::now();
    //     for (Slice *v : deleteItem )
    //     {
    //         cbfArr[i]->counting_bloom_remove(*v);
    //     }
    //     auto t2 = steady_clock::now();       
    //     double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size());  
    //     sumDeleteTime += tempV;
    //     minDeleteTime = tempV < minDeleteTime ? tempV : minDeleteTime;
    //     maxDeleteTime = tempV > maxDeleteTime ? tempV : maxDeleteTime;
    //     cbfArr[i]->printHaveInsertNum();
    //     deleteItem.clear();
    // }  
    // std::cout<<"Average Deletion Latency (Time(ns)/Key) : "<<(sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2)<<endl;
    
    // double sumFNR_afterDelete=0;
    // for (size_t i = 0; i < datasetNum; i++)
    // {
    //     size_t fnr_c = 0;
    //     auto t1 = steady_clock::now();
    //     for (Slice *v : dlSet[i].pos_keys_ )
    //     {
    //         if(!cbfArr[i]->counting_bloom_check(*v)){
    //             fnr_c++;
    //         } 
    //     }
    
    //     sumFNR_afterDelete += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    // }
    // std::cout<<"FNR(After delete): "<< (double)sumFNR_afterDelete / datasetNum << "%" <<endl; 
    for (size_t i = 0; i < datasetNum; i++)
    {
        delete cbfArr[i];
        delete adpPerfectCF[i];
    }
}

void TestPerfectCFCAIDA_adp_potential(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio) 
{
    size_t nfuncs;
    size_t bitPerCounter = 4;
    countingBloom::CountingBloom * cbfArr[datasetNum];
    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;
    cuckoofilter::PerfectCF<size_t, 32> * perfectCF[datasetNum];
    vector<size_t>* potential[5];
    size_t potentialNum = 128;
    size_t thr = 2;
    size_t cntBits = pow(2, std::ceil(std::log2(thr)));
    for (size_t i = 0; i < datasetNum; i++)
    {

        // perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(dlSet[i].pos_keys_.size()*vulNegKeyRatio*2);
        if (REM_SIZE == 4) {
            perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(8192);//REM_SIZE 4
            // perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(256);//REM_SIZE 4
        } else if (REM_SIZE == 8) {
            perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(1024);//REM_SIZE 8
            // perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(1024);//REM_SIZE 8
        } else if (REM_SIZE == 16) {
            perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(16);//REM_SIZE 16
        } 
        potential[i] = new vector<size_t>;
        for (size_t j = 0; j < potentialNum; j++)
        {
            // size_t temp = 1ULL << 32;
            potential[i]->push_back(0);
        }
        size_t totalBits = pow(2, std::ceil(std::log2(dlSet[i].pos_keys_.size()))) * bits_per_key_sets[i];
        std::cout<<"totalBits: "<<totalBits<<endl;
        double new_bits_per_key_sets = (totalBits - perfectCF[i]->SizeInBytes() * 8.0 - potentialNum * (32 + cntBits)) / dlSet[i].pos_keys_.size();
        std::cout<<"new_bits_per_key_sets: "<<new_bits_per_key_sets<<endl;
        nfuncs = new_bits_per_key_sets * log(2);
        cbfArr[i] = new countingBloom::CountingBloom(totalBits - perfectCF[i]->SizeInBytes() * 8.0 - potentialNum * (32 + cntBits), bitPerCounter, nfuncs);

    }

    for (size_t i = 0; i < datasetNum; i++)
    {
        auto t1 = steady_clock::now(); 
        for (Slice *v : dlSet[i].pos_keys_ )
        {    
            cbfArr[i]->counting_bloom_add( *v);
        }
        // 创建一个最大堆，元素是 <size_t ipaddr, size_t cnt>
        // std::priority_queue<std::pair<size_t, size_t>, std::vector<std::pair<size_t, size_t>>, CompareCnt> maxHeap;
        // size_t isOneNum = dlSet[i].pos_keys_.size() * vulNegKeyRatio;
        // for (int j=0; j<isOneNum; j++)
        // {
        //     Slice *v = dlSet[i].neg_keys_[j];
        //     // s[i]->insert(XXH3_128bits(&v->str[0], v->str.size()).low64);
        //     perfectCF[i]->Add32(v->ipaddr);
        // }
        // std::cout<<perfectCF[i]->Info()<<endl;
        // std::cout<<"set size: "<<s[i]->size()<<endl;

        auto t2 = steady_clock::now();
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;
    }

    std::cout<<"Average Insertion Latency (Time(ns)/Key): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";

    double sumFNR = 0;
    double sumQueryTime = 0;
    double minQueryTime = MAXFLOAT;
    double maxQueryTime = 0;
    double sumWeightedFPR = 0;
    double sumFPR = 0;
    double max_weightedFpr = 0;
    double min_weightedFpr = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0;
        double wfpr_c = 0;
        size_t fpr_Num = 0;
        double total_cost_ = 0;
        std::cout<<"neg_keys size: "<<dlSet[i].neg_keys_.size()<<endl;
        for (Slice *v : dlSet[i].neg_keys_ )
        {
            total_cost_ += v->cost;
        }

        auto t1 = steady_clock::now();
        for (Slice *v : dlSet[i].pos_keys_ )
        {
            if(cbfArr[i]->counting_bloom_check(*v) && perfectCF[i]->Contain32(v->ipaddr) != cuckoofilter::Ok)
            {
                bool isPotential = false;
                bool isEmpty = false;
                size_t emptyIdx = -1;
                size_t smallestNum = 1 << cntBits;
                size_t smallestIdx = -1;
                for (size_t j = 0; j < potentialNum; j++)
                {
                    if (v->ipaddr == potential[i]->at(j) & 0xFFFFFFFF)
                    {
                        if (potential[i]->at(j) >> 32 == thr - 1) {
                            potential[i]->at(j) = 0;
                            perfectCF[i]->Add32(v->ipaddr);
                        } else {
                            potential[i]->at(j) = ((((potential[i]->at(j) >> 32) + 1)) << 32) + v->ipaddr;
                        }
                        isPotential = true;
                        break;
                    }
                    if (potential[i]->at(j) == 0)
                    {
                        isEmpty = true;
                        emptyIdx = j;
                    }
                    if ((potential[i]->at(j) >> 32) < smallestNum)
                    {
                        smallestNum = potential[i]->at(j) >> 32;
                        smallestIdx = j;
                    }
                }
                if (isPotential) {
                    fnr_c++;
                    continue;
                }
                if (dlSet[i].pos_keys_ipaddr_set.find(v->ipaddr) == dlSet[i].pos_keys_ipaddr_set.end())
                {
                    fnr_c++;
                    if (isEmpty) {
                        potential[i]->at(emptyIdx) += v->ipaddr;
                        if ((potential[i]->at(emptyIdx) >> 32) == thr) {
                            perfectCF[i]->Add32(potential[i]->at(emptyIdx) & 0xFFFFFFFF);
                            potential[i]->at(emptyIdx) = 0;
                        }
                    } else {
                        size_t cntNum = potential[i]->at(smallestIdx) >> 32;
                        size_t randNum = rand() % (cntNum + 1);
                        if (randNum < 1) {
                            potential[i]->at(smallestIdx) = ((cntNum + 1) << 32) + v->ipaddr;
                        } else potential[i]->at(smallestIdx) = ((cntNum + 1) << 32) + (potential[i]->at(smallestIdx) & 0xFFFFFFFF);
                        if ((potential[i]->at(smallestIdx) >> 32) == thr) {
                            perfectCF[i]->Add32(potential[i]->at(smallestIdx) & 0xFFFFFFFF);
                            potential[i]->at(smallestIdx) = 0;
                        }
                    }
                    // if (perfectCF[i]->LoadFactor() < 0.95) {
                    //     perfectCF[i]->Add32(v->ipaddr);
                    //     cnt2++;
                    // } else cnt++;
                }
            }

        }
        size_t cnt = 0;
        size_t cnt2 = 0;
        for (Slice *v : dlSet[i].query_keys )
        {
            if(cbfArr[i]->counting_bloom_check(*v) && perfectCF[i]->Contain32(v->ipaddr) != cuckoofilter::Ok)
            {
                bool isPotential = false;
                bool isEmpty = false;
                size_t emptyIdx = -1;
                size_t smallestNum = 1 << cntBits;
                size_t smallestIdx = -1;
                for (size_t j = 0; j < potentialNum; j++)
                {
                    if (v->ipaddr == (potential[i]->at(j) & 0xFFFFFFFF))
                    {
                        if (potential[i]->at(j) >> 32 == thr - 1) {
                            potential[i]->at(j) = 0;
                            perfectCF[i]->Add32(v->ipaddr);
                        } else {
                            potential[i]->at(j) = ((((potential[i]->at(j) >> 32) + 1)) << 32) + v->ipaddr;
                        }
                        isPotential = true;
                        break;
                    }
                    if (potential[i]->at(j) == 0)
                    {
                        isEmpty = true;
                        emptyIdx = j;
                    }
                    if ((potential[i]->at(j) >> 32) < smallestNum)
                    {
                        smallestNum = potential[i]->at(j) >> 32;
                        smallestIdx = j;
                    }
                }
                // std::cout<<"smallestNum: "<<smallestNum<<endl;
                // std::cout<<"smallestIdx: "<<smallestIdx<<endl;
                if (isPotential) continue;
                wfpr_c ++;
                fpr_Num ++;
                if (dlSet[i].pos_keys_ipaddr_set.find(v->ipaddr) == dlSet[i].pos_keys_ipaddr_set.end())
                {
                    
                    if (isEmpty) {
                        // std::cout<<"emptyIdx: "<<emptyIdx<<endl;
                        potential[i]->at(emptyIdx) += v->ipaddr + (1ULL << 32);
                        if ((potential[i]->at(emptyIdx) >> 32) == thr) {
                            perfectCF[i]->Add32(potential[i]->at(emptyIdx) & 0xFFFFFFFF);
                            potential[i]->at(emptyIdx) = 0;
                        } 
                    } else {
                        cnt++;
                        // std::cout<<"smallestIdx: "<<smallestIdx<<endl;
                        size_t randNum = rand() % (smallestNum + 1);
                        if (randNum < 1) {
                            cnt2++;
                            // std::cout<<"cnt2: "<<cnt2<<endl;
                            potential[i]->at(smallestIdx) = ((smallestNum + 1) << 32) + v->ipaddr;
                        } else potential[i]->at(smallestIdx) = ((smallestNum + 1) << 32) + (potential[i]->at(smallestIdx) & 0xFFFFFFFF);
                        if ((potential[i]->at(smallestIdx) >> 32) == thr) {
                            perfectCF[i]->Add32(potential[i]->at(smallestIdx) & 0xFFFFFFFF);
                            potential[i]->at(smallestIdx) = 0;
                        }

                    }
                    // if (perfectCF[i]->LoadFactor() < 0.95) {
                    //     perfectCF[i]->Add32(v->ipaddr);
                    //     cnt2++;
                    // } else cnt++;
                    // for (size_t j = 0; j < potentialNum; j++)
                    // {
                    //     size_t cnt = potential[i]->at(j) >> 32;
                    //     // std::cout<<"cnt: "<<cnt<<endl;
                    //     size_t ipaddr = potential[i]->at(j) & 0xFFFFFFFF;
                    //     // std::cout<<"ipaddr: "<<ipaddr<<endl;
                    


                    // }
                }
            }

        }
        std::cout<<"cnt: "<<cnt<<endl;
        std::cout<<"cnt2: "<<cnt2<<endl;
        std::cout<<"fpr_Num: "<<fpr_Num<<endl;
        std::cout<<"load factor: "<<perfectCF[i]->LoadFactor()<<endl;
        auto t2 = steady_clock::now();   
        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size()+dlSet[i].query_keys.size());
        sumQueryTime += tempV;
        minQueryTime = tempV < minQueryTime ? tempV : minQueryTime;
        maxQueryTime = tempV > maxQueryTime ? tempV : maxQueryTime;
    
        sumFPR += (double)(fpr_Num*100) / (double)dlSet[i].query_keys.size();

        sumWeightedFPR += (double)(wfpr_c*100) / dlSet[i].query_keys.size();
        if((double)(wfpr_c*100) /  dlSet[i].query_keys.size() > max_weightedFpr){
            max_weightedFpr = (double)(wfpr_c*100) /  dlSet[i].query_keys.size();
        }
        if((double)(wfpr_c*100) /  dlSet[i].query_keys.size() < min_weightedFpr){
            min_weightedFpr = (double)(wfpr_c*100) /  dlSet[i].query_keys.size();
        }
        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    }
    
    std::cout<<"Average Query Latency (Time(ns)/Key): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"sumFPR: "<< sumFPR <<endl;
    std::cout<<"datasetNum: "   << datasetNum <<endl;
    std::cout<<"sumFNR: "<< sumFNR <<endl;
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;
    wFPR[4].push_back((double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2));
    insert_latency[4].push_back((sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2));
    query_latency[4].push_back((sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2));   
    // std::cout<<"~~~~~~~CBF Start Deleting~~~~~~~\n";
    double sumDeleteTime = 0;
    double minDeleteTime = MAXFLOAT;
    double maxDeleteTime = 0;
    for (size_t i = 0; i < datasetNum; i++)
    {
        //random shuffle
        vector<Slice *> deleteItem(dlSet[i].pos_keys_.size());
        for(size_t j=0; j<dlSet[i].pos_keys_.size(); j++)
        {
            deleteItem[j] = dlSet[i].pos_keys_[j]; 
        }
        random_shuffle(deleteItem.begin(),deleteItem.end());
        
        auto t1 = steady_clock::now();
        for (Slice *v : deleteItem )
        {
            cbfArr[i]->counting_bloom_remove(*v);
        }
        auto t2 = steady_clock::now();       
        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size());  
        sumDeleteTime += tempV;
        minDeleteTime = tempV < minDeleteTime ? tempV : minDeleteTime;
        maxDeleteTime = tempV > maxDeleteTime ? tempV : maxDeleteTime;
        cbfArr[i]->printHaveInsertNum();
    }  
    std::cout<<"Average Deletion Latency (Time(ns)/Key) : "<<(sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2)<<endl;
    
    double sumFNR_afterDelete=0;
    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0;
        auto t1 = steady_clock::now();
        for (Slice *v : dlSet[i].pos_keys_ )
        {
            if(!cbfArr[i]->counting_bloom_check(*v)){
                fnr_c++;
            } 
        }
    
        sumFNR_afterDelete += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    }
    std::cout<<"FNR(After delete): "<< (double)sumFNR_afterDelete / datasetNum << "%" <<endl; 
}

void TestPerfectCFCAIDA(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio) 
{
    size_t nfuncs;
    size_t bitPerCounter = 4;
    countingBloom::CountingBloom * cbfArr[5];
    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;
    cuckoofilter::PerfectCF<size_t, 32> * perfectCF[5];
    for (size_t i = 0; i < datasetNum; i++)
    {
        // bits_per_key_sets[i] *= 0.9;
        perfectCF[i] = new cuckoofilter::PerfectCF<size_t, 32>(dlSet[i].pos_keys_.size()*vulNegKeyRatio);
        double alpha = 1 - perfectCF[i]->SizeInBytes() * 8.0 / (bits_per_key_sets[0] * dlSet[0].pos_keys_.size());
        std::cout<<"alpha: "<<alpha<<endl;
        nfuncs = bits_per_key_sets[i] * alpha/ 4 * log(2);
        cbfArr[i] = new countingBloom::CountingBloom(bits_per_key_sets[0] *alpha* dlSet[0].pos_keys_.size() /bitPerCounter, bitPerCounter, nfuncs);
    }

    for (size_t i = 0; i < datasetNum; i++)
    {
        auto t1 = steady_clock::now(); 
        for (Slice *v : dlSet[i].pos_keys_ )
        {    
            cbfArr[i]->counting_bloom_add( *v);
        }
        size_t isOneNum = dlSet[i].pos_keys_.size() * vulNegKeyRatio;
        for (int j=0; j<isOneNum; j++)
        {
            Slice *v = dlSet[i].neg_keys_[j];
            // s[i]->insert(XXH3_128bits(&v->str[0], v->str.size()).low64);
            perfectCF[i]->Add32(v->ipaddr);
        }
        // std::cout<<perfectCF[i]->Info()<<endl;
        // std::cout<<"set size: "<<s[i]->size()<<endl;

        auto t2 = steady_clock::now();
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;
    }

    std::cout<<"Average Insertion Latency (Time(ns)/Key): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";

    double sumFNR = 0;
    double sumQueryTime = 0;
    double minQueryTime = MAXFLOAT;
    double maxQueryTime = 0;
    double sumWeightedFPR = 0;
    double sumFPR = 0;
    double max_weightedFpr = 0;
    double min_weightedFpr = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0;
        double wfpr_c = 0;
        size_t fpr_Num = 0;
        double total_cost_ = 0;
        std::cout<<"neg_keys size: "<<dlSet[i].neg_keys_.size()<<endl;
        for (Slice *v : dlSet[i].neg_keys_ )
        {
            total_cost_ += v->cost;
        }

        auto t1 = steady_clock::now();
        for (Slice *v : dlSet[i].pos_keys_ )
        {
            if (!cbfArr[i]->counting_bloom_check(*v) || perfectCF[i]->Contain32(v->ipaddr)== cuckoofilter::Ok)
            {
                fnr_c++;
            }

        }
        
        for (Slice *v : dlSet[i].neg_keys_ )
        {
            if(cbfArr[i]->counting_bloom_check(*v) && perfectCF[i]->Contain32(v->ipaddr)!= cuckoofilter::Ok)
            {
                wfpr_c += v->cost;
                fpr_Num ++;
            }

        }
        std::cout<<"fpr_Num: "<<fpr_Num<<endl;
        auto t2 = steady_clock::now();   
        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size()+dlSet[i].neg_keys_.size());
        sumQueryTime += tempV;
        minQueryTime = tempV < minQueryTime ? tempV : minQueryTime;
        maxQueryTime = tempV > maxQueryTime ? tempV : maxQueryTime;
    
        sumFPR += (double)(fpr_Num*100) / (double)dlSet[i].neg_keys_.size();

        sumWeightedFPR += (double)(wfpr_c*100) / total_cost_;
        if((double)(wfpr_c*100) / total_cost_ > max_weightedFpr){
            max_weightedFpr = (double)(wfpr_c*100) / total_cost_;
        }
        if((double)(wfpr_c*100) / total_cost_ < min_weightedFpr){
            min_weightedFpr = (double)(wfpr_c*100) / total_cost_;
        }
        sumFNR += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    }
    
    std::cout<<"Average Query Latency (Time(ns)/Key): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"sumFPR: "<< sumFPR <<endl;
    std::cout<<"datasetNum: "   << datasetNum <<endl;
    std::cout<<"sumFNR: "<< sumFNR <<endl;
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;
    // wFPR[0].push_back((double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2));
    // insert_latency[0].push_back((sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2));
    // query_latency[0].push_back((sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2));   
    std::cout<<"~~~~~~~CBF Start Deleting~~~~~~~\n";
    double sumDeleteTime = 0;
    double minDeleteTime = MAXFLOAT;
    double maxDeleteTime = 0;
    for (size_t i = 0; i < datasetNum; i++)
    {
        //random shuffle
        vector<Slice *> deleteItem(dlSet[i].pos_keys_.size());
        for(size_t j=0; j<dlSet[i].pos_keys_.size(); j++)
        {
            deleteItem[j] = dlSet[i].pos_keys_[j]; 
        }
        random_shuffle(deleteItem.begin(),deleteItem.end());
        
        auto t1 = steady_clock::now();
        for (Slice *v : deleteItem )
        {
            cbfArr[i]->counting_bloom_remove(*v);
        }
        auto t2 = steady_clock::now();       
        double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size());  
        sumDeleteTime += tempV;
        minDeleteTime = tempV < minDeleteTime ? tempV : minDeleteTime;
        maxDeleteTime = tempV > maxDeleteTime ? tempV : maxDeleteTime;
        cbfArr[i]->printHaveInsertNum();
    }  
    std::cout<<"Average Deletion Latency (Time(ns)/Key) : "<<(sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2)<<endl;
    
    double sumFNR_afterDelete=0;
    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t fnr_c = 0;
        auto t1 = steady_clock::now();
        for (Slice *v : dlSet[i].pos_keys_ )
        {
            if(!cbfArr[i]->counting_bloom_check(*v)){
                fnr_c++;
            } 
        }
    
        sumFNR_afterDelete += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    }
    std::cout<<"FNR(After delete): "<< (double)sumFNR_afterDelete / datasetNum << "%" <<endl; 
}

bool compare_func(Slice * key1, Slice * key2){
    return key1->cost > key2->cost;
}