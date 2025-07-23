#include <time.h>
#include <vector>
#include <set>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <cassert>
#include <cstdlib>

#include "util/dataloader.h"
#include "nonlearnedfilter/SSCF/SeesawCF.h"
#include "nonlearnedfilter/countingBloom.h"
#include "nonlearnedfilter/wcbf.h"
#include "nonlearnedfilter/stackedfilter/StackedFilter.h"
#include "nonlearnedfilter/cuckoofilter/src/perfectCF.h"


#define SHALLA_NAME "shalla_"
#define YCSBT_NAME "ycsbt_"
#define HOSEFIRE_NAME "hosefire_"

using namespace std;
using namespace std::chrono;
using cuckoofilter::PerfectCF;

void TestSSCF(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio);
void TestCountingBloom(dataloader * dlSet, double * bits_per_key_sets);
void TestWeightedCountingBloom(dataloader * dlSet,  double * bits_per_key_sets);
void TestStackedFilter(dataloader * dlSet,  double * bits_per_key_sets, double vulNegKeyRatio);
// void TestSimple(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio);

void TestPerfectCF(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio);
void TestPerfectCFCAIDA(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio);

// void TestTelescopingFilter(dataloader * dlSet);
// void TestTelescopingFilterCAIDA(dataloader * dlSet);

void ConvertDataset(dataloader *dlSet, std::vector<std::vector<StringElement>> &vec_positives, std::vector<std::vector<StringElement>> &vec_negatives, std::vector<std::vector<double>> &vec_cdf);
bool compare_func(Slice * key1, Slice * key2);
bool compare_WcbfKey_func(weightedCountingBloom::WCBFKey &key1, weightedCountingBloom::WCBFKey &key2);


vector<vector<double>> wFPR = vector<vector<double>>(5);
vector<vector<double>> insert_latency = vector<vector<double>>(5);
vector<vector<double>> query_latency = vector<vector<double>>(5);
vector<vector<double>> delete_latency = vector<vector<double>>(5);
vector<double> load_factor = vector<double>();


double vulNegKeyRatio =0.05;
const size_t datasetNum = 5; 
// 0->CAIDA; 1->Shalla; 2->YCSB 
vector<int> dataSets = {1, 0};
// vector<double> skewnessSet = {0.25, 0.5, 0.75, 0.99, 1.0, 1.5, 2.0, 2.5};
vector<vector<double>> skewness = {{1.0}, {1.0, 1.5, 2.0, 2.5}, {}};
void display() {
    vector<StringElement> filters = {
        StringElement("RF"),
        StringElement("SSCF"),
        StringElement("CBF"),
        StringElement("WCBF"),
        StringElement("SF")
    };
    // display wFPR
    std::cout << "--------------wFPR testing--------------" << endl;
    for (size_t i = 0; i < wFPR.size(); i++)
    {
        std::cout << filters[i].get_value() << ": ";
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
        std::cout << filters[i].get_value() << ": ";
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
        std::cout << filters[i].get_value() << ": ";
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
            dataloader * dlSet = new dataloader[5];
            for (size_t i = 0; i < datasetNum; i++)
            {
                switch (dataSet)
                {
                case 0:
                    // load CAIDA
                    if(!dlSet[i].loadZipf_semi(to_string(i),s))
                    {
                        std::cout<< "error!\n";
                    }
                    break;
                case 1:
                    // load shalla
                    if(!dlSet[i].loadShalla_semi(to_string(i), s))
                    {
                        std::cout<< "error!\n";
                    }
                    break;
                case 2: 
                    // load ycsb
                    if(!dlSet[i].loadYCSB_semi(to_string(i),s))
                    {
                        std::cout<< "error!\n";
                    }
                    break;
                }
            }
            for(size_t i = 0; i<datasetNum; i++){
                if(!is_sorted(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func)){
                    sort(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func);
                }
            }

            double *bits_per_key_sets = new double[5];
            for (int round = 0; round < 1 ; round++){
                // for (int bpk = 12; bpk <= 36; bpk += 1) {
                for (int bpk = 20; bpk <= 20; bpk += 1) {
                    std::cout << "--------------bits per key: " << bpk << "--------------" << endl;
                    for (size_t i = 0; i < datasetNum; i++) {
                        bits_per_key_sets[i] = bpk;
                    }
                    if (dataSet == 1) {
                        std::cout << "--------------PerfectCFShalla testing--------------" << endl;
                        TestPerfectCF(dlSet, bits_per_key_sets, vulNegKeyRatio);
                    } else {
                        std::cout << "--------------PerfectCFCAIDA testing--------------" << endl;
                        TestPerfectCFCAIDA(dlSet, bits_per_key_sets, vulNegKeyRatio);
                    }

                    std::cout << "--------------SSCF testing--------------" << endl;
                    TestSSCF(dlSet, bits_per_key_sets, vulNegKeyRatio);

                    std::cout << "--------------countingBloom testing--------------" << endl;
                    TestCountingBloom(dlSet, bits_per_key_sets);

                    std::cout << "--------------WCBF testing--------------" << endl;
                    TestWeightedCountingBloom(dlSet, bits_per_key_sets);

                    std::cout << "--------------stackedFilter testing--------------" << endl;
                    TestStackedFilter(dlSet, bits_per_key_sets, vulNegKeyRatio);   
                    
                }
            }
            delete[] dlSet;
            delete[] bits_per_key_sets;
            std::cout.clear();
            std::cout << "--------------display--------------" << endl;
            std::cout << "skewness: "  << s <<endl;
            display();
        }    
    }


}


void TestSSCF(dataloader * dlSet,  double * bits_per_key_sets, double vulNegKeyRatio)
{
    SeesawCF::SeesawCF *sscfArr[datasetNum];

    for (size_t i = 0; i < datasetNum; i++)
    {
        if(!is_sorted(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func)){
            std::cout<<"error!\n";
        }
        sscfArr[i] = new SeesawCF::SeesawCF(vulNegKeyRatio, bits_per_key_sets[i], dlSet[i].pos_keys_.size(), dlSet[i].neg_keys_.size(), dlSet[i].neg_keys_ , true);
    }
    
    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        auto t1 = steady_clock::now();

        for (Slice *v : dlSet[i].pos_keys_ )
        {          
            sscfArr[i]->Add( *v);
        }

        auto t2 = steady_clock::now();
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;
    }
    
    std::cout<<"Average Insertion Latency (Time/Key ns): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";

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
        size_t j=0;
        for (j=0; j<dlSet[i].neg_keys_.size(); j++)
        {
            total_cost_ += dlSet[i].neg_keys_[j]->cost;
        }        

        auto t1 = steady_clock::now();

        for (j=0; j<dlSet[i].pos_keys_.size(); j++)
        {
            if (!sscfArr[i]->Contain(dlSet[i].pos_keys_[j]))
            {
                fnr_c++;
            } 
        }
            
        for (j=0; j<dlSet[i].neg_keys_.size(); j++)
        {
            if (sscfArr[i]->Contain(dlSet[i].neg_keys_[j]))
            {
                wfpr_c += dlSet[i].neg_keys_[j]->cost;
                fpr_Num ++;
            }
        }
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
    
    std::cout<<"Average Query Latency (Time/Key ns): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;
    wFPR[1].push_back((double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2));
    insert_latency[1].push_back((sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2));
    query_latency[1].push_back((sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2));
    // std::cout<<"~~~~~~~SSCF Start Deleting~~~~~~~\n";
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
    //     for( Slice * v : deleteItem)
    //     {
    //         sscfArr[i]->Deletion(v);
    //     }
    //     auto t2 = steady_clock::now();   

    //     double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size());  
    //     sumDeleteTime += tempV;
    //     minDeleteTime = tempV < minDeleteTime ? tempV : minDeleteTime;
    //     maxDeleteTime = tempV > maxDeleteTime ? tempV : maxDeleteTime;
    //     sscfArr[i]->PrintCounterAboveOne();
    // }    
    // std::cout<<"Average Deletion Latency (Time/Key ns): "<<(sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2)<<endl;
    // delete_latency[1].push_back((sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2));
    // double sumFNR_afterDelete = 0;
    // for (size_t i = 0; i < datasetNum; i++)
    // {
        
    //     size_t fnr_c = 0;
    //     for(size_t j=0; j<dlSet[i].pos_keys_.size(); j++)
    //     {
    //         if(!sscfArr[i]->Contain(dlSet[i].pos_keys_[j]))
    //         {
    //             fnr_c++;
    //         } 
    //     }
    //     sumFNR_afterDelete += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();
    // }    
    // std::cout<<"FNR(After delete): "<< (double)sumFNR_afterDelete / datasetNum << "%" <<endl;

    for (size_t i = 0; i < datasetNum; i++) delete sscfArr[i] ;
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
        isOneNum = isOneNum > dlSet[i].neg_keys_.size() ? dlSet[i].neg_keys_.size() : isOneNum;
        for (int j=0; j<isOneNum; j++)
        {
            Slice *v = dlSet[i].neg_keys_[j];
            // s[i]->insert(XXH3_128bits(&v->str[0], v->str.size()).low64);
            perfectCF[i]->Add32(v->ipaddr);
            
        }
        std::cout<<perfectCF[i]->Info()<<endl;
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
    wFPR[0].push_back((double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2));
    insert_latency[0].push_back((sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2));
    query_latency[0].push_back((sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2));   
    // std::cout<<"~~~~~~~CBF Start Deleting~~~~~~~\n";
    for (size_t i = 0; i < datasetNum; i++)
    {
        delete cbfArr[i];
        delete perfectCF[i];
    }
}

void TestPerfectCF(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio) 
{
    size_t nfuncs;
    size_t bitPerCounter = 4;
    countingBloom::CountingBloom * cbfArr[5];
    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;
    cuckoofilter::PerfectCF<size_t> * perfectCF[5];
    for (size_t i = 0; i < datasetNum; i++)
    {
        // bits_per_key_sets[i] *= 0.9;
        perfectCF[i] = new cuckoofilter::PerfectCF<size_t>(dlSet[i].pos_keys_.size()*vulNegKeyRatio);
        double alpha = 1 - perfectCF[i]->SizeInBytes() * 8.0 / (bits_per_key_sets[0] * dlSet[0].pos_keys_.size());
        std::cout<<"alpha: "<<alpha<<endl;
        nfuncs = bits_per_key_sets[i] * alpha/ 4 * log(2);
        cbfArr[i] = new countingBloom::CountingBloom(bits_per_key_sets[0] *alpha* dlSet[0].pos_keys_.size() /bitPerCounter, bitPerCounter, nfuncs);
    }
    std::cout<<"start adding---------------------"<<endl;
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
            size_t temp=XXH3_128bits(&v->str[0], v->str.size()).low64;
            perfectCF[i]->Add64(temp);
        }
        std::cout<<perfectCF[i]->Info()<<endl;
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
            if (!cbfArr[i]->counting_bloom_check(*v) || perfectCF[i]->Contain64(XXH3_128bits(&v->str[0], v->str.size()).low64)== cuckoofilter::Ok) {
                fnr_c++;
            }

        }
        
        for (Slice *v : dlSet[i].neg_keys_ )
        {
            if(cbfArr[i]->counting_bloom_check(*v) && perfectCF[i]->Contain64(XXH3_128bits(&v->str[0], v->str.size()).low64)!= cuckoofilter::Ok)
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
    wFPR[0].push_back((double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2));
    std::cout<<"~~~~~~~CBF Start Deleting~~~~~~~\n";
    insert_latency[0].push_back((sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2));
    query_latency[0].push_back((sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2));
    // load_factor.push_back(perfectCF[0]->LoadFactor());
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
    // }  
    // std::cout<<"Average Deletion Latency (Time(ns)/Key) : "<<(sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2)<<endl;
    // delete_latency[0].push_back((sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2));
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




void TestCountingBloom(dataloader * dlSet, double * bits_per_key_sets){

    size_t nfuncs;
    size_t bitPerCounter = 4;
    countingBloom::CountingBloom * cbfArr[5];
    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        nfuncs = bits_per_key_sets[i] / 4 * log(2);
        cbfArr[i] = new countingBloom::CountingBloom(bits_per_key_sets[0] * dlSet[0].pos_keys_.size() /bitPerCounter, bitPerCounter, nfuncs);
    }

    for (size_t i = 0; i < datasetNum; i++)
    {
        auto t1 = steady_clock::now(); 

        for (Slice *v : dlSet[i].pos_keys_ )
        {    
            cbfArr[i]->counting_bloom_add( *v);
        }   

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

        for (Slice *v : dlSet[i].neg_keys_ )
        {
            total_cost_ += v->cost;
        }

        auto t1 = steady_clock::now();
        for (Slice *v : dlSet[i].pos_keys_ )
        {
            if (!cbfArr[i]->counting_bloom_check(*v))
            {
                fnr_c++;
            } 
        }
        
        for (Slice *v : dlSet[i].neg_keys_ )
        {
            if(cbfArr[i]->counting_bloom_check(*v))
            {
                wfpr_c += v->cost;
                fpr_Num ++;
            }  
        }
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
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;
    wFPR[2].push_back((double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2));
    std::cout<<"~~~~~~~CBF Start Deleting~~~~~~~\n";
    insert_latency[2].push_back((sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2));
    query_latency[2].push_back((sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2));
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
    // }  
    // std::cout<<"Average Deletion Latency (Time(ns)/Key) : "<<(sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2)<<endl;
    // delete_latency[2].push_back((sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2));
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
    }
}

void TestWeightedCountingBloom(dataloader * dlSet, double * bits_per_key_sets)
{
    size_t bits_Per_Counter = 4;
    weightedCountingBloom::WCBF *wcbfArr[datasetNum];
    for (size_t i = 0; i < datasetNum; i++)
    {
        wcbfArr[i] =  new weightedCountingBloom::WCBF(bits_per_key_sets[i], bits_Per_Counter, dlSet[i].pos_keys_, dlSet[i].neg_keys_);
    }

    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t pos_n = dlSet[i].pos_keys_.size();
        size_t neg_n = dlSet[i].neg_keys_.size();
        weightedCountingBloom::WCBFKey* wbf_pos_keys_ = new weightedCountingBloom::WCBFKey[dlSet[i].pos_keys_.size()];
        weightedCountingBloom::WCBFKey* wbf_neg_keys_ = new weightedCountingBloom::WCBFKey[dlSet[i].neg_keys_.size()];
        
        for(int j=0; j<pos_n; ++j){
             wbf_pos_keys_[j].data_ = dlSet[i].pos_keys_[j];
        }
        for(int j=0; j<neg_n; ++j){
            wbf_neg_keys_[j].data_ = dlSet[i].neg_keys_[j];
        }
      
        if(!is_sorted(wbf_neg_keys_, wbf_neg_keys_+neg_n, compare_WcbfKey_func)){
            std::sort(wbf_neg_keys_, wbf_neg_keys_+neg_n, compare_WcbfKey_func);
        }            

        auto t1 = steady_clock::now();
        wcbfArr[i]->AddAll_latency(wbf_pos_keys_, wbf_neg_keys_, pos_n, neg_n);
        auto t2 = steady_clock::now();
        
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;
        delete [] wbf_pos_keys_;
        delete [] wbf_neg_keys_;
    }
    std::cout<<"Average Insertion (Time(ns)/Key): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";

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
        double wfpr_c = 0;
        size_t fnr_c = 0;
        double total_cost_ = 0.0;
        size_t fpr_Num = 0;
        for (int j=0; j<dlSet[i].neg_keys_.size(); j++)
        {
            total_cost_ += dlSet[i].neg_keys_[j]->cost;
        }

        auto t1 = steady_clock::now();

        for (int j=0; j<dlSet[i].pos_keys_.size(); j++)
        {
            if(!wcbfArr[i]->Contain(*dlSet[i].pos_keys_[j]))
            {
                fnr_c++;
            } 
        }
            
        for (int j=0; j<dlSet[i].neg_keys_.size(); j++)
        {
            if (wcbfArr[i]->Contain(*dlSet[i].neg_keys_[j]))
            {
                wfpr_c += dlSet[i].neg_keys_[j]->cost;
                fpr_Num++;
            }
        }
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

    std::cout<<"Average Query Latency (Time(ns)/key): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR: "<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;
    wFPR[3].push_back((double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2));
    std::cout<<"~~~~~~~WCBF Start Deleting~~~~~~~\n";
    insert_latency[3].push_back((sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2));
    query_latency[3].push_back((sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2));
    // double sumDeleteTime = 0;
    // double minDeleteTime = MAXFLOAT;
    // double maxDeleteTime = 0;
    // double sum_DeleteTime=0;
    // for (size_t i = 0; i < datasetNum; i++)
    // {
    //     vector<Slice *> deleteItem(dlSet[i].pos_keys_.size());
    //     for(size_t j=0; j<dlSet[i].pos_keys_.size(); j++)
    //     {
    //         deleteItem[j] = dlSet[i].pos_keys_[j]; 
    //     }
    //     random_shuffle(deleteItem.begin(),deleteItem.end());        
 
    //     auto t1 = steady_clock::now();

    //     wcbfArr[i]->DeleteAll(deleteItem);

    //     auto t2 = steady_clock::now();

    //     double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size());  
    //     sumDeleteTime += tempV;
    //     minDeleteTime = tempV < minDeleteTime ? tempV : minDeleteTime;
    //     maxDeleteTime = tempV > maxDeleteTime ? tempV : maxDeleteTime;
    //     wcbfArr[i]->PrintHaveInsertNum();
    // }
    // std::cout<<"Average Deletion Latency (Time(ns)/Key) : "<<(sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2)<<endl;
    // delete_latency[3].push_back((sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2));
    // double sumFNR_afterDelete = 0;
    // for (size_t i = 0; i < datasetNum; i++)
    // {
    //     int fnr_c = 0;
    //     for(int j=0; j<dlSet[i].pos_keys_.size(); j++)
    //     {
    //         if(!wcbfArr[i]->Contain(*dlSet[i].pos_keys_[j]))
    //         {
    //             fnr_c++;
    //         } 
    //     }
    //     sumFNR_afterDelete += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();   
    // }
    // std::cout<<"FNR(After delete): "<< (double)sumFNR_afterDelete / datasetNum << "%" <<endl;
    for (size_t i = 0; i < datasetNum; i++)
    {
        delete wcbfArr[i];
    }
}

void TestStackedFilter(dataloader * dlSet, double * bits_per_key_sets, double vulNegKeyRatio)
{     
    std::vector<std::vector<StringElement>> vec_positives;
    std::vector<std::vector<StringElement>> vec_negatives;
    std::vector<std::vector<double>> vec_cdf;
    StackedFilter<CountingBloomFilterLayer, StringElement> * filterArr[datasetNum];
    size_t bit_Per_Counter =4;
    vec_positives.resize(datasetNum);
    vec_negatives.resize(datasetNum);
    vec_cdf.resize(datasetNum);

    for (size_t i = 0; i < datasetNum; i++)
    {
        vec_negatives[i].resize(dlSet[i].neg_keys_.size());
        vec_positives[i].resize(dlSet[i].pos_keys_.size());
        vec_cdf[i].resize(dlSet[i].neg_keys_.size());
    }
    
    ConvertDataset(dlSet, vec_positives, vec_negatives, vec_cdf);

    for (size_t i = 0; i < datasetNum; i++)
    {
        filterArr[i] =  new StackedFilter<CountingBloomFilterLayer, StringElement>(bits_per_key_sets[i] * vec_positives[i].size(), vec_positives[i], vec_negatives[i], vec_cdf[i], bit_Per_Counter, vulNegKeyRatio);
    }

    double sumInsertTime = 0;
    double maxInsertTime = 0;
    double minInsertTime = MAXFLOAT;

    for (size_t i = 0; i < datasetNum; i++)
    {
        size_t haveAddCount= 0;
        size_t bits_per_item_ite = 50;
        size_t sumSize_temp = bits_per_key_sets[i] * vec_positives[i].size();

        auto t1 = steady_clock::now();

        for(StringElement &s : vec_positives[i])
        {
            filterArr[i]->InsertPositiveElement(s);
        }

        auto t2 = steady_clock::now();
        double tempV = 1000000000.0 *  duration<double>(t2 - t1).count()/ (double)(dlSet[i].pos_keys_.size()) ;
        sumInsertTime += tempV;
        minInsertTime = tempV < minInsertTime ? tempV : minInsertTime;
        maxInsertTime = tempV > maxInsertTime ? tempV : maxInsertTime;
    }
    std::cout<<"Average Insertion (Time(ns)/Key): "<<(sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2)<<"\n";

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
        size_t fnr_c = 0 ;
        size_t fpr_Num = 0;
        double total_cost_ = 0;
        double wfpr_c = 0;
        for (int j=0; j<vec_negatives[i].size(); ++j)
        {
            total_cost_ += vec_cdf[i][j];
        }

        auto t1 = steady_clock::now();
        for (int j=0; j<vec_positives[i].size(); ++j)
        {
            if (!filterArr[i]->LookupElement(vec_positives[i][j]))
            {
                fnr_c++;
            }
        }

        for (int j=0; j<vec_negatives[i].size(); ++j)
        {
            if (filterArr[i]->LookupElement(vec_negatives[i][j]))
            {
                wfpr_c +=  vec_cdf[i][j];
                fpr_Num++;
            }
        }
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

    std::cout<<"Average Query Latency (Time(ns)/key): "<< (sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2)<<endl;
    std::cout<<"Average Weighted FPR:"<< (double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2) << "%" <<endl;
    std::cout<<"Average FPR: "<< (double)sumFPR / datasetNum << "%" <<endl;
    std::cout<<"Average FNR (Before delete): "<< (double)sumFNR / datasetNum << "%" <<endl;
    wFPR[4].push_back((double)(sumWeightedFPR - max_weightedFpr - min_weightedFpr) / (datasetNum-2));
    std::cout<<"~~~~~~~SF Start Deleting~~~~~~~\n";
    insert_latency[4].push_back((sumInsertTime - minInsertTime - maxInsertTime) / (datasetNum - 2));
    query_latency[4].push_back((sumQueryTime - minQueryTime - maxQueryTime) / (datasetNum - 2));
    // double sumDeleteTime = 0;
    // double minDeleteTime = MAXFLOAT;
    // double maxDeleteTime = 0;
    // for (size_t i = 0; i < datasetNum; i++)
    // {

    //     random_shuffle(vec_positives[i].begin(),vec_positives[i].end());
    //     auto t1 = steady_clock::now();
    //     for(int j=0; j<vec_positives[i].size(); ++j){
    //         filterArr[i]->DeleteElement(vec_positives[i][j]);
    //     }
    //     auto t2 = steady_clock::now();
    //     double tempV = 1000000000.0 * duration<double>(t2 - t1).count() / (double)(dlSet[i].pos_keys_.size());  
    //     sumDeleteTime += tempV;
    //     minDeleteTime = tempV < minDeleteTime ? tempV : minDeleteTime;
    //     maxDeleteTime = tempV > maxDeleteTime ? tempV : maxDeleteTime;
    // }   
    // std::cout<<"~~~~~~~SF Finish Deleting~~~~~~~\n";  
    // std::cout<<"Average Deletion Latency (Time(ns)/Key) : "<<(sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2)<<endl;  
    // delete_latency[4].push_back((sumDeleteTime - minDeleteTime - maxDeleteTime) / (datasetNum - 2));
    // double sumFNR_afterDelete = 0;
    // for (size_t i = 0; i < datasetNum; i++)
    // {
    //     size_t fnr_c = 0 ;
    //     for (int j=0; j<vec_positives[i].size(); ++j){
    //         if (!filterArr[i]->LookupElement(vec_positives[i][j]))
    //         {
    //             fnr_c++;
    //         }
    //     }  
    //     sumFNR_afterDelete += (double)(fnr_c*100) / (double)dlSet[i].pos_keys_.size();        
    // } 
    // std::cout<<"FNR(After delete): "<< (double)sumFNR_afterDelete / datasetNum << "%" <<endl;    
    for (size_t i = 0; i < datasetNum; i++)
    {
        delete filterArr[i];
    }
}

void ConvertDataset(dataloader *dlSet, std::vector<std::vector<StringElement>> &vec_positives, std::vector<std::vector<StringElement>> &vec_negatives, std::vector<std::vector<double>> &vec_cdf){
    
    for (size_t i = 0; i < datasetNum; i++)
    {

        if(!is_sorted(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func)){
            sort(dlSet[i].neg_keys_.begin(), dlSet[i].neg_keys_.end(), compare_func);
        }
 
        for(int j=0; j<dlSet[i].pos_keys_.size(); ++j){
            vec_positives[i][j] = dlSet[i].pos_keys_[j]->str;
        }
        for(int j=0; j<dlSet[i].neg_keys_.size(); ++j){
            vec_negatives[i][j] = dlSet[i].neg_keys_[j]->str;
        }
        for(int j=0; j<dlSet[i].neg_keys_.size(); ++j){
            vec_cdf[i][j] = dlSet[i].neg_keys_[j]->cost;
        }        
    }
}

bool compare_func(Slice * key1, Slice * key2){
    return key1->cost > key2->cost;
}
bool compare_WcbfKey_func(weightedCountingBloom::WCBFKey &key1, weightedCountingBloom::WCBFKey &key2){
    return key1.data_->cost > key2.data_->cost;
}
