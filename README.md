# Reconfirm Filter

# About this repo
This repo contains all the source codes of RF together with other comparison algorithms, which are detailed in the following table.

|Algorithm| Description|
|----|----|
|TAF|D. J. Lee, S. McCauley, S. Singh, and M. Stein, “Telescoping Filter: A Practical Adaptive Filter,” in ESA, 2021. [[Implemention]](/nonlearnedfilter/telescoping-filter)|
|SSCF|M. Li, D. Chen, H. Dai, R. Xie, S. Luo, R. Gu, T. Yang, and G. Chen,“Seesaw counting filter: A dynamic filtering framework for vulnerable negative keys,” IEEE Transactions on Knowledge and Data Engineering, vol. 35, no. 12, pp. 12987–13001, 2023. [[Implemention]](/nonlearnedfilter/SSCF/SeesawCF.h)|
|Counting Bloom filter|L. Fan, P. Cao, J. Almeida, and A. Broder, “Summary cache: a scalable wide-area web cache sharing protocol,” IEEE/ACM Transactions on Networking, vol. 8, no. 3, pp. 281–293, 2000. [[Implementation]](/nonlearnedfilter/countingBloom.h)|
|WCBF|J. Bruck, J. Gao, and A. Jiang, “Weighted bloom filter,” in 2006 IEEE International Symposium on Information Theory, 2006. [[Implementation]](/nonlearnedfilter/wcbf.h)|
|SF|K. Deeds, B. Hentschel, and S. Idreos, “Stacked filters: learning to filter by structure,” Proc. VLDB Endow., vol. 14, no. 4, p. 600–612, 2020. [[Implementation]](/nonlearnedfilter/stackedfilter/)|

# Requirement 
   1. cmake@3+
   2. make
   3. recommended compiler: g++ or clang
# How to build
```bash
mkdir -p build && cd build && cmake .. && make
```
# How to run the benchmark
```bash
./experiment_semi_adaptive
./experiment_full_adaptive
```


