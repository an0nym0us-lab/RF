# Reconfirm Filter

## About this Repo

This repository contains all the source code for RF (Reconfirm Filter) along with several comparison algorithms, as detailed in the following table:

| Algorithm              | Description |
|------------------------|-------------|
| **TAF**                | D. J. Lee, S. McCauley, S. Singh, and M. Stein, “Telescoping Filter: A Practical Adaptive Filter,” in ESA, 2021. [[Implementation]](/nonlearnedfilter/telescoping-filter) |
| **SSCF**               | M. Li, D. Chen, H. Dai, R. Xie, S. Luo, R. Gu, T. Yang, and G. Chen, “Seesaw counting filter: A dynamic filtering framework for vulnerable negative keys,” IEEE Transactions on Knowledge and Data Engineering, vol. 35, no. 12, pp. 12987–13001, 2023. [[Implementation]](/nonlearnedfilter/SSCF/SeesawCF.h) |
| **Counting Bloom filter** | L. Fan, P. Cao, J. Almeida, and A. Broder, “Summary cache: a scalable wide-area web cache sharing protocol,” IEEE/ACM Transactions on Networking, vol. 8, no. 3, pp. 281–293, 2000. [[Implementation]](/nonlearnedfilter/countingBloom.h) |
| **WCBF**               | J. Bruck, J. Gao, and A. Jiang, “Weighted bloom filter,” in 2006 IEEE International Symposium on Information Theory, 2006. [[Implementation]](/nonlearnedfilter/wcbf.h) |
| **SF**                 | K. Deeds, B. Hentschel, and S. Idreos, “Stacked filters: learning to filter by structure,” Proc. VLDB Endow., vol. 14, no. 4, p. 600–612, 2020. [[Implementation]](/nonlearnedfilter/stackedfilter/) |

---

## Dataset

The dataset used in our experiments is available at:  
**[https://github.com/an0nym0us-lab/RF/releases/tag/v1.0](https://github.com/an0nym0us-lab/RF/releases/tag/v1.0)**

> **Note:**  
> The release includes the large-scale dataset `dataset_caida`, which may be too large for direct inclusion in the repository. Please download it from the release page if needed.

---

## Project Structure Note

This version contains the previously-submodule projects (**cuckoofilter** and **telescoping-filter**) as regular directories, rather than git submodules. This change is made for ease of packaging and to facilitate anonymous distribution.

---

## Requirements

1. CMake >= 3.0
2. make
3. Recommended compiler: g++ or clang

---

## How to Build

```bash
mkdir -p build && cd build && cmake .. && make
```

---

## How to Run the Benchmark

```bash
./experiment_semi_adaptive
./experiment_full_adaptive
```

---

If you have any questions or issues, please feel free to open an issue on GitHub.