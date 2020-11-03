#include <iostream>
#include <cstdlib>
#include <vector>

#include <fstream>
#include <iterator>
#include <iomanip>
#include <string>
#include <algorithm>

#include <jellyfish/jellyfish.hpp>
#include <gzstream.h>

using jellyfish::mer_dna;

struct mercount {
  uint64_t mer;
  uint64_t count;
};
bool operator<(const mercount& m1, const mercount& m2) {
  return m1.mer < m2.mer;
}
struct hashinfo {
  std::vector<mercount> mers;
  double                norm;
  hashinfo() : norm(0.0) {
    mers.reserve(100000);
  }
};

std::vector<hashinfo> readKmerCounts(std::vector<std::string> argv) {
  int argc = argv.size();
  std::vector<hashinfo> res(argc);

#pragma omp parallel for
  for(int i = 0; i < argc; ++i) {
    auto& mers = res[i].mers;
    igzstream in(argv[i].c_str());
    if(!in.good()) {
      #pragma omp critical
      std::cerr << "Error openinig file '" << argv[i] << '\'' << std::endl;
      exit(1);
    }

    mer_dna  mer;
    uint64_t value;
    double   norm = 0.0;
    while(true) {
      in >> mer >> value;
      if(in.eof())
        break;
      if(!in.good()) {
        #pragma omp critical
        std::cerr << "Error reading file '" << argv[i] << '\'' << std::endl;
        exit(1);
      }
      mers.push_back({ mer.get_bits(0, 2*mer_dna::k()), value });
      norm += value * value;
    }

    std::sort(mers.begin(), mers.end());
    res[i].norm = std::sqrt(norm);
  }

  return res;
}

std::vector<hashinfo> readOneKmerCounts(std::string argv) {
  std::vector<hashinfo> res(1);

  auto& mers = res[0].mers;
  igzstream in(argv.c_str());
  if(!in.good()) {
    std::cerr << "Error openinig file '" << argv << '\'' << std::endl;
    exit(1);
  }

  mer_dna  mer;
  uint64_t value;
  double   norm = 0.0;
  while(true) {
    in >> mer >> value;
    if(in.eof())
      break;
    if(!in.good()) {
      std::cerr << "Error reading file '" << argv << '\'' << std::endl;
      exit(1);
    }
    mers.push_back({ mer.get_bits(0, 2*mer_dna::k()), value });
    norm += value * value;
  }

  std::sort(mers.begin(), mers.end());
  res[0].norm = std::sqrt(norm);

  return res;
}

std::vector<std::string> readDatasetsFile(char *filename) {
  std::vector<std::string> list_datasets;

  std::ifstream datasets_file(filename);
  if (!datasets_file.is_open()) {
    std::cerr << "Error openinig file '" << filename << '\'' << std::endl;
    exit(1);
  }

  std::string dataset;
  while (datasets_file >> dataset) {
    list_datasets.push_back(dataset);
  }
  datasets_file.close();

  return list_datasets;
}

double computeSimilarity(const hashinfo& mers1, const hashinfo& mers2) {
  auto it1 = mers1.mers.begin(), it2 = mers2.mers.begin();
  const auto end1 = mers1.mers.end(), end2 = mers2.mers.end();

  double product = 0.0;
  while(it1 != end1 && it2 != end2) {
    if(*it1 < *it2) {
      ++it1;
    } else if(*it2 < *it1) {
      ++it2;
    } else {
      product += it1->count * it2->count;
      ++it1;
      ++it2;
    }
  }

  return product / (mers1.norm * mers2.norm);
}

int main(int argc, char *argv[]) {
  if(argc < 5) {
    std::cerr << "Usage: " << argv[0] << " klen full_set_datasets_file rep_set_datasets_file q" << std::endl;
    exit(1);
  }

  const int     klen = std::atoi(argv[1]);
  if(klen <= 0) {
    std::cerr << "Invalid k-mer length '" << klen << '\'' << std::endl;
    exit(1);
  }

  std::vector<std::string> full_set_datasets = readDatasetsFile(argv[2]);
  std::vector<std::string> rep_set_datasets = readDatasetsFile(argv[3]);

  const double q = std::atof(argv[4]);
  const int K = (int) ((1.0-q)*full_set_datasets.size());

  jellyfish::mer_dna::k(klen); // Set k-mer length for Jellyfish
  std::vector<hashinfo> counts = readKmerCounts(rep_set_datasets);

  std::vector<double> min_distances(full_set_datasets.size());

#pragma omp parallel for 
  for(size_t k = 0; k < min_distances.size(); ++k) {
    std::vector<hashinfo> one_count = readOneKmerCounts(full_set_datasets[k]);
    std::vector<double> distances(counts.size());
    for(size_t i = 0; i < counts.size(); ++i) {
      distances[i] = 1.0 - computeSimilarity(one_count[0], counts[i]);
    }
    min_distances[k] = *std::min_element(distances.begin(), distances.end());
  }

  std::sort(min_distances.begin(), min_distances.end());
  double d_H = min_distances[min_distances.size()-1];
  double d_HK = min_distances[K-1];

  std::cout << "K = " << K << std::endl;
  std::cout << std::fixed << std::setprecision(9);
  std::cout << "Hausdorff distance: d_H = " << d_H << std::endl;
  std::cout << "Partial Hausdorff distance (with q = " << q << "): d_HK = " << d_HK << std::endl;

  return 0;
}
