#include <iostream>
#include <cstdlib>
#include <vector>

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

std::vector<hashinfo> readKmerCounts(int argc, char* argv[]) {
  std::vector<hashinfo> res(argc);

#pragma omp parallel for
  for(int i = 0; i < argc; ++i) {
    auto& mers = res[i].mers;
    igzstream in(argv[i]);
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

// Position in a vector representing a lower triangular matrix (with
// diagonal). Pre-condition: i >= j
inline size_t triangle(size_t i, size_t j = 0) {
  return i * (i + 1) / 2 + j;
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
  if(argc < 3) {
    std::cerr << "Usage: " << argv[0] << " klen file.gz..." << std::endl;
    exit(1);
  }

  const int     klen = std::atoi(argv[1]);
  if(klen <= 0) {
    std::cerr << "Invalid k-mer length '" << klen << '\'' << std::endl;
    exit(1);
  }

  jellyfish::mer_dna::k(klen); // Set k-mer length for Jellyfish
  std::vector<hashinfo> counts = readKmerCounts(argc - 2, argv + 2);

  std::vector<double> matrix(triangle(counts.size() + 1));
  matrix[0] = 1.0;
  for(size_t i = 1; i < counts.size(); ++i) {
    matrix[triangle(i, i)] = 1.0;
    #pragma omp parallel for
    for(size_t j = 0; j < i; ++j) {
      const double sim = computeSimilarity(counts[i], counts[j]);
      matrix[triangle(i, j)] = sim;
    }
  }

  auto elt = matrix.cbegin();
  for(size_t i = 0; i < counts.size(); ++i) {
    for(size_t j = 0; j <= i; ++j, ++elt)
      std::cout << *elt << ' ';
    std::cout << '\n';
  }

  return 0;
}
