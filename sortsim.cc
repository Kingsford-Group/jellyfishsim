#include <iostream>
#include <cstdlib>
#include <vector>

#include <fstream>
#include <iterator>
#include <iomanip>
#include <string>

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

// Position in a vector representing a lower triangular matrix (with
// diagonal). Pre-condition: i >= j
inline size_t triangle(size_t i, size_t j = 0) {
  return i * (i + 1) / 2 + j;
}

inline void invTriangle(size_t index, size_t& i, size_t& j) {
  i = std::floor((std::sqrt(8*index + 1) - 1.0) / 2);
  j = index - triangle(i);
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
    std::cerr << "Usage: " << argv[0] << " klen list_datasets_file" << std::endl;
    exit(1);
  }

  const int     klen = std::atoi(argv[1]);
  if(klen <= 0) {
    std::cerr << "Invalid k-mer length '" << klen << '\'' << std::endl;
    exit(1);
  }

  std::vector<std::string> list_datasets = readDatasetsFile(argv[2]);

  jellyfish::mer_dna::k(klen); // Set k-mer length for Jellyfish
  std::vector<hashinfo> counts = readKmerCounts(list_datasets);

  std::vector<double> matrix(triangle(counts.size()));

#pragma omp parallel for
  for(size_t k = 0; k < matrix.size(); ++k) {
    size_t i, j;
    invTriangle(k, i, j);
    if(i != j)
      matrix[k] = computeSimilarity(counts[i], counts[j]);
    else
      matrix[k] = 1.0;
  }

  std::vector<std::vector<double> > sim_matrix( counts.size(), std::vector<double> (counts.size()));

  auto elt = matrix.cbegin();
  for(size_t i = 0; i < counts.size(); ++i) {
    for(size_t j = 0; j <= i; ++j, ++elt) {
      sim_matrix[i][j] = *elt;
      sim_matrix[j][i] = *elt;
    }
  }

  std::ofstream output_file("similarity_matrix");
  output_file << std::fixed << std::setprecision(9);
  std::ostream_iterator<double> output_iterator(output_file, " ");
  for(const auto& vt : sim_matrix) {
    std::copy(vt.cbegin(), vt.cend(), output_iterator);
    output_file << '\n';
  }

  return 0;
}
