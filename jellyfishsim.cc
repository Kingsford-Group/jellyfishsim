#include <iostream>
#include <cstdlib>
#include <vector>

#include <jellyfish/jellyfish.hpp>
#include <gzstream.h>


typedef jellyfish::large_hash::unbounded_array<jellyfish::mer_dna> merhash;

// Read every file in argv (files containing kmer/count pairs, assumed
// gziped), and store each in its own mer hash
std::vector<merhash> readKmerCounts(size_t size, int argc, char* argv[]) {
  std::vector<merhash> res;
  jellyfish::RectangularBinaryMatrix matrix(jellyfish::ceilLog2(size), jellyfish::mer_dna::k());

  matrix.randomize_pseudo_inverse();
  res.reserve(argc);
  for(int i = 0; i < argc; ++i) {
    res.emplace_back(size, 2 * jellyfish::mer_dna::k(), 6, 100,
                     jellyfish::RectangularBinaryMatrix(matrix)); // somehow we need to copy matrix
  }

#pragma omp parallel for
  for(int i = 0; i < argc; ++i) {
    auto& hash = res[i];
    igzstream in(argv[i]);
    if(!in.good()) {
      std::cerr << "Error openinig file '" << argv[i] << '\'' << std::endl;
      exit(1);
    }

    jellyfish::mer_dna mer;
    uint64_t           value;
    while(true) {
      in >> mer >> value;
      if(in.eof())
        break;
      if(!in.good()) {
        std::cerr << "Error reading file '" << argv[i] << '\'' << std::endl;
        exit(1);
      }
      hash.add(mer, value);
    }
  }

  return res;
}

// Compare elements in hash point to by iterator
template<typename IT>
int comp(IT& it1, IT& it2) {
  size_t pos1 = it1.pos(), pos2 = it2.pos();

  if(pos1 < pos2) {
    return -1;
  } else if(pos1 > pos2) {
    return 1;
  } else if(it1.key() < it2.key()) {
    return -1;
  } else if(it1.key() < it2.key()) {
    return 1;
  } else {
    return 0;
  }
}

// Compute cosine-similarity from two mer hashes built using the same
// hash matrix (hence follow same k-mer order).
double computeSimilarity(const merhash& h1, const merhash& h2) {
  auto it1 = h1.begin(), it2 = h2.begin();
  const auto end1 = h1.end(), end2 = h2.end();

  double norm1 = 0.0, norm2 = 0.0, product = 0.0;
  while(it1 != end1 && it2 != end2) {
    switch(comp(it1, it2)) {
    case 0:
      product += it1.val() * it2.val();
      norm1   += it1.val() * it1.val();
      norm2   += it2.val() * it2.val();
      ++it1;
      ++it2;
      break;
    case -1:
      norm1 += it1.val() * it1.val();
      ++it1;
      break;
    case 1:
      norm2 += it2.val() * it2.val();
      ++it2;
      break;
    }
  }

  return product / (std::sqrt(norm1 * norm2));
}

int main(int argc, char *argv[]) {
  if(argc < 4) {
    std::cerr << "Usage: " << argv[0] << " klen size file.gz..." << std::endl;
    exit(1);
  }

  const int     klen = std::atoi(argv[1]);
  const ssize_t size = std::atol(argv[2]);
  if(klen <= 0) {
    std::cerr << "Invalid k-mer length '" << klen << '\'' << std::endl;
    exit(1);
  }
  if(size <= 0) {
    std::cerr << "Invalid size '" << size << '\'' << std::endl;
    exit(1);
  }

  jellyfish::mer_dna::k(klen); // Set k-mer length for Jellyfish
  std::vector<merhash> counts = readKmerCounts(size, argc - 3, argv + 3);

  std::vector<std::vector<double>> matrix(counts.size());
  for(size_t i = 0; i < counts.size(); ++i) {
    matrix[i].resize(counts.size());
    matrix[i][i] = 1.0;
  }

  for(size_t i = 0; i < counts.size() - 1; ++i) {
    #pragma omp parallel for
    for(size_t j = i + 1; j < counts.size(); ++j) {
      const auto sim = computeSimilarity(counts[i], counts[j]);
      matrix[i][j] = matrix[j][i] = sim;
    }
  }

  for(auto& row : matrix) {
    for(double elt : row)
      std::cout << elt << ' ';
    std::cout << '\n';
  }

  return 0;
}
