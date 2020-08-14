#include <iostream>
#include <cstdlib>
#include <vector>

#include <jellyfish/jellyfish.hpp>
#include <gzstream.h>


typedef jellyfish::large_hash::unbounded_array<jellyfish::mer_dna> merhash;
typedef jellyfish::large_hash::region_iterator_base<merhash> hashiterator;

struct hashinfo {
  merhash hash;
  double  norm;
  // somehow we need to copy matrix
  hashinfo(size_t size, int value_len, jellyfish::RectangularBinaryMatrix matrix)
    : hash(size, 2 * jellyfish::mer_dna::k(), value_len, 100, std::move(matrix))
    , norm(0.0)
  { }
};

// Read every file in argv (files containing kmer/count pairs, assumed
// gziped), and store each in its own mer hash
std::vector<hashinfo> readKmerCounts(size_t size, int argc, char* argv[]) {
  std::vector<hashinfo>              res;
  jellyfish::RectangularBinaryMatrix matrix(jellyfish::ceilLog2(size), 2*jellyfish::mer_dna::k());
  matrix.randomize_pseudo_inverse();
  res.reserve(argc);
  for(int i = 0; i < argc; ++i)
    res.emplace_back(size, 6, matrix);

#pragma omp parallel for
  for(int i = 0; i < argc; ++i) {
    auto& hash = res[i].hash;
    igzstream in(argv[i]);
    if(!in.good()) {
      #pragma omp critical
      std::cerr << "Error openinig file '" << argv[i] << '\'' << std::endl;
      exit(1);
    }

    jellyfish::mer_dna mer;
    uint64_t           value;
    double             norm = 0.0;
    size_t             nb = 0;
    while(true) {
      in >> mer >> value;
      if(in.eof())
        break;
      if(!in.good()) {
        #pragma omp critical
        std::cerr << "Error reading file '" << argv[i] << '\'' << std::endl;
        exit(1);
      }
      if(!hash.add(mer, value)) {
        #pragma omp critical
        std::cerr << "Failed to insert kmer " << size << ' ' << mer << ' ' << value
                  << ' ' << nb << " from '" << argv[i] << "'\n";
        exit(1);
      }
      norm += value * value;
      ++nb;
    }
    res[i].norm = std::sqrt(norm);
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
  } else if(it1.key() > it2.key()) {
    return 1;
  } else {
    return 0;
  }
}

// Compute cosine-similarity from two mer hashes built using the same
// hash matrix (hence follow same k-mer order).
double computeSimilarity(const hashinfo& h1, const hashinfo& h2) {
  hashiterator it1(&h1.hash, 0, h1.hash.size());
  hashiterator it2(&h2.hash, 0, h2.hash.size());

  double product = 0.0;
  bool   done    = !it1.next() || !it2.next(); // Pull first element

  while(!done) {
    std::cerr << "> " << it1.key() << ' ' << it2.key() << ' ' << it1.val() << ' ' << it2.val() << ' ' << it1.id() << ':' << it1.oid() << ' ' << it2.id() << ':' << it2.oid() << '\n';
    switch(comp(it1, it2)) {
    case 0:
      product += it1.val() * it2.val();
      std::cerr << it1.key() << ' ' << it1.val() << ' ' << it2.val() << ' ' << it1.id() << ':' << it1.oid() << ' ' << it2.id() << ':' << it2.oid() << '\n';
      done = !it1.next() || !it2.next();
      break;
    case -1:
      // std::cerr << product << ' ' <<  it1.pos() << ' ' << it2.pos() << ' ' << it1.key() << " < " << it2.key() << std::endl;
      std::cerr << it1.key() << ' ' << it1.val() << ' ' << 0 << ' ' << it1.id() << ':' << it1.oid() << ' ' << "-:-" << '\n';;
      done = !it1.next();
      break;
    case 1:
      // std::cerr << product << ' ' <<  it1.pos() << ' ' << it2.pos() << ' ' << it1.key() << " > " << it2.key() << std::endl;
      std::cerr << it2.key() << ' ' << 0 << ' ' << it2.val() << ' ' << "-:-" << ' ' << it2.id() << ':' << it2.oid() << '\n';
      done = !it2.next();
      break;
    }
  }
  // #pragma omp critical
  //   std::cerr << product << ' ' << h1.norm << ' ' << h2.norm << '\n';

  return product / (h1.norm * h2.norm);
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
  std::vector<hashinfo> counts = readKmerCounts(size, argc - 3, argv + 3);

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
