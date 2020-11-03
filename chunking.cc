#include <iostream>
#include <cstdlib>
#include <vector>

#include <fstream>
#include <iterator>
#include <iomanip>
#include <string>
#include <algorithm>
#include <experimental/random>

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

int gen_rand_int(int range) {
  std::experimental::reseed(19951228);
  int random_number = std::experimental::randint(0, range-1);

  return random_number;
}

int arg_max(std::vector<double> v) {
  int argmax = std::distance(v.begin(), std::max_element(v.begin(), v.end()));

  return argmax;
}

int arg_min(std::vector<double> v) {
  int argmin = std::distance(v.begin(), std::min_element(v.begin(), v.end()));

  return argmin;
}

std::vector<double> pairwise_min(std::vector<double> v1, std::vector<double> v2) {
  std::vector<double> result;

  std::transform(v1.begin(), v1.end(), v2.begin(), std::back_inserter(result), [](double a, double b) {return std::min(a, b);});

  return result;
}

std::vector<double> slicing_vector(std::vector<double> in_vector, std::vector<int> s_indices) {
  std::vector<double> result;

  std::transform(s_indices.begin(), s_indices.end(), std::back_inserter(result), [&in_vector](int idx) {return in_vector[idx];});

  return result;
}

std::vector<int> getUnfullChunks(int chunk_size, std::vector<std::vector<std::string>> list_chunks) {
  std::vector<int> result;
  int idx = 0;

  std::for_each(list_chunks.begin(), list_chunks.end(), [&](std::vector<std::string> chunk) {if (chunk.size() < chunk_size) result.push_back(idx); ++idx;});

  return result;
}

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

std::vector<double> getDistancesToSeed(int seed, std::vector<hashinfo> rand_set_counts) {
  std::vector<double> dists_to_seed(rand_set_counts.size());

#pragma omp parallel for
  for(size_t k = 0; k < rand_set_counts.size(); ++k) {
    dists_to_seed[k] = 1.0 - computeSimilarity(rand_set_counts[seed], rand_set_counts[k]);
  }

  return dists_to_seed;
}

int main(int argc, char *argv[]) {
  if(argc < 6) {
    std::cerr << "Usage: " << argv[0] << " klen full_set_datasets_file rand_set_datasets_file num_chunks chunk_size" << std::endl;
    exit(1);
  }

  const int     klen = std::atoi(argv[1]);
  if(klen <= 0) {
    std::cerr << "Invalid k-mer length '" << klen << '\'' << std::endl;
    exit(1);
  }

  std::vector<std::string> full_set_datasets = readDatasetsFile(argv[2]);
  std::vector<std::string> rand_set_datasets = readDatasetsFile(argv[3]);

  const int K = std::atoi(argv[4]);
  const int chunk_size = std::atoi(argv[5]);

  jellyfish::mer_dna::k(klen); // Set k-mer length for Jellyfish
  std::vector<hashinfo> rand_set_counts = readKmerCounts(rand_set_datasets);

  std::vector<int> seeds(K);
  std::vector<std::string> seeds_datasets(K);
  std::vector<std::vector<std::string>> list_chunks(K, std::vector<std::string>());

  auto pend = full_set_datasets.end();

  // Farthest point sampling for selecting seeds of chunks, using rand set
  seeds[0] = gen_rand_int(rand_set_datasets.size());
  std::vector<double> crt_dists_to_seeds = getDistancesToSeed(seeds[0], rand_set_counts);
  seeds_datasets[0] = rand_set_datasets[seeds[0]];
  pend = std::remove(full_set_datasets.begin(), pend, seeds_datasets[0]);
  list_chunks[0].push_back(seeds_datasets[0]);

  for(int i = 1; i < K; ++i) {
    seeds[i] = arg_max(crt_dists_to_seeds);
    crt_dists_to_seeds = pairwise_min(crt_dists_to_seeds, getDistancesToSeed(seeds[i], rand_set_counts));
    seeds_datasets[i] =  rand_set_datasets[seeds[i]];
    pend = std::remove(full_set_datasets.begin(), pend, seeds_datasets[i]);
    list_chunks[i].push_back(seeds_datasets[i]);
  }

  // Create full-set datasets vs. seeds distances 2D array
  std::vector<std::vector<double> > distances_2d_array( full_set_datasets.size()-K, std::vector<double> (K));

#pragma omp parallel for
  for(size_t k = 0; k < full_set_datasets.size()-K; ++k) {
    std::vector<hashinfo> one_count = readOneKmerCounts(full_set_datasets[k]);
    std::vector<double> distances(K);
    for(size_t i = 0; i < K; ++i) {
      distances[i] = 1.0 - computeSimilarity(one_count[0], rand_set_counts[seeds[i]]);
    }
    distances_2d_array[k] = distances;
  }

  // Assign full-set datasets to its closest seed
  std::vector<int> unfull_chunks;
  std::vector<double> crt_row;
  int seed_idx;
  for(size_t i = 0; i < full_set_datasets.size()-K; ++i) {
    unfull_chunks = getUnfullChunks(chunk_size, list_chunks);
    crt_row = slicing_vector(distances_2d_array[i], unfull_chunks);
    seed_idx = arg_min(crt_row);
    list_chunks[unfull_chunks[seed_idx]].push_back(full_set_datasets[i]);
  }

#pragma omp parallel for
  for(int i = 0; i < K; ++i) {
    std::ofstream chunk_file("fullset_datasets_chunk_" + std::to_string(i));
    std::ostream_iterator<std::string> chunk_output_iterator(chunk_file, "\n");
    std::copy(list_chunks[i].begin(), list_chunks[i].end(), chunk_output_iterator);
  }

  std::cout << "Number of seeds = " << K << std::endl;
  std::cout << "Seeds indices in the rand set:" << std::endl;
  for (auto sd : seeds) {
    std::cout << sd << ", ";
  }
  std::cout << std::endl;

  /*
  std::ofstream seeds_file("seeds_datasets");
  std::ostream_iterator<std::string> seeds_output_iterator(seeds_file, "\n");
  std::copy(seeds_datasets.begin(), seeds_datasets.end(), seeds_output_iterator);

  std::ofstream fullset_noseeds_file("fullset_datasets_no_seeds");
  std::ostream_iterator<std::string> fullset_noseeds_output_iterator(fullset_noseeds_file, "\n");
  std::copy(full_set_datasets.begin(), pend, fullset_noseeds_output_iterator);

  std::ofstream output_file("distances_2d_array");
  output_file << std::fixed << std::setprecision(9);
  std::ostream_iterator<double> output_iterator(output_file, " ");
  for(const auto& vt : distances_2d_array) {
    std::copy(vt.cbegin(), vt.cend(), output_iterator);
    output_file << '\n';
  }
  */

  return 0;
}
