/*
 * Run a simulated annealing (SA) or Evolutionary Strategies (ES) optimization
 * on the generalized integer One-Max problem from Rothlauf's book:
 * "Representations for Genetic and Evolutionary Algorithms", 2nd ed., Sec. 5.4.2.
 *
 * In main(), you can control all simulation parameters, including which
 * representation to use, which fitness function, and SA/ES generations.
 * Prerequisite: Intel TBB library (libtbb-dev on debian distributions).
 * If you don't have TBB, use the commented-out loop in main().
 *
 * Compile with:
   g++ -Wall -Wextra -pedantic -O3 -march=native -std=c++17 onemax.cc  -ltbb -o onemax
 *
 * author: Eitan Frachtenberg
 */

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
using namespace tbb;

std::random_device randdev;
using bits_t = std::vector<bool>;

/*
 * An Organism lets you construct a random bit sequence, mutate, and compute
 * its fitness.
 */

class Organism {
 public:
  using fitness_fun_t = std::function<double (const bits_t&)>;

  // Construct a random sequence of 'len' bits, with a given fitness function
  // and mutation probability.
  Organism(size_t len, fitness_fun_t fit, double p_m)
  : bits_(len), fitness_(fit), p_m_(p_m)
  {
    std::default_random_engine reng(randdev());
    std::uniform_real_distribution dist(0., 1.);

    for (size_t i = 0; i < len; ++i) {
      if (dist(reng) < 0.5) {
        flip(i);
      }
    }
  }

  double fitness() const { return fitness_(bits_); }  // Compute fitness

  void flip(size_t idx) { bits_[idx] = bits_[idx] xor 1; }   // Flip a single bit

  // Mutate all bits with probability p_m_
  void mutate_all()
  {
    std::default_random_engine reng(randdev());
    std::uniform_real_distribution dist(0., 1.);

    for (size_t i = 0; i < bits_.size(); ++i) {
      if (dist(reng) < p_m_) {
        flip(i);
      }
    }
  }

  friend std::ostream& operator<<(std::ostream&, const Organism&);

 private:
  bits_t bits_;
  fitness_fun_t fitness_;
  double p_m_;
};

std::ostream&
operator<<(std::ostream& os, const Organism& o)
{
  for (auto b : o.bits_) {
    os << (b? "1" : "0");
  }
  return os;
}


/////////////////////////////////////////////////////////////////////////////
// A Sim class runs a single-organism simulated annealing or (1+1)-Es with a
// specified initial temperature, temperature adjustment factor, and a number
// of genotype in Organism units of a given length
class Sim {
 public:
  Sim(size_t units, size_t len, Organism::fitness_fun_t f, double p_m,
      double temp = 50, double t_adjust = 0.995)
  : genotype_(), temp_(temp), tadj_(t_adjust), eng_(randdev())
  , prob_dist_(0., 1.), org_dist_(0, units - 1), bit_dist_(0, len - 1)
  {
    for (size_t i = 0; i < units; ++i) {
      genotype_.push_back(Organism(len, f, p_m));
    }
  }

  // Run a single generation of simulated annealing: pick a random organism and flip a
  // random bit in it. If it improves fitness (or if it draws a "success" in
  // a random Boltzmann distribution), replace the organism with the new one.
  void SA_generation()
  {
    const auto org = org_dist_(eng_);
    const auto bit = bit_dist_(eng_);
    auto neworg = genotype_[org];

    const auto f0 = neworg.fitness();
    neworg.flip(bit);
    const auto f1 = neworg.fitness();

    if (f1 > f0 || prob_dist_(eng_) < exp((f1 - f0) / temp_)) {
      genotype_[org] = neworg;
    }

    temp_ *= tadj_;
  }

  // Run a single generation of (1+1)-ES: pick a random organism and flip a
  // random bit in it. If it improves fitness (or if it draws a "success" in
  // a random Boltzmann distribution), replace the organism with the new one.
  void ES_generation()
  {
    const auto org = org_dist_(eng_);
    auto neworg = genotype_[org];

    const auto f0 = neworg.fitness();
    neworg.mutate_all();
    const auto f1 = neworg.fitness();

    if (f1 > f0) {
       genotype_[org] = neworg;
    }
  }

  // Count how many organisms have optimal fitness
  unsigned num_optimal(double optimum) const
  {
    return std::count_if(genotype_.cbegin(), genotype_.cend(),
        [=](const Organism& o) { return o.fitness() == optimum; });
  }

  // Sum up individual organisms' fitness into one fitness:
  double fitness() const
  {
    return std::accumulate(genotype_.cbegin(), genotype_.cend(), 0,
        [](double sum, const Organism& o) { return sum + o.fitness(); });
  }

  friend std::ostream& operator<<(std::ostream&, const Sim&);

 private:
  std::vector<Organism> genotype_;
  double temp_;
  const double tadj_;
  std::default_random_engine eng_;
  std::uniform_real_distribution<double> prob_dist_;
  std::uniform_int_distribution<size_t> org_dist_, bit_dist_;
};

std::ostream&
operator<<(std::ostream& os, const Sim& sim)
{
  for (const auto& o : sim.genotype_) {
    os << "\t" << o << "\tFitness: " << sim.fitness();
  }
  return os;
}

/////////////////////////////////////////////////////////////////////////////
// Collection of representation encodings. A representation is
// a mapping from a bit vector (genotype) to an integer value (phenotype).

using phenotype_t = uint64_t;
using rep_t = std::function<phenotype_t (const bits_t&)>;

// Standard binary encoding: phenotype and genotype are identical:
phenotype_t
std_binary_rep(const bits_t& bits)
{
  phenotype_t ret = 0;
  int p = 1;
  for (auto bit = bits.crbegin(); bit != bits.crend(); bit++, p++) {
    ret += phenotype_t(*bit) << (p - 1);
  }
  assert(ret < (phenotype_t(1) << bits.size()));
  return ret;
}

// Binary-reflected gray encoding
phenotype_t
brg_rep(const bits_t& bits)
{
  auto it = bits.cbegin();
  phenotype_t ret = phenotype_t(*it);

  for (it++; it != bits.cend(); ++it) {
    const auto prev_bit = ret & 0x1;
    ret <<= 1;
    ret |= (*it)? prev_bit ^ 1 : prev_bit;
  }

  assert(ret < (phenotype_t(1) << bits.size()));
  return ret;
}

// Representation from an explicit mapping of bits to values.
// The mapping is given as an array, where the value in the n-th location
// is the mapped value from the n-th bitstring (using standrard binary ordering)
phenotype_t
explicit_rep(const bits_t& bits, const std::vector<phenotype_t>& mapping)
{
  const phenotype_t loc = std_binary_rep(bits);
  return mapping[loc];
}

// A few example mappings for len=3 bits (assume a=4)
const std::vector<phenotype_t> one_maxima =   { 5, 4, 1, 6, 7, 3, 0, 2 };
const std::vector<phenotype_t> two_maxima =   { 7, 2, 0, 5, 1, 6, 4, 3 };
const std::vector<phenotype_t> three_maxima = { 0, 5, 4, 7, 1, 3, 6, 2 };
const std::vector<phenotype_t> four_maxima =  { 5, 7, 6, 4, 1, 3, 2, 0 };
const std::vector<phenotype_t> different_four_maxima =  { 3, 7, 0, 2, 1, 4, 5, 6 };

// "Worst" representation for len=5, a=15
const std::vector<phenotype_t> five_worst =
  { 4, 30, 29, 13, 24, 8, 2, 18, 21, 15, 10, 25, 14, 31, 17, 1,
    28, 9, 3, 27, 7, 20, 16, 5, 0, 23, 26, 6, 19, 12, 11, 22 };


// Yet another "Worst" representation for len=5, a=15
const std::vector<phenotype_t> five_ubl =
  { 24, 1, 4, 19, 15, 16, 21, 13, 9, 26, 18, 0, 23, 12, 6, 22,
    3, 28, 20, 14, 30, 7, 5, 27, 29, 10, 8, 31, 2, 17, 25, 11 };

// "Non-greedy Gray encoding" for len=5:
const std::vector<phenotype_t> five_ngg =
  { 0, 1, 19, 2, 31, 28, 20, 3, 23, 26, 24, 25, 22, 27, 21, 4,
    13, 14, 18, 15, 30, 29, 17, 16,12, 9, 11, 10, 7, 8, 6, 5 };

/////////////////////////////////////////////////////////////////////////////
// Fitness functions for the one-max problem: given an 'a' value and a
// representation, compute the phenotypical value of the input bits given the
// representation, and calculate a linear scaling fitness that maximizes at
// the 'a' value.
double
onemax(const phenotype_t a, const rep_t& rep, const bits_t& bits)
{
  const auto phenotype = rep(bits);
  const auto maxfit = (1 << bits.size()) - 1;
  assert(double(a) <= maxfit);
  return maxfit - abs(double(phenotype) - a);
}

double
count_ones(const phenotype_t, const rep_t&, const bits_t& bits)
{
  return std::accumulate(bits.cbegin(), bits.cend(), 0);
}

/////////////////////////////////////////////////////////////////////////////
void usage()
{
  std::cerr << "Try running with the following integer arguments: a p g e\n";
  std::cerr << "a:\tThe binary value to strive to (default: 31)\n";
  std::cerr << "p:\tPopulation size, how many bitstrings are concatenated\n";
  std::cerr << "g:\tNumber of generations (fitness evaluations) to run for\n";
  std::cerr << "e:\tNumber of experiments to run concurrently\n";
}

/////////////////////////////////////////////////////////////////////////////
// Simulation main loop
// First, simulation parameters are chosen, including which representation
// to interpret the bit-string with.
// Algorithm: Loop over number of generations. In each generation, mutate
// each organism (there are 'experiments' of them), and decide whether to use
// the mutated offspring instead of the parent organism for the next gen.
// The decisions on how to mutate and when to replace a parent are based on
// the specific GEA chosen (SA / ES).
// Results are saved per generation, aggregated over all experiments, and
// reported on a generation-by-generation basis.
//
int main(int argc, char* argv[])
{
  const size_t len = 5;  // How many bits per organism?
  int a = (1 << len) - 1;   // Value to maximize to
  int popsize = 1;    // How many organisms per SA?
  unsigned generations = 2000; // How many fitness evaluations to run for?
  unsigned experiments = 100000;  // How many different SAs to average over?

  if (argc == 1) {
    usage();
  }
  if (argc > 1) {
    a = atoi(argv[1]);
  }
  if (argc > 2) {
    popsize = atoi(argv[2]);
  }
  if (argc > 3) {
    generations = atoi(argv[3]);
  }
  if (argc > 4) {
    experiments = atoi(argv[4]);
  }

  //const auto rep = [](const bits_t& bits){ return explicit_rep(bits, five_ngg); };
  //const auto rep = [](const bits_t& bits){ return explicit_rep(bits, five_ubl); };
  //const auto rep = [](const bits_t& bits){ return brg_rep(bits); };
  const auto rep = [](const bits_t& bits){ return std_binary_rep(bits); };

  const auto fit = [=](const bits_t& bits) { return onemax(a, rep, bits); };
  const auto maxfit = (1 << len) - 1;
  //const auto fit = [=](const bits_t& bits) { return count_ones(a, rep, bits); };
  //const auto maxfit = len;

  std::vector<Sim> sims;
  std::vector<unsigned> gens(experiments, generations);

  assert(len > 0);
  for(size_t i = 0; i < experiments; ++i) {
    sims.push_back(Sim(popsize, len, fit, 1. / len));
  }


  std::cout << "# Generation\tratio_optimal\tmean_fitness\n";

  // Main loops: generations and experiments
  for (unsigned g = 1; g <= generations; ++g) {
    std::atomic<unsigned> opt_count = 0;
    std::atomic<uint64_t> sum_fitness = 0;

    parallel_for(size_t(0), sims.size(), [&](size_t i) {
      opt_count += sims[i].num_optimal(maxfit);
      sum_fitness += sims[i].fitness();
      if (sims[i].fitness() == maxfit && g < gens[i]) {
        gens[i] = g;
      }
      sims[i].ES_generation();
    });

/* Sequential version of inner loop, if TBB is missing:
     for (auto& sim : sims) {
       opt_count += sim.num_optimal(maxfit);
      sum_fitness += sim.fitness();
       sim.generation();
     }
*/

    std::cout << g << "\t";
    std::cout << double(opt_count) / (experiments * popsize) << "\t";
    std::cout << double(sum_fitness) / (experiments * popsize) << "\n";
  }

  std::vector<unsigned> completed;
  std::copy_if(gens.cbegin(), gens.cend(), std::back_inserter(completed),
      [generations](unsigned g){ return g < generations; });

  std::cerr << "Mean generation to optimal solution: ";
  std::cerr << std::accumulate(completed.cbegin(), completed.cend(), 0) / double(completed.size());
  std::cerr << std::endl;
  return 0;
}
