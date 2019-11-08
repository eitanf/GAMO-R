/*
 * This program lets you run a simulated annealing (SA) optimization on the
 * integer One-Max problem from Rothlauf's book "Representations for Genetic
 * and Evolutionary Algorithms", 2nd ed., Sec. 5.4.2.
 */

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

std::random_device randdev;
using bits_t = std::vector<bool>;

/*
 * An Organism lets you construct a random bit sequence, mutate, and compute
 * its fitness.
 */

class Organism {
 public:
  using fitness_fun_t = std::function<double (const bits_t&)>;

  // Construct a random sequence of 'len' bits
  Organism(size_t len, fitness_fun_t fit)
  : bits_(len), fitness_(fit)
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

  friend std::ostream& operator<<(std::ostream&, const Organism&);

 private:
  bits_t bits_;
  fitness_fun_t fitness_;
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
// A SA class runs a single-organism simulated annealing simulation with a
// specified initial temperature, temperature adjustment factor, and a number
// of genotype in Organism units of a given length
class SA {
 public:
  SA(size_t units, size_t len, Organism::fitness_fun_t f, double temp = 50, double t_adjust = 0.995)
  : genotype_(), temp_(temp), tadj_(t_adjust), eng_(randdev())
  , prob_dist_(0., 1.), org_dist_(0, units-1), bit_dist_(0, len - 1)
  {
    for (size_t i = 0; i < units; ++i) {
      genotype_.push_back(Organism(len, f));
    }
  }

  // Run a single generation of simulation: pick a random organism and flip a
  // random bit in it. If it improves fitness (or if it draws a "success" in
  // a random Boltzmann distribution), replace the organism with the new one.
  void generation()
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

  friend std::ostream& operator<<(std::ostream&, const SA&);

 private:
  std::vector<Organism> genotype_;
  double temp_;
  const double tadj_;
  std::default_random_engine eng_;
  std::uniform_real_distribution<double> prob_dist_;
  std::uniform_int_distribution<size_t> org_dist_, bit_dist_;
};

std::ostream&
operator<<(std::ostream& os, const SA& sa)
{
  for (const auto& o : sa.genotype_) {
    os << "\t" << o << "\tFitness: " << sa.fitness();
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
  assert(ret < (1LL << bits.size()));
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

  assert(ret < (1LL << bits.size()));
  return ret;
}


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
  assert(a <= maxfit);
  return maxfit - abs(double(phenotype) - a);
}

/////////////////////////////////////////////////////////////////////////////
void usage()
{
  std::cerr << "Try running with the following integer arguments: a p g e\n";
  std::cerr << "a:\tThe binary value to strive to (default: 31)\n";
  std::cerr << "p:\tPopulation size, how many bitstrings are concatenated\n";
  std::cerr << "g:\tNumber of generations (fitness evaluations) to run for\n";
  std::cerr << "e:\tNumber of experiments to run concurrently\n";

int main(int argc, char* argv[])
{
  const size_t len = 5;  // How many bits per organism?
  int a = (1 << len) - 1;   // Value to maximize to
  int popsize = 1;    // How many organisms per SA?
  unsigned generations = 20; // How many fitness evaluations to run for?
  unsigned experiments = 1000;  // How many different SAs to average over?

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

  const auto fit = [a](const bits_t& bits) { return onemax(a, std_binary_rep, bits); };
  const auto maxfit = (1 << len) - 1;

  std::vector<SA> sims;
  std::vector<unsigned> opt_counts(1000);

  for(size_t i = 0; i < experiments; ++i) {
    sims.push_back(SA(popsize, len, fit));
  }


  std::cout << "# Generation  ratio_optimal\n";
  for (int g = 0; g < generations; ++g) {
    unsigned opt_count = 0;
    for (auto& sim : sims) {
      opt_count += sim.num_optimal(maxfit);
      sim.generation();
    }

    std::cout << g << "\t";
    std::cout << double(opt_count) / (experiments * popsize) << "\n";
  }
  return 0;
}
