/*
 * This tool lets you compute the locality of a bit-to-integer mapping,
 * as defined by Rothlauf's "Representations for Genetic and Evolutionary
 * Algorithms", 2nd ed., p. 77, eq. 3.23.
 */

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <numeric>

// In this implementation, all bit strings are represented as a simple
// integer (with the typical binary representation), up to 256 bits.
using bits_t = uint8_t;
using num_t = uint64_t;  // A natural number, what we're representing.

constexpr bits_t BITS = 3;   // How many bits long is the representation?
constexpr num_t N = 1ULL << BITS;   // No. of distinct values in representation.

// A representation is simply a mapping from bit string to an integer value
// in the range 0..N-1. There are N such values, so we map each bit string
// (in the normal binary enumeration order) to a value.
using rep_t = std::array<num_t, N>;

// A list of all neighbors (single-bit mutations) for a given bitstring
using neighbors_t = std::array<bits_t, BITS>;


///////////////////
// Helper functions:

// Check whether a set of values is actually a representation, namely,
// a permutation of the values 0..N-1.
bool
is_representation(const rep_t& rep)
{
  rep_t reference;
  std::iota(reference.begin(), reference.end(), 0);
  return std::is_permutation(rep.begin(), rep.end(), reference.begin());
}

// Return an array of all the bit strings that are single-bit neighbors
// of a given bit string.
neighbors_t
bit_neighbors(bits_t bits)
{
  neighbors_t ret;
  for (bits_t i = 0; i < BITS; i++) {
    ret[i] = bits ^ (1 << i);
  }
  return ret;
}


// Main utility function: for a given representation, compute its locality.
uint64_t
locality(const rep_t& rep)
{
  assert(is_representation(rep));

  // Compute a mapping from any bitstring to its neighbors:
  std::array<neighbors_t, N> n_map;
  for (num_t i = 0; i < N; ++i) {
    n_map[i] = bit_neighbors(i);
  }

  // Now, we can sum up the integer value differences between all neighbors
  uint64_t sum = 0;

  for (num_t i = 0; i < N; ++i) { // Loop over all bit strings
    const auto phenotype1 = rep[i];
    for (bits_t n = 0; n < BITS; ++n) { // Loop over all neighbors of string
      const auto phenotype2 = rep[n_map[i][n]];
      if (phenotype1 > phenotype2) {
        sum += phenotype1 - phenotype2 - 1;
      }
    }
  }
  // We only sum up half the cases, because of symmetry, so the total
  // sum must be doubled.
  return sum * 2;
}

int main()
{
  auto x = bit_neighbors(6);
  rep_t bin = { 0, 1, 2, 3, 4, 5, 6, 7 };
  rep_t brg = { 0, 1, 3, 2, 7, 6, 4, 5 };
  rep_t ngg = { 0, 7, 1, 2, 5, 6, 4, 3 };
  rep_t worst = { 0, 5, 6, 3, 7, 1, 2, 4 };

  std::cout << "locality of binary: " << locality(bin) << "\n";
  std::cout << "locality of binary reflected gray: " << locality(brg) << "\n";
  std::cout << "locality of non-greedy gray: " << locality(ngg) << "\n";
  std::cout << "locality of worst: " << locality(worst) << "\n";
  return 0;
}

