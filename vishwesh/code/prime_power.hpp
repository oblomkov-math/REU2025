#ifndef PRIME_POWER_HPP
#define PRIME_POWER_HPP

#include <cstdint>
#include <cmath>

/**
 * @brief Check whether a given integer n >= 2 is a prime power (i.e., n = p^k for some prime p and integer k >= 1).
 * @param n The integer to check.
 * @return true if n is a prime power, false otherwise.
 */
inline bool isPrimePower(std::uint64_t n) {
    if (n < 2) {
        return false;
    }

    // Handle factor 2 efficiently
    if ((n & 1) == 0) {
        while ((n & 1) == 0) {
            n >>= 1;
        }
        return (n == 1);
    }

    // Handle factor 3
    if (n % 3 == 0) {
        while (n % 3 == 0) {
            n /= 3;
        }
        return (n == 1);
    }

    // Check possible prime factors of form 6k +/- 1 up to sqrt(n)
    std::uint64_t limit = static_cast<std::uint64_t>(std::sqrt(static_cast<long double>(n)));
    for (std::uint64_t i = 5; i <= limit; i += 6) {
        // Check i
        if (n % i == 0) {
            std::uint64_t p = i;
            while (n % p == 0) {
                n /= p;
            }
            return (n == 1);
        }
        // Check i + 2
        std::uint64_t q = i + 2;
        if (n % q == 0) {
            while (n % q == 0) {
                n /= q;
            }
            return (n == 1);
        }
    }

    // If no small factors found, n is prime => n = prime^1
    return true;
}

#endif // PRIME_POWER_HPP

