# mrtab

An iteration threshold table generator for the Miller-Rabin random probable prime search.
#

In a random, probable prime search, the Miller-Rabin strong probable prime (or *compositeness*) test requires a certain number of iterations to achieve an acceptably low probability of failure - that is, the probability that a composite will be incorrectly declared prime. The **mrtab** utility generates a **k-bit** threshold table for the number of iterations required to achieve a specified upper-bound for the error probability - denoted in the literature as **p(k, t)**.

These concepts are central to the generation of *industrial strength primes*, as used in public key cryptography. A minimal iteration count, based on number-theoretic results, is of particular interest for efficient key generation in constrained environments, such as low-power devices and embedded systems.
