# unitary-construction-and-decomposition

The code 'Unitary_generation.m' has an in-built Quantum Inference Protocol that takes in a bit-string and builds quantum memory states based on L-length past (L is also an input parameter). In order to build a unitary operator that is of acceptable and implementable dimensions, we opt to merge quantum memory states if the conditional futures of the states are within delta-similarity. The parameter delta has a default value but one has the option to implement any value of delta.

Once merged, generating the unitary operator can begin. The methodology follows that in the paper Phys. Rev. Lett. 120, 240502 (2018).

In order to implement the unitary operator on quantum devices, one has to decompose the unitary operator into implementable elementary quantum gates such as controlled rotations, NOTs, phase gates, etc. Our unitary operator is numerically not unitary because multiplying its Hermitian conjugate only gives the identity matrix after accounting for rounding errors; this is a signature of having a finite bit-string that affects the probability amplitudes of the quantum memory states. As such, we write an algorithm called 'csd_gsvd.m' that decomposes a unitary operator based on the cosine-sine decomposition (CSD) and the generalised singular value decomposition (GSVD) for our unitary operator.

We provide an example of how to use the code(s) is given in the file 'Script_How_to_use_unitary_construction_decomposition.m'.

This set of codes is made up of two parts with an example:
(1) Constructing the unitary operator.
(2) Decomposing the unitary operator.
(3) Script with an example.

The following article shall be cited should this code be used: http://arxiv.org/abs/2105.06448.










