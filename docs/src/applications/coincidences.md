```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```
# Coincidence rates

## Outline

- [Introduction](#introduction)
- [Two-photon coincidence rate with delay](#two-photon-coincidence-rate-with-delay)
- [Three photons with one delayed input](#three-photons-with-one-delayed-input)

## Introduction

This page studies coincidence rates for photons in linear optical networks. The discussion below introduces the experimental setting, the overlap parameter that captures delay-induced distinguishability, and the shorthand used throughout the examples.

The basic experiment is as follows:
1. A device produces a given number of photons in specific spatial modes on demand, or at least *heralds the preparation of such a state* (for instance, with multiplexed SPDC), so that the photons occupy well-defined time windows.
2. The photons can be delayed relative to one another (e.g. by having a heralding photon traverse a longer distance).
3. We pass the photons through a linear optical system and measure coincident clicks at specific output ports.

The relative delay controls photon distinguishability and thereby affects the observed coincidence rates. We model this dependence with a frequency-varying photon profile, which corresponds roughly to a time profile after a Fourier transform. Each spatial mode $k$ contains a continuum of frequency modes $\omega$, with annihilation operators labelled $\hat a_k(\omega)$. An operator that creates a photon with a Gaussian profile in mode $k$, centered at time $\tau$, is

```math
\hat{A}^\dagger_k(\tau) \coloneqq
\int d\omega\,e^{i\omega\tau}\varphi(\omega)\hat{a}^\dagger_k(\omega),
```
As a simple concrete example, $\hat{A}^\dagger$ creates a photon with a Gaussian spectrum when $\varphi$ has the following frequency profile:

```math
|\varphi(\omega)|^2 \coloneqq \frac{\exp({-(\omega-\omega_0)^2/(2\sigma^2)})}{\sqrt{2\pi}\sigma},
```

For two photons with the same spectral profile but centered at different times $\tau, \tau'$, the wavefunction overlap $\langle 0\vert \hat A (\tau) \hat A^\dagger(\tau') \vert 0\rangle$ is

```math
\zeta \coloneqq \exp({-\frac12\sigma^2(\tau-\tau')^2}).
```

The next two sections calculate coincidence rates as a function of this overlap. For convenience, we define the squared overlap

$$z \coloneqq \vert \zeta\vert^2.$$

To track which photon goes where, we introduce the shorthand

```math
R(\text{inputs} \to \text{outputs}, z)
```

Here, $R$ denotes the coincidence probability for fixed input and output occupation patterns. One photon, specified in the text, is delayed and has squared overlap $z$ with the other photons.

For example, $R(23 \to 13; z)$ means that one photon enters input 2, one enters input 3, and one is detected at each of outputs 1 and 3. Similarly, $R(234 \to \alpha\alpha4; z)$ means that one photon enters each of inputs 2, 3, and 4, while two photons are detected at output $\alpha$ and one at output 4.

Repeated output labels encode multiplicity. Thus, $R(234\to pp4;z)$ means that we observe two photons at port $p$.

## Two-photon coincidence rate with delay

For two photons entering input modes $k,l$ and detected at outputs $m,n$, the coincidence rate at squared overlap $z$ is

```math
R(kl \to mn;z)
=\frac{1}{2}(1+z)\left|\mathrm{Per}(U_{kl\to mn})\right|^2
+\frac{1}{2}(1-z)\left|\mathrm{Det}(U_{kl\to mn})\right|^2,
```

where $U_{kl \to mn}$ is constructed from $U$ by selecting the specified columns and rows (see [the background page](../background/group_functions.md#Bosons:-the-permanent)).

If the delay is zero, the photons are completely indistinguishable ($z=1$), and the result reduces to a standard permanent:

```math
\tau=0 \Rightarrow R=\left|\mathrm{Per}(U_{kl\to mn})\right|^2,
```

With infinite delay, the photons are completely distinguishable ($z=0$), and the rate is an equal mixture of the permanent and determinant:

```math
\tau\to\infty \Rightarrow
R=\frac{1}{2}\left|\mathrm{Per}(U_{kl\to mn})\right|^2
+\frac{1}{2}\left|\mathrm{Det}(U_{kl\to mn})\right|^2.
```


We can verify these equations with our library. Tracking a continuum of frequency modes is infeasible, but only the overlap matters. We can therefore simulate the behavior by introducing internal states in addition to the spatial modes.

We use two spatial modes to illustrate the method. The calculations use a basis in which the occupation numbers $(a,a',b,b')$ have the following meanings: $a$ counts "on-time" photons in spatial mode 1, $a'$ counts "late" (out-of-time) photons in spatial mode 1, $b$ counts "on-time" photons in spatial mode 2, and $b'$ counts "late" photons in spatial mode 2.

In this extended space, the linear optical transformation on the spatial modes is extended by the identity on the internal modes because the transformations are constant on coherence timescales:

```math
V \coloneqq U_{\text{spatial}} \otimes \operatorname{id}_{\text{internal}}
```

We therefore extend a beamsplitter and define a group function in the extended space:
```@repl coincidences2
using GroupFunctions
using LinearAlgebra: kron
two_photon_spatial_bs = su2_block(2, 1, (1.5, 0.9, 0.2)); # 2×2 spatial beam splitter
two_photon_extended_u = kron(two_photon_spatial_bs, [1.0 0.0; 0.0 1.0]);

two_photon_group_vals, two_photon_occ_pats = group_function([2, 0, 0, 0], two_photon_extended_u);
```

We next define the initial state:


```@repl coincidences2
overlap_amp = 0.37; # amplitude overlap

two_photon_overlap_sq = overlap_amp^2; # squared overlap
#initial state is spanned by these two:
input_occ_on_time = [1, 0, 1, 0]; # one photon in sp. mode 1, one in sp. mode 2, both "on time"
input_occ_delayed = [1, 0, 0, 1]; # one photon "on time" in sp. mode 1, one "late" in 2
two_photon_occ_idx(occupation) = findfirst(p -> occupation_number(p) == occupation, two_photon_occ_pats);
two_photon_input_vec = zeros(ComplexF64, length(two_photon_occ_pats));
two_photon_input_vec[two_photon_occ_idx(input_occ_on_time)] = overlap_amp;
two_photon_input_vec[two_photon_occ_idx(input_occ_delayed)] = sqrt(1 - overlap_amp^2);

```

We calculate the coincidence rate by summing over four possible output states. Each contains one photon in spatial mode 1 and one in spatial mode 2, but the states have different internal ("on-time"/"late") indices:


```@repl coincidences2
output_occupation_patterns = [
    [1,0,1,0],
    [1,0,0,1],
    [0,1,1,0],
    [0,1,0,1]
];
```

A direct verification shows that the calculation in the extended space matches the hand calculation:

```@repl coincidences2
two_photon_output_amp(occupation) = (two_photon_group_vals * two_photon_input_vec)[two_photon_occ_idx(occupation)]
extended_space_rate = sum(abs2(two_photon_output_amp(occupation)) for occupation in output_occupation_patterns)

perm2(matrix) = matrix[1, 1] * matrix[2, 2] + matrix[1, 2] * matrix[2, 1]
det2(matrix) = matrix[1, 1] * matrix[2, 2] - matrix[1, 2] * matrix[2, 1]

calc_rate = (1 + two_photon_overlap_sq) / 2 * abs2(perm2(two_photon_spatial_bs)) +
            (1 - two_photon_overlap_sq) / 2 * abs2(det2(two_photon_spatial_bs))
extended_space_rate, calc_rate, extended_space_rate ≈ calc_rate
```



## Three photons with one delayed input

Consider photons entering inputs $2,3,4$ of a four-mode interferometer, with the photon in mode 4 delayed relative to the others so that it has squared overlap $z$. For the output pattern `pp4`, the coincidence rate is

```math
R(234\to pp4;z)
=\frac{1}{3}(1+2z)\left|\mathrm{Per}(U_{234\to pp4})\right|^2
+\frac{2}{3}(1-z)\left|\mathrm{Imm}^{(2,1)}(U_{234\to pp4})\right|^2.
```


At zero delay, $z=1$, and the photons are fully indistinguishable; only the permanent term survives in the `pp4` rate.

As before, we can compute the `pp4` rate in two ways: directly from the permanent and immanant of the spatial submatrix, or from a single permanent in an extended space obtained by extending each spatial mode into early/late internal sub-modes. The two calculations agree. Below, we calculate the rates for $p=2$.

```@repl coincidences
using GroupFunctions
using LinearAlgebra: kron

#auxiliary functions
perm3(matrix) = matrix[1,1]*matrix[2,2]*matrix[3,3] + matrix[1,1]*matrix[2,3]*matrix[3,2] +
                matrix[1,2]*matrix[2,1]*matrix[3,3] + matrix[1,2]*matrix[2,3]*matrix[3,1] +
                matrix[1,3]*matrix[2,1]*matrix[3,2] + matrix[1,3]*matrix[2,2]*matrix[3,1];
imm21(matrix) = 2*matrix[1,1]*matrix[2,2]*matrix[3,3] - matrix[1,2]*matrix[2,3]*matrix[3,1] - matrix[1,3]*matrix[2,1]*matrix[3,2];


# linear optical transformation in 4 modes
four_mode_spatial_u = su2_block(4,1,(0.5,0.9,0.2)) * su2_block(4,2,(0.4,0.8,0.1)) *
    su2_block(4,3,(0.7,0.3,0.6)) * su2_block(4,2,(0.2,0.5,0.9));

overlap_amp = 0.37; # overlap

three_photon_overlap_sq = overlap_amp^2; # squared overlap
submatrix = four_mode_spatial_u[[2,2,4],[2,3,4]]; # we go from 2,3,4 to 2,2,4
#direct calculation from the formula above:
calc_rate = (1/2) * ((1/3) * (1 + 2 * three_photon_overlap_sq) * abs2(perm3(submatrix)) +
                     (2/3) * (1 - three_photon_overlap_sq) * abs2(imm21(submatrix)))

#… or extend each spatial mode into (early, late); one permanent in 8 modes
three_photon_extended_u = kron(four_mode_spatial_u, [1.0 0.0; 0.0 1.0]);
three_photon_group_vals, three_photon_occ_pats = group_function([3,0,0,0,0,0,0,0], three_photon_extended_u);

three_photon_occ_idx(occupation) = findfirst(p -> occupation_number(p) == occupation, three_photon_occ_pats);
three_photon_input_vec = zeros(ComplexF64, length(three_photon_occ_pats));
three_photon_input_vec[three_photon_occ_idx([0,0,1,0,1,0,1,0])] = overlap_amp; # photon in mode 4 is late
three_photon_input_vec[three_photon_occ_idx([0,0,1,0,1,0,0,1])] = sqrt(1 - overlap_amp^2); # photon in mode 4 is on time

#all observations of three photons in spatial modes 2 (two photons) and 4 (one photon),
#including all "late"/"on time" combinations:
output_occupation_patterns = [[0,0,2,0,0,0,1,0], [0,0,2,0,0,0,0,1],
                               [0,0,1,1,0,0,1,0], [0,0,1,1,0,0,0,1],
                               [0,0,0,2,0,0,1,0], [0,0,0,2,0,0,0,1]];
three_photon_output_amp(occupation) = (three_photon_group_vals * three_photon_input_vec)[three_photon_occ_idx(occupation)];
extended_space_rate = sum(abs2(three_photon_output_amp(occupation)) for occupation in output_occupation_patterns)

extended_space_rate ≈ calc_rate
```
