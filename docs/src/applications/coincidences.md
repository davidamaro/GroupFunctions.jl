```@meta
CurrentModule = GroupFunctions
CollapsedDocStrings = true
DocTestSetup = GroupFunctions.doctestsetup()
```
# Coincidence rates

The contents of this page can be thought of as a mathematical description and formalisation of the following experiment:
1. We have a device that can output given number of photons in specific spatial modes at will, or at least *herald preparation of such a state* (for instance with multiplexed SPDC), so that the photons are contained in well-defined time windows.
2. The photons can be delayed with respect to each other (e.g. by having the heralding photon traverse a longer distance). 
3. We pass the photons through linear optical system, and measure coincidences of observing clicks in specific output ports.

Distinguishability of the photons, controlled by the relative delay, affects the observed coincidences rates. This can be mathematically modelled by a frequency-varying photon profile: this is roughly the same as a time profile after a Fourier transform. Each spatial mode $k$ contains a continuum of frequency modes $\omega$ (with annilation operators labelled $\hat a_k(\omega)$), and an operator creating a photon with a Gaussian profile, in mode $k$, centered at time $\tau$, is

```math
\hat{A}^\dagger_k(\tau)=
\int d\omega\,e^{i\omega\tau}\varphi(\omega)\hat{a}^\dagger_k(\omega),
```
As a concrete simple example, a photon with Gaussian spectrum is created by such $\hat{A}^\dagger$ with the following frequency profile $\varphi$:

```math
|\varphi(\omega)|^2=\frac{\exp({-(\omega-\omega_0)^2/(2\sigma^2)})}{\sqrt{2\pi}\sigma},
```

If we have two photons with the same spectral profile, but centered at different times $\tau, \tau'$, the overlap of their wavefunctions $\langle 0\vert \hat A (\tau) \hat A^\dagger(\tau') \vert 0\rangle$ is

```math
\zeta=\exp({-\frac12\sigma^2(\tau-\tau')^2}).
```

In the next two sections we will calculate coincidence rates as a function of the overlap.
 It will be useful to define the squared overlap 

$$z=\vert \zeta\vert^2.$$

 To keep track of the indices "which photon goes where", the following shorthand is introduced:

```math
R(\text{inputs} \to \text{outputs}, z)
```

meaning: coincidence probability for a fixed input occupation pattern, output occupation pattern; one (specified in the text) photon is delayed, with squared overlap $z$ with respect to the other ones.

For example, $R(23 \to 13; z)$ means "one photon enters input 2 and one enters input 3; detect one at output 1 and one at output 3". Similarly, $R(234 \to \alpha\alpha4; z)$ reads "one photon enters each of inputs 2,3,4; detect two photons at output $\alpha$ and one at output 4".

Repeated output labels encode multiplicity, $R(234\to pp4;z)$ meaning that we observe two photons at port $p$.

## Two-photon coincidence rate with delay

For two photons entering input modes $k,l$ and detected at outputs $m,n$, with squared overlap $z$:

```math
R(kl \to mn;z)
=\frac{1}{2}(1+z)\left|\mathrm{Per}(U_{kl\to mn})\right|^2
+\frac{1}{2}(1-z)\left|\mathrm{Det}(U_{kl\to mn})\right|^2,
```

where $U_{kl \to mn}$ is a matrix constructed from $U$ by picking the specified columns and rows (see [the background page](../background/group_functions.md#Bosons:-the-permanent)).

If the delay is zero, photons are completely indistinguishable ($z=1$) and the result collapses to a standard permanent:

```math
\tau=0 \Rightarrow R=\left|\mathrm{Per}(U_{kl\to mn})\right|^2,
```

With infinite delay, the photons are completely distinguishable ($z=0$) and the rate is an equal mixture of permanent and determinant: 

```math
\tau\to\infty \Rightarrow
R=\frac{1}{2}\left|\mathrm{Per}(U_{kl\to mn})\right|^2
+\frac{1}{2}\left|\mathrm{Det}(U_{kl\to mn})\right|^2.
```


We can verify these equations with our library. Naturally, it is infeasible to track a continuum of frequency modes. But the only thing that matters is the overlap, and we can simulate the behavior by introducing internal states in addition to the spatial modes. 

We restrict to two spatial modes to exemplify an usage. The calculationss are done in the basis where occupation numbers $(a,a',b,b')$ denote $a=$photons "on time" in spatial mode 1, $a'$=photons "late" (out of time) in spatial mode 1; $b$=photons "on time" in spatial mode 2, $b'$=photon "late" in spatial mode 2.

In this space, linear optical transformation in spatial modes get inflated (transformations are constant on coherence timescales, hence the identity):

```math
V = U_{\text{spatial}} \otimes \operatorname{id}_{\text{internal}} 
```

So we inflate a beamsplitter and define a group function in the enlarged space:
```@repl coincidences2
using GroupFunctions
using LinearAlgebra: kron
u2 = su2_block(2,1,(1.5,0.9,0.2));          # 2×2 spatial beam splitter
V2=kron(u2, [1.0 0.0; 0.0 1.0]);

vals, pats = group_function([2,0,0,0], V2);
```

Now we define the initial state:


```@repl coincidences2
ζ = 0.37; # amplitude overlap

z = ζ^2; # squared overlap
#initial state is spanned by these two:
occ_a = [1,0,1,0]; # one photon in sp. mode 1, one in sp. mode 2, both "on time"
occ_b = [1,0,0,1]; # one photon "on time" in sp. mode 1, one "late" in 2
idx(occ) = findfirst(p -> occupation_number(p) == occ, pats);
in_vec = zeros(ComplexF64, length(pats));
in_vec[idx(occ_a)] = ζ ;
in_vec[idx(occ_b)] = sqrt(1-ζ^2);

```

The coincidence rate are calculated with four possible states; all are "one photon in spatial mode 1, one in spatial mode 2", but they have different internal ("on time"/"late") indices:


```@repl coincidences2
out_occs = [
    [1,0,1,0],
    [1,0,0,1],
    [0,1,1,0],
    [0,1,0,1]
];
```

Direct verification shows that the calculations using inflated space and hand-calculated one match:

```@repl coincidences2
amp(occ) = (vals * in_vec)[idx(occ)]
R_infl = sum(abs2(amp(o)) for o in out_occs)

permanent2(M) = M[1,1]*M[2,2]+M[1,2]*M[2,1]
det2(M) = M[1,1]*M[2,2]-M[1,2]*M[2,1]

R_calc = (1+z)/2*abs2(permanent2(u2)) + (1-z)/2*abs2(det2(u2))
R_infl, R_calc, R_infl ≈ R_calc
```



## Three photons with one delayed input

For inputs $2,3,4$ in a four-mode interferometer, with the photon in mode 4 delayed with respect to the other ones (so that it has squared overlap $z$), and output pattern `pp4`:

```math
R(234\to pp4;z)
=\frac{1}{3}(1+2z})\left|\mathrm{Per}(U_{234\to pp4})\right|^2
+\frac{2}{3}(1-z)\left|\mathrm{Imm}^{(2,1)}(U_{234\to pp4})\right|^2.
```


At zero delay, $z=1$ and photons are fully indistinguishable; only the permanent term survives in the `pp4` rate.

Just as before, the `pp4` rate can be computed two ways: directly from the permanent and immanant of the spatial submatrix, or by inflating each spatial mode into early/late internal sub-modes and using a single permanent in the enlarged space. They agree; below we calculate the rates for $p=2$.

```@repl coincidences
using GroupFunctions
using LinearAlgebra: kron

#auxiliary functions
permanent3(M) = M[1,1]*M[2,2]*M[3,3] + M[1,1]*M[2,3]*M[3,2] +
                M[1,2]*M[2,1]*M[3,3] + M[1,2]*M[2,3]*M[3,1] +
                M[1,3]*M[2,1]*M[3,2] + M[1,3]*M[2,2]*M[3,1];
imm21(M) = 2*M[1,1]*M[2,2]*M[3,3] - M[1,2]*M[2,3]*M[3,1] - M[1,3]*M[2,1]*M[3,2];


# linear optical transformation in 4 modes
u = su2_block(4,1,(0.5,0.9,0.2)) * su2_block(4,2,(0.4,0.8,0.1)) *
    su2_block(4,3,(0.7,0.3,0.6)) * su2_block(4,2,(0.2,0.5,0.9));

ζ = 0.37;

z = ζ^2;
M = u[[2,2,4],[2,3,4]]; # we go from 2,3,4 to 2,2,4
#direct calculation from the formula above:
R_calc = (1/2) * ((1/3)*(1+2z)*abs2(permanent3(M)) + (2/3)*(1-z)*abs2(imm21(M)))

#… or inflate each spatial mode into (early, late); one permanent in 8 modes
V = kron(u, [1.0 0.0; 0.0 1.0]);
vals, pats = group_function([3,0,0,0,0,0,0,0], V);

idx(occ) = findfirst(p -> occupation_number(p) == occ, pats);
in_vec = zeros(ComplexF64, length(pats));
in_vec[idx([0,0,1,0,1,0,1,0])] = ζ; # photon in mode 4 is late
in_vec[idx([0,0,1,0,1,0,0,1])] = sqrt(1-ζ^2); #photon in mode 4 is on time

#all observations of three photons in spatial modes 2 (two photons) and 4 (one photon),
#including all "late"/"on time" combinations:
out_occs = [[0,0,2,0,0,0,1,0], [0,0,2,0,0,0,0,1],
            [0,0,1,1,0,0,1,0], [0,0,1,1,0,0,0,1],
            [0,0,0,2,0,0,1,0], [0,0,0,2,0,0,0,1]]; 
amp(occ) = (vals * in_vec)[idx(occ)];
R_infl = sum(abs2(amp(o)) for o in out_occs)

R_infl ≈ R_calc
```

