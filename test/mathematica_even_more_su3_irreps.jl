function mathematica_even_more_su3_irrep_cases()
    return [
        (
            [2, 0, 0],
            "u12_mix",
            ComplexF64[
                ComplexF64(0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0);
                ComplexF64(-0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(1, 0)
            ],
            ComplexF64[
                ComplexF64(1, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0.5, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(-0.70710678118654757, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0.5, 0)
            ],
        ),
        (
            [2, 0, 0],
            "u23_mix",
            ComplexF64[
                ComplexF64(1, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0);
                ComplexF64(0, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0)
            ],
            ComplexF64[
                ComplexF64(0.5, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0, 0);
                ComplexF64(0.5, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(1, 0)
            ],
        ),
        (
            [2, 0, 0],
            "u_full",
            ComplexF64[
                ComplexF64(0.70710678118654757, 0) ComplexF64(0.5, 0) ComplexF64(0.5, 0);
                ComplexF64(-0.5, 0) ComplexF64(-0.14644660940672624, 0) ComplexF64(0.85355339059327373, 0);
                ComplexF64(0.5, 0) ComplexF64(-0.85355339059327373, 0) ComplexF64(0.14644660940672624, 0)
            ],
            ComplexF64[
                ComplexF64(0.021446609406726238, 0) ComplexF64(-0.17677669529663689, 0) ComplexF64(0.10355339059327376, 0) ComplexF64(0.72855339059327373, 0) ComplexF64(-0.60355339059327373, 0) ComplexF64(0.25, 0);
                ComplexF64(0.17677669529663689, 0) ComplexF64(-0.75, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(0.17677669529663689, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(-0.35355339059327379, 0);
                ComplexF64(0.10355339059327376, 0) ComplexF64(-0.35355339059327379, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(-0.60355339059327373, 0) ComplexF64(-0.35355339059327379, 0) ComplexF64(0.5, 0);
                ComplexF64(0.72855339059327373, 0) ComplexF64(-0.17677669529663689, 0) ComplexF64(-0.60355339059327373, 0) ComplexF64(0.021446609406726238, 0) ComplexF64(0.10355339059327376, 0) ComplexF64(0.25, 0);
                ComplexF64(0.60355339059327373, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(-0.10355339059327376, 0) ComplexF64(-0.35355339059327379, 0) ComplexF64(-0.5, 0);
                ComplexF64(0.25, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(0.5, 0) ComplexF64(0.25, 0) ComplexF64(0.5, 0) ComplexF64(0.5, 0)
            ],
        ),
        (
            [2, 2, 0],
            "u12_mix",
            ComplexF64[
                ComplexF64(0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0);
                ComplexF64(-0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(1, 0)
            ],
            ComplexF64[
                ComplexF64(0.5, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0.5, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(1, 0)
            ],
        ),
        (
            [2, 2, 0],
            "u23_mix",
            ComplexF64[
                ComplexF64(1, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0);
                ComplexF64(0, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0)
            ],
            ComplexF64[
                ComplexF64(1, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0.5, 0);
                ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.70710678118654757, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0.5, 0)
            ],
        ),
        (
            [2, 2, 0],
            "u_full",
            ComplexF64[
                ComplexF64(0.70710678118654757, 0) ComplexF64(0.5, 0) ComplexF64(0.5, 0);
                ComplexF64(-0.5, 0) ComplexF64(-0.14644660940672624, 0) ComplexF64(0.85355339059327373, 0);
                ComplexF64(0.5, 0) ComplexF64(-0.85355339059327373, 0) ComplexF64(0.14644660940672624, 0)
            ],
            ComplexF64[
                ComplexF64(0.5, 0) ComplexF64(-0.5, 0) ComplexF64(0.25, 0) ComplexF64(0.5, 0) ComplexF64(-0.35355339059327379, 0) ComplexF64(0.25, 0);
                ComplexF64(0.5, 0) ComplexF64(-0.35355339059327379, 0) ComplexF64(0.10355339059327376, 0) ComplexF64(-0.35355339059327379, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(-0.60355339059327373, 0);
                ComplexF64(0.25, 0) ComplexF64(-0.10355339059327376, 0) ComplexF64(0.021446609406726238, 0) ComplexF64(-0.60355339059327373, 0) ComplexF64(0.17677669529663689, 0) ComplexF64(0.72855339059327373, 0);
                ComplexF64(0.5, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(-0.60355339059327373, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(0.10355339059327376, 0);
                ComplexF64(0.35355339059327379, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(-0.17677669529663689, 0) ComplexF64(-0.35355339059327379, 0) ComplexF64(-0.75, 0) ComplexF64(-0.17677669529663689, 0);
                ComplexF64(0.25, 0) ComplexF64(0.60355339059327373, 0) ComplexF64(0.72855339059327373, 0) ComplexF64(0.10355339059327376, 0) ComplexF64(0.17677669529663689, 0) ComplexF64(0.021446609406726238, 0)
            ],
        ),
        (
            [3, 1, 0],
            "u12_mix",
            ComplexF64[
                ComplexF64(0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0);
                ComplexF64(-0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(1, 0)
            ],
            ComplexF64[
                ComplexF64(0.70710678118654757, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(1, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(-0.61237243569579447, 0) ComplexF64(0.61237243569579447, 0) ComplexF64(-0.35355339059327379, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.61237243569579447, 0) ComplexF64(-0.35355339059327379, 0) ComplexF64(-0.35355339059327379, 0) ComplexF64(0.61237243569579447, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.61237243569579447, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(-0.35355339059327379, 0) ComplexF64(-0.61237243569579447, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(0.61237243569579447, 0) ComplexF64(0.61237243569579447, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0.5, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(-0.70710678118654757, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0.5, 0)
            ],
        ),
        (
            [3, 1, 0],
            "u23_mix",
            ComplexF64[
                ComplexF64(1, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0);
                ComplexF64(0, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0.70710678118654757, 0)
            ],
            ComplexF64[
                ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(-0.5, 0) ComplexF64(0, 0) ComplexF64(-0.35355339059327379, 0) ComplexF64(0, 0) ComplexF64(0.57735026918962573, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.20412414523193151, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.35355339059327379, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.5, 0) ComplexF64(0, 0) ComplexF64(-0.40824829046386302, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.28867513459481287, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(-0.5, 0) ComplexF64(0, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.61237243569579447, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(-0.57735026918962573, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.40824829046386302, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0.57735026918962573, 0) ComplexF64(0.40824829046386302, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.23570226039551584, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.33333333333333331, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.57735026918962573, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.57735026918962573, 0) ComplexF64(0, 0) ComplexF64(0.33333333333333331, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.47140452079103168, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.57735026918962573, 0) ComplexF64(0, 0);
                ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0.20412414523193151, 0) ComplexF64(-0.28867513459481287, 0) ComplexF64(0, 0) ComplexF64(0.61237243569579447, 0) ComplexF64(0, 0) ComplexF64(-0.33333333333333331, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.58925565098878963, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.20412414523193151, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.40824829046386302, 0) ComplexF64(0, 0) ComplexF64(-0.47140452079103168, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.66666666666666663, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.40824829046386302, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(-0.70710678118654757, 0);
                ComplexF64(0, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(0, 0) ComplexF64(0.57735026918962573, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.20412414523193151, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.35355339059327379, 0) ComplexF64(0, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(0, 0) ComplexF64(0.57735026918962573, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.40824829046386302, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.5, 0) ComplexF64(0, 0);
                ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0) ComplexF64(0, 0) ComplexF64(0, 0) ComplexF64(0.70710678118654757, 0)
            ],
        ),
        (
            [3, 1, 0],
            "u_full",
            ComplexF64[
                ComplexF64(0.70710678118654757, 0) ComplexF64(0.5, 0) ComplexF64(0.5, 0);
                ComplexF64(-0.5, 0) ComplexF64(-0.14644660940672624, 0) ComplexF64(0.85355339059327373, 0);
                ComplexF64(0.5, 0) ComplexF64(-0.85355339059327373, 0) ComplexF64(0.14644660940672624, 0)
            ],
            ComplexF64[
                ComplexF64(0.015165042944955322, 0) ComplexF64(-0.010723304703363119, 0) ComplexF64(0.015165042944955322, 0) ComplexF64(-0.125, 0) ComplexF64(0.11427669529663688, 0) ComplexF64(-0.051776695296636879, 0) ComplexF64(-0.10206207261596575, 0) ComplexF64(0.059786577934525069, 0) ComplexF64(0.51516504294495535, 0) ComplexF64(-0.55877696061835858, 0) ComplexF64(0.34846171252933794, 0) ComplexF64(-0.125, 0) ComplexF64(0.36427669529663687, 0) ComplexF64(-0.30177669529663687, 0) ComplexF64(0.125, 0);
                ComplexF64(0.010723304703363119, 0) ComplexF64(-0.0031407832308854582, 0) ComplexF64(-0.02588834764831844, 0) ComplexF64(-0.088388347648318447, 0) ComplexF64(0.054917478527522337, 0) ComplexF64(-0.015165042944955322, 0) ComplexF64(0.17423085626466897, 0) ComplexF64(-0.10206207261596575, 0) ComplexF64(0.36427669529663687, 0) ComplexF64(-0.30799954989171524, 0) ComplexF64(0.14433756729740643, 0) ComplexF64(-0.03661165235168156, 0) ComplexF64(-0.62185921676911449, 0) ComplexF64(0.51516504294495535, 0) ComplexF64(-0.21338834764831843, 0);
                ComplexF64(0.015165042944955322, 0) ComplexF64(0.02588834764831844, 0) ComplexF64(-0.14016504294495533, 0) ComplexF64(-0.125, 0) ComplexF64(-0.099111652351681553, 0) ComplexF64(0.125, 0) ComplexF64(0.45052378514530372, 0) ComplexF64(-0.26391072316645658, 0) ComplexF64(0.51516504294495535, 0) ComplexF64(0.15928421178103772, 0) ComplexF64(-0.49279927982674437, 0) ComplexF64(0.30177669529663687, 0) ComplexF64(0.15088834764831843, 0) ComplexF64(-0.125, 0) ComplexF64(0.051776695296636879, 0);
                ComplexF64(0.125, 0) ComplexF64(-0.088388347648318447, 0) ComplexF64(0.125, 0) ComplexF64(-0.5303300858899106, 0) ComplexF64(0.44194173824159222, 0) ComplexF64(-0.17677669529663689, 0) ComplexF64(-0.4330127018922193, 0) ComplexF64(0.20412414523193151, 0) ComplexF64(0.125, 0) ComplexF64(0.15309310892394862, 0) ComplexF64(-0.28867513459481287, 0) ComplexF64(0.17677669529663689, 0) ComplexF64(0.088388347648318447, 0) ComplexF64(0.17677669529663689, 0) ComplexF64(-0.17677669529663689, 0);
                ComplexF64(0.11427669529663688, 0) ComplexF64(-0.054917478527522337, 0) ComplexF64(-0.099111652351681553, 0) ComplexF64(-0.44194173824159222, 0) ComplexF64(0.35669417382415924, 0) ComplexF64(-0.16161165235168157, 0) ComplexF64(0.37835500149660051, 0) ComplexF64(-0.10206207261596575, 0) ComplexF64(-0.23927669529663689, 0) ComplexF64(0.070355451604885239, 0) ComplexF64(0.14433756729740643, 0) ComplexF64(-0.14016504294495533, 0) ComplexF64(-0.32008252147247768, 0) ComplexF64(-0.33838834764831843, 0) ComplexF64(0.3901650429449553, 0);
                ComplexF64(0.051776695296636879, 0) ComplexF64(-0.015165042944955322, 0) ComplexF64(-0.125, 0) ComplexF64(-0.17677669529663689, 0) ComplexF64(0.16161165235168157, 0) ComplexF64(-0.051776695296636879, 0) ComplexF64(0.34846171252933794, 0) ComplexF64(-0.34846171252933794, 0) ComplexF64(-0.30177669529663687, 0) ComplexF64(-0.093306530989423569, 0) ComplexF64(0.18661306197884714, 0) ComplexF64(-0.073223304703363121, 0) ComplexF64(0.51516504294495535, 0) ComplexF64(0.30177669529663687, 0) ComplexF64(-0.42677669529663687, 0);
                ComplexF64(0.10206207261596575, 0) ComplexF64(0.17423085626466897, 0) ComplexF64(-0.45052378514530372, 0) ComplexF64(-0.4330127018922193, 0) ComplexF64(-0.37835500149660051, 0) ComplexF64(0.34846171252933794, 0) ComplexF64(-0.09763107293781749, 0) ComplexF64(0.23570226039551584, 0) ComplexF64(0.10206207261596575, 0) ComplexF64(0.26725889843221229, 0) ComplexF64(0.16666666666666666, 0) ComplexF64(-0.34846171252933794, 0) ComplexF64(0.029893288967262534, 0) ComplexF64(0.059786577934525069, 0) ComplexF64(-0.059786577934525069, 0);
                ComplexF64(0.059786577934525069, 0) ComplexF64(0.10206207261596575, 0) ComplexF64(-0.26391072316645658, 0) ComplexF64(-0.20412414523193151, 0) ComplexF64(-0.10206207261596575, 0) ComplexF64(0.34846171252933794, 0) ComplexF64(-0.23570226039551584, 0) ComplexF64(-0.16666666666666666, 0) ComplexF64(-0.34846171252933794, 0) ComplexF64(-0.51011002862997024, 0) ComplexF64(-0.11785113019775792, 0) ComplexF64(0.49279927982674437, 0) ComplexF64(-0.10206207261596575, 0) ComplexF64(-0.059786577934525069, 0) ComplexF64(0.084550989362881371, 0);
                ComplexF64(0.51516504294495535, 0) ComplexF64(-0.36427669529663687, 0) ComplexF64(0.51516504294495535, 0) ComplexF64(-0.125, 0) ComplexF64(-0.23927669529663689, 0) ComplexF64(0.30177669529663687, 0) ComplexF64(-0.10206207261596575, 0) ComplexF64(-0.34846171252933794, 0) ComplexF64(0.015165042944955322, 0) ComplexF64(0.053595475077435992, 0) ComplexF64(0.059786577934525069, 0) ComplexF64(-0.125, 0) ComplexF64(0.010723304703363119, 0) ComplexF64(0.051776695296636879, 0) ComplexF64(0.125, 0);
                ComplexF64(0.55877696061835858, 0) ComplexF64(-0.30799954989171524, 0) ComplexF64(-0.15928421178103772, 0) ComplexF64(0.15309310892394862, 0) ComplexF64(-0.070355451604885239, 0) ComplexF64(-0.093306530989423569, 0) ComplexF64(0.26725889843221229, 0) ComplexF64(0.51011002862997024, 0) ComplexF64(-0.053595475077435992, 0) ComplexF64(-0.1188980579413864, 0) ComplexF64(-0.014297739604484159, 0) ComplexF64(0.18298639789121116, 0) ComplexF64(-0.052844368351800862, 0) ComplexF64(-0.19536860360538932, 0) ComplexF64(-0.32732396518861762, 0);
                ComplexF64(0.34846171252933794, 0) ComplexF64(-0.14433756729740643, 0) ComplexF64(-0.49279927982674437, 0) ComplexF64(0.28867513459481287, 0) ComplexF64(0.14433756729740643, 0) ComplexF64(-0.18661306197884714, 0) ComplexF64(-0.16666666666666666, 0) ComplexF64(-0.11785113019775792, 0) ComplexF64(0.059786577934525069, 0) ComplexF64(0.014297739604484159, 0) ComplexF64(-0.083333333333333329, 0) ComplexF64(-0.084550989362881371, 0) ComplexF64(0.14433756729740643, 0) ComplexF64(0.39073720721077865, 0) ComplexF64(0.49279927982674437, 0);
                ComplexF64(0.125, 0) ComplexF64(-0.03661165235168156, 0) ComplexF64(-0.30177669529663687, 0) ComplexF64(0.17677669529663689, 0) ComplexF64(0.14016504294495533, 0) ComplexF64(-0.073223304703363121, 0) ComplexF64(-0.34846171252933794, 0) ComplexF64(-0.49279927982674437, 0) ComplexF64(0.125, 0) ComplexF64(0.18298639789121116, 0) ComplexF64(0.084550989362881371, 0) ComplexF64(-0.073223304703363121, 0) ComplexF64(-0.21338834764831843, 0) ComplexF64(-0.42677669529663687, 0) ComplexF64(-0.42677669529663687, 0);
                ComplexF64(0.36427669529663687, 0) ComplexF64(0.62185921676911449, 0) ComplexF64(0.15088834764831843, 0) ComplexF64(-0.088388347648318447, 0) ComplexF64(-0.32008252147247768, 0) ComplexF64(-0.51516504294495535, 0) ComplexF64(-0.029893288967262534, 0) ComplexF64(-0.10206207261596575, 0) ComplexF64(0.010723304703363119, 0) ComplexF64(0.052844368351800862, 0) ComplexF64(0.14433756729740643, 0) ComplexF64(0.21338834764831843, 0) ComplexF64(0.0031407832308854582, 0) ComplexF64(0.015165042944955322, 0) ComplexF64(0.03661165235168156, 0);
                ComplexF64(0.30177669529663687, 0) ComplexF64(0.51516504294495535, 0) ComplexF64(0.125, 0) ComplexF64(0.17677669529663689, 0) ComplexF64(0.33838834764831843, 0) ComplexF64(0.30177669529663687, 0) ComplexF64(0.059786577934525069, 0) ComplexF64(0.059786577934525069, 0) ComplexF64(-0.051776695296636879, 0) ComplexF64(-0.19536860360538932, 0) ComplexF64(-0.39073720721077865, 0) ComplexF64(-0.42677669529663687, 0) ComplexF64(-0.015165042944955322, 0) ComplexF64(-0.051776695296636879, 0) ComplexF64(-0.073223304703363121, 0);
                ComplexF64(0.125, 0) ComplexF64(0.21338834764831843, 0) ComplexF64(0.051776695296636879, 0) ComplexF64(0.17677669529663689, 0) ComplexF64(0.3901650429449553, 0) ComplexF64(0.42677669529663687, 0) ComplexF64(0.059786577934525069, 0) ComplexF64(0.084550989362881371, 0) ComplexF64(0.125, 0) ComplexF64(0.32732396518861762, 0) ComplexF64(0.49279927982674437, 0) ComplexF64(0.42677669529663687, 0) ComplexF64(0.03661165235168156, 0) ComplexF64(0.073223304703363121, 0) ComplexF64(0.073223304703363121, 0)
            ],
        ),
    ]
end

@testset "Mathematica references for more SU(3) irreps" begin
    for (irrep, case_name, input_matrix, expected_representation) in mathematica_even_more_su3_irrep_cases()
        @testset "$irrep $case_name" begin
            actual_representation, _ = group_function(irrep, input_matrix)
            @test isapprox(actual_representation, expected_representation; atol=1e-12, rtol=1e-12)
        end
    end
end
