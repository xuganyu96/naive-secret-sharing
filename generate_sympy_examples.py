from sympy import Symbol, Poly
import random
import math


def sample_polynomial(indeterminate, degree, modulus):
    """Sample a polynomial up to the specified degree, with coefficients in the
    specified modulus
    """
    expr = 0
    for p in range(degree + 1):
        coeff = random.randint(0, modulus - 1)
        term = coeff * (indeterminate**p)
        expr += term
    return Poly(expr, modulus=modulus)


def hex_encode_coeffs(poly, gf_bits):
    """Convert the coefficients of polynomial into a left-padded hexadecimal
    string. Coeff of higher-power term will be encoded at the beginning of the
    string
    For example:
    in GF(2 ** 128), the polynomial (x ** 16 + 1) will be encoded as:
    0x0000_0000_0000_0000_0000_0000_0001_0001
    """
    digits = math.ceil(gf_bits / 4)
    encoding_num = 0
    for p in range(gf_bits):
        if poly.as_expr().coeff(Symbol("x"), p):
            encoding_num += 2**p
    hex_str = f"{encoding_num:X}"
    pad = digits - len(hex_str)
    hex_str = "0" * pad + hex_str
    return hex_str


def convert_to_limbs(hex_str, word_bits) -> list[str]:
    """Split the hexadecimal string according to word size in bits"""
    assert (word_bits % 4) == 0, "Word bits must be multiples of 4"
    word_hex_size = word_bits // 4
    i = 0
    limbs = []
    while i + word_hex_size <= len(hex_str):
        limb = hex_str[i : i + word_hex_size]
        limbs.append(f"0x{limb}")
        i += word_hex_size
    return limbs


def generate_extfield_widening_mul(gf_bits: int):
    lhs = sample_polynomial(Symbol("x"), gf_bits - 1, 2)
    rhs = sample_polynomial(Symbol("x"), gf_bits - 1, 2)
    prod = lhs * rhs

    lhs_limbs = convert_to_limbs(hex_encode_coeffs(lhs, gf_bits), 16)
    lhs_str = ", ".join(lhs_limbs)
    rhs_limbs = convert_to_limbs(hex_encode_coeffs(rhs, gf_bits), 16)
    rhs_str = ", ".join(rhs_limbs)
    prod_limbs = convert_to_limbs(hex_encode_coeffs(prod, gf_bits * 2), 16)
    prod_high_tokens = prod_limbs[: len(prod_limbs) // 2]
    prod_low_tokens = prod_limbs[len(prod_limbs) // 2 :]
    prod_high_str = ", ".join(prod_high_tokens)
    prod_low_str = ", ".join(prod_low_tokens)

    print(
        f"""
    let lhs = GF_2_{gf_bits}::from_limbs([
            {lhs_str}
    ]);
    let rhs = GF_2_{gf_bits}::from_limbs([
            {rhs_str}
    ]);
    let prod = (
        GF_2_{gf_bits}::from_limbs([
            {prod_high_str}
        ]),
        GF_2_{gf_bits}::from_limbs([
            {prod_low_str}
        ]),
    );
    assert_eq!(lhs.widening_gf_mul(&rhs), prod);"""
    )


def generate_short_euclidean_division(gf_bits: int):
    lhs = sample_polynomial(Symbol("x"), gf_bits - 1, 2)
    rhs = sample_polynomial(Symbol("x"), (gf_bits - 1) // 2, 2)
    rem = lhs % rhs
    quot = (lhs - rem) // rhs
    assert (rhs * quot + rem) == lhs

    lhs_limbs = convert_to_limbs(hex_encode_coeffs(lhs, gf_bits), 16)
    lhs_str = ", ".join(lhs_limbs)
    rhs_limbs = convert_to_limbs(hex_encode_coeffs(rhs, gf_bits), 16)
    rhs_str = ", ".join(rhs_limbs)
    quot_limbs = convert_to_limbs(hex_encode_coeffs(quot, gf_bits), 16)
    quot_str = ", ".join(quot_limbs)
    rem_limbs = convert_to_limbs(hex_encode_coeffs(rem, gf_bits), 16)
    rem_str = ", ".join(rem_limbs)

    print(
        f"""        assert_eq!(
            GF_2_128::from_limbs([
                {lhs_str}
            ]).euclidean_div(&GF_2_128::from_limbs([
                {rhs_str}
            ])),
            (GF_2_128::from_limbs([
                {quot_str}
            ]), GF_2_128::from_limbs([
                {rem_str}
            ]))
        );"""
    )


def random_gf2_128_modmul():
    x = Symbol("x")
    gf_bits, wordsize = 128, 16
    lhs = sample_polynomial(x, gf_bits - 1, 2)
    rhs = sample_polynomial(x, gf_bits - 1, 2)
    gf2_128_modulus = Poly(x**128 + x**7 + x**2 + 1, modulus=2)
    lhs_limbs = convert_to_limbs(hex_encode_coeffs(lhs, gf_bits), wordsize)
    lhs_str = ", ".join(lhs_limbs)
    rhs_limbs = convert_to_limbs(hex_encode_coeffs(rhs, gf_bits), wordsize)
    rhs_str = ", ".join(rhs_limbs)
    modulus_limbs = convert_to_limbs(
        hex_encode_coeffs(gf2_128_modulus, 2 * gf_bits), wordsize
    )
    # modulus_high_limbs = modulus_limbs[:len(modulus_limbs) // 2]
    modulus_low_limbs = modulus_limbs[len(modulus_limbs) // 2 :]
    modulus_low_str = ", ".join(modulus_low_limbs)
    rem = lhs * rhs % gf2_128_modulus
    rem_limbs = convert_to_limbs(hex_encode_coeffs(rem, 128), wordsize)
    rem_str = ", ".join(rem_limbs)
    # print(
    #     f"""        let lhs = F2_128::from_limbs([
    #         {lhs_str}
    #     ]);
    #     let rhs = F2_128::from_limbs([
    #         {rhs_str}
    #     ]);
    #     let modulus = WideF2X::<8>::from_f2x(
    #         F2_128::ONE,
    #         F2_128::from_limbs([
    #             {modulus_low_str}
    #         ]),
    #     );
    #     let rem = F2_128::from_limbs([
    #         {rem_str}
    #     ]);
    #     assert_eq!(lhs.modmul(&rhs, &modulus), rem);"""
    # )
    print(
        f"""
        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            {lhs_str}
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            {rhs_str}
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            {rem_str}
        ]));
        assert_eq!(lhs.mul(&rhs), rem);"""
    )


if __name__ == "__main__":
    # generate_extfield_widening_mul(128)
    # generate_short_euclidean_division(128)
    random_gf2_128_modmul()
