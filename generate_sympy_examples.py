from sympy import Symbol, Poly
import random
import math

# Generate random test cases for widening multiplication in GF(2 ** 128)
gf_bits = 128

# Generate random LHS and RHS, then compute product
lhs = Poly(
    sum([random.randint(0, 1) * Symbol("x") ** p for p in range(gf_bits)]), 
    modulus=2
)
rhs = Poly(
    sum([random.randint(0, 1) * Symbol("x") ** p for p in range(gf_bits)]), 
    modulus=2
)
prod = lhs * rhs

# Encode coefficients into binaries
def hex_encode_coeffs(poly, bits):
    digits = math.ceil(bits / 4)
    encoding_num = 0
    for p in range(bits):
        if poly.as_expr().coeff(Symbol("x") ** p):
            encoding_num += 2 ** p
    hex_str = f"{encoding_num:X}"
    pad = digits - len(hex_str)
    hex_str = "0" * pad + hex_str
    return hex_str

def convert_to_limbs(hex_str, word_bits):
    word_hex_size = word_bits // 4
    i = 0
    limbs = []
    while (i + word_hex_size <= len(hex_str)):
        limb = hex_str[i:i+word_hex_size]
        limbs.append(f"0x{limb}")
        i += word_hex_size
    return ", ".join(limbs)


print("LHS: ", convert_to_limbs(hex_encode_coeffs(lhs, gf_bits), 16))
print("RHS: ", convert_to_limbs(hex_encode_coeffs(rhs, gf_bits), 16))
print("Prod: ", convert_to_limbs(hex_encode_coeffs(prod, gf_bits * 2), 16))
