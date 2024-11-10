// #![no_std]

pub mod f2x;

use crate::f2x::{F2x, WideF2x};

/// An element of the extension field GF(2^128)
#[derive(Debug, PartialEq, Eq)]
pub struct GF2_128 {
    pub poly: F2x<{ Self::LIMBS }>,
}

impl GF2_128 {
    pub const LIMBS: usize = 8;
    /// NOTE: Mathematically all finite fields of the same order are isomorphic, so the choice of
    /// the modulus polynomial should not have mattered. However, there might be a need to use a
    /// specific modulus, so future builds might allow users to specify custom modulus
    /// This is (x^128 + x^77 + x^35 + x^11 + 1)
    pub const MODULUS: WideF2x<{ Self::LIMBS }> = WideF2x::from_f2x(
        F2x::<{ Self::LIMBS }>::ONE,
        F2x::<{ Self::LIMBS }>::from_limbs([
            0x0000, 0x0000, 0x0000, 0x2000, 0x0000, 0x0008, 0x0000, 0x0801,
        ]),
    );

    pub const ONE: Self = Self::from_poly(F2x::<{ Self::LIMBS }>::ONE);

    pub const fn from_poly(poly: F2x<{ Self::LIMBS }>) -> Self {
        Self { poly }
    }

    pub fn add(&self, rhs: &Self) -> Self {
        Self::from_poly(self.poly.add(&rhs.poly))
    }

    pub fn sub(&self, rhs: &Self) -> Self {
        self.add(rhs)
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        let prod = self.poly.widening_mul(&rhs.poly);
        let (_, rem) = prod.div_rem(&Self::MODULUS);
        Self::from_poly(rem.truncate())
    }

    pub fn inv(&self) -> Option<Self> {
        let inverse = self.poly.modinv(&Self::MODULUS);
        inverse.map_or(None, |poly| Some(Self::from_poly(poly)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gf2_128_mul() {
        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x1254, 0x4198, 0x8DA7, 0x29BD, 0xECF1, 0x64DE, 0xFBA7, 0xB692,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x7D89, 0xD76A, 0x644E, 0x3A1C, 0x047C, 0xB60A, 0x1B98, 0x30F0,
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x38B8, 0xB93F, 0xE30A, 0xE55E, 0xFC24, 0x2F4D, 0x5A16, 0xBD14,
        ]));
        assert_eq!(lhs.mul(&rhs), rem);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x5E17, 0xA183, 0x0FB7, 0xC7A9, 0xD3D3, 0xD1A4, 0x8962, 0xC3F5,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x6ECC, 0x895B, 0x76CC, 0x855E, 0x9C14, 0xEF6F, 0x587A, 0x4A04,
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xC42C, 0xBD01, 0xE0B0, 0x693D, 0x5E49, 0x6D28, 0xA69F, 0x701A,
        ]));
        assert_eq!(lhs.mul(&rhs), rem);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x6C70, 0x49B5, 0x97A4, 0x65D1, 0x2370, 0x8DBE, 0x127F, 0x5EFB,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x3387, 0xF3D0, 0xBD53, 0xADF3, 0x2994, 0x3B7A, 0xB2A8, 0x974A,
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x4402, 0xF342, 0x347A, 0xA630, 0x0B5D, 0x31BB, 0xBED8, 0x76A2,
        ]));
        assert_eq!(lhs.mul(&rhs), rem);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x98AA, 0xD35C, 0x02F5, 0x612C, 0x67A1, 0x9B5C, 0xF1FE, 0x98C7,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x4F13, 0x9428, 0x6D75, 0x61C6, 0x2CE7, 0xA102, 0xB546, 0xC183,
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x09D9, 0x9D26, 0x5665, 0x2DAE, 0x2A22, 0x9928, 0xC29C, 0xD153,
        ]));
        assert_eq!(lhs.mul(&rhs), rem);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x1EA8, 0xE4CA, 0x33E9, 0x7C78, 0xD1E4, 0x903F, 0x6F70, 0x20F8,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xEA22, 0x6D4C, 0xC5D4, 0x0C82, 0xF584, 0x3185, 0x36F4, 0x12B2,
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x8AF2, 0x5D11, 0xB56A, 0x68C0, 0xC3DD, 0xD9A1, 0xAC6E, 0x930B,
        ]));
        assert_eq!(lhs.mul(&rhs), rem);
    }

    #[test]
    fn random_gf2_128_inv() {
        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xA5A0, 0xD9E5, 0xCD78, 0xB647, 0x3A5E, 0xBCC6, 0x3BAB, 0xAAF4,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xA843, 0xF793, 0x72D1, 0xF7C3, 0xFAD0, 0x2824, 0xFC83, 0x72CE,
        ]));
        assert_eq!(lhs.inv().unwrap(), rhs);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xDD55, 0xB556, 0x029D, 0x9FB2, 0xCCFA, 0x34DB, 0x51C3, 0x1A72,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x44B8, 0x88EF, 0xF7F0, 0xAC9B, 0x7EE4, 0x0871, 0x19CF, 0x0C09,
        ]));
        assert_eq!(lhs.inv().unwrap(), rhs);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x8369, 0x1C2E, 0x68E5, 0x9EAA, 0x0B5A, 0xAF97, 0x359B, 0x49B0,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xEC62, 0x19D1, 0x056C, 0x8568, 0xCD7A, 0x0482, 0xFD99, 0xD508,
        ]));
        assert_eq!(lhs.inv().unwrap(), rhs);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x5E01, 0x689F, 0xB0F9, 0x5C2E, 0xA122, 0x9596, 0xE179, 0x74D6,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xD1FE, 0x3C26, 0x487E, 0x7952, 0x2AF6, 0x3FC2, 0x8CC6, 0x2CF6,
        ]));
        assert_eq!(lhs.inv().unwrap(), rhs);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x192A, 0x0A05, 0x3E08, 0x6DF8, 0x05B1, 0x3EA8, 0x6575, 0x6A4A,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x6A8E, 0x5EAB, 0x61DB, 0xD68F, 0x7DF2, 0xA04A, 0x2374, 0x46C1,
        ]));
        assert_eq!(lhs.inv().unwrap(), rhs);
    }
}
