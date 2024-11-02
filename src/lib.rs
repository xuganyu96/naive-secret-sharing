#![no_std]

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
    pub const MODULUS: WideF2x<{ Self::LIMBS }> = WideF2x::from_f2x(
        F2x::<{ Self::LIMBS }>::ONE,
        F2x::<{ Self::LIMBS }>::from_limbs([
            0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0085,
        ]),
    );

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
        let (_, rem) = prod.euclidean_div(&Self::MODULUS);
        Self::from_poly(rem.truncate())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gf2_128_mul() {
        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x77DD, 0x8AD8, 0xA825, 0x22FE, 0xBC18, 0xCF93, 0xE41B, 0x9F89,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xEB77, 0xFBFF, 0xF38A, 0x7CBD, 0xA326, 0x6CFB, 0x87AF, 0xF50A,
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xCBC8, 0xC9F4, 0x2C12, 0x3174, 0x7303, 0xF74F, 0x1F31, 0x2D9A,
        ]));
        assert_eq!(lhs.mul(&rhs), rem);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x00C6, 0x7CCF, 0xB07B, 0xE43D, 0xAD74, 0x0A18, 0x7DD8, 0xFC71,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x2F22, 0x4F9D, 0x999F, 0xEB7E, 0xC730, 0x8E79, 0xCDCA, 0x4263,
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x82D8, 0x05F2, 0xCCD1, 0x6307, 0xA462, 0xBF8F, 0xF655, 0xD97C,
        ]));
        assert_eq!(lhs.mul(&rhs), rem);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xACA8, 0x9F05, 0xB940, 0xAC95, 0x1D79, 0x6B88, 0x92C1, 0x2955,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x25C6, 0x08D5, 0xF50B, 0x41F9, 0x9E25, 0x32B5, 0x275D, 0x5D46,
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x3B02, 0xBF37, 0xCF8E, 0x6721, 0xA98B, 0x6D13, 0xCCD3, 0x4B0C,
        ]));
        assert_eq!(lhs.mul(&rhs), rem);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x4BF7, 0xEED9, 0x4B00, 0x1514, 0x3261, 0xE3B5, 0xD1E5, 0x9DE3,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xA022, 0x4E61, 0x64FF, 0x7188, 0x7042, 0x8983, 0x84E4, 0xE0B4,
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x1EBF, 0x97E8, 0x6F39, 0xDFEE, 0xE039, 0x17FF, 0x77D1, 0x993C,
        ]));
        assert_eq!(lhs.mul(&rhs), rem);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x01C9, 0xD525, 0x19C7, 0x7CBC, 0x4BC7, 0x8B00, 0x21C0, 0x4196,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x137B, 0x3993, 0x92E9, 0x81D8, 0xD12C, 0xADE0, 0x71F1, 0x1823,
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x42D1, 0xFBB2, 0x9373, 0xF45F, 0x54D9, 0x4D63, 0xC805, 0x7E78,
        ]));
        assert_eq!(lhs.mul(&rhs), rem);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x42F4, 0x54EC, 0xEAC0, 0x22E7, 0x5944, 0x2826, 0x88F5, 0x2888,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x88B9, 0x7413, 0x4B20, 0xE5EB, 0xB70D, 0xEE15, 0x6A6C, 0x1ED6,
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xF7F7, 0x0EE1, 0x8587, 0x03F8, 0xBF07, 0xE219, 0x9D25, 0x0C3D,
        ]));
        assert_eq!(lhs.mul(&rhs), rem);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xC131, 0x5E79, 0x87B5, 0xCDF6, 0x6BE4, 0x7D8F, 0x7EB5, 0xBE04,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xC490, 0xDC8E, 0x6521, 0xBAA3, 0xE9A1, 0x0756, 0x8B17, 0xDE92,
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x79E8, 0x2EB4, 0x76BE, 0xA5EF, 0xED83, 0xC5DD, 0xBB1A, 0x6EAA,
        ]));
        assert_eq!(lhs.mul(&rhs), rem);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x353A, 0x5218, 0x3E32, 0xA1DB, 0x1D1E, 0x2F91, 0x1307, 0x35CA,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x2089, 0x198C, 0xEFF2, 0x8C4E, 0x41C9, 0xE1AE, 0x72AE, 0xA2A9,
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xD4B3, 0x0D1B, 0x2646, 0x7A32, 0x8B5D, 0x3455, 0x0341, 0xEB0A,
        ]));
        assert_eq!(lhs.mul(&rhs), rem);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0xCCB4, 0xDB7C, 0x7F81, 0xCBA1, 0xAF89, 0x8090, 0xC43A, 0xDE76,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x7E48, 0x41F2, 0xF3AA, 0x620C, 0x4B52, 0x5B44, 0x3FC6, 0x0464,
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x5E9F, 0x91A2, 0xB4FB, 0x8857, 0x5C7A, 0x823F, 0x128E, 0xFFF3,
        ]));
        assert_eq!(lhs.mul(&rhs), rem);

        let lhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x8498, 0x7F6D, 0x5FFF, 0x225D, 0x21B2, 0xDBFF, 0x38D4, 0x63E8,
        ]));
        let rhs = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x94C2, 0xFB93, 0xB3ED, 0xB84C, 0x4A8B, 0xC091, 0xBC00, 0xEFCB,
        ]));
        let rem = GF2_128::from_poly(F2x::<8>::from_limbs([
            0x7037, 0x730C, 0x7BC6, 0x4447, 0xE67F, 0x9235, 0x55F5, 0xBB30,
        ]));
        assert_eq!(lhs.mul(&rhs), rem);
    }
}
