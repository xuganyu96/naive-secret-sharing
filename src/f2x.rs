//! The polynomial ring over F2[x], efficiently encoded in an integer array

pub type Word = u16;

/// Carryless multiplciation of words
/// e.g. mul(0b1111, 0b1111) = 15 * 15 = 225 = 0b11100001
///     clmul(0b1111, 0b1111) = 0b1010101
/// TODO: this is not constant time!
pub fn widening_clmul(a: Word, b: Word) -> (Word, Word) {
    let mut prod: (Word, Word) = (0, 0);

    for i in 0..(Word::BITS) {
        if ((1 << i) & b) != 0 {
            // Need to "widening left shift" a by i positions
            let high_bits = a.checked_shr(Word::BITS - i).unwrap_or(0);
            let low_bits = a.checked_shl(i).unwrap_or(0);
            prod = (prod.0 ^ high_bits, prod.1 ^ low_bits);
        }
    }

    return prod;
}

/// The degree of a polynomial is the highest power of term with non-zero coefficient
/// The degree of (x ** 4 + 1) is 4
/// The degree of 1 is 0, the degree of 0 is minus infinity
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum Degree {
    NonNegative(usize),
    NegativeInfinity,
}

impl PartialOrd for Degree {
    /// Comparison for degree; it will always return Some(_)
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        match (self, other) {
            (Self::NegativeInfinity, Self::NegativeInfinity) => Some(core::cmp::Ordering::Equal),
            (Self::NegativeInfinity, Self::NonNegative(_)) => Some(core::cmp::Ordering::Less),
            (Self::NonNegative(_), Self::NegativeInfinity) => Some(core::cmp::Ordering::Greater),
            (Self::NonNegative(a), Self::NonNegative(b)) => a.partial_cmp(b),
        }
    }
}

impl Ord for Degree {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.partial_cmp(other)
            .expect("Degree comparison should always be well-defined")
    }
}

/// An element of the polynomial ring F2[x]
///
/// Each bit encodes the a coefficient. The most significant bit encodes the coefficient of the
/// highest power term. F2x::<L> can encode a polynomial with degree up to (Word::BITS * L - 1).
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct F2x<const L: usize> {
    limbs: [Word; L],
}

impl<const L: usize> F2x<L> {
    pub const ZERO: Self = Self::zero();
    pub const ONE: Self = Self::one();
    pub const BITS: usize = (Word::BITS as usize) * L;

    pub const fn from_limbs(limbs: [Word; L]) -> Self {
        Self { limbs }
    }

    /// Get a referene to the limb at the specified location
    pub fn get_limb(&self, i: usize) -> Option<&Word> {
        self.limbs.get(i)
    }

    /// Get a mutable referene to the limb at the specified location
    pub fn get_mut_limb(&mut self, i: usize) -> Option<&mut Word> {
        self.limbs.get_mut(i)
    }

    /// The 0 polynomial
    pub const fn zero() -> Self {
        Self::from_limbs([0; L])
    }

    /// Return true if self is the 0 polynomial
    pub fn is_zero(&self) -> bool {
        self.limbs == Self::ZERO.limbs
    }

    /// The 1 polynomial
    pub const fn one() -> Self {
        let mut limbs = [0; L];
        limbs[L - 1] = 1;

        Self::from_limbs(limbs)
    }

    /// The number of leading zeros, counting from higher power terms
    /// e.g. In GF(2^128), the polynomial "1" has 127 leading zeros
    pub fn leading_zeros(&self) -> usize {
        let mut count = 0;
        let mut i = 0;

        while i < L && self.limbs[i] == 0 {
            count += Word::BITS as usize;
            i += 1;
        }
        if i < L {
            count += self.limbs[i].leading_zeros() as usize;
        }

        return count;
    }

    /// The degree of this polynomial
    pub fn degree(&self) -> Degree {
        if self.is_zero() {
            return Degree::NegativeInfinity;
        }

        return Degree::NonNegative(Self::BITS - self.leading_zeros() - 1);
    }

    /// Equivalent to applying the bitflip operator "!"
    pub fn not(&self) -> Self {
        let mut output = Self::ZERO;

        self.limbs
            .iter()
            .enumerate()
            .for_each(|(i, limb)| output.limbs[i] = !limb);

        return output;
    }

    /// Addition in GF(2^m) is a simple XOR and will never overflow
    pub fn add(&self, other: &Self) -> Self {
        let mut limbs = [0; L];

        for i in 0..L {
            // No need for bound check; guaranteed to be within bounds.
            limbs[i] = self.limbs[i] ^ other.limbs[i];
        }

        return Self::from_limbs(limbs);
    }

    /// Subtraction is identical to addition in GF(2^m) because -1 = 1
    pub fn sub(&self, other: &Self) -> Self {
        self.add(other)
    }

    /// School book multiplication with L^2 steps
    pub fn widening_mul(&self, other: &Self) -> WideF2x<L> {
        let (mut high, mut low) = (Self::ZERO, Self::ZERO);
        for i in 0..L {
            for j in 0..L {
                let (high_limb, low_limb) =
                    widening_clmul(self.limbs[L - i - 1], other.limbs[L - j - 1]);
                if (i + j) < L {
                    low.limbs[L - (i + j) - 1] ^= low_limb;
                } else {
                    high.limbs[L - (i + j - L) - 1] ^= low_limb;
                }
                if (i + j + 1) < L {
                    low.limbs[L - (i + j + 1) - 1] ^= high_limb;
                } else {
                    high.limbs[L - (i + j + 1 - L) - 1] ^= high_limb;
                }
            }
        }

        WideF2x::<L>::from_f2x(high, low)
    }

    /// Attempt to left shift (e.g. 0xFFFF.overflowing_shl(4) = 0xFFF0)
    /// If the shift will cause overflow, this method will panic. This is consistent with integer
    /// arithmetics
    pub fn shl(&self, rhs: usize) -> Self {
        let mut shifted = Self::ZERO;

        if rhs >= Self::BITS {
            panic!("attempt to shift left with overflow");
        }

        let limb_offset = rhs / (Word::BITS as usize);
        let limb_fraction = rhs % (Word::BITS as usize);
        self.limbs.iter().enumerate().for_each(|(i, limb)| {
            if limb_offset <= i {
                let near_loc = i - limb_offset;
                let near_limb = limb << limb_fraction;
                shifted.limbs[near_loc] ^= near_limb;
            }
            if limb_offset + 1 <= i {
                let far_loc = i - limb_offset - 1;
                let far_limb = if limb_fraction == 0 {
                    0
                } else {
                    limb >> (Word::BITS as usize - limb_fraction)
                };
                shifted.limbs[far_loc] ^= far_limb;
            }
        });

        shifted
    }

    /// Attempt to shift right by the specified number of bits
    /// e.g. 0xFFFF.overflowing_shr(4) = 0x0FFF
    /// Attempt to shift with overflow will cause panic
    pub fn shr(&self, rhs: usize) -> Self {
        let mut shifted = Self::ZERO;

        if rhs >= Self::BITS {
            panic!("attempt to shift right with overflow");
        }

        let limb_offset = rhs / Word::BITS as usize;
        let limb_fraction = rhs % Word::BITS as usize;
        self.limbs.iter().enumerate().for_each(|(i, limb)| {
            if (i + limb_offset) < L {
                let near_loc = i + limb_offset;
                let near_limb = limb >> limb_fraction;
                shifted.limbs[near_loc] ^= near_limb;
            }
            if (i + limb_offset + 1) < L {
                let far_loc = i + limb_offset + 1;
                let far_limb = if limb_fraction == 0 {
                    0
                } else {
                    limb << (Word::BITS as usize - limb_fraction)
                };
                shifted.limbs[far_loc] ^= far_limb;
            }
        });

        return shifted;
    }

    /// Euclidean long division, will panic if rhs is zero
    pub fn euclidean_div(&self, rhs: &Self) -> (Self, Self) {
        if rhs.is_zero() {
            panic!("attempt to divide by zero");
        }
        let mut quot = Self::ZERO;
        let mut rem = self.clone();

        while rem.degree() >= rhs.degree() {
            // Both rem and rhs are guaranteed to be NonNegative
            let rem_degree: usize = match rem.degree() {
                Degree::NonNegative(degree) => degree,
                _ => panic!("Remainder is unexpectedly zero"),
            };
            let rhs_degree: usize = match rhs.degree() {
                Degree::NonNegative(degree) => degree,
                _ => panic!("Divisor is unexpectedly zero"),
            };
            let degree_diff = rem_degree - rhs_degree;
            quot = quot.add(&Self::ONE.shl(degree_diff));
            rem = rem.sub(&rhs.shl(degree_diff));
        }

        return (quot, rem);
    }

    /// Multiplication followed by modulus reduction. This function assumes that the remainder will
    /// fit into a single F2x<L>. If not, then
    pub fn modmul(&self, rhs: &Self, modulus: &WideF2x<L>) -> Self {
        if self.is_zero() || rhs.is_zero() {
            return Self::ZERO;
        }
        // The degree of modulus must be no more than there are bits; for example, the modulus in
        // GF(2 ** 128) must have degree 128 or lower
        if modulus.degree() > Degree::NonNegative(Self::BITS) {
            panic!("modulus degree is too high");
        }
        let prod = self.widening_mul(rhs);
        let (_, rem) = prod.euclidean_div(modulus);

        rem.truncate()
    }
}

/// Wide F2[x] is useful for performing reduction after widening multiplication
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WideF2x<const L: usize> {
    /// The higher power terms
    high: F2x<L>,

    /// The lower power terms
    low: F2x<L>,
}

impl<const L: usize> WideF2x<L> {
    const BITS: usize = 2 * (Word::BITS as usize) * L;
    const ZERO: Self = Self::zero();
    const ONE: Self = Self::one();

    pub const fn from_f2x(high: F2x<L>, low: F2x<L>) -> Self {
        Self { high, low }
    }

    /// Get the limb at the specified location if it exists
    pub fn get_limb(&self, loc: usize) -> Option<&Word> {
        if loc < L {
            return self.high.limbs.get(loc);
        } else if loc < 2 * L {
            return self.low.limbs.get(loc - L);
        }
        return None;
    }

    /// Get a mutable reference to the limb at the specified location
    pub fn get_mut_limb(&mut self, loc: usize) -> Option<&mut Word> {
        if loc < L {
            return self.high.limbs.get_mut(loc);
        } else if loc < 2 * L {
            return self.low.limbs.get_mut(loc - L);
        }
        return None;
    }

    pub const fn zero() -> Self {
        Self {
            high: F2x::<L>::ZERO,
            low: F2x::<L>::ZERO,
        }
    }

    pub const fn one() -> Self {
        Self {
            high: F2x::<L>::ZERO,
            low: F2x::<L>::ONE,
        }
    }

    pub fn is_zero(&self) -> bool {
        self.high.is_zero() && self.low.is_zero()
    }

    pub fn leading_zeros(&self) -> usize {
        if self.high.is_zero() {
            self.high.leading_zeros() + self.low.leading_zeros()
        } else {
            self.high.leading_zeros()
        }
    }

    pub fn degree(&self) -> Degree {
        if self.is_zero() {
            return Degree::NegativeInfinity;
        }
        return Degree::NonNegative(Self::BITS - self.leading_zeros() - 1);
    }

    pub fn add(&self, rhs: &Self) -> Self {
        Self::from_f2x(self.high.add(&rhs.high), self.low.add(&rhs.low))
    }

    pub fn sub(&self, rhs: &Self) -> Self {
        self.add(rhs)
    }

    /// Left shift by the specified amount. Panic if rhs is greater than or equal to Self::BITS
    pub fn shl(&self, rhs: usize) -> Self {
        let mut shifted = Self::ZERO;

        if rhs >= Self::BITS {
            panic!("attempt to shift left with overflow");
        }

        let limb_offset = rhs / (Word::BITS as usize);
        let limb_fraction = rhs % (Word::BITS as usize);
        for i in 0..(2 * L) {
            let limb = self.get_limb(i).expect("unexpected out-of-bound");
            if limb_offset <= i {
                let near_loc = i - limb_offset;
                let near_limb = (*limb) << limb_fraction;
                let shifted_limb = shifted
                    .get_mut_limb(near_loc)
                    .expect("unexpected out-of-bound");
                *shifted_limb ^= near_limb;
            }
            if limb_offset + 1 <= i {
                let far_loc = i - limb_offset - 1;
                let far_limb = if limb_fraction == 0 {
                    0
                } else {
                    (*limb) >> (Word::BITS as usize - limb_fraction)
                };
                let shifted_limb = shifted
                    .get_mut_limb(far_loc)
                    .expect("unexpected out-of-bound");
                *shifted_limb ^= far_limb;
            }
        }

        shifted
    }

    pub fn shr(&self, rhs: usize) -> Self {
        let mut shifted = Self::ZERO;

        if rhs >= Self::BITS {
            panic!("attempt to shift right with overflow");
        }

        let limb_offset = rhs / Word::BITS as usize;
        let limb_fraction = rhs % Word::BITS as usize;
        for i in 0..(2 * L) {
            let limb = self.get_limb(i).expect("unexpected out-of-bound");
            if (i + limb_offset) < 2 * L {
                let near_loc = i + limb_offset;
                let near_limb = *limb >> limb_fraction;
                let shifted_limb = shifted
                    .get_mut_limb(near_loc)
                    .expect("unexpected out-of-bound");
                *shifted_limb ^= near_limb;
            }
            if (i + limb_offset + 1) < 2 * L {
                let far_loc = i + limb_offset + 1;
                let far_limb = if limb_fraction == 0 {
                    0
                } else {
                    *limb << (Word::BITS as usize - limb_fraction)
                };
                let shifted_limb = shifted
                    .get_mut_limb(far_loc)
                    .expect("unexpected out-of-bound");
                *shifted_limb ^= far_limb;
            }
        }

        return shifted;
    }

    /// Euclidean division, returning (quotient, remainder)
    /// Will panic if divisor is zero
    pub fn euclidean_div(&self, rhs: &Self) -> (Self, Self) {
        if rhs.is_zero() {
            panic!("attempt to divide by zero");
        }
        let mut quot = Self::ZERO;
        let mut rem = self.clone();

        while rem.degree() >= rhs.degree() {
            // Both rem and rhs are guaranteed to be NonNegative
            let rem_degree: usize = match rem.degree() {
                Degree::NonNegative(degree) => degree,
                _ => panic!("Remainder is unexpectedly zero"),
            };
            let rhs_degree: usize = match rhs.degree() {
                Degree::NonNegative(degree) => degree,
                _ => panic!("Divisor is unexpectedly zero"),
            };
            let degree_diff = rem_degree - rhs_degree;
            quot = quot.add(&Self::ONE.shl(degree_diff));
            rem = rem.sub(&rhs.shl(degree_diff));
        }

        return (quot, rem);
    }

    pub fn truncate(&self) -> F2x<L> {
        if !self.high.is_zero() {
            panic!("high-order limbs are not zeros");
        }
        self.low.clone()
    }
}

// TODO: what if I need to implement GF(2^12), such as in classic McEliece
#[cfg(test)]
mod tests {
    use super::*;

    type F2_128 = F2x<8>;

    #[test]
    fn test_widening_clmul() {
        assert_eq!(widening_clmul(15, 15), (0, 0b1010101));
        assert_eq!(widening_clmul(0xFFFF, 0xFFFF), (0x5555, 0x5555));
        assert_eq!(widening_clmul(0xE223, 0x672F), (0x267B, 0xB291));
        assert_eq!(widening_clmul(0, 0), (0, 0));
        assert_eq!(widening_clmul(0, 1), (0, 0));
        assert_eq!(widening_clmul(1, 0), (0, 0));
    }

    #[test]
    fn test_extfield_widening_mul() {
        assert_eq!(
            F2_128::ZERO.not().widening_mul(&F2_128::ZERO),
            WideF2x::from_f2x(F2_128::ZERO, F2_128::ZERO)
        );
        // 0xFFFF * 0xFFFF = 0x5555,0x5555
        let fives = F2_128::from_limbs([0x5555; 8]);
        assert_eq!(
            F2_128::ZERO.not().widening_mul(&F2_128::ZERO.not()),
            WideF2x::from_f2x(fives, fives)
        );

        // Random cases generated by SymPy
        let lhs = F2_128::from_limbs([
            0x3DCC, 0x5CE2, 0x8A9D, 0x3FE3, 0x5309, 0x07F3, 0xC9FD, 0x43B6,
        ]);
        let rhs = F2_128::from_limbs([
            0x8370, 0x7DA9, 0x108D, 0xF5B7, 0x30C9, 0xAEB8, 0x719A, 0xEDB5,
        ]);
        let prod = WideF2x::from_f2x(
            F2_128::from_limbs([
                0x1EAB, 0x66E7, 0x4160, 0x869E, 0xA3A7, 0x038E, 0x03AB, 0x25BF,
            ]),
            F2_128::from_limbs([
                0x77C9, 0xE332, 0x2107, 0x2707, 0x8AFD, 0x8E14, 0xE779, 0x45CE,
            ]),
        );
        assert_eq!(lhs.widening_mul(&rhs), prod);

        let lhs = F2_128::from_limbs([
            0x102D, 0x2BD4, 0x66AC, 0xBCB1, 0xF7C7, 0x5FE9, 0xBBC2, 0x335D,
        ]);
        let rhs = F2_128::from_limbs([
            0xEB90, 0xC40B, 0xFD14, 0xE019, 0xDFC5, 0xE087, 0x23EF, 0xA19F,
        ]);
        let prod = WideF2x::from_f2x(
            F2_128::from_limbs([
                0x0EA0, 0x7C5D, 0xFBA7, 0x0792, 0x1B33, 0x323D, 0xE533, 0x6BF7,
            ]),
            F2_128::from_limbs([
                0x19B6, 0xA88E, 0x62E5, 0xBB2E, 0x06AF, 0xAB14, 0x6A88, 0xE42B,
            ]),
        );
        assert_eq!(lhs.widening_mul(&rhs), prod);

        let lhs = F2_128::from_limbs([
            0xA95D, 0x01B0, 0xD0A6, 0x81A9, 0x92A5, 0xA216, 0xC971, 0x961A,
        ]);
        let rhs = F2_128::from_limbs([
            0x0817, 0x2EE5, 0xB309, 0x150F, 0x2BF1, 0x5A62, 0x2197, 0xB1C8,
        ]);
        let prod = WideF2x::from_f2x(
            F2_128::from_limbs([
                0x0543, 0x30B1, 0x8B03, 0x4A8E, 0x43F7, 0x29DB, 0x10A6, 0xBCB3,
            ]),
            F2_128::from_limbs([
                0x7D9B, 0x512D, 0x94F5, 0x12B3, 0x8D64, 0xE68F, 0xEAFD, 0xC150,
            ]),
        );
        assert_eq!(lhs.widening_mul(&rhs), prod);
    }

    #[test]
    fn test_f2x_leading_zeros() {
        assert_eq!(F2_128::ZERO.leading_zeros(), 128);
        assert_eq!(F2_128::ONE.leading_zeros(), 127);
        assert_eq!(F2_128::ZERO.not().leading_zeros(), 0);
        assert_eq!(
            F2_128::from_limbs([0x0000, 0x01B0, 0xD0A6, 0x81A9, 0x92A5, 0xA216, 0xC971, 0x961A,])
                .leading_zeros(),
            23
        );
    }

    #[test]
    fn test_f2x_degree() {
        assert_eq!(F2_128::ZERO.degree(), Degree::NegativeInfinity);
        assert_eq!(F2_128::ONE.degree(), Degree::NonNegative(0));
        assert_eq!(F2_128::ZERO.not().degree(), Degree::NonNegative(127));
        assert_eq!(
            F2_128::from_limbs([0x0000, 0x01B0, 0xD0A6, 0x81A9, 0x92A5, 0xA216, 0xC971, 0x961A,])
                .degree(),
            Degree::NonNegative(104)
        );
    }

    #[test]
    fn test_f2x_shl() {
        let expected_poly = F2_128::from_limbs([
            0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0002,
        ]);
        assert_eq!(F2_128::ONE.shl(1), expected_poly);
        assert_eq!(
            F2_128::from_limbs([0x0000, 0x01B0, 0xD0A6, 0x81A9, 0x92A5, 0xA216, 0xC971, 0x961A,])
                .shl(4),
            F2_128::from_limbs([0x0000, 0x1B0D, 0x0A68, 0x1A99, 0x2A5A, 0x216C, 0x9719, 0x61A0]),
        );
        assert_eq!(
            F2_128::from_limbs([0x0000, 0x01B0, 0xD0A6, 0x81A9, 0x92A5, 0xA216, 0xC971, 0x961A,])
                .shl(16),
            F2_128::from_limbs([0x01B0, 0xD0A6, 0x81A9, 0x92A5, 0xA216, 0xC971, 0x961A, 0x0000]),
        );
    }

    #[test]
    fn test_f2x_shr() {
        assert_eq!(
            F2_128::from_limbs([0x0000, 0x01B0, 0xD0A6, 0x81A9, 0x92A5, 0xA216, 0xC971, 0x961A,])
                .shr(4),
            F2_128::from_limbs([0x0000, 0x001B, 0x0D0A, 0x681A, 0x992A, 0x5A21, 0x6C97, 0x1961,]),
        );
        assert_eq!(
            F2_128::from_limbs([0x0000, 0x01B0, 0xD0A6, 0x81A9, 0x92A5, 0xA216, 0xC971, 0x961A,])
                .shr(16),
            F2_128::from_limbs([0x0000, 0x0000, 0x01B0, 0xD0A6, 0x81A9, 0x92A5, 0xA216, 0xC971,]),
        );
    }

    #[test]
    fn test_f2x_euclidean_div() {
        assert_eq!(
            F2_128::ZERO.euclidean_div(&F2_128::ONE),
            (F2_128::ZERO, F2_128::ZERO)
        );
        assert_eq!(
            F2_128::ONE.euclidean_div(&F2_128::ONE),
            (F2_128::ONE, F2_128::ZERO)
        );
        assert_eq!(
            F2_128::ZERO.not().euclidean_div(&F2_128::ONE),
            (F2_128::ZERO.not(), F2_128::ZERO)
        );

        // Random test cases
        assert_eq!(
            F2_128::from_limbs([0x6F88, 0x6586, 0x5701, 0x2619, 0x964B, 0x2BCD, 0x0A5E, 0xAD0C])
                .euclidean_div(&F2_128::from_limbs([
                    0x1A02, 0x3B04, 0xCBB4, 0x3729, 0x5B9B, 0x4756, 0x35A0, 0x530B
                ])),
            (
                F2_128::from_limbs([
                    0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0004
                ]),
                F2_128::from_limbs([
                    0x0780, 0x8995, 0x79D1, 0xFABC, 0xF826, 0x3695, 0xDCDF, 0xE120
                ])
            )
        );
        assert_eq!(
            F2_128::from_limbs([0x7AB3, 0xC7AC, 0x5F7B, 0xC3DC, 0x29AD, 0x137D, 0xAAD0, 0x7920])
                .euclidean_div(&F2_128::from_limbs([
                    0x4374, 0xFB75, 0x7A4D, 0x1BBC, 0xD872, 0xF253, 0x9CE8, 0x3F10
                ])),
            (
                F2_128::from_limbs([
                    0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001
                ]),
                F2_128::from_limbs([
                    0x39C7, 0x3CD9, 0x2536, 0xD860, 0xF1DF, 0xE12E, 0x3638, 0x4630
                ])
            )
        );
        assert_eq!(
            F2_128::from_limbs([0x3C84, 0xBCBF, 0xFAB8, 0xAC77, 0xDBC8, 0xA478, 0x1D91, 0xBC64])
                .euclidean_div(&F2_128::from_limbs([
                    0x0000, 0x0000, 0x0000, 0x0000, 0xF943, 0xE449, 0xAF54, 0xEA20
                ])),
            (
                F2_128::from_limbs([
                    0x0000, 0x0000, 0x0000, 0x0000, 0x474D, 0x1926, 0x9A35, 0x429A
                ]),
                F2_128::from_limbs([
                    0x0000, 0x0000, 0x0000, 0x0000, 0x4C23, 0xF643, 0x1158, 0xCB24
                ])
            )
        );
    }

    #[test]
    fn test_widef2x_degree() {
        let lhs = WideF2x::from_f2x(F2_128::ONE, F2_128::ZERO);
        assert_eq!(lhs.degree(), Degree::NonNegative(128));
        let lhs = WideF2x::from_f2x(F2_128::ZERO, F2_128::ONE);
        assert_eq!(lhs.degree(), Degree::NonNegative(0));
        let lhs = WideF2x::from_f2x(F2_128::ZERO, F2_128::ZERO);
        assert_eq!(lhs.degree(), Degree::NegativeInfinity);
    }

    #[test]
    fn test_widef2x_shl() {
        assert_eq!(
            WideF2x::from_f2x(F2_128::ZERO, F2_128::ONE.shl(127)).shl(1),
            WideF2x::from_f2x(F2_128::ONE, F2_128::ZERO)
        );
    }

    #[test]
    fn test_widef2x_shr() {
        assert_eq!(
            WideF2x::from_f2x(F2_128::ONE, F2_128::ZERO).shr(1),
            WideF2x::from_f2x(F2_128::ZERO, F2_128::ONE.shl(127)),
        );
    }

    #[test]
    fn test_widef2x_euclidean_div() {
        let lhs = F2_128::from_limbs([
            0x3DCC, 0x5CE2, 0x8A9D, 0x3FE3, 0x5309, 0x07F3, 0xC9FD, 0x43B6,
        ])
        .widening_mul(&F2_128::from_limbs([
            0x8370, 0x7DA9, 0x108D, 0xF5B7, 0x30C9, 0xAEB8, 0x719A, 0xEDB5,
        ]));
        let quot = WideF2x::from_f2x(
            F2_128::ZERO,
            F2_128::from_limbs([
                0x1EAB, 0x66E7, 0x4160, 0x869E, 0xA3A7, 0x038E, 0x03AB, 0x25BF,
            ]),
        );
        let rem = WideF2x::from_f2x(
            F2_128::ZERO,
            F2_128::from_limbs([
                0x77C9, 0xE332, 0x2107, 0x2707, 0x8AFD, 0x8E14, 0xE779, 0x45CE,
            ]),
        );
        let rhs = WideF2x::from_f2x(F2_128::ONE, F2_128::ZERO);
        assert_eq!(lhs.euclidean_div(&rhs), (quot, rem));
    }

    #[test]
    fn test_f2x_modmul() {
        let lhs = F2_128::from_limbs([
            0x1A02, 0xC5D5, 0x3035, 0x2794, 0xBAE0, 0x9A71, 0x0B95, 0xD81A,
        ]);
        let rhs = F2_128::from_limbs([
            0xC9BE, 0xD58C, 0x706B, 0x3D2B, 0xF2FE, 0x1C1B, 0xAFAA, 0x1F84,
        ]);
        let modulus = WideF2x::from_f2x(
            F2_128::ONE,
            F2_128::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0085,
            ]),
        );
        let rem = F2_128::from_limbs([
            0x2CD6, 0x23A4, 0x9F7D, 0xEA15, 0x62DD, 0x69E1, 0x63D9, 0xB940,
        ]);

        assert_eq!(lhs.modmul(&rhs, &modulus), rem);

        let lhs = F2_128::from_limbs([
            0x292C, 0xFAE0, 0x70BB, 0xC697, 0xEF2B, 0x325B, 0x3CE3, 0x1AEF,
        ]);
        let rhs = F2_128::from_limbs([
            0x212F, 0x8E8A, 0xCF1E, 0x5621, 0x12C1, 0xADE4, 0x326B, 0x3D5F,
        ]);
        let modulus = WideF2x::from_f2x(
            F2_128::ONE,
            F2_128::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0085,
            ]),
        );
        let rem = F2_128::from_limbs([
            0x248A, 0xC6E8, 0xEF4F, 0xD895, 0x5C44, 0xC8A9, 0x7F71, 0x3518,
        ]);

        assert_eq!(lhs.modmul(&rhs, &modulus), rem);
    }
}
