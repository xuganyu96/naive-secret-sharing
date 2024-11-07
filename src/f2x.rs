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

    /// The degree of this polynomial. See `Degree` for the mathematical definition.
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

    /// Apply bitwise XOR
    pub fn xor(&self, other: &Self) -> Self {
        let mut limbs = [0; L];

        for i in 0..L {
            // No need for bound check; guaranteed to be within bounds.
            limbs[i] = self.limbs[i] ^ other.limbs[i];
        }

        return Self::from_limbs(limbs);
    }

    /// Addition in GF(2^m) is a simple XOR and will never overflow
    pub fn add(&self, other: &Self) -> Self {
        self.xor(other)
    }

    /// Subtraction is identical to addition in GF(2^m) because -1 = 1
    pub fn sub(&self, other: &Self) -> Self {
        self.xor(other)
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

    /// Use [Euclid's algorithm](https://en.wikipedia.org/wiki/Euclidean_algorithm) to compute the
    /// highest-degree polynomial that divides both lhs and rhs.
    pub fn gcd(lhs: &Self, rhs: &Self) -> Self {
        let (mut a, mut b): (Self, Self) = (lhs.clone(), rhs.clone());

        while !b.is_zero() {
            let (_, rem) = a.euclidean_div(&b);
            (a, b) = (b, rem);
        }

        a
    }

    /// Use [Extended Euclid's algorithm](https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm)
    /// to compute (s, t, d) such that s * lhs + t * rhs = d and d is the GCD between lhs and rhs
    pub fn xgcd(lhs: &Self, rhs: &Self) -> (Self, Self, Self) {
        let (mut rr, mut r): (Self, Self) = (lhs.clone(), rhs.clone());
        let (mut ss, mut s): (Self, Self) = (Self::ONE, Self::ZERO);
        let (mut tt, mut t): (Self, Self) = (Self::ZERO, Self::ONE);

        while !r.is_zero() {
            let (quot, rem) = rr.euclidean_div(&r);
            (rr, r) = (r, rem);
            (ss, s) = (s, ss.sub(&quot.widening_mul(&s).truncate()));
            (tt, t) = (t, tt.sub(&quot.widening_mul(&t).truncate()));
        }

        (ss, tt, rr)
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

    /// Use [Euclid's algorithm](https://en.wikipedia.org/wiki/Euclidean_algorithm) to compute the
    /// highest-degree polynomial that divides both lhs and rhs.
    pub fn gcd(lhs: &Self, rhs: &Self) -> Self {
        let (mut a, mut b): (Self, Self) = (lhs.clone(), rhs.clone());

        while !b.is_zero() {
            let (_, rem) = a.euclidean_div(&b);
            (a, b) = (b, rem);
        }

        a
    }

    // TODO: we don't have a non-widening multiplication
    pub fn xgcd(lhs: &Self, rhs: &Self) -> (Self, Self, Self) {
        todo!("first implement non-widening multiplication");
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

    /// Randomly generated GCD from Python
    #[test]
    fn random_f2x_gcd() {
        assert_eq!(
            F2_128::gcd(
                &F2_128::from_limbs([
                    0x88D3, 0xC663, 0x66EC, 0x77EC, 0xF526, 0x6510, 0x5C19, 0x6517
                ]),
                &F2_128::from_limbs([
                    0x23E6, 0x86AB, 0x9E60, 0x2E74, 0x6BE2, 0x87A8, 0x3D6B, 0x651D
                ]),
            ),
            F2_128::from_limbs([0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x000D])
        );

        assert_eq!(
            F2_128::gcd(
                &F2_128::from_limbs([
                    0xF0D9, 0x2027, 0xBC65, 0x0F0D, 0xE5EB, 0xD78D, 0x967D, 0x489C
                ]),
                &F2_128::from_limbs([
                    0x8F6A, 0xD588, 0xC423, 0xA33D, 0x6819, 0x6C4E, 0x676C, 0xE06F
                ]),
            ),
            F2_128::from_limbs([0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001])
        );

        assert_eq!(
            F2_128::gcd(
                &F2_128::from_limbs([
                    0xB3D5, 0x9549, 0xCF5D, 0x897A, 0xC79E, 0xC029, 0xB84B, 0xAC6A
                ]),
                &F2_128::from_limbs([
                    0xA495, 0xDB48, 0x0E91, 0x2ED4, 0x6B75, 0xF726, 0xBC16, 0xF311
                ]),
            ),
            F2_128::from_limbs([0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0007])
        );

        assert_eq!(
            F2_128::gcd(
                &F2_128::from_limbs([
                    0xC6D2, 0xA56B, 0xFEBC, 0x0BC0, 0xFAB8, 0x54AC, 0xF72B, 0x8EE6
                ]),
                &F2_128::from_limbs([
                    0x1585, 0xD7E5, 0x319B, 0x34E5, 0x0628, 0xD2B8, 0x84CB, 0x9DB3
                ]),
            ),
            F2_128::from_limbs([0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x000B])
        );

        assert_eq!(
            F2_128::gcd(
                &F2_128::from_limbs([
                    0xA2D7, 0xB6DF, 0xECE2, 0xF5F9, 0x9358, 0x001F, 0xC58B, 0x6998
                ]),
                &F2_128::from_limbs([
                    0x53B3, 0xDE89, 0x8495, 0x65F2, 0x2746, 0x2197, 0xF043, 0x0E20
                ]),
            ),
            F2_128::from_limbs([0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0308])
        );

        assert_eq!(
            F2_128::gcd(
                &F2_128::from_limbs([
                    0xC01A, 0x9409, 0xDF5B, 0x988E, 0x427B, 0xA208, 0x6DEE, 0x5082
                ]),
                &F2_128::from_limbs([
                    0x80C7, 0xDA18, 0x255E, 0xC295, 0x4796, 0x263A, 0x0263, 0x1603
                ]),
            ),
            F2_128::from_limbs([0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001])
        );

        assert_eq!(
            F2_128::gcd(
                &F2_128::from_limbs([
                    0x1BA5, 0x59C7, 0x51A8, 0x8633, 0x6996, 0xDDE2, 0xFD11, 0xBC10
                ]),
                &F2_128::from_limbs([
                    0x3873, 0x9182, 0xEAC5, 0x205D, 0x7D38, 0x04AE, 0x972C, 0x906C
                ]),
            ),
            F2_128::from_limbs([0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0004])
        );

        assert_eq!(
            F2_128::gcd(
                &F2_128::from_limbs([
                    0x478F, 0x4705, 0x220D, 0xAC34, 0x171F, 0xAF89, 0x15AC, 0x1700
                ]),
                &F2_128::from_limbs([
                    0xC9C4, 0xCE2D, 0x7266, 0xC78B, 0x3DBB, 0x589C, 0xDCDF, 0x8929
                ]),
            ),
            F2_128::from_limbs([0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001])
        );

        assert_eq!(
            F2_128::gcd(
                &F2_128::from_limbs([
                    0x6480, 0x2AA3, 0x2F8A, 0x8B81, 0x1371, 0x75ED, 0x7B98, 0x17D9
                ]),
                &F2_128::from_limbs([
                    0xF899, 0x9ABB, 0xD6F7, 0x2998, 0xB5BF, 0xBB6D, 0xF0EC, 0xABA9
                ]),
            ),
            F2_128::from_limbs([0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001])
        );

        assert_eq!(
            F2_128::gcd(
                &F2_128::from_limbs([
                    0xD641, 0x00B5, 0xC808, 0x3961, 0x2CD7, 0xD722, 0xD9C9, 0x3C09
                ]),
                &F2_128::from_limbs([
                    0x2E55, 0xB6A1, 0x3FB4, 0x2C56, 0x5E1C, 0xEE4D, 0xA7A2, 0x6130
                ]),
            ),
            F2_128::from_limbs([0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001])
        );
    }

    #[test]
    fn random_f2x_xgcd() {
        let lhs = F2_128::from_limbs([
            0x4746, 0x4DE9, 0xB867, 0x3F15, 0xDB58, 0xBFE4, 0xAE50, 0x9F35,
        ]);
        let rhs = F2_128::from_limbs([
            0xCD9E, 0x5830, 0x5EDB, 0x11F2, 0x9232, 0x419A, 0xB9E3, 0xA334,
        ]);
        let expected_s = F2_128::from_limbs([
            0x0DA5, 0x30C5, 0x52A7, 0x9083, 0x5354, 0x1B1F, 0xCA7C, 0xE235,
        ]);
        let expected_t = F2_128::from_limbs([
            0x04A4, 0xF401, 0x528E, 0xB496, 0xE555, 0x39AE, 0xB7A9, 0xE274,
        ]);
        let expected_d = F2_128::from_limbs([
            0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001,
        ]);
        let (s, t, divisor) = F2_128::xgcd(&lhs, &rhs);
        assert_eq!((s, t, divisor), (expected_s, expected_t, expected_d));

        let lhs = F2_128::from_limbs([
            0xDC28, 0x1585, 0x36FC, 0x9CE3, 0xE075, 0x7C5B, 0x0B0F, 0xEE28,
        ]);
        let rhs = F2_128::from_limbs([
            0x8456, 0xF362, 0x3C86, 0x0D2F, 0x2A71, 0x2C57, 0x42EA, 0xBEF3,
        ]);
        let expected_s = F2_128::from_limbs([
            0x7760, 0xB406, 0xE230, 0x3311, 0x2E48, 0x971D, 0xE703, 0xCA4E,
        ]);
        let expected_t = F2_128::from_limbs([
            0x44E5, 0xDF01, 0xAFCB, 0x32F7, 0x610A, 0x6AA8, 0x1E8E, 0xA15F,
        ]);
        let expected_d = F2_128::from_limbs([
            0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001,
        ]);
        let (s, t, divisor) = F2_128::xgcd(&lhs, &rhs);
        assert_eq!((s, t, divisor), (expected_s, expected_t, expected_d));

        let lhs = F2_128::from_limbs([
            0x0A0B, 0xA3AD, 0x178E, 0x4E43, 0xA028, 0x0178, 0xE8C1, 0xBA6E,
        ]);
        let rhs = F2_128::from_limbs([
            0xC3FA, 0x7569, 0x4509, 0xC88A, 0xA16B, 0x0121, 0x2A00, 0xD254,
        ]);
        let expected_s = F2_128::from_limbs([
            0x13A6, 0xE946, 0xFF13, 0x6807, 0x329C, 0x936D, 0xDEF4, 0x8027,
        ]);
        let expected_t = F2_128::from_limbs([
            0x01A1, 0x86C5, 0x07BB, 0x927E, 0x9837, 0x8143, 0xB62B, 0xAFFA,
        ]);
        let expected_d = F2_128::from_limbs([
            0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0002,
        ]);
        let (s, t, divisor) = F2_128::xgcd(&lhs, &rhs);
        assert_eq!((s, t, divisor), (expected_s, expected_t, expected_d));

        let lhs = F2_128::from_limbs([
            0xAE40, 0x001C, 0x1F77, 0x6239, 0x1B8D, 0xDDBD, 0x4010, 0x571B,
        ]);
        let rhs = F2_128::from_limbs([
            0x3F65, 0x16D7, 0x80C3, 0x4D2F, 0xE702, 0xA2AC, 0xF894, 0xF8B2,
        ]);
        let expected_s = F2_128::from_limbs([
            0x0ECE, 0xC430, 0xA308, 0x7B8A, 0x61B4, 0x9431, 0xE44F, 0x9D61,
        ]);
        let expected_t = F2_128::from_limbs([
            0x2C21, 0xF0E7, 0xEA61, 0x023D, 0x05AC, 0x45D7, 0x6AB9, 0x09FC,
        ]);
        let expected_d = F2_128::from_limbs([
            0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0003,
        ]);
        let (s, t, divisor) = F2_128::xgcd(&lhs, &rhs);
        assert_eq!((s, t, divisor), (expected_s, expected_t, expected_d));

        let lhs = F2_128::from_limbs([
            0x70C3, 0x2E36, 0xC1F8, 0x0858, 0x9668, 0x4B8B, 0x7DC7, 0x8C0B,
        ]);
        let rhs = F2_128::from_limbs([
            0x8B91, 0x441A, 0x1C82, 0x83F5, 0x0BC4, 0x499A, 0x9BAE, 0xEC2E,
        ]);
        let expected_s = F2_128::from_limbs([
            0x1AB7, 0x6742, 0x08B8, 0x4488, 0x159F, 0xC45A, 0x8A1F, 0x07BF,
        ]);
        let expected_t = F2_128::from_limbs([
            0x080E, 0xE174, 0xA328, 0x94B8, 0x8266, 0x4093, 0xC73B, 0x9BD9,
        ]);
        let expected_d = F2_128::from_limbs([
            0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0007,
        ]);
        let (s, t, divisor) = F2_128::xgcd(&lhs, &rhs);
        assert_eq!((s, t, divisor), (expected_s, expected_t, expected_d));

        let lhs = F2_128::from_limbs([
            0x94DF, 0xEE78, 0x4A39, 0x9065, 0x4ECE, 0xE5CC, 0x57D4, 0x5C90,
        ]);
        let rhs = F2_128::from_limbs([
            0xE79F, 0xE6A0, 0xAA2C, 0x5817, 0x68A9, 0x2F9F, 0x5660, 0x73EF,
        ]);
        let expected_s = F2_128::from_limbs([
            0x277D, 0x4C20, 0x7788, 0xBB9F, 0x5912, 0x8E66, 0xD790, 0x5957,
        ]);
        let expected_t = F2_128::from_limbs([
            0x34F2, 0xE96F, 0x5AD0, 0x5CC2, 0x37FB, 0xFA75, 0xFE74, 0x1F43,
        ]);
        let expected_d = F2_128::from_limbs([
            0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001,
        ]);
        let (s, t, divisor) = F2_128::xgcd(&lhs, &rhs);
        assert_eq!((s, t, divisor), (expected_s, expected_t, expected_d));

        let lhs = F2_128::from_limbs([
            0x95DF, 0x4A39, 0x5AEA, 0xD3E9, 0x3B31, 0xCF9E, 0xA082, 0x0C02,
        ]);
        let rhs = F2_128::from_limbs([
            0x60F3, 0x0B74, 0xB7AB, 0x8B48, 0x9D1F, 0x693B, 0xEFE8, 0xAFAC,
        ]);
        let expected_s = F2_128::from_limbs([
            0x1A1D, 0xBED7, 0x7BD9, 0x99BE, 0x5491, 0x182E, 0x7AB1, 0x3DF1,
        ]);
        let expected_t = F2_128::from_limbs([
            0x223E, 0x1D1E, 0x8229, 0xF29F, 0xF020, 0x2121, 0x2CEB, 0x1B68,
        ]);
        let expected_d = F2_128::from_limbs([
            0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0002,
        ]);
        let (s, t, divisor) = F2_128::xgcd(&lhs, &rhs);
        assert_eq!((s, t, divisor), (expected_s, expected_t, expected_d));

        let lhs = F2_128::from_limbs([
            0xB2AC, 0x5425, 0xD230, 0x1A5B, 0x847A, 0xF133, 0x37BD, 0xC23B,
        ]);
        let rhs = F2_128::from_limbs([
            0x5FD9, 0xF185, 0x1193, 0x6CD9, 0x0277, 0xEDAA, 0xD520, 0xCDBE,
        ]);
        let expected_s = F2_128::from_limbs([
            0x31CA, 0x7C63, 0x8B6E, 0x9028, 0xF895, 0x9009, 0xEA7D, 0xDB1D,
        ]);
        let expected_t = F2_128::from_limbs([
            0x6784, 0x1600, 0x6CD1, 0xDC98, 0x9808, 0x1058, 0x18C9, 0x7801,
        ]);
        let expected_d = F2_128::from_limbs([
            0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001,
        ]);
        let (s, t, divisor) = F2_128::xgcd(&lhs, &rhs);
        assert_eq!((s, t, divisor), (expected_s, expected_t, expected_d));

        let lhs = F2_128::from_limbs([
            0x4253, 0x5575, 0x05C3, 0x2E2F, 0x47BB, 0xDD26, 0x156C, 0x9397,
        ]);
        let rhs = F2_128::from_limbs([
            0x65E3, 0xA49E, 0x4049, 0xA7E7, 0xA96D, 0xC55C, 0xB7F2, 0x975F,
        ]);
        let expected_s = F2_128::from_limbs([
            0x01D3, 0x574C, 0xE524, 0x0E33, 0xA6C1, 0xB763, 0x65BD, 0xCD90,
        ]);
        let expected_t = F2_128::from_limbs([
            0x0176, 0xAB05, 0x0265, 0x4A3F, 0x6D58, 0xEF03, 0xA981, 0xE445,
        ]);
        let expected_d = F2_128::from_limbs([
            0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0013,
        ]);
        let (s, t, divisor) = F2_128::xgcd(&lhs, &rhs);
        assert_eq!((s, t, divisor), (expected_s, expected_t, expected_d));

        let lhs = F2_128::from_limbs([
            0x2E46, 0x1AB1, 0x8719, 0xB670, 0x49D4, 0x3EE3, 0xD927, 0x019E,
        ]);
        let rhs = F2_128::from_limbs([
            0x3FA7, 0xAA63, 0xE426, 0x84EA, 0x03BB, 0xB6F5, 0x8E43, 0xF797,
        ]);
        let expected_s = F2_128::from_limbs([
            0x1F46, 0xF302, 0x24E0, 0x27EA, 0xFADB, 0x5FEA, 0x37B8, 0xA278,
        ]);
        let expected_t = F2_128::from_limbs([
            0x17DD, 0x9394, 0xD239, 0x14F8, 0x0868, 0xD45D, 0x9347, 0x217B,
        ]);
        let expected_d = F2_128::from_limbs([
            0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001,
        ]);
        let (s, t, divisor) = F2_128::xgcd(&lhs, &rhs);
        assert_eq!((s, t, divisor), (expected_s, expected_t, expected_d));
    }

    #[test]
    fn random_widef2x_gcd() {
        let lhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0xA84F, 0xED17, 0xFF3E, 0x9B60, 0x2EA8, 0x3480, 0xC4E6, 0x6981,
            ]),
            F2x::<8>::from_limbs([
                0x8174, 0x57C6, 0x54ED, 0x250C, 0xC8A6, 0x3D3B, 0x6F64, 0x7A27,
            ]),
        );
        let rhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x7139, 0x9134, 0x3CAA, 0x16DD, 0xD47E, 0x83CF, 0x7006, 0xF923,
            ]),
            F2x::<8>::from_limbs([
                0x203A, 0xEF7E, 0x0745, 0x6DC9, 0x003B, 0xF271, 0xF403, 0xCF9F,
            ]),
        );
        let gcd = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
            ]),
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001,
            ]),
        );
        assert_eq!(WideF2x::<8>::gcd(&lhs, &rhs), gcd);

        let lhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0xEC10, 0x1284, 0x0415, 0x8B0E, 0xD842, 0x5F1C, 0xCCC3, 0x2903,
            ]),
            F2x::<8>::from_limbs([
                0xBC0E, 0xFBFC, 0x0430, 0x7C73, 0x40AB, 0x5717, 0x3F59, 0xCBA5,
            ]),
        );
        let rhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x5C85, 0xE20D, 0x1371, 0x9C06, 0x02C2, 0x7CD7, 0xCF41, 0xD5BB,
            ]),
            F2x::<8>::from_limbs([
                0xEFC2, 0x1B46, 0xCBCD, 0x1887, 0xCC3E, 0x08D7, 0xBC57, 0x429D,
            ]),
        );
        let gcd = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
            ]),
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001,
            ]),
        );
        assert_eq!(WideF2x::<8>::gcd(&lhs, &rhs), gcd);

        let lhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0xACC6, 0x7E82, 0x1082, 0x4A1A, 0xE365, 0x2B0B, 0x5092, 0x64C8,
            ]),
            F2x::<8>::from_limbs([
                0x9A5D, 0x8B28, 0xAFC0, 0x78A5, 0x995A, 0x724B, 0x73D2, 0xC1C0,
            ]),
        );
        let rhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x9BDF, 0x3946, 0x97B7, 0x9D84, 0x9DBF, 0x79C9, 0x9317, 0x25B9,
            ]),
            F2x::<8>::from_limbs([
                0xC8CB, 0x2F7F, 0xD8D6, 0xFF02, 0xD357, 0x83D6, 0x2A59, 0xCBED,
            ]),
        );
        let gcd = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
            ]),
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001,
            ]),
        );
        assert_eq!(WideF2x::<8>::gcd(&lhs, &rhs), gcd);

        let lhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x2A28, 0x37B7, 0xB61E, 0x8DE7, 0x53C3, 0xAA1F, 0x293B, 0xE8AA,
            ]),
            F2x::<8>::from_limbs([
                0xBA21, 0x04EE, 0x9E59, 0x90EA, 0x0DB1, 0xA340, 0x6193, 0xD7CD,
            ]),
        );
        let rhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0xE463, 0x29E4, 0xCF47, 0x6F8F, 0x1C3F, 0x06C2, 0x133B, 0xB016,
            ]),
            F2x::<8>::from_limbs([
                0x2AF6, 0xD239, 0xC37D, 0xC04F, 0x12BB, 0x7879, 0x1DBA, 0x1D84,
            ]),
        );
        let gcd = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
            ]),
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0003,
            ]),
        );
        assert_eq!(WideF2x::<8>::gcd(&lhs, &rhs), gcd);

        let lhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x15AB, 0x501E, 0xB2A5, 0x382C, 0x98A4, 0x42F9, 0xC3BC, 0x9003,
            ]),
            F2x::<8>::from_limbs([
                0x1749, 0x6784, 0xE98B, 0x35A8, 0xC618, 0x1915, 0x625C, 0x8ACF,
            ]),
        );
        let rhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0xBB17, 0x7D58, 0x7663, 0x6B7B, 0x05EE, 0x0271, 0x3316, 0x0462,
            ]),
            F2x::<8>::from_limbs([
                0x9E0A, 0x2111, 0xA529, 0x1CE7, 0x4E60, 0x8467, 0xBA0C, 0x01AE,
            ]),
        );
        let gcd = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
            ]),
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001,
            ]),
        );
        assert_eq!(WideF2x::<8>::gcd(&lhs, &rhs), gcd);

        let lhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x58B4, 0x0CF0, 0x4912, 0x6384, 0xFAC8, 0x9F45, 0x6DA5, 0x7FF5,
            ]),
            F2x::<8>::from_limbs([
                0x88B5, 0xE58C, 0x412A, 0x54B3, 0xE38D, 0xF939, 0xBAEE, 0x9242,
            ]),
        );
        let rhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x0AAD, 0xBF23, 0x3AAB, 0x5B91, 0x27A5, 0xE9CD, 0xD58C, 0x4D15,
            ]),
            F2x::<8>::from_limbs([
                0xE926, 0xFCF8, 0x7674, 0x2312, 0xDCA7, 0xE614, 0x1191, 0x4412,
            ]),
        );
        let gcd = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
            ]),
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0002,
            ]),
        );
        assert_eq!(WideF2x::<8>::gcd(&lhs, &rhs), gcd);

        let lhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x8A8D, 0x0230, 0x0EE1, 0x50AB, 0x474C, 0x9828, 0xF8A8, 0x41AE,
            ]),
            F2x::<8>::from_limbs([
                0x0C45, 0x2A82, 0x9E81, 0x0FD4, 0x98D9, 0x8CF1, 0x9889, 0xD500,
            ]),
        );
        let rhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0xE0BB, 0x71DD, 0xEB35, 0xD698, 0xEF9D, 0x86DB, 0xABAC, 0x8FCF,
            ]),
            F2x::<8>::from_limbs([
                0x6C0D, 0x748C, 0x6EE4, 0x64C8, 0xE6D1, 0xD8D8, 0xC366, 0x3F66,
            ]),
        );
        let gcd = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
            ]),
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0002,
            ]),
        );
        assert_eq!(WideF2x::<8>::gcd(&lhs, &rhs), gcd);

        let lhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x1AA4, 0x210A, 0xA6E5, 0x8D2E, 0x33B9, 0x803D, 0x1E66, 0x033D,
            ]),
            F2x::<8>::from_limbs([
                0x7413, 0x3B1E, 0x288D, 0x644C, 0x305F, 0xFA3D, 0x72BD, 0x5066,
            ]),
        );
        let rhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0xE651, 0x233D, 0xF98D, 0xDDDB, 0x718F, 0x4E9B, 0x0BFC, 0x00CC,
            ]),
            F2x::<8>::from_limbs([
                0xFBEB, 0x8F7B, 0x4B61, 0x16EF, 0x5585, 0x195B, 0xA4B2, 0x5693,
            ]),
        );
        let gcd = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
            ]),
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0003,
            ]),
        );
        assert_eq!(WideF2x::<8>::gcd(&lhs, &rhs), gcd);

        let lhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0xCDDF, 0x1FC0, 0xA77B, 0x92ED, 0x195C, 0xBD3E, 0x1A50, 0x715F,
            ]),
            F2x::<8>::from_limbs([
                0x0D74, 0x3B00, 0x6B22, 0x4722, 0x92EE, 0xA275, 0x08BD, 0xD36B,
            ]),
        );
        let rhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0xAB9D, 0x6D37, 0xC731, 0xDD1F, 0x6BEC, 0x0DD6, 0xC930, 0x8821,
            ]),
            F2x::<8>::from_limbs([
                0x858E, 0x9FA6, 0x1BE4, 0x231D, 0xE4C3, 0x10E8, 0xF6AE, 0x71BC,
            ]),
        );
        let gcd = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
            ]),
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0001,
            ]),
        );
        assert_eq!(WideF2x::<8>::gcd(&lhs, &rhs), gcd);

        let lhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0xD88A, 0xF734, 0x5A39, 0xA5F3, 0x7F9C, 0xB4D1, 0x06BF, 0x452A,
            ]),
            F2x::<8>::from_limbs([
                0xE550, 0xC758, 0x63DC, 0x431F, 0x0D26, 0xF96B, 0x5C4B, 0xC30E,
            ]),
        );
        let rhs = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x1D9D, 0x91C4, 0x8EA4, 0xBBDA, 0x8EF7, 0x9270, 0x50D5, 0x6712,
            ]),
            F2x::<8>::from_limbs([
                0xA228, 0xBD5A, 0x6662, 0x28F3, 0x6322, 0x9D24, 0x4E86, 0x34C6,
            ]),
        );
        let gcd = WideF2x::<8>::from_f2x(
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
            ]),
            F2x::<8>::from_limbs([
                0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0002,
            ]),
        );
        assert_eq!(WideF2x::<8>::gcd(&lhs, &rhs), gcd);
    }
}
