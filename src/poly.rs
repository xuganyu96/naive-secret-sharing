//! General polynomials
use crate::f2x::Degree;
use crate::galoisfields::{FieldArithmetic, GF2p256};
use rand::Rng;

/// The canonical representation of a polynomial using its coefficients
/// Coefficients are organized in little-endian order: the value at lower index encodes the
/// coefficient of a lower-power term; coeffs[0] encodes the constant term
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Poly<E: FieldArithmetic> {
    pub coeffs: Vec<E>,
}

impl<E: FieldArithmetic> Poly<E> {
    pub fn from_coeffs(coeffs: Vec<E>) -> Self {
        Self { coeffs }
    }

    pub fn bytes(&self) -> usize {
        self.capacity() * E::bytes()
    }

    pub fn serialize(&self, dst: &mut [u8]) {
        if self.bytes() != dst.len() {
            panic!("buffer size is incorrect");
        }
        self.coeffs.iter().enumerate().for_each(|(i, coeff)| {
            let start = E::bytes() * i;
            let stop = start + E::bytes();
            coeff.write_be_bytes(&mut dst[start..stop]);
        })
    }

    pub fn deserialize(src: &[u8], capacity: usize) -> Self {
        if capacity * E::bytes() != src.len() {
            panic!("buffer size is incorrect");
        }
        let mut output = Self::zero_with_capacity(capacity);

        output.coeffs.iter_mut().enumerate().for_each(|(i, coeff)| {
            let start = E::bytes() * i;
            let stop = start + E::bytes();
            *coeff = E::from_be_bytes(&src[start..stop]);
        });

        output
    }

    /// The length of the coefficient vector.
    ///
    /// The capacity of a polynomial is determined when it is created and will not change. Each
    /// polynomial will only iteract with other polynomials of the same length, or the iteraction
    /// will panic.
    pub fn capacity(&self) -> usize {
        self.coeffs.len()
    }

    /// Return a zero polynomial with the specified number of coefficients
    pub fn zero_with_capacity(capacity: usize) -> Self {
        let coeffs = vec![E::zero(); capacity];
        Self { coeffs }
    }

    /// Fill self with random elements from the finite field
    pub fn fill_random(&mut self) {
        self.coeffs.iter_mut().for_each(|coeff| {
            *coeff = E::random();
        });
    }

    /// Fill self with random elements according to the input RNG
    pub fn fill_random_with_rng(&mut self, rng: &mut impl Rng) {
        self.coeffs.iter_mut().for_each(|coeff| {
            *coeff = E::random_with_rng(rng);
        });
    }

    /// Return true if self encodes a zero polynomial. An empty vector encodes a zero polynomial
    pub fn is_zero(&self) -> bool {
        for coeff in self.coeffs.iter() {
            if !coeff.is_zero() {
                return false;
            }
        }
        return true;
    }

    /// The degree of a polynomial is the power of the highest-power term with a non-zero
    /// coefficient. The degree of 0 is minus infinity.
    pub fn degree(&self) -> Degree {
        if self.is_zero() {
            return Degree::NegativeInfinity;
        }
        let degree = self
            .coeffs
            .iter()
            .enumerate()
            .filter_map(|(i, coeff)| {
                if coeff.is_zero() {
                    return None;
                } else {
                    return Some(i);
                }
            })
            .max()
            .expect("Non-zero polynomial has no non-zero coefficient");
        return Degree::NonNegative(degree);
    }

    pub fn add(&self, rhs: &Self) -> Self {
        if self.capacity() != rhs.capacity() {
            panic!("Polynomial capacities do not match");
        }
        let sum_coeffs = self
            .coeffs
            .iter()
            .zip(rhs.coeffs.iter())
            .map(|(a, b)| a.modadd(b))
            .collect::<Vec<E>>();

        Self::from_coeffs(sum_coeffs)
    }

    pub fn sub(&self, rhs: &Self) -> Self {
        if self.capacity() != rhs.capacity() {
            panic!("Polynomial capacities do not match");
        }
        let sum_coeffs = self
            .coeffs
            .iter()
            .zip(rhs.coeffs.iter())
            .map(|(a, b)| a.modmul(b))
            .collect::<Vec<E>>();

        Self::from_coeffs(sum_coeffs)
    }

    /// School-book multiplication. Will panic if the product will overflow.
    pub fn mul(&self, rhs: &Self) -> Self {
        if self.capacity() != rhs.capacity() {
            panic!("Polynomial capacities do not match");
        }
        let mut prod = Self::zero_with_capacity(self.capacity());

        for (i, a) in self.coeffs.iter().enumerate() {
            for (j, b) in rhs.coeffs.iter().enumerate() {
                let c = a.modmul(b);
                if !c.is_zero() {
                    if (i + j) >= self.capacity() {
                        panic!("attempt to multiply polynomials with overflow");
                    } else {
                        prod.coeffs[i + j] = prod.coeffs[i + j].modadd(&c);
                    }
                }
            }
        }

        prod
    }

    /// Multiplication by a single scaler. Avoids unnecessary convolution.
    pub fn mul_coeff(&self, rhs: &E) -> Self {
        let coeffs = self
            .coeffs
            .iter()
            .map(|coeff| coeff.modmul(rhs))
            .collect::<Vec<E>>();

        Self::from_coeffs(coeffs)
    }

    /// Evaluate the polynomial at the given points
    pub fn evaluate(&self, at: &E) -> E {
        let mut acc = E::zero();
        let mut indeterminate = E::one();

        for i in 0..self.capacity() {
            acc = acc.modadd(&self.coeffs[i].modmul(&indeterminate));
            indeterminate = indeterminate.modmul(at);
        }

        return acc;
    }

    /// Compute the Lagrange polynomial on the points (x, f(x)).
    /// Points are assumed to be distinct. The behavior of interpoate is undefined if there are
    /// duplicate points.
    pub fn interpolate(points: &[(E, E)], capacity: usize) -> Self {
        if points.len() > capacity {
            panic!("Cannot interpolate more points than capacity");
        }
        let mut lagrange = Self::zero_with_capacity(capacity);

        for (alpha_i, r) in points {
            let mut basis = Self::zero_with_capacity(capacity);
            basis.coeffs[0] = E::one();

            for (alpha_j, _) in points {
                if alpha_j != alpha_i {
                    let mut factor = Self::zero_with_capacity(capacity);
                    // factor is (x - alpha_j)/(alpha_i - alpha_j)
                    factor.coeffs[0] = *alpha_j;
                    factor.coeffs[1] = E::one();
                    factor = factor.mul_coeff(&alpha_i.modsub(alpha_j).modinv().unwrap());
                    basis = basis.mul(&factor);
                }
            }

            lagrange = lagrange.add(&basis.mul_coeff(r));
        }

        lagrange
    }
}

pub type Poly256 = Poly<GF2p256>;
/// A point contains two values (x, f(x))
pub type Poly256Point = (GF2p256, GF2p256);

#[cfg(test)]
mod tests {
    use super::{Poly, Poly256};
    use crate::f2x::Degree;
    use crate::galoisfields::{FieldArithmetic, GF2p256, F3329};

    type Poly3329 = Poly<F3329>;

    #[test]
    fn sanity() {
        let cap = 10;
        let mut poly = Poly256::zero_with_capacity(cap);
        poly.fill_random();
        let mut one = Poly256::zero_with_capacity(cap);
        one.coeffs[0] = GF2p256::ONE;
        assert_eq!(poly.add(&poly), Poly256::zero_with_capacity(cap));
        assert!(Poly256::zero_with_capacity(cap)
            .evaluate(&GF2p256::random())
            .is_zero());
        assert_eq!(poly.mul(&one), poly);
    }

    #[test]
    fn poly3329_eval() {
        let mut poly = Poly3329::zero_with_capacity(3);
        poly.coeffs[0] = F3329::from(536);
        poly.coeffs[1] = F3329::from(49);
        poly.coeffs[2] = F3329::from(1873);
        assert_eq!(poly.evaluate(&F3329::from(1)), F3329::from(2458));
        assert_eq!(poly.evaluate(&F3329::from(2)), F3329::from(1468));
        assert_eq!(poly.evaluate(&F3329::from(3)), F3329::from(895));
    }

    #[test]
    fn interpolate_random_polynomial() {
        let cap = 3;
        let mut poly = Poly256::zero_with_capacity(cap);
        poly.fill_random();
        let points = (0..cap)
            .map(|_| {
                // Does not check distinctness of field ordering, but the chance of collision is
                // cryptographically small
                let alpha = GF2p256::random();
                let r = poly.evaluate(&alpha);
                (alpha, r)
            })
            .collect::<Vec<(GF2p256, GF2p256)>>();
        let interpolate = Poly256::interpolate(&points, cap);
        assert_eq!(interpolate.degree(), Degree::NonNegative(cap - 1));
        assert_eq!(poly, Poly256::interpolate(&points, cap));
    }

    #[test]
    fn random_poly256_serde() {
        let cap = 10;
        let mut lhs = Poly256::zero_with_capacity(cap);
        lhs.fill_random();
        let mut buf = vec![0u8; lhs.bytes()];
        lhs.serialize(&mut buf);
        let rhs = Poly256::deserialize(&buf, cap);
        assert_eq!(lhs, rhs);
    }
}
