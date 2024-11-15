//! Shamir's secret sharing is a cryptographic protocol in which some secret value is split into n
//! shares (n is called redundancy) such that any t of those n shares (t is called threshold) is
//! enough to recover the original secret, but fewer than t shares will not reveal any information
//! about the secret.
//!
//! # How does it work?
//! Thanks to [Lagrange's interpolation](https://en.wikipedia.org/wiki/Lagrange_polynomial) we know
//! that a degree (t-1) polynomial is uniquely determined by its evaluation at t distinct points.
//! In this secret sharing scheme, the secret is the degree (t-1) polynomial. We evaluate the
//! polynomial at n distinct points, then distribute each point as a single share. When t shares
//! are assembled, the polynomial can be recovered. When only (t-1) shares are assembled, there are
//! as many possible polynomials as there are unique values the polynomial could have evaluate to,
//! so as long as we choose a cryptographically large field to build the polynomial on, the scheme
//! will be secure (it is in fact information theoretically secure).
pub mod f2x;
pub mod gf2;

use f2x::Degree;
use gf2::{FieldArithmetic, GF2p256};

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

    /// Compute the Lagrange polynomial.
    /// Points are assumed to be distinct. The behavior of interpoate is undefined if there are
    /// duplicate points.
    pub fn interpoate(points: &[(E, E)], capacity: usize) -> Self {
        if points.len() > capacity {
            panic!("Cannot interpolate more points than capacity");
        }
        let mut lagrange = Self::zero_with_capacity(capacity);

        for (alpha_i, r) in points {
            let mut basis = Self::zero_with_capacity(capacity);
            basis.coeffs[0] = E::one();

            for (alpha_j, _) in points {
                if alpha_j == alpha_i {
                    continue;
                }
                let mut factor = Self::zero_with_capacity(capacity);
                // factor is (x - alpha_j)/(alpha_i - alpha_j)
                factor.coeffs[0] = E::zero().modsub(alpha_j);
                factor.coeffs[1] = E::one();
                factor.mul_coeff(&alpha_i.modsub(alpha_j));
                basis = basis.mul(&factor);
            }

            lagrange = lagrange.add(&basis.mul_coeff(r));
        }

        lagrange
    }
}

pub type Poly256 = Poly<GF2p256>;

#[cfg(test)]
mod tests {
    use super::*;

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
}
