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
use gf2::FieldArithmetic;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Poly<E: FieldArithmetic> {
    coeffs: Vec<E>,
}

impl<E: FieldArithmetic> Poly<E> {
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
            .unwrap(); // because self is not zero, this is guaranteed to have something
        return Degree::NonNegative(degree);
    }

    pub fn add(&self, _rhs: &Self) -> Self {
        todo!();
    }

    pub fn sub(&self, _rhs: &Self) -> Self {
        todo!();
    }

    pub fn mul(&self, _rhs: &Self) -> Self {
        todo!();
    }

    pub fn evaluate(&self, _at: E) -> E {
        todo!();
    }

    pub fn interpoate(_points: &[(E, E)]) -> Self {
        todo!();
    }
}
