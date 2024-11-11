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

use gf2::GF2p128;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Poly {
    /// Coefficients of the polynomial are laid out in little-endian order: lower vector index
    /// encodes coefficient for lower-power term
    coeffs: Vec<GF2p128>,
}

impl Poly {
    /// Return a zero polynomial with the specified number of coefficients
    pub fn zero(threshold: usize) -> Self {
        let coeffs = vec![GF2p128::ZERO; threshold];
        Self { coeffs }
    }

    pub fn add(&self, rhs: &Self) -> Self {
        todo!();
    }

    pub fn sub(&self, rhs: &Self) -> Self {
        todo!();
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        todo!();
    }

    pub fn evaluate(&self, at: GF2p128) -> GF2p128 {
        todo!();
    }

    pub fn interpoate(points: &[(GF2p128, GF2p128)]) -> Self {
        todo!();
    }
}
