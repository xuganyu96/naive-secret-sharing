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

use crate::poly::Poly256;
use galoisfields::GF2p256;
use rand::Rng;
mod f2x;
mod galoisfields;
pub mod poly;

pub type SecretSharingResult<T> = Result<T, SecretSharingError>;

#[derive(Debug)]
pub enum SecretSharingError {}

impl core::fmt::Display for SecretSharingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl core::error::Error for SecretSharingError {}

pub trait SecretSharing: Sized {
    type Poly;
    type Field;
    type Shard;

    /// Sample a random polynomial, then hash the polynomial into a 256-bit key
    /// Threshold indicates the minimum number of distinct shares needed to recover the entire
    /// polynomial, which is equal to the number of coefficients.
    fn init_with_rng(rng: &mut impl Rng, threshold: usize) -> Self;
    fn init(threshold: usize) -> Self;
    fn update(&mut self, msg: &[u8]);
    fn finalize(&mut self);
    fn split_once_with_rng(&self, rng: &mut impl Rng) -> Self::Shard;
    fn split_once(&self) -> Self::Shard;
    fn split_with_rng(&self, n: usize, rng: &mut impl Rng) -> Vec<Self::Shard>;
    fn split(&self, n: usize) -> Vec<Self::Shard>;
    fn assemble(shards: &[Self::Shard]) -> SecretSharingResult<Self>;
}

/// secret sharing using polynomials in GF(2^256)[x]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SecretSharing256 {
    secret_poly: Poly256,
}

#[derive(Debug)]
pub struct SecretSharing256Shard;

impl SecretSharing for SecretSharing256 {
    type Poly = Poly256;
    type Field = GF2p256;
    type Shard = SecretSharing256Shard;

    fn init_with_rng(rng: &mut impl Rng, threshold: usize) -> Self {
        let mut secret_poly = Poly256::zero_with_capacity(threshold);
        secret_poly.fill_random_with_rng(rng);

        Self { secret_poly }
    }
    fn init(threshold: usize) -> Self {
        todo!()
    }
    fn update(&mut self, msg: &[u8]) {
        todo!()
    }
    fn finalize(&mut self) {
        todo!()
    }
    fn split_once_with_rng(&self, rng: &mut impl Rng) -> Self::Shard {
        todo!()
    }
    fn split_once(&self) -> Self::Shard {
        todo!()
    }
    fn split_with_rng(&self, n: usize, rng: &mut impl Rng) -> Vec<Self::Shard> {
        todo!()
    }
    fn split(&self, n: usize) -> Vec<Self::Shard> {
        todo!()
    }
    fn assemble(shards: &[Self::Shard]) -> SecretSharingResult<Self> {
        todo!()
    }
}
