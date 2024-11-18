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
use aes_gcm::{
    aead::{Aead, AeadCore, KeyInit},
    Aes256Gcm, Key,
};
use galoisfields::{FieldArithmetic, GF2p256};
use rand::{rngs::OsRng, CryptoRng, Rng};
use sha3::{Digest, Sha3_256};
mod f2x;
pub mod galoisfields;
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

/// secret sharing using polynomials in GF(2^256)[x]
pub struct SecretSharing256 {
    secret_poly: Poly256,
    // TODO: Aes256Gcm does not support Debug trait; maybe I need to implement Debug for this
    // struct manually
    cipher: Aes256Gcm,
    /// AES-256-GCM uses 96-bit nonce
    nonce: [u8; 12],
    /// A ciphertext of length 0 indicates there is no ciphertext yet
    ciphertext: Vec<u8>,
    /// The collection of random points
    shards: Vec<SecretSharingShard<GF2p256>>,
}

#[derive(Debug)]
pub struct SecretSharingShard<E> {
    x: E,
    poly_at_x: E,
}

impl<E: FieldArithmetic> SecretSharingShard<E> {
    pub fn from_vals(x: E, poly_at_x: E) -> Self {
        Self { x, poly_at_x }
    }
}

impl SecretSharing256 {
    /// Generate a random polynomial and hash it into an AES key. The nonce will be randomly
    /// generated and will be included in the ciphertext.
    pub fn init_with_rng(rng: &mut (impl Rng + CryptoRng), threshold: usize) -> Self {
        let mut secret_poly = Poly256::zero_with_capacity(threshold);
        secret_poly.fill_random_with_rng(rng);

        let mut hasher: Sha3_256 = Digest::new();
        secret_poly.update_hasher(&mut hasher);
        let result = hasher.finalize(); // GenericArray<u8, OutputSize>
        let key: Key<Aes256Gcm> = result.into();
        let cipher = Aes256Gcm::new(&key);
        let nonce: [u8; 12] = Aes256Gcm::generate_nonce(rng).into();

        Self {
            secret_poly,
            cipher,
            nonce,
            ciphertext: vec![],
            shards: vec![],
        }
    }

    /// Use the default OsRng to initialize a secret sharing session
    pub fn init(threshold: usize) -> Self {
        Self::init_with_rng(&mut OsRng, threshold)
    }

    /// Feed a segment of the msg into the cipher
    pub fn encrypt(&mut self, msg: &[u8]) -> Result<(), aes_gcm::Error> {
        let ciphertext = self.cipher.encrypt(&self.nonce.into(), msg)?;
        self.ciphertext = ciphertext;
        Ok(())
    }

    /// Check if `val` has already been sampled before
    pub fn contains_point(&self, val: &GF2p256) -> bool {
        self.shards.iter().any(|shard| shard.x == *val)
    }

    /// Sample the specified number of points. This method will check to ensure all points are
    /// distinct, which might incur performance penalty
    pub fn safe_split_with_rng(&mut self, n: usize, rng: &mut (impl Rng + CryptoRng)) {
        (0..n).for_each(|_| 'rejsample: loop {
            let x = GF2p256::random_with_rng(rng);
            if !self.contains_point(&x) {
                self.shards.push(SecretSharingShard::from_vals(
                    x,
                    self.secret_poly.evaluate(&x),
                ));
                break 'rejsample;
            }
        });
    }

    /// Sample the specified number of points. This method will not check that all points are
    /// distinct and is thus faster. However, the probability of collision is cryptographically
    /// small with a large field such as GF(2^256)
    pub fn fast_split_with_rng(&mut self, n: usize, rng: &mut (impl Rng + CryptoRng)) {
        for _ in 0..n {
            let x = GF2p256::random_with_rng(rng);
            self.shards.push(SecretSharingShard::from_vals(
                x,
                self.secret_poly.evaluate(&x),
            ));
        }
    }

    /// Sample the specified number of points using the default OsRng. This method will check all
    /// points are distinct
    pub fn safe_split(&mut self, n: usize) {
        self.safe_split_with_rng(n, &mut OsRng);
    }

    /// Sample the specified number of points using the default OsRng. This method will not check
    /// all points are distinct, though the probability of collision is cryptographically small
    pub fn fast_split(&mut self, n: usize) {
        self.fast_split_with_rng(n, &mut OsRng);
    }

    /// Decrypt the ciphertext using a symmetric key derived from hashing the secret polynomial
    /// interpolating the input shards
    /// Caller is responsible for supplying the correct number of shards
    pub fn decrypt(
        ciphertext: &[u8],
        nonce: &[u8],
        shards: &[SecretSharingShard<GF2p256>],
    ) -> Result<Vec<u8>, aes_gcm::Error> {
        let points = shards
            .iter()
            .map(|shard| (shard.x.clone(), shard.poly_at_x.clone()))
            .collect::<Vec<(GF2p256, GF2p256)>>();
        let recovered_poly = Poly256::interpolate(&points, shards.len());
        let mut hasher: Sha3_256 = Digest::new();
        recovered_poly.update_hasher(&mut hasher);
        let result = hasher.finalize(); // GenericArray<u8, OutputSize>
        let key: Key<Aes256Gcm> = result.into();
        let cipher = Aes256Gcm::new(&key);
        cipher.decrypt(nonce.into(), ciphertext)
    }
}
