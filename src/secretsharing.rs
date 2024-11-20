//! Top level secret sharing data structures and routines
use crate::galoisfields::{FieldArithmetic, GF2p256};
use crate::poly::{Poly256, Poly256Point};
use aes_gcm::{
    aead::{Aead, AeadCore, KeyInit},
    Aes256Gcm, Key,
};
use base64::{prelude::BASE64_STANDARD, Engine};
use rand::{rngs::OsRng, CryptoRng, Rng};
use sha3::{Digest, Sha3_256};

pub type SecretSharingResult<T> = Result<T, SecretSharingError>;

#[derive(Debug)]
pub enum SecretSharingError {
    NoCiphertext,
    TooFewShards {
        expect: usize,
        has: usize,
    },
    InvalidShares(String),
    InvalidNonce(String),
    EmptySharesInput,
    AesGcmError(aes_gcm::Error),
    /// some input base64 encoding is invalid
    Base64DecodingError,
}

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
    pub nonce: [u8; 12],
    /// A ciphertext of length 0 indicates there is no ciphertext yet
    pub ciphertext: Vec<u8>,
    /// The collection of random points
    pub shards: Vec<Poly256Point>,
}

/// Convenient struct for serde
#[derive(Debug)]
pub struct SecretShare {
    threshold: usize,
    nonce: String,
    ciphertext: String,
    secret_x: String,
    secret_fx: String,
}

impl SecretShare {
    pub fn to_string(&self) -> String {
        todo!()
    }
    pub fn from_string(_s: &str) -> Result<Self, SecretSharingError> {
        todo!()
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
        println!("Polynomial hashes into: {:?}", &result);
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

    /// The minimum number of shares needed to recover the secret polynomial
    pub fn threshold(&self) -> usize {
        self.secret_poly.capacity()
    }

    /// Use the default OsRng to initialize a secret sharing session
    pub fn init(threshold: usize) -> Self {
        Self::init_with_rng(&mut OsRng, threshold)
    }

    /// Feed a segment of the msg into the cipher
    pub fn encrypt(&mut self, msg: &[u8]) -> Result<(), SecretSharingError> {
        let ciphertext = self
            .cipher
            .encrypt(&self.nonce.into(), msg)
            .map_err(|aes_gcm_error| SecretSharingError::AesGcmError(aes_gcm_error))?;
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
                self.shards
                    .push(Poly256Point::from_vals(x, self.secret_poly.evaluate(&x)));
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
            self.shards
                .push(Poly256Point::from_vals(x, self.secret_poly.evaluate(&x)));
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
        shards: &[Poly256Point],
    ) -> Result<Vec<u8>, SecretSharingError> {
        let points = shards
            .iter()
            .map(|shard| (shard.x.clone(), shard.fx.clone()))
            .collect::<Vec<(GF2p256, GF2p256)>>();
        let recovered_poly = Poly256::interpolate(&points, shards.len());
        let mut hasher: Sha3_256 = Digest::new();
        recovered_poly.update_hasher(&mut hasher);
        let result = hasher.finalize(); // GenericArray<u8, OutputSize>
        println!("Polynomial hashes into: {:?}", &result);
        let key: Key<Aes256Gcm> = result.into();
        let cipher = Aes256Gcm::new(&key);

        match cipher.decrypt(nonce.into(), ciphertext) {
            Ok(decryption) => Ok(decryption),
            Err(e) => Err(SecretSharingError::AesGcmError(e)),
        }
    }

    /// Return a list of string, where each string contains (threshold, nonce, ciphertext, shard),
    /// which can then be written to a file for further storage.
    /// If there is no ciphertext or if the number of shards is less than the threshold, then
    /// return appropriate error.
    pub fn stringify_shards(&self) -> Result<Vec<SecretShare>, SecretSharingError> {
        if self.ciphertext.len() == 0 {
            return Err(SecretSharingError::NoCiphertext);
        }
        if self.shards.len() < self.threshold() {
            return Err(SecretSharingError::TooFewShards {
                expect: self.threshold(),
                has: self.shards.len(),
            });
        }
        let nonce_str = BASE64_STANDARD.encode(&self.nonce);
        let ciphertext_str = BASE64_STANDARD.encode(&self.ciphertext);
        let shares = self
            .shards
            .iter()
            .map(|shard| {
                let mut buf = [0u8; GF2p256::BYTES];
                shard.x.write_be_bytes(&mut buf);
                let secret_x_str = BASE64_STANDARD.encode(&buf);
                shard.fx.write_be_bytes(&mut buf);
                let secret_fx_str = BASE64_STANDARD.encode(&buf);

                SecretShare {
                    threshold: self.threshold(),
                    nonce: nonce_str.clone(),
                    ciphertext: ciphertext_str.clone(),
                    secret_x: secret_x_str,
                    secret_fx: secret_fx_str,
                }
            })
            .collect::<Vec<SecretShare>>();

        Ok(shares)
    }

    /// Attempt to decrypt using the input set of shares
    /// The shares need to have identical threshold, nonce, ciphertext, and distinct secret_x, or
    /// they will be considered invalid input
    pub fn decrypt_from_secret_shares(
        shares: &[SecretShare],
    ) -> Result<Vec<u8>, SecretSharingError> {
        if shares.len() == 0 {
            return Err(SecretSharingError::EmptySharesInput);
        }
        for i in 1..(shares.len()) {
            let share = &shares[i];
            if share.nonce != shares[0].nonce
                || share.threshold != shares[0].threshold
                || share.ciphertext != shares[0].ciphertext
            {
                return Err(SecretSharingError::InvalidShares(
                    "Nonce, threshold, or ciphertext is inconsistent".into(),
                ));
            }
        }
        let nonce = BASE64_STANDARD
            .decode(&shares[0].nonce)
            .map_err(|_b64_err| SecretSharingError::Base64DecodingError)?;
        if nonce.len() != aes_gcm::Aes256Gcm::generate_nonce(&mut OsRng).len() {
            return Err(SecretSharingError::InvalidNonce(
                "expect nonce to be 96-bit".into(),
            ));
        }
        let ciphertext = BASE64_STANDARD
            .decode(&shares[0].ciphertext)
            .map_err(|_b64_err| SecretSharingError::Base64DecodingError)?;
        let threshold = shares[0].threshold;
        if shares.len() < threshold {
            return Err(SecretSharingError::TooFewShards {
                expect: threshold,
                has: shares.len(),
            });
        }
        let top_shares = match shares.get(..threshold) {
            Some(slice) => slice,
            None => {
                return Err(SecretSharingError::TooFewShards {
                    expect: threshold,
                    has: shares.len(),
                });
            }
        };
        // This is the coolest shit!
        // An Iterator<Item = Result<T, E>> can be collected into a Result<Vec<T>, E>
        let points = top_shares
            .iter()
            .map(|share| {
                let x = GF2p256::from_be_bytes(
                    &BASE64_STANDARD
                        .decode(&share.secret_x)
                        .map_err(|_b64err| SecretSharingError::Base64DecodingError)?,
                );
                let fx = GF2p256::from_be_bytes(
                    &BASE64_STANDARD
                        .decode(&share.secret_fx)
                        .map_err(|_b64err| SecretSharingError::Base64DecodingError)?,
                );
                let point = Poly256Point::from_vals(x, fx);
                Ok(point)
            })
            .collect::<Result<Vec<Poly256Point>, SecretSharingError>>()?;

        Self::decrypt(&ciphertext, &nonce, &points)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn random_secret_sharing_256() -> Result<(), SecretSharingError> {
        let threshold = 3;
        let redundancy = 10;
        let secret_msg = b"Hello, world";
        let mut secret_sharing = SecretSharing256::init(threshold);
        secret_sharing.safe_split(redundancy);
        secret_sharing.encrypt(secret_msg)?;
        let decryption = SecretSharing256::decrypt(
            &secret_sharing.ciphertext,
            &secret_sharing.nonce,
            &secret_sharing.shards[..threshold],
        )?;
        assert_eq!(decryption, secret_msg);

        Ok(())
    }

    #[test]
    fn random_secret_sharing_256_stringified() -> Result<(), SecretSharingError> {
        let threshold = 3;
        let redundancy = 10;
        let secret_msg = b"Hello, world";
        let mut secret_sharing = SecretSharing256::init(threshold);
        secret_sharing.safe_split(redundancy);
        secret_sharing.encrypt(secret_msg)?;
        let shares = secret_sharing.stringify_shards()?;
        let decryption = SecretSharing256::decrypt_from_secret_shares(&shares)?;
        assert_eq!(decryption, secret_msg);

        Ok(())
    }
}
