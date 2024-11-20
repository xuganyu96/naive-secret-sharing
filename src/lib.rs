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
mod f2x;
pub mod galoisfields;
pub mod poly;
pub mod secretsharing;
