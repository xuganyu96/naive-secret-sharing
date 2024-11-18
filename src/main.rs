use secretsharing::SecretSharing256;

fn main() {
    let threshold = 3;
    let redundancy = 10;
    let secret_msg = b"Hello, world";
    let mut secret_sharing = SecretSharing256::init(threshold);
    secret_sharing.safe_split(redundancy);
    secret_sharing.encrypt(secret_msg).unwrap();
    let _decryption = SecretSharing256::decrypt(
        &secret_sharing.ciphertext,
        &secret_sharing.nonce,
        &secret_sharing.shards[..threshold],
    )
    .unwrap();
}
