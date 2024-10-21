# Secret Sharing
Say I have a secret key (think master password, SSH secret key, signing key, decryption key), and I want to prevent it from getting lost if I lose my device or forget my master password. One way to recover is by sharing it with other people, but sharing my master password with other people directly is too risky.

What are my security assumptions:
- My friends will not collude
- My friends may be forgetful
- My friends may unintentionally leak my secret

A naive secret sharing scheme would be as follows:
- Generate a random 256-bit key $k$ (AES-256-GCM or ChaCha20-Poly1305) and encrypt my secret using this key (we call it the master key). Publish the ciphertext
- I want to distribute my key to n trustees such that any t of them can assemble the master key. This means I need to break up the key $n \choose t$ times.
- For each of the $n \choose t$ possible parties of $t$ trustees:
    - Generate $t - 1$ random masks $s_1, s_2, \ldots, s_{t-1}$
    - The last mask is $s_t = k \oplus s_1 \oplus s_2 \oplus \ldots \oplus s_{t-1}$
- Each trustee should receive $n-1 \choose t-1$ masks.

This naive secret sharing can be made even more general. It doesn't have to be symmetric key; it can just be some secret seeds. 

## Usage
The main function of secret sharing is to generate a user-specified number of bits of entropy, then split it across the user-specified number of trustees at the user-specified number of threshold.

```bash
# Generate n files in the specified folder named "1.share", "2.share", etc.
secretshare generate -n <trustees> -t <threshold> -b <entropybytes> -o <folder>

# Given files, attempt to recover the secret
secretshare assemble [file1] [file2] ...
```
