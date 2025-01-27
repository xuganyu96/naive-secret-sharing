# Next:
This project is now minimally viable. I'm moving on.

# Shamir's secret sharing
Shamir's secret sharing is actually quite simple. Its security is based on the fact that any degree- $(t-1)$ polynomial can be uniquely determined by its evaluations at $t$ points using the Lagrange Interpolation. The procedure is as follows:

- Determine some base field $K$ to operate on. This field needs to be cryptographically large, such as $\mathbb{F}_{2^{128}}$
- Determine the number of shares $n$ and the threshold $t$. This means that we will distribute shares of the secret to $n$ parties, among which any $t$ of them can be used to recover the secret
- The secret is a degree- $(t-1)$ polynomial $f \in K[x]$, which we will define using its canonical form, which contains $t$ field elements
- Choose $n$ distinct field elements $s_1, s_2, \ldots, s_n \in K$ and evaluate the secret polynomial at those n points $r_i = f(s_i)$. Distribute $(s_i, r_i)$ to each of the $n$ trustees
- Any $t$ of trustees can pool their points and uniquely determine the secret polynomial $f$ using Lagrange Interpolation

Let field size be $k = \log_2 \vert K \vert$, denote the number of shares by $n$ and the threshold by $t$, then:
- The secret polynomial needs $tk$ bits
- Each share is $2k$ bits
- Evaluating a single points requires $t$ field multiplication, so generating all shares requires $nt$ field multiplications
- Lagrange interpolation takes how many field multiplication?

# Getting started
```bash
# Split "Hello, world", each share is a file in the output ZIP archive
echo "Hello, world" | ./target/debug/shamir split -t 2 -n 3 -o shamir.zip
unzip shamir.zip
# Assemble the shares back into the original secret
./target/debug/shamir assemble 0.txt 1.txt 2.txt
# clean-up
rm shamir.zip 0.txt 1.txt 2.txt
```
