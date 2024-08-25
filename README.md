# on-proving-pairing

The pairing is costly. `On Proving Pairing` has figured out a solution: prove and verify for pairing.

By precompute miller lines and avoid the `final_exponentiation`, can reduce the cost of paring.


## How does it work
### pairing
1. compute f with miller_loop
2. compute final_f with `final_exponentiation`

## On Proving Pairing
1. setup for pvk: precompute lines and find_c.
2. prove: generate pairing proof `final_f`
3. verify: check the pairing proof with hint.

> Note: Only support Bn254 for now.


## How to use this
Add dependency on Cargo.toml:
```toml
on-proving-pairings = {git="https://github.com/SuccinctPaul/on-proving-pairings.git", package = "on-proving-pairings"}

# As some structs and funcitons are private in arkworks-rs v0.4.0. So needs to use the modified one.
[patch.crates-io]
ark-ff = { git = "https://github.com/SuccinctPaul/arkworks-algebra.git",  branch = "v0.4.2"}
ark-ec = { git = "https://github.com/SuccinctPaul/arkworks-algebra.git",  branch = "v0.4.2"}
ark-serialize = { git = "https://github.com/SuccinctPaul/arkworks-algebra.git",  branch = "v0.4.2"}
ark-poly = { git = "https://github.com/SuccinctPaul/arkworks-algebra.git",  branch = "v0.4.2"}
```

## Examples
* [Groth16 Verifier](./groth16-verifier)
* [Fflonk Verifier](https://github.com/SuccinctPaul/ark-fflonk-verifier)


## Reference
* [On Proving Pairings](https://eprint.iacr.org/2024/640)
* [The Realm of the Pairings](https://eprint.iacr.org/2013/722)
* [A remark on the computation of cube roots in finite fields](https://eprint.iacr.org/2009/457)
* [How we implemented the BN254 Ate pairing in lambdaworks](https://blog.lambdaclass.com/how-we-implemented-the-bn254-ate-pairing-in-lambdaworks/)
* [See More Implements References](REFERENCES.md)