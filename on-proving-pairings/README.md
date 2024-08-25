# on-proving-pairing

The pairing is costly. `On Proving Pairing` has figured out a solution: prove and verify for pairing.

By precompute miller lines and avoid the `final_exponentiation`, can reduce the cost of paring.


## How does it works
### pairing
1. compute f with miller_loop
2. compute final_f with `final_exponentiation`

## On Proving Pairing
1. setup for pvk: precompute lines and find_c.
2. prove: generate pairing proof `final_f`
3. verify: check the pairing proof with hint.

>> Note: Only support Bn254 for now.


## Reference
* [On Proving Pairings](https://eprint.iacr.org/2024/640)
* [Optimal Pairings](https://eprint.iacr.org/2008/096)
* [A remark on the computation of cube roots in finite fields](https://eprint.iacr.org/2009/457)