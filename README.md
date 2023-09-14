# TFHECDSA (TFHE ECDSA)

This repo is a tutorial/experiment/technical-report on implementing Fully Homomorphic Encrypted Elliptic Curve Digital Signature Algorithm (ECDSA) for using TFHE. I'm try to give some intro to each key concept. feel free to skip ahead.

## YOLO run

if you just want to run the experiment now.

- [install rust](https://www.rust-lang.org/tools/install)
- clone this repo.
- update `Cargo.toml` to fit your machine. see [TFHE Supported platforms](https://docs.zama.ai/tfhe-rs/getting-started/installation#supported-platforms) for more info.
- (optional) open `tmux` or equivalence if you're on a remote server.
- run `cargo run --release --example ecdsa`.
- now you have `~1-2 days` to read this doc if you're on a 64-cores machine.

## Logging

There's 3 level of logging built-in in the crates: `Info`, `Debug`, `Traces`. Default value is configured to `Debug`.

- Info - High level operation logging: all ECDSA operation, all group operation, scalar mul bit operation.
- Debug - Low level operation logging: all field operation.
- Trace - Auxilary, loop, and super low level logging: inverse mod bit, pow mod bit, selector, modular reduction,etc.

### Custom Logging Level

Use `RUST_LOG="level"` and then cargo command to force Rust to use that level of logging. For example,

```bash
RUST_LOG="info" cargo test --release -- --nocapture correct_inverse
```

### Tweak Loggin

We use `logging_timer` to automatically log our function runtime.

To remove/add timer, find these 2 patterns of code.

- Function macros - `#[time("info", "Group Projective Add")]` decorated on top of function will print out the timer when the function ends
- Inline macros
  - `let cal_bits_tmr = timer!(Level::Debug; "Calculating bits");` - `timer!` will print out timer when the timer object was dropped
  - `let _tmr = stimer!(Level::Info; "Scalar Mul", "Bits {:?}", _ic);` - `stimer!` will print out timer when the timer object was initialized and dropped
