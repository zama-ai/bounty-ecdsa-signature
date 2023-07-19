# TFHECDSA (TFHE ECDSA)
This repo is a tutorial/experiment/technical-report on implementing Fully Homomorphic Encrypted Elliptic Curve Digital Signature Algorithm (ECDSA) for using TFHE. I'm try to give some intro to each key concept. feel free to skip ahead. 
table of content
- Intro to tech
- Design constraint
- Technical detail

## Intro
### (Fully) Homomorphic Encryption (FHE)
`tl;dr` is homomorphic encryption allow us to do math on encrypted data. informally you can think about it as
```
E(a) + E(b) = E(a + b)it will ended up the same
or 
E(f(a)) = f(E(a))
```

the fully homomorphic encryption basically mean you can do both `+` and `x` with the encrypted data. long story short, if you can do both you can compute anything. (might be super slow tho)

### TFHE 
it's one of the leading FHE lib right now. you can read more about it [here](https://zama.ai).

### Elliptic Curve Digital Signature Algorithm (ECDSA)
so it's digital signature is basically digital version of signature that you use to sign things and allow others to verify authenticity of the things you signed. if you're in blockchain this is what happened when you sign a transaction. this algorithm also use behind the scene of a lot of things in your everyday life. 

Elliptic Curve is a weird curve lol
$$y^{2}=x^{3}+ax+b$$
I'll talk about it's property later. 
for now just know that this curve is the key building block for the signature algorithm that we working on.

### Finite Field Arithmetic
tl;dr just imagine you have limited number to work with like do math with your fingures or do math on a clock. turn out this is super useful in cryptography and super slow to compute.

## Design constraint
before we dive into how we build this. let's give some context about constraint in working with (T)FHE.

### speed
bitops < add, sub << mul <<< sub

### control flow
you can't do treditional control flow on encrypted data.
