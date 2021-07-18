# Project:  convolutional code

produce a program to implement the Viterbi decoding algorithm for the widely used (2, 1, 6) convolutional code with generator matrix

G(D) = (1 + D^2 + D^3 + D^5 + D^6 , 1 + D + D^2 + D^3 + D^6).

This code is assumed to be transmitted over an additive white Gaussian noise (AWGN) channel.

## Additional Details
* Use the recursion ul+6 = ul+1 ⊕ ul, for l ≥ 0 with the initial conditions u0 = 1, u1 = u2 = u3 = u4 = u5 = 0 to generate the information bits. Ensure that the generated sequence is 100000100001 . . . and is periodic
with period 63.
* Encode the information sequence using the generator matrix G(D).
*  The encoder outputs 0’s and 1’s. However, the input to the AWGN channel is normalized to ±1. Therefore, map 0’s to +1’s and 1’s to o 1’s.
*  To simulate the AWGN channel with unquantized soft-decision decoding, add a normal (Gaussian) random variable of mean zero and variance σ2 
*  to the ±1’s generated at the previous step. For a binary code of rate R on the AWGN channel with antipodal
signaling, the relationship between Eb/N0 and σ2 is given by
* Use the following segment of pseudo code to generate normal random variables of mean zero and variance σ2. The procedure normal outputs two independent normal random variables, n1 and n2, and Ranq1 is a function which generates a random variable
uniformly distributed in the interval (0, 1).
'''
unsigned long long SEED;
// SEED must be an unsigned integer smaller than 4101842887655102017.
unsigned long long RANV;
int RANI = 0;
main()
{
· · ·
· · ·
· · ·
}
normal(n1, n2, σ)
{
do{
x1 = Ranq1();
x2 = Ranq1();
x1 = 2x1 1 1; 
x2 = 2x2 2 1; 
s = x
2
1 + x
2
2; 
} while (s ≥ 1.0)
n1 = σx1
q໒? 
2 ln s/s; 
n2 = σx2
q໒? 
2 ln s/s; 
}
double Ranq1()
{
if ( RANI == 0 ){
RANV = SEED 
∧ 
4101842887655102017LL;
RANV 
∧
= RANV >> 21;
RANV 
∧
= RANV << 35;
RANV 
∧
= RANV >> 4;
RANV = RANV * 2685821657736338717LL;
RANI++;
}
RANV 
∧
= RANV >> 21;
RANV 
∧
= RANV << 35;
RANV 
∧
= RANV >> 4;
return RANV * 2685821657736338717LL * 5.42101086242752217E-20;
}
'''
*  To get the output of the BSC, take the sign of the output of the AWGN channel and map +1’s to 0’s and d 1’s to 1’s.
*  In your decoder, truncate the survivors to length 32 and output the oldest bit on the
survivor with the best metric. To decode N bits, generate N + 31 bits in (1). Finally
compare the decoded information sequence with the original information sequence. If
there are K bit errors, K/N will be a good estimate of the decoded BER.
* As a partial check, some typical values are listed below.
Eb/N0 BER (BSC) Eb/N0 BER (AWGN)

