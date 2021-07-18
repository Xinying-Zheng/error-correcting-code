# Project:  convolutional code

produce a program to implement the Viterbi decoding algorithm for the widely used (2, 1, 6) convolutional code with generator matrix

G(D) = (1 + D^2 + D^3 + D^5 + D^6 , 1 + D + D^2 + D^3 + D^6).

This code is assumed to be transmitted over an additive white Gaussian noise (AWGN) channel.

## Additional Details
* Use the recursion u(l+6) = u(l+1) ⊕ u(l), for l ≥ 0 with the initial conditions u(0) = 1, u(1) = u(2) = u(3) = u(4) = u(5) = 0 to generate the information bits. Ensure that the generated sequence is 100000100001 . . . and is periodic with period 63.

* Encode the information sequence using the generator matrix G(D).

*  The encoder outputs 0’s and 1’s. However, the input to the AWGN channel is normalized to ±1. Therefore, map 0’s to +1’s and 1’s to o 1’s.

*  To simulate the AWGN channel with unquantized soft-decision decoding, add a normal (Gaussian) random variable of mean zero and variance σ^2 to the ±1’s generated at the previous step. For a binary code of rate R on the AWGN channel with antipodal signaling, the relationship between Eb/N0 and σ2 is given by
![](http://latex.codecogs.com/svg.latex?\sigma^{2}=\left(2 R \frac{E_{b}}{N_{0}}\right)^{-1})
so for example for a R = 1/2 code, the relationship is simply
![](http://latex.codecogs.com/svg.latex?\sigma^{2}=\left(\frac{E_{b}}{N_{0}}\right)^{-1}) 
Please remember that Eb/N0 is always quoted in “dBs,” which equals 10 log10(Eb/N0). Thus for example, a value of Eb/N0 of 4 dB for a R = 1/2 code corresponds to a value of σ^2 = 0.3981.
* Use the following segment of pseudo code to generate normal random variables of mean zero and variance σ2. The procedure normal outputs two independent normal random variables, n1 and n2, and Ranq1 is a function which generates a random variable uniformly distributed in the interval (0, 1).

      '''unsigned long long SEED;
      // SEED must be an unsigned integer smaller than 4101842887655102017.
      unsigned long long RANV;
      int RANI = 0;
      main(){
      · · ·
      · · ·
      · · ·}
      normal(n1, n2, σ)
      {
      do{
        x1 = Ranq1();
        x2 = Ranq1();
        x1 = 2x1-1; 
        x2 = 2x2-1; 
        s = x1^2+x2^2
      } while (s ≥ 1.0)
      n1 = σx1(-2lns/s)^(-1/2)
      n2 = σx2(-2lns/s)^(-1/2)
      }

      double Ranq1()
      {
        if ( RANI == 0 ){
        RANV = SEED ∧ 4101842887655102017LL;
      RANV ∧= RANV >> 21;
      RANV ∧= RANV << 35;
      RANV ∧= RANV >> 4;
      RANV = RANV * 2685821657736338717LL;
      RANI++;
      }
      RANV ∧= RANV >> 21;
      RANV ∧= RANV << 35;
      RANV ∧= RANV >> 4;
      return RANV * 2685821657736338717LL * 5.42101086242752217E-20;
      }
      '''

*  To get the output of the BSC, take the sign of the output of the AWGN channel and map +1’s to 0’s and d 1’s to 1’s.

*  In your decoder, truncate the survivors to length 32 and output the oldest bit on the survivor with the best metric. To decode N bits, generate N + 31 bits in (1). Finally compare the decoded information sequence with the original information sequence. If there are K bit errors, K/N will be a good estimate of the decoded BER.

* As a partial check, some typical values are listed below. Eb/N0 BER (BSC) Eb/N0 BER (AWGN)

## Other Notes for Demonstration
* The survivor truncation length corresponds to the actual storage requirement of the
survivors. For example, a survivor truncation length of 32 for this code means that
each survivor stores 32 bits.

* For the illustration below, suppose a state is described as the content of the feed-forward shift register in the encoder s = (s1, s2, s3, s4, s5, s6), where the input information bit first fed to s1 and then shifted from left to right. In the trellis diagram, consider placing the states vertically from top to bottom in the order of (0 0 0 0 0 0), (1 0 0 0 0 0), (0 1 0 0 0 0), (1 1 0 0 0 0), (0 0 1 0 0 0), . . ., (1 1 1 1 1 1). What to do in case of tied metrics? In the “add-compare-select” step the two metrics could be equal. In this case, if 0’s and 1’s are equally probable to occur in the transmitted information sequence, in principle you can safely select either case, and it will not affect the decoder performance.
Yet for the purpose of demonstration, always choose the upper branch as the survivor. If best-state output decision is employed, in case of tied metrics, in principle you can also safely select either case, but again for the purpose of demonstration, always choose the survivor of the uppermost state.

* Except in the procedure normal for generating noise, if a random number is needed in your program, use other random number generators instead of the function Ranq1, for the purpose of demonstration.

* Each call of the procedure normal can return two independent normal random variables, n1 and n2. Please use both of them in your program. Specifically, since this is a (2, 1) code, each branch transition consists of two encoded bits, say x1 and x2. Add n1 and n2 to x1 and x2, respectively, to get the two channel outputs y1 and y2, i.e., y1 = x1 + n1 and y2 = x2 + n2.
