EECS475PS4
==========

heyoh@umich.edu
hning@umich.edu
milewa@umich.edu

I'm not sure if any of this actually works, so it's probably good to look it over hah - hning.

Command Line Arguments
-------------
Takes optional command line arguments. By default, runs all algorithms.

A value of 0 zero disables the given algorithm. Partial inputs are accepted. 

For example the following disables montgomery reduction algorithms:
> ./modex 1 1 0 1 0

* Arg0 [1/0]: Run mpz_powm
* Arg1 [1/0]: Run S&M
* Arg2 [1/0]: Run S&M w/ Montgomery Reductions
* Arg3 [1/0]: Run CRT 
* Arg4 [1/0]: Run CRT w/ Montgomery Reductions
 
Sample Output
-------------

x: 13
c: 1023
a: 881
b: 883
Info           |Results
---------------|---------------
Base size      |13 bits
Exponent size  |10 bits
Prime size     |1000 bits
Result: 401718
S&M            |1ms
Result: 401718
S&M w/ Mont    |5ms
Result: 401718
CRT            |0ms
Result: 401718
CRT w/ Mont    |4ms

x: 13
c: 1024
a: 881
b: 883
Info           |Results
---------------|---------------
Base size      |13 bits
Exponent size  |10 bits
Prime size     |1000 bits
Result: 554796
S&M            |0ms
Result: 554796
S&M w/ Mont    |4ms
Result: 554796
CRT            |0ms
Result: 554796
CRT w/ Mont    |5ms

x: 13
c: 1025
a: 881
b: 883
Info           |Results
---------------|---------------
Base size      |13 bits
Exponent size  |10 bits
Prime size     |1000 bits
Result: 211041
S&M            |0ms
Result: 211041
S&M w/ Mont    |4ms
Result: 211041
CRT            |0ms
Result: 211041
CRT w/ Mont    |4ms