EECS475PS4
==========

heyoh@umich.edu
hning@umich.edu
milewa@umich.edu

Command Line Arguments
-------------
Takes optional command line arguments

By default, runs all algorithms.

A value of 0 zero disables the given algorithm.

> ./modex 0 0 0 0

* Arg0 [1/0]: Run S&M
* Arg1 [1/0]: Run S&M w/ Montgomery Reductions
* Arg2 [1/0]: Run CRT 
* Arg3 [1/0]: Run CRT w/ Montgomery Reductions
 
Sample Output
-------------

Info           |Results
---------------|---------------
Base size      |100000
Prime size     |500
Generate Prime |35ms
mpz_powm       |151ms
S&M            |5280ms
S&M w/ Mont    |938ms
CRT            |0ms
CRT w/ Mont    |0ms

Info           |Results
---------------|---------------
Base size      |100000
Prime size     |1000
Generate Prime |123ms
mpz_powm       |254ms
S&M            |8671ms
S&M w/ Mont    |2155ms
CRT            |1ms
CRT w/ Mont    |1ms

Info           |Results
---------------|---------------
Base size      |200000
Prime size     |500
Generate Prime |47ms
mpz_powm       |139ms
S&M            |19151ms
S&M w/ Mont    |2267ms
CRT            |0ms
CRT w/ Mont    |0ms

Info           |Results
---------------|---------------
Base size      |200000
Prime size     |1000
Generate Prime |237ms
mpz_powm       |508ms
S&M            |33983ms
S&M w/ Mont    |4405ms
CRT            |2ms
CRT w/ Mont    |2ms

