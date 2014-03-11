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

Info           |Results
---------------|---------------
Base size      |50000 bits
Prime size     |1000 bits
Generate Nums  |179ms
mpz_powm       |281ms
S&M            |2397ms
S&M w/ Mont    |1048ms
CRT            |51ms
CRT w/ Mont    |38ms

Info           |Results
---------------|---------------
Base size      |50000 bits
Prime size     |2000 bits
Generate Nums  |1618ms
mpz_powm       |482ms
S&M            |3980ms
S&M w/ Mont    |2694ms
CRT            |189ms
CRT w/ Mont    |148ms

Info           |Results
---------------|---------------
Base size      |100000 bits
Prime size     |1000 bits
Generate Nums  |262ms
mpz_powm       |258ms
S&M            |8746ms
S&M w/ Mont    |2310ms
CRT            |103ms
CRT w/ Mont    |39ms

Info           |Results
---------------|---------------
Base size      |100000 bits
Prime size     |2000 bits
Generate Nums  |603ms
mpz_powm       |937ms
S&M            |14635ms
S&M w/ Mont    |5312ms
CRT            |368ms
CRT w/ Mont    |152ms

