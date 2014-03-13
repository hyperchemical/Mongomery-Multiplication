// +, -, *, /, >> (bit shift), &, |, ^ (bitwise and, or, and xor), <, <=, ==, !=, >=, >.

// We will also provide an inverse function, so you are welcome to use gmp's inverse for now.

/*
#Compiler Settings
PATH := /usr/um/gcc-4.7.0/bin:$(PATH)
LD_LIBRARY_PATH := /usr/um/gcc-4.7.0./lib64
LD_RUN_PATH := /usr/um/gcc-4.7.0/lib64

g++ -lgmpxx -lgmp -std=c++11 -O3 modex.cpp -o modex

Linux caen-vnc14.engin.umich.edu 2.6.32-431.3.1.el6.x86_64 #1 SMP 
Fri Dec 13 06:58:20 EST 2013 x86_64 x86_64 x86_64 GNU/Linux
*/

#include <string> 
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <chrono>
#include <map>
#include <utility>
#include <vector>
#include <stdexcept> 
#include <cassert>
#include <iomanip>

#include "uberzahl.h"

using namespace std;

typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::milliseconds milliseconds;

class Timer{
	public: 
		Timer(){
			t0 = Clock::now();
		}

		void begin(){
			t0 = Clock::now();
		}

		void end(){
			Clock::time_point t1 = Clock::now();
			ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
		}

		uberzahl diff(){
			return ms.count();
		}
	private:
		Clock::time_point t0;
		milliseconds ms;
};

//Global caches for montgomery reductions
//Takes in M and T, returns m'
map<pair<uberzahl,uberzahl>, uberzahl> t_m_cache;
//Takes in M, returns Minv and M' 
map<uberzahl, pair<uberzahl, uberzahl>> m_cache;

//Generates a b bit prime number
//uberzahl generate_prime(long b, gmp_randstate_t& rand_state);
uberzahl square_and_multiply_exp(uberzahl x, uberzahl c, 
	uberzahl a, uberzahl b);
uberzahl montgomery_crt_exp(uberzahl x, uberzahl c, uberzahl a, uberzahl b);
uberzahl chinese_remainder_exp(uberzahl x, uberzahl c,
	uberzahl a, uberzahl b);
uberzahl montgomery_reduction(uberzahl T, long exponent, uberzahl r,  uberzahl M);
uberzahl square_and_multiply_montgomery_exp(uberzahl x, uberzahl c, uberzahl a, uberzahl b);
void assert_all_equal(vector<uberzahl>& vec);

// uberzahl chinese_remainder_exp_sm(uberzahl x, uberzahl c,
// 	uberzahl a, uberzahl b);


int main(int argc, char *argv[])
{
	vector<long> base_sizes = {100, 200};
	vector<long> exponent_sizes = {100, 200};
	vector<long> prime_sizes = {1000, 2000};

	Timer timer;
	//Types of exponentiation
	vector<string> function_names = {"S&M", "S&M w/ Mont","CRT", "CRT w/ Mont"};
	vector<bool> active(function_names.size()+1, true);
	vector<uberzahl (*)(uberzahl, uberzahl, uberzahl, uberzahl)> functions
		= {square_and_multiply_exp, square_and_multiply_montgomery_exp, chinese_remainder_exp,
			montgomery_crt_exp};

	//Deactivate algorithms as needed
	for(int i = 1; i < argc; i++){
		if(!atoi(argv[i])){
			active[i-1] = false;
		}
	}
	srand(time(NULL));

	Clock::time_point t0 = Clock::now();
	Clock::time_point t1;
	milliseconds ms;

	//t0 = Clock::now();
	uberzahl x, c, n, a, b, result;
	int line_width = 15;
	cout << left;

	//Run algorithms for all pairs of sizes
	for(int i = 0; i < base_sizes.size(); i++){
		for(int j = 0; j < prime_sizes.size(); j++){
			vector<uberzahl> values;
			long num_size = base_sizes[i];
			long prime_size = prime_sizes[j];
			//Generate random x, c, a, b

			timer.begin();
			x = x.random(base_sizes[i]);//x.random(10);
			c = x.random(exponent_sizes[i]);//c.random(10);
			a = 881;
			b = 883;
			n = a*b;
			timer.end();

			cout << "x: " << x << endl;
			cout << "c: " << c << endl;
			// cout << "a: " << a << endl;
			// cout << "b: " << b << endl;

			//Table output
			cout << setw(line_width) << "Info" << "|Results\n";
			cout << setfill('-') << setw(line_width) << "-" <<  "|" << setw(line_width) << "-" << endl;
			cout << setfill(' ');
			cout << setw(line_width) << "Base size" << "|" << num_size << " bits" << endl;
			cout << setw(line_width) << "Prime size" << "|" << prime_size << " bits" << endl;
			cout << setw(line_width) << "Generate Nums" << "|" << timer.diff() << "ms\n";

			if(active[0]){
				// //Generate correct answer
				// timer.begin();
				// uberzahl x3 = x;
				// for(uberzahl i = 1; i < c; i = i + 1){
				// 	x3 = (x3 * x3) % n;
				// }
				// cout << x3 << endl;
				// //mpz_powm(result, x, c, n);
				// timer.end();
				// values.push_back(x3);
				// cout << setw(line_width) << "mpz_powm" << "|" << timer.diff() << "ms\n";
			}

			for(int i = 0; i < functions.size(); i++){
				if(!active[i+1]) continue;
				timer.begin();
				values.push_back(functions[i](x,c,a,b));
				timer.end();
				ms = std::chrono::duration_cast<milliseconds>(t1 - t0); //timing
				cout << setw(line_width) << function_names[i] << "|" << timer.diff() << "ms\n";
			}

			//Verification of correctness
			assert_all_equal(values);
			cout << endl;
		}
	}
}

// uberzahl generate_prime(long b, gmp_randstate_t& rand_state){
// 	uberzahl prime;
// 	mpz_urandomb(prime, rand_state, b);
// 	mpz_nextprime(prime, prime);
// 	return prime;
// }

uberzahl square_and_multiply_exp(uberzahl x, uberzahl c, 
	uberzahl a, uberzahl b){

	uberzahl n;
	n = a*b;

	uberzahl z("1");
	unsigned int length = c.bitLength() - 1; //bit_length(c);

	while(1){
		z = (z*z) % n;
		if(c.bit(length) == 1)
		{
			z = (z*x) % n;
		}
		if(length == 0) break;
		length = length - 1;
	}

	return z;
}


uberzahl chinese_remainder_exp(uberzahl x, uberzahl c,
	uberzahl a, uberzahl b){
	uberzahl dp, dq, t, m1, m2, u, one("1");
	// dp = mod_by_div(c, (a-1));
	// dq = mod_by_div(c, (b-1));
	dp = c % (a-1);
	dq = c % (b-1);
	//Invert b
	t = b.inverse(a);
	m1 = square_and_multiply_exp(x, dp, a, one) % a;
	m2 = square_and_multiply_exp(x, dq, b, one) % a;
	// u1 = mod_by_div((m1-m2), a);
	// u2 = mod_by_div((u1*t), a);


	if(m2 > m1){
		m1 = m1 + a;
	}

	u = (m1-m2)*t % a;
	return m2 + u*b;
}

uberzahl montgomery_crt_exp(uberzahl x, uberzahl c, uberzahl a, uberzahl b){
	uberzahl n, r, rinv, one;
	one = 1;
	n = a*b;
	r = 2;
	long exponent = 1;
	for(;r < a; exponent++){
		r = r*2;
	}

	uberzahl dp, dq, t, m1, m2, u;
	dp = c % (a-1);
	dq = c % (b-1);
	//Invert b
	t = b.inverse(a);
	m1 = square_and_multiply_montgomery_exp(x, dp, a, one) % a;
	m2 = square_and_multiply_montgomery_exp(x, dq, b, one) % a;


	if(m2 > m1){
		m1 = m1 + a;
	}

	u = (m1-m2)*t % a;
	return m2 + u*b;
}

uberzahl montgomery_reduction(uberzahl T, long exponent, uberzahl M, uberzahl r){
	uberzahl m, Minv, Mprime, t;
	if(m_cache.find(M) == m_cache.end()){
		Minv = M.inverse(r);
		Mprime = r - Minv; //(negone*Minv) % r;


		m_cache[M] = {Minv, Mprime};
	}
	else{
		Minv = m_cache[M].first;
		Mprime = m_cache[M].second;
	}

	//Not cached
	if(t_m_cache.find({M,T}) == t_m_cache.end()){
		//m = mod_by_div((T*Mprime), r);
		m = (T*Mprime) % r;
		//Right shifts instead of mod
		//mpz_fdiv_r_2exp(m, m, exponent);
		t_m_cache[{M,T}] = m;
	}
	else { //Cached
		m = t_m_cache[{M,T}];
	}
	
	t = (T+m*M)/r;
	return t >= M ? t - M : t;
}


uberzahl square_and_multiply_montgomery_exp(uberzahl x, uberzahl c, uberzahl a, uberzahl b){
	uberzahl n, r, rinv, one;
	n = a*b;
	one = 1;
	r = 2;
	long exponent = 1;
	for(;r < n; exponent++){
		r = r*2;
	}

	//rinv = mod_by_div(square_and_multiply_exp(r, (a-1)*(b-1)-1, a, b), n);
	rinv = r.inverse(n);
	uberzahl z = 1;
	//z = mod_by_div((z*r), n); //Create residues
	//x = mod_by_div((x*r), n); 
	z = (z*r) % n;
	x = (x*r) % n;
	unsigned int length = c.bitLength() - 1; //bit_length(c);
	while(1){
		//z = (z*z*rinv) % n;
		z = z * z;
		z = montgomery_reduction(z, exponent, n, r);
		if(c.bit(length) == 1)
		{
			z = z * x;
			z = montgomery_reduction(z, exponent, n, r);
		}
		if(length == 0) break;
		length = length - 1;
	}

	return montgomery_reduction(z, exponent, n, r);
}

void assert_all_equal(vector<uberzahl>& vec)
{
	for(int i = 0; i < vec.size(); i++){
		for(int j = i; j < vec.size(); j++)
		{
			if(vec[i]!=vec[j])
			{
				cout << "Array Index: " << j << endl;
				throw runtime_error("ERROR: Assert Equal Failed");
			}
		}
	}
}