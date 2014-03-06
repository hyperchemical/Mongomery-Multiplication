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
#include <stdlib.h>
#include <gmp.h>
#include <iostream>
#include <gmpxx.h>
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

		long long diff(){
			return ms.count();
		}
	private:
		Clock::time_point t0;
		milliseconds ms;
};

//Global caches for montgomery reductions
//Takes in M and T, returns m'
map<pair<mpz_class,mpz_class>, mpz_class> t_m_cache;
//Takes in M, returns Minv and M' [SIGNIFICANT SPEEDUP]
map<mpz_class, pair<mpz_class, mpz_class>> m_cache;

//Generates a b bit prime number
mpz_class generate_prime(long b, gmp_randstate_t& rand_state);
mpz_class square_and_multiply_exp(mpz_class x, mpz_class c, 
	mpz_class a, mpz_class b);
mpz_class montgomery_crt_exp(mpz_class x, mpz_class c, mpz_class a, mpz_class b);
mpz_class chinese_remainder_exp(mpz_class x, mpz_class c,
	mpz_class a, mpz_class b);
mpz_class montgomery_reduction(mpz_class T, long exponent, mpz_class& r, mpz_class& M);
mpz_class montgomery_exp(mpz_class x, mpz_class c, mpz_class a, mpz_class b);
void assert_all_equal(vector<mpz_class>& vec);


int main(int argc, char *argv[])
{
	vector<long> base_sizes = {500000, 1000000};
	vector<long> prime_sizes = {1000, 2000};

	Timer timer;
	//Types of exponentiation
	vector<string> function_names = {"S&M", "S&M w/ Mont","CRT", "CRT w/ Mont"};
	vector<bool> active(function_names.size()+1, true);
	vector<mpz_class (*)(mpz_class, mpz_class, mpz_class, mpz_class)> functions
		= {square_and_multiply_exp, montgomery_exp, chinese_remainder_exp,
			montgomery_crt_exp};

	//Deactivate algorithms as needed
	for(int i = 1; i < argc; i++){
		if(!atoi(argv[i])){
			active[i-1] = false;
		}
	}

	gmp_randstate_t rand_state;
	gmp_randinit_default(rand_state);
	srand(time(NULL));

	Clock::time_point t0 = Clock::now();
	Clock::time_point t1;
	milliseconds ms;

	//t0 = Clock::now();
	mpz_class x, c, n, a, b, result;
	int line_width = 15;
	cout << left;

	//Run algorithms for all pairs of sizes
	for(int i = 0; i < base_sizes.size(); i++){
		for(int j = 0; j < prime_sizes.size(); j++){
			vector<mpz_class> values;
			long num_size = base_sizes[i];
			long prime_size = prime_sizes[j];
			//Generate random x, c, a, b

			timer.begin();
			mpz_urandomb(x.get_mpz_t(), rand_state, num_size);
			mpz_urandomb(c.get_mpz_t(), rand_state, num_size);
			a = generate_prime(prime_size, rand_state);
			b = generate_prime(prime_size, rand_state);
			n = a*b;
			timer.end();

			cout << setw(line_width) << "Info" << "|Results\n";
			cout << setfill('-') << setw(line_width) << "-" <<  "|" << setw(line_width) << "-" << endl;
			cout << setfill(' ');
			cout << setw(line_width) << "Base size" << "|" << num_size << " bits" << endl;
			cout << setw(line_width) << "Prime size" << "|" << prime_size << " bits" << endl;
			cout << setw(line_width) << "Generate Nums" << "|" << timer.diff() << "ms\n";

			if(active[0]){
				//Generate correct answer
				timer.begin();
				mpz_powm(result.get_mpz_t(), x.get_mpz_t(), c.get_mpz_t(), n.get_mpz_t());
				timer.end();
				values.push_back(result);
				cout << setw(line_width) << "mpz_powm" << "|" << timer.diff() << "ms\n";
			}

			for(int i = 0; i < functions.size(); i++){
				if(!active[i+1]) continue;
				timer.begin();
				values.push_back(functions[i](x,c,a,b));
				//cout << function_names[i] << " " << values.back() << endl;
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

mpz_class generate_prime(long b, gmp_randstate_t& rand_state){
	mpz_class prime;
	// while(1){
	// 	prime = "1";
	// 	for(int i = 0; i < b-1; i++){
	// 		prime *= 2;
	// 		if(rand()%2 == 0)
	// 			mpz_setbit(prime.get_mpz_t(), 0);
	// 	}

	// 	if(mpz_probab_prime_p(prime.get_mpz_t(), 25) >= 1)
	// 		return prime;
	// }
	mpz_urandomb(prime.get_mpz_t(), rand_state, b);
	mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
	return prime;
}

long bit_length(mpz_class x){
	long num_ones = mpz_popcount(x.get_mpz_t());
	long num_bits = 0;
	for(long i = 0; num_ones > 0; i++)
	{
		if(mpz_tstbit(x.get_mpz_t(), i) == 1){
			num_ones--;
			num_bits = i+1;
		}
	}
	return num_bits;
}

mpz_class square_and_multiply_exp(mpz_class x, mpz_class c, 
	mpz_class a, mpz_class b){

	mpz_class n;
	n = a*b;

	mpz_class z;
	z = "1";
	long length = bit_length(c);
	for(long i = length-1; i >= 0; i--){
		z = (z*z) % n;
		if(mpz_tstbit(c.get_mpz_t(),i) == 1)
		{
			z = (z*x) % n;	
		}
	}

	return z;
}

mpz_class chinese_remainder_exp(mpz_class x, mpz_class c,
	mpz_class a, mpz_class b){
	mpz_class dp, dq, t, m1, m2, u1, u2;
	dp = c % (a-1);
	dq = c % (b-1);
	mpz_invert(t.get_mpz_t(), b.get_mpz_t(), a.get_mpz_t());
	mpz_powm(m1.get_mpz_t(), x.get_mpz_t(), dp.get_mpz_t(), a.get_mpz_t());
	mpz_powm(m2.get_mpz_t(), x.get_mpz_t(), dq.get_mpz_t(), b.get_mpz_t());
	u1 = (m1-m2);
	mpz_mod(u1.get_mpz_t(), u1.get_mpz_t(), a.get_mpz_t());
	u2 = (u1*t);
	mpz_mod(u2.get_mpz_t(), u2.get_mpz_t(), a.get_mpz_t());
	return m2 + u2*b;
}

mpz_class montgomery_crt_exp(mpz_class x, mpz_class c, mpz_class a, mpz_class b){
	mpz_class n, r, rinv;
	n = a*b;
	r = "2";
	long exponent = 1;
	for(;r < a; exponent++){
		r = r*2;
	}
	mpz_invert(rinv.get_mpz_t(), r.get_mpz_t(), n.get_mpz_t());

	mpz_class dp, dq, t, m1, m2, u1, u2;
	dp = c % (a-1);
	dq = c % (b-1);
	mpz_invert(t.get_mpz_t(), b.get_mpz_t(), a.get_mpz_t());
	mpz_powm(m1.get_mpz_t(), x.get_mpz_t(), dp.get_mpz_t(), a.get_mpz_t());
	mpz_powm(m2.get_mpz_t(), x.get_mpz_t(), dq.get_mpz_t(), b.get_mpz_t());

	u1 = (m1-m2)*r;
	mpz_mod(u1.get_mpz_t(), u1.get_mpz_t(), a.get_mpz_t());
	t = (t*r);
	mpz_mod(t.get_mpz_t(), t.get_mpz_t(), a.get_mpz_t());
	u2 = montgomery_reduction((u1*t), exponent, a, r);
	u2 = montgomery_reduction(u2, exponent, a, r);
	return m2 + u2*b;
}

mpz_class montgomery_reduction(mpz_class T, long exponent, mpz_class& M, mpz_class & r){
	mpz_class m, Minv, Mprime, t;

	if(m_cache.find(M) == m_cache.end()){
		mpz_invert(Minv.get_mpz_t(), M.get_mpz_t(), r.get_mpz_t());
		Mprime = (-1 * (Minv));
		//Right shifts instead of mod
		mpz_fdiv_r_2exp(Mprime.get_mpz_t(), Mprime.get_mpz_t(), exponent);
		//mpz_mod(Mprime.get_mpz_t(), Mprime.get_mpz_t(), r.get_mpz_t());
		m_cache[M] = {Minv, Mprime};
	}
	else{
		Minv = m_cache[M].first;
		Mprime = m_cache[M].second;
	}

	//Not cached
	if(t_m_cache.find({M,T}) == t_m_cache.end()){
		m = (T*Mprime);
		//Right shifts instead of mod
		mpz_fdiv_r_2exp(m.get_mpz_t(), m.get_mpz_t(), exponent);
		//mpz_mod(m.get_mpz_t(), m.get_mpz_t(), r.get_mpz_t());
		t_m_cache[{M,T}] = m;
	}
	else { //Cached
		m = t_m_cache[{M,T}];
	}
	
	t = (T+m*M)/r;
	return t >= M ? t - M : t;
}


mpz_class montgomery_exp(mpz_class x, mpz_class c, mpz_class a, mpz_class b){
	mpz_class n, r, rinv;
	n = a*b;
	r = "2";
	long exponent = 1;
	for(;r < n; exponent++){
		r = r*2;
	}
	mpz_invert(rinv.get_mpz_t(), r.get_mpz_t(), n.get_mpz_t());

	mpz_class z;
	z = "1";
	z = (z*r) % n; //Create residues
	x = (x*r) % n; 
	long length = bit_length(c);
	for(long i = length-1; i >= 0; i--){
		z = montgomery_reduction((z*z),exponent, n, r);;
		if(mpz_tstbit(c.get_mpz_t(),i) == 1)
		{
			z = montgomery_reduction((z*x),exponent,n, r);	
		}
	}

	return montgomery_reduction(z,exponent,n, r);
}

void assert_all_equal(vector<mpz_class>& vec)
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