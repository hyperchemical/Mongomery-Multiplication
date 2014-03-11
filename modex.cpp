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
map<pair<long long,long long>, long long> t_m_cache;
//Takes in M, returns Minv and M' 
map<long long, pair<long long, long long>> m_cache;

//Generates a b bit prime number
//long long generate_prime(long b, gmp_randstate_t& rand_state);
long long square_and_multiply_exp(long long x, long long c, 
	long long a, long long b);
long long montgomery_crt_exp(long long x, long long c, long long a, long long b);
long long chinese_remainder_exp(long long x, long long c,
	long long a, long long b);
long long montgomery_reduction(long long T, long exponent, const long long& r, const long long& M);
long long square_and_multiply_montgomery_exp(long long x, long long c, long long a, long long b);
void assert_all_equal(vector<long long>& vec);

// long long chinese_remainder_exp_sm(long long x, long long c,
// 	long long a, long long b);

long long mod_by_div(long long x, long long n){
	if(x > 0){
		while(x > n){
			x -= n;
		}
	}
	else{
		while(x < 0){
			x += n;
		}
	}

	return x;
}

int main(int argc, char *argv[])
{
	vector<long> base_sizes = {50000, 100000};
	vector<long> prime_sizes = {1000, 2000};

	Timer timer;
	//Types of exponentiation
	vector<string> function_names = {"S&M", "S&M w/ Mont","CRT", "CRT w/ Mont"};
	vector<bool> active(function_names.size()+1, true);
	vector<long long (*)(long long, long long, long long, long long)> functions
		= {square_and_multiply_exp, square_and_multiply_montgomery_exp, chinese_remainder_exp,
			montgomery_crt_exp};

	//Deactivate algorithms as needed
	for(int i = 1; i < argc; i++){
		if(!atoi(argv[i])){
			active[i-1] = false;
		}
	}

	// gmp_randstate_t rand_state;
	// gmp_randinit_default(rand_state);
	srand(time(NULL));

	Clock::time_point t0 = Clock::now();
	Clock::time_point t1;
	milliseconds ms;

	//t0 = Clock::now();
	long long x, c, n, a, b, result;
	int line_width = 15;
	cout << left;

	//Run algorithms for all pairs of sizes
	for(int i = 0; i < base_sizes.size(); i++){
		for(int j = 0; j < prime_sizes.size(); j++){
			vector<long long> values;
			long num_size = base_sizes[i];
			long prime_size = prime_sizes[j];
			//Generate random x, c, a, b

			timer.begin();
			// mpz_urandomb(x, rand_state, num_size);
			// mpz_urandomb(c, rand_state, num_size);
			//a = generate_prime(prime_size, rand_state);
			//b = generate_prime(prime_size, rand_state);
			x = 3;
			c = 2;
			a = 5;
			b = 7;
			n = a*b;
			timer.end();

			// cout << "x: " << x << endl;
			// cout << "c: " << c << endl;
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
				//Generate correct answer
				timer.begin();
				double x1 = x + 0.0;
				double x2 = c + 0.0;
				long long x3 = pow(x1,x2);
				result = mod_by_div(x3,n); //x3%n;
				cout << result << endl;
				//mpz_powm(result, x, c, n);
				timer.end();
				values.push_back(result);
				cout << setw(line_width) << "mpz_powm" << "|" << timer.diff() << "ms\n";
			}

			for(int i = 0; i < functions.size(); i++){
				if(!active[i+1]) continue;
				timer.begin();
				values.push_back(functions[i](x,c,a,b));
				cout << values.back() << endl;
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

// long long generate_prime(long b, gmp_randstate_t& rand_state){
// 	long long prime;
// 	mpz_urandomb(prime, rand_state, b);
// 	mpz_nextprime(prime, prime);
// 	return prime;
// }

long long bit_length(long long x){
	// long num_ones = mpz_popcount(x);
	// long num_bits = 0;
	// for(long i = 0; num_ones > 0; i++)
	// {
	// 	if(mpz_tstbit(x, i) == 1){
	// 		num_ones--;
	// 		num_bits = i+1;
	// 	}
	// }
	// return num_bits;
	long long num_bits = 0;
	while(x > 0){
		x = x >> 1;
		num_bits++;
	}

	return num_bits;
}

long long square_and_multiply_exp(long long x, long long c, 
	long long a, long long b){

	long long n;
	n = a*b;

	long long z;
	z = 1;
	long long length = bit_length(c); //bit_length(c);
	for(long long i = length-1; i >= 0; i--){
		z = mod_by_div((z*z), n);
		//if(mpz_tstbit(c,i) == 1)
		if((c >> i) & 1)
		{
			z = mod_by_div((z*x), n);
		}
	}

	return z;
}


long long chinese_remainder_exp(long long x, long long c,
	long long a, long long b){
	long long dp, dq, t, m1, m2, u1, u2, one;
	one = 1;
	dp = mod_by_div(c, (a-1));
	dq = mod_by_div(c, (b-1));
	//Invert b
	t = square_and_multiply_exp(b, a-2, a, one);
	m1 = square_and_multiply_exp(x, dp, a, one);
	m2 = square_and_multiply_exp(x, dq, b, one);
	u1 = mod_by_div((m1-m2), a);
	u2 = mod_by_div((u1*t), a);
	return m2 + u2*b;
}

long long montgomery_crt_exp(long long x, long long c, long long a, long long b){
	long long n, r, rinv, one;
	one = 1;
	n = a*b;
	r = 2;
	long exponent = 1;
	for(;r < a; exponent++){
		r = r*2;
	}

	long long dp, dq, t, m1, m2, u1, u2;
	dp = mod_by_div(c, (a-1));
	dq = mod_by_div(c, (b-1));
	//Invert b
	t = square_and_multiply_montgomery_exp(b, a-2, a, one);
	m1 = square_and_multiply_montgomery_exp(x, dp, a, one);
	m2 = square_and_multiply_montgomery_exp(x, dq, b, one);

	u1 = mod_by_div((m1-m2), a);
	u2 = mod_by_div((u1*t), a);
	u2 = montgomery_reduction((u1*t), exponent, a, r);
	u2 = montgomery_reduction(u2, exponent, a, r);
	return m2 + u2*b;
}

long long montgomery_reduction(long long T, long exponent, const long long& M, const long long & r){
	long long m, Minv, Mprime, t, one;
	one = 1;

	if(m_cache.find(M) == m_cache.end()){
		Minv = mod_by_div(square_and_multiply_exp(M, r - r/2 - 1, r, one), r);
		//mpz_fdiv_r_2exp(Mprime, Mprime, exponent);
		Mprime = mod_by_div((-1 * (Minv)), r);

		//Right shifts instead of mod
		//mpz_fdiv_r_2exp(Mprime, Mprime, exponent);
		m_cache[M] = {Minv, Mprime};
	}
	else{
		Minv = m_cache[M].first;
		Mprime = m_cache[M].second;
	}

	//Not cached
	if(t_m_cache.find({M,T}) == t_m_cache.end()){
		m = mod_by_div((T*Mprime), r);
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


long long square_and_multiply_montgomery_exp(long long x, long long c, long long a, long long b){
	long long n, r, rinv, one;
	cout << x << " " << c << " " << a << " " << b << endl;
	n = a*b;
	one = 1;
	r = 2;
	long exponent = 1;
	for(;r < n; exponent++){
		r = r*2;
	}

	rinv = mod_by_div(square_and_multiply_exp(r, (a-1)*(b-1)-1, a, b), n);
	cout << r << " Rinv " << rinv << " " << n <<  endl;

	long long z = 1;
	z = mod_by_div((z*r), n); //Create residues
	x = mod_by_div((x*r), n); 
	long long length = bit_length(c);
	for(long long i = length-1; i >= 0; i--){
		//z = montgomery_reduction((z*z),exponent, n, r);
		cout << "z: " << z * z * rinv << endl;
		z = mod_by_div(z*z*rinv, n); 
		//if(mpz_tstbit(c,i) == 1)
		if((c >> i) && 1)
		{
			//z = montgomery_reduction((z*x),exponent, n, r);	
			z = mod_by_div(z*x*rinv,n);
		}
	}

	return mod_by_div(z*rinv,n);//montgomery_reduction(z,exponent,n, r);
}

void assert_all_equal(vector<long long>& vec)
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