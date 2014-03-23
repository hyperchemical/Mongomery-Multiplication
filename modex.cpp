/*
EECS 475
Team: Lord Helix
Team #: 26

Heyung Joo Oh (heyoh)
Haoran Ning (hning)
Michael Wang (milewa)

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

//Timer Infrastructure
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

		long diff(){
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
//t cache for b^-1 mod a
map<pair<uberzahl, uberzahl>, uberzahl> binv_cache;

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


int main(int argc, char *argv[])
{
	long base_size = 128;
	long n_size = 256;
	int num_times_same_n = 10;

	Timer timer;
	//Types of exponentiation
	vector<long> n_sizes = {16, 32, 48, 64, 80, 96, 112, 128, 256};
	vector<string> function_names = {"S&M", "S&M w/ Mont","CRT", "CRT w/ Mont"};
	vector<double> total_times(function_names.size());
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

	uberzahl x, c, n, a, b, result;
	int line_width = 15;
	cout << left;

	for(int j = 0; j < n_sizes.size(); j++){
		t_m_cache.clear();
		m_cache.clear();
		binv_cache.clear();
		for(auto i : total_times){
			i = 0;
		}

		a = a.random(n_sizes[j]/2);
		a = nextprime(a, 50);
		b = b.random(n_sizes[j]/2);
		b = nextprime(b, 50);
		n = a*b;

		for(int i = 0; i < num_times_same_n; i++){
			x = x.random(base_size);
			c = c.random(n_sizes[j]);

			vector<uberzahl> values;

			for(int i = 0; i < functions.size(); i++){
				if(!active[i+1]) continue;
				timer.begin();
				values.push_back(functions[i](x,c,a,b));
				timer.end();
				long time_diff = timer.diff();
				total_times[i] += time_diff;
			}

			//Verification of correctness
			assert_all_equal(values);

		}

		cout << "Average Times with " << num_times_same_n << " trials:\n";
		cout << "N size: " << n_sizes[j] << endl;
		for(int i = 0; i < functions.size(); i++){
			cout << function_names[i] << ": " << total_times[i]/num_times_same_n << " ms\n";
		}
		cout << endl;
	}
}


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
	dp = c % (a-1);
	dq = c % (b-1);

	t = b.inverse(a);

	m1 = square_and_multiply_exp(x, dp, a, one);
	m2 = square_and_multiply_exp(x, dq, b, one);

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


	uberzahl dp, dq, t, m1, m2, u;
	dp = c % (a-1);
	dq = c % (b-1);
	//Invert b
	if(binv_cache.find({b,a}) == binv_cache.end())
	{
		t = b.inverse(a);
		binv_cache[{b,a}] = t;
	} else {
		t = binv_cache[{b,a}];
	}

	m1 = square_and_multiply_montgomery_exp(x, dp, a, one); 
	m2 = square_and_multiply_montgomery_exp(x, dq, b, one); 

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
		Mprime = r - Minv; 


		m_cache[M] = {Minv, Mprime};
	}
	else{
		Minv = m_cache[M].first;
		Mprime = m_cache[M].second;
	}

	//Not cached
	if(t_m_cache.find({M,T}) == t_m_cache.end()){

		m = (T*Mprime) & (r-1);

		t_m_cache[{M,T}] = m;
	}
	else { //Cached
		m = t_m_cache[{M,T}];
	}
	
	//t = (T+m*M)/r;
	t = (T+m*M);
	t = t >> (exponent); 
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
	//rinv = r.inverse(n);
	uberzahl z = 1;
	z = (z*r) % n;
	x = (x*r) % n;
	unsigned int length = c.bitLength() - 1; //bit_length(c);
	while(1){
		// z = z * z;
		z = montgomery_reduction(z*z, exponent, n, r);
		if(c.bit(length) == 1)
		{
			//z = z * x;
			z = montgomery_reduction(z*x, exponent, n, r);
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