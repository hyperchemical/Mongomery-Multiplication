#include <string>
#include <sys/time.h>  
#include <stdlib.h>
#include <gmp.h>
#include <iostream>
#include <gmpxx.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <chrono>
#include <map>

typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::milliseconds milliseconds;
using namespace std;

int count = 0;
int count2 = 0;
map<mpz_class, mpz_class> t_m_cache;
mpz_class Mprime;

mpz_class generate_prime(int b){
	mpz_class prime;
	while(1){
		prime = "1";
		for(int i = 0; i < b-1; i++){
			prime *= 2;
			if(rand()%2 == 0)
				mpz_setbit(prime.get_mpz_t(), 0);
		}

		if(mpz_probab_prime_p(prime.get_mpz_t(), 25) >= 1)
			return prime;
	}
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
	//cout << x << " num_bits: " << num_bits << endl;
	return num_bits;
}

mpz_class square_and_multiply_exp(mpz_class x, mpz_class c, 
	mpz_class n){

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
	while(u1 < 0){
		u1 += a;
	}
	u1 = u1 % a;
	u2 = (u1*t) % a;
	return m2 + u2*b;
}

mpz_class montgomery_crt_exp(mpz_class x, mpz_class c, mpz_class a, mpz_class b){
	mpz_class n, r, rinv;
	n = a*b;
	r = "2";
	while(r < n){
		r = r*2;
	}
	mpz_invert(rinv.get_mpz_t(), r.get_mpz_t(), n.get_mpz_t());

	mpz_class dp, dq, t, m1, m2, u1, u2;
	dp = c % (a-1);
	dq = c % (b-1);
	mpz_invert(t.get_mpz_t(), b.get_mpz_t(), a.get_mpz_t());
	mpz_powm(m1.get_mpz_t(), x.get_mpz_t(), dp.get_mpz_t(), a.get_mpz_t());
	mpz_powm(m2.get_mpz_t(), x.get_mpz_t(), dq.get_mpz_t(), b.get_mpz_t());
	u1 = (m1-m2);
	while(u1 < 0){
		u1 += a;
	}
	u1 = u1 % a;
	u2 = (u1*t) % a;
	return m2 + u2*b;
}

mpz_class montgomery_reduction(mpz_class T, mpz_class& r, mpz_class& M){
	mpz_class m, Minv, t;
	mpz_invert(Minv.get_mpz_t(), M.get_mpz_t(), r.get_mpz_t());

	count2++;
	if(Mprime == -1){
		Mprime = (-1 * (Minv%r));
		while(Mprime < 0){
			Mprime = Mprime + r;
		}
		Mprime = Mprime % r;
	}

	T = T%r;
	if(t_m_cache.find(T) == t_m_cache.end()){
		m = (T*Mprime) % r;
		t_m_cache[T] = m;
	}
	else {
		m = t_m_cache[T];
		count++;
	}
	
	t = (T+m*M)/r;
	return t >= M ? t - M : t;
}


mpz_class montgomery_exp(mpz_class x, mpz_class c, mpz_class a, mpz_class b){
	mpz_class n, r, rinv;
	n = a*b;

	r = "2";
	while(r < n){
		r = r*2;
	}
	mpz_invert(rinv.get_mpz_t(), r.get_mpz_t(), n.get_mpz_t());

	mpz_class z;
	z = "1";
	z = (z*r) % n;
	x = (x*r) % n;
	long length = bit_length(c);
	for(long i = length-1; i >= 0; i--){
		z = montgomery_reduction((z*z),r, n);
		if(mpz_tstbit(c.get_mpz_t(),i) == 1)
		{
			z = montgomery_reduction((z*x),r,n);	
		}
	}

	return montgomery_reduction(z,r,n);
	//return montgomery_reduction(z, r, n);
}


int main()
{
	Mprime = "-1";
	srand(time(NULL));

	Clock::time_point t0 = Clock::now();
	Clock::time_point t1;
	milliseconds ms;

	t0 = Clock::now();
	mpz_class x, c, n, a, b, result;
	x = generate_prime(1000);
	c = generate_prime(1000);
	a = generate_prime(300);
	b = generate_prime(300);
	n = a*b;
	t1 = Clock::now();
	ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
	cout << "Generate Prime: " <<  ms.count() << "ms\n";

	t0 = Clock::now();
	mpz_powm(result.get_mpz_t(), x.get_mpz_t(), c.get_mpz_t(), n.get_mpz_t());
	cout << "mpz_powm: " << result << endl;
	t1 = Clock::now();
	ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
	cout << "Time: " << ms.count() << "ms\n";

	t0 = Clock::now();
	cout << "S&M: " << square_and_multiply_exp(x, c, n) << endl;
	t1 = Clock::now();
	ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
	cout << "Time: " << ms.count() << "ms\n";

	t0 = Clock::now();
	cout << "CRT: " <<  chinese_remainder_exp(x,c,a,b) << endl;
	t1 = Clock::now();
	ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
	cout << "Time: " << ms.count() << "ms\n";

	t0 = Clock::now();
	cout << "S&M w/ Montgomery: " <<  montgomery_exp(x,c,a,b) << endl;
	t1 = Clock::now();
	ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
	cout << "Time: " << ms.count() << "ms\n";

	cout << count << endl;
	cout << count2 << endl;
	cout << (0.0+count)/count2 << endl;
	// t0 = Clock::now();
	// cout << "CRT w/ Montgomery: " <<  montgomery_crt_exp(x,c,a,b) << endl;
	// t1 = Clock::now();
	// ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
	// cout << "Time: " << ms.count() << "ms\n";
}