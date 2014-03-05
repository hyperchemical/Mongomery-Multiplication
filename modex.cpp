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
#include <vector>
#include <stdexcept> 

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

void assert_all_equal(vector<mpz_class>& vec)
{
	for(int i = 0; i < vec.size(); i++){
		for(int j = i; j < vec.size(); j++)
		{
			if(vec[i]!=vec[j])
			{
				throw runtime_error("ERROR: Assert Equal Failed");
			}
		}
	}
}

int main()
{
	vector<string> function_names = {"S&M", "CRT", "S&M w/ Montgomery"};
	vector<mpz_class (*)(mpz_class, mpz_class, mpz_class, mpz_class)> functions
		= {square_and_multiply_exp, chinese_remainder_exp, montgomery_exp};
	vector<mpz_class> values;

	Mprime = "-1";
	srand(time(NULL));


	Clock::time_point t0 = Clock::now();
	Clock::time_point t1;
	milliseconds ms;

	t0 = Clock::now();
	mpz_class x, c, n, a, b, result;
	x = "1198149879571789438280899179054834842764382882709690800846302623909981933857960623903110405022053847724692601166911024953741917161287312694381224661764272148719173353701787070197506564165415050566378224397189522616401305529822509494871407827542307560920013284182619715357187579547729841287118072854134618923926232453293472017981796136430152429351051237402639294447950148719448597229833401403703285811230637539608443920424962853861861070612627261727458227848218588380026797993061905975890910591813615484632501667223344057718152803309542158891328789507738222946519390620355456258609472150146854368643870933081765786285901356077609120335719652457795816120622368839794307848372934612612714317602700217815610874463286425213465028961645857286284653119163799751763037509103110729418441641873973321577694109943899159848403056280877964306131529409571892842171058552572325189595634476637871926763075755939993146011";
	c = "1056727158623602207244671395916435520435082370915576644965062078868130345228438443824468207124383595755531818235058013694654070963931934093961965780353571720379766648260892966451516596772127387096402982376058972991538704667031228758971022402543943586174493584363331529631588104468792817679521820624811448338058879831601582015261359811376536828272804966173150814262097853827343922095299192515322908057565463037620249476781082512025178659398610120130363679179149799172722684503613505123313101878116246628001662193757973205518739251699588845015000713607405078476078887666670107890545530664437240293383821849425744170397510815369258904407897452926143468494150837451623912189777174687964990564621691874190501022152064978790032918152697869690667232196410399328138416478193324575637757359218185069006403281915856147986572212662598593516473959100173560063214823112104571953647857814646712866314838175451461506687";
	a = generate_prime(300);
	b = generate_prime(300);
	n = a*b;
	t1 = Clock::now();
	ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
	cout << "Generate Prime: " <<  ms.count() << "ms\n";

	t0 = Clock::now();
	mpz_powm(result.get_mpz_t(), x.get_mpz_t(), c.get_mpz_t(), n.get_mpz_t());
	values.push_back(result);
	cout << "mpz_powm: " << result << endl;
	t1 = Clock::now();
	ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
	cout << "Time: " << ms.count() << "ms\n";

	for(int i = 0; i < functions.size(); i++){
		t0 = Clock::now();
		result = functions[i](x,c,a,b);
		cout << function_names[i] << " " << result << endl;
		values.push_back(result);
		t1 = Clock::now();
		ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
		cout << "Time: " << ms.count() << "ms\n";
	}

	assert_all_equal(values);

	// t0 = Clock::now();
	// cout << "S&M: " << square_and_multiply_exp(x, c, a, b) << endl;
	// t1 = Clock::now();
	// ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
	// cout << "Time: " << ms.count() << "ms\n";

	// t0 = Clock::now();
	// cout << "CRT: " <<  chinese_remainder_exp(x,c,a,b) << endl;
	// t1 = Clock::now();
	// ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
	// cout << "Time: " << ms.count() << "ms\n";

	// t0 = Clock::now();
	// cout << "S&M w/ Montgomery: " <<  montgomery_exp(x,c,a,b) << endl;
	// t1 = Clock::now();
	// ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
	// cout << "Time: " << ms.count() << "ms\n";

	cout << count << endl;
	cout << count2 << endl;
	cout << (0.0+count)/count2 << endl;
	// t0 = Clock::now();
	// cout << "CRT w/ Montgomery: " <<  montgomery_crt_exp(x,c,a,b) << endl;
	// t1 = Clock::now();
	// ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
	// cout << "Time: " << ms.count() << "ms\n";
}