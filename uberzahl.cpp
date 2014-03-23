#include<iostream>
#include<string>
#include<vector>
#include<assert.h>
#include<cstdlib>

#include"uberzahl.h"
#define maxBits 32
#define mask 0xffffffffL
#define smallType unsigned int
#define mediumType unsigned long long
#define largeType unsigned long long

uberzahl::uberzahl ( void )
  : string_value("0")
{
  convert_to_numeric();
}

uberzahl::uberzahl ( largeType number ){
  value_vector.push_back(number);
  while ( number != 0 ){
    number = number >> maxBits;
    value_vector.push_back(number);
  }
}

uberzahl::uberzahl ( const char* number ){
  string_value = number;
  convert_to_numeric();
}
/*
uberzahl::uberzahl(const mpz_class& number){
	mpz_class scale = mask+1;
	mpz_class current = number;
	while(current>0) {
		mpz_class t = current%scale;
		value_vector.push_back(t.get_ui());
		current/=scale;
	}
}*/

uberzahl::uberzahl ( const uberzahl& number ){
  string_value = number.string_value;
  for ( size_t i=0; i < number.value_vector.size(); ++i )
    value_vector.push_back( number.value_vector[i] );
  clean_bits();
}

uberzahl::~uberzahl ( void )
{ ;;; }

const uberzahl& uberzahl::operator = ( const uberzahl& number )
{
  if ( this == &number ) return *this;
  string_value = number.string_value;
  value_vector.clear();
  for ( size_t i=0; i < number.value_vector.size(); ++i )
    value_vector.push_back( number.value_vector[i] );
  clean_bits();
  return *this;
}

void uberzahl::clean_bits ( void ){
  // remove leading 0's
  while ( value_vector.size() > 1 && !value_vector.back() )
    value_vector.pop_back();
}


uberzahl uberzahl::operator << ( smallType shift ) const
{
  // binary left shift
  uberzahl retval;
  retval.value_vector.pop_back();
  smallType largeshift = shift / maxBits;
  smallType smallshift = shift % maxBits;

  for ( size_t i=0; i < largeshift + value_vector.size() + 1; ++i )
    retval.value_vector.push_back(0);

  for ( size_t i=0; i < value_vector.size(); ++i )
    retval.value_vector[i+largeshift] = (value_vector[i] << smallshift)&mask;
  for ( size_t i=0; i < value_vector.size(); ++i ){
    mediumType workspace = value_vector[i];
    workspace = workspace >> ( maxBits - smallshift );
    retval.value_vector[i+largeshift+1] += workspace;
  }

  retval.clean_bits();
  return retval;
}

uberzahl uberzahl::operator >> ( smallType shift ) const
{
  // binary right shift
  uberzahl retval;
  smallType largeshift = shift / maxBits;
  smallType smallshift = shift % maxBits;
  for ( int i=0; i < (int)value_vector.size() - (int)largeshift - 1; ++i )
    retval.value_vector.push_back(0);

  for ( int i=0; i < (int)value_vector.size() - (int)largeshift; ++i )
    retval.value_vector[i] = value_vector[i + largeshift] >> smallshift;
  for ( int i=0; i < (int)value_vector.size() - (int)largeshift - 1; ++i ){
    mediumType workspace = value_vector[i + largeshift + 1];
    workspace = workspace << ( maxBits - smallshift );
    retval.value_vector[i] += workspace&mask;
  }

  retval.clean_bits();
  return retval;
}

uberzahl uberzahl::operator + ( const uberzahl& input ) const
{
  uberzahl x = *this;
  uberzahl y = input;
  mediumType workbench = 0;
  uberzahl retval = "0";
  retval.value_vector.clear();

  // pad extra zeros onto the left of the smaller
  while ( x.value_vector.size() != y.value_vector.size() )
    if ( x.value_vector.size() > y.value_vector.size() )
      y.value_vector.push_back(0);
    else
      x.value_vector.push_back(0);

  // perform addition operation
  for ( size_t i = 0; i < x.value_vector.size(); ++i ){
    workbench = workbench + x.value_vector[i] + y.value_vector[i];
    retval.value_vector.push_back(workbench&mask);
    workbench = workbench >> maxBits;
  }

  // add carry bit
  retval.value_vector.push_back(workbench);
  retval.clean_bits();
  return retval;
}

uberzahl uberzahl::operator - ( const uberzahl& input ) const
{
  uberzahl x = *this;
  uberzahl y = input;
  mediumType workbench = 0;
  uberzahl retval = "0";

  // constraint that left side !< right side
  if ( x < y ) return retval;
  retval.value_vector.clear();

  // pad extra zeros onto the left of the smaller
  while ( x.value_vector.size() != y.value_vector.size() )
    if ( x.value_vector.size() > y.value_vector.size() )
      y.value_vector.push_back(0);
    else
      x.value_vector.push_back(0);

  // perform subtraction
  for ( size_t i = 0; i < x.value_vector.size(); ++i ){
    workbench = -workbench + x.value_vector[i] - y.value_vector[i];
    retval.value_vector.push_back(workbench&mask);
    workbench = workbench >> maxBits;
    if ( workbench ) workbench = 1;
  }

  retval.clean_bits();
  return retval;
}

uberzahl uberzahl::operator * ( const uberzahl& y ) const
{
  size_t n = value_vector.size() - 1;
  size_t t = y.value_vector.size() - 1;
  uberzahl retval = "0";
  retval.value_vector.clear();

  smallType carry = 0;
  largeType workbench = 0;

  for ( size_t i = 0; i <= n + t + 1; ++i )
    retval.value_vector.push_back(0);
  for ( size_t i = 0; i <= t; ++i ){
    carry = 0;
    for ( size_t j = 0; j <= n; ++ j ){
      workbench = retval.value_vector[i+j] + ((largeType) value_vector[j])*y.value_vector[i] + carry;
      retval.value_vector[i+j] = workbench & ((1ULL<<maxBits)-1);
      carry = workbench >> maxBits;
    }
	retval.value_vector[n+i+1] += carry;
  }
  retval.clean_bits();
  return retval;
}

uberzahl uberzahl::operator / ( const uberzahl& number ) const
{
	uberzahl x = *this;
	uberzahl y = number;
	uberzahl q = 0ULL;
	assert( y != "0" ); // y can not be 0 in our division algorithm
	if ( x < y ) return q; // return 0 since y > x
	x.clean_bits();
  y.clean_bits();
	size_t n = x.value_vector.size() - 1;
	size_t t = y.value_vector.size() - 1;

	// initialize q to the correct size
	for ( size_t i = 0; i < n - t; ++i )
		q.value_vector.push_back(0);
	y = y << (maxBits*(n-t+1));
	for(int i=0;i<maxBits*(n-t+1);i++) {
		y = y>>1;
		q = q<<1;
		if(x>=y) {
			x = x-y;
			q = q+1;
		}
	}
	return q;
/*
  uberzahl x = *this;
  uberzahl y = number;
  uberzahl q = "0";

  assert( y != "0" ); // y can not be 0 in our division algorithm
  if ( x < y ) return q; // return 0 since y > x

  x.clean_bits();
  x.value_vector.push_back(0);
  y.clean_bits();

  size_t n = x.value_vector.size() - 1;
  size_t t = y.value_vector.size() - 1;

  // step 1 -- initialize q to the correct size
  for ( size_t i = 0; i < n - t - 1; ++i )
    q.value_vector.push_back(0);

  // step 2 -- begin long division on first 2^16-digit
  for ( size_t i=n; i > t; --i ){
//    x.value_vector.push_back(0);
	if(i+1>x.value_vector.size())
		continue;
    largeType workbench = x.value_vector[i];
    workbench = workbench << maxBits;
    workbench = workbench + x.value_vector[i-1];
    // workbench = x[i]B + x[i-1]

    q.value_vector[i-t-1] = workbench / y.value_vector[t];
    while ( x < (y << ((i-t-1)*maxBits))*q.value_vector[i-t-1] )
      q.value_vector[i-t-1] = q.value_vector[i-t-1] - 1;
    x = x - (y << ((i-t-1)*maxBits))*q.value_vector[i-t-1];
  }

  q.clean_bits();
  return q;*/
}


uberzahl uberzahl::operator % ( const uberzahl& number ) const
{
  uberzahl retval = *this - number*( *this / number );
  retval.clean_bits();
  return retval;
}

uberzahl uberzahl::operator / (smallType divisor) const
{
  uberzahl retval = *this;
  largeType current = 0;
  for(int i=value_vector.size()-1;i>=0;i--) {
    current <<= maxBits;
    current+=value_vector[i];
    retval.value_vector[i] = current/divisor;
    current %=divisor;
  }
  retval.clean_bits();
  return retval;
}

smallType uberzahl::operator % (smallType modulus) const
{
  largeType retval = 0;
  largeType coefficient = 1;
  for(int i=0;i<value_vector.size();i++) {
    retval+=coefficient*value_vector[i];
    retval%=modulus;
    coefficient<<=maxBits;
    coefficient%=modulus;
  }
  return retval;
}

// convert the stored numeric_value into a string
std::string uberzahl::convert_to_string ( void ) const
{
  if(*this == "0")
    return "0";

  uberzahl temp = *this;
  std::string reversed = "";
  while(temp>"0") {
    reversed+=temp%10+'0';
    temp=temp/10;
  }

  std::string retval = "";
  for(int i=reversed.size()-1;i>=0;i--)
    retval+=reversed[i];

  return retval;
}

// take the string_value and convert it into a numeric_value
void uberzahl::convert_to_numeric ( void ){
  std::string workbench = string_value;
  value_vector.clear();
  smallType numeric_value = 0;
  size_t bits = 0;

  while ( workbench.length() ){
    ++bits;

    size_t len = workbench.length();
    numeric_value = numeric_value >> 1;
    if ( workbench[len-1] % 2 == 1 ){
      numeric_value = numeric_value | ( 1 << (maxBits-1) );
      workbench[len-1] = workbench[len-1] ^ 1;
    }

    for ( size_t i = workbench.length(); i > 0; --i ){
      // constant consistancy check! if this algorithm fails I want to know
      assert( workbench[i-1] >= '0' && workbench[i-1] <= '9' );

      if ( workbench[i-1] % 2 ){ // odd decimal digit
        workbench[i] += 0x05; // adjust the lower term up
      }
      workbench[i-1] -= '0'; // cut off 0x30 ascii
      workbench[i-1] = workbench[i-1] >> 1; // divide by 2
      workbench[i-1] += '0'; // add back the ascii
    }

    while ( workbench[0] == '0' )
      workbench = workbench.substr(1,len);

    // we have run out of room in numeric_value
    if ( bits % maxBits == 0 ){
      value_vector.push_back(numeric_value);
      numeric_value = 0;
    }
  }

  // take care of the remaining bits
  numeric_value = numeric_value >> (-bits % maxBits);
  value_vector.push_back(numeric_value);
}

std::ostream& operator << ( std::ostream& ost, const uberzahl& number ){
  /*  ost << "string : " << number.convert_to_string() << std::endl << "base-2**" << maxBits << " : ";
      for ( size_t i = 0; i < number.value_vector.size(); ++ i )
      ost << number.value_vector[i] << ' ';
      ost << std::endl;
      */
  ost << number.convert_to_string();
  return ost;
}

// Comparator operators
// Removed the padding necessity to allow them to be 
// true const& passes
bool uberzahl::operator <= (const uberzahl& rhs) const
{
  size_t lhs_size = value_vector.size();
  size_t rhs_size = rhs.value_vector.size();
  if ( lhs_size < rhs_size )
    return rhs >= *this;

  for ( size_t i=rhs_size; i < lhs_size; ++i )
    if ( value_vector[i] > 0 )
      return false;

  for ( size_t i=rhs_size; i > 0; --i )
    if ( value_vector[i-1] > rhs.value_vector[i-1] )
      return false;
    else if ( value_vector[i-1] < rhs.value_vector[i-1] )
      return true;

  return true; 
}

bool uberzahl::operator >= (const uberzahl& rhs) const
{
  size_t lhs_size = value_vector.size();
  size_t rhs_size = rhs.value_vector.size();
  if ( lhs_size < rhs_size )
    return rhs <= *this;

  for ( size_t i=rhs_size; i < lhs_size; ++i )
    if ( value_vector[i] > 0 )
      return true;

  for ( size_t i=rhs_size; i > 0; --i )
    if ( value_vector[i-1] < rhs.value_vector[i-1] )
      return false;
    else if ( value_vector[i-1] > rhs.value_vector[i-1] )
      return true;

  return true; 
}

bool uberzahl::operator < ( const uberzahl& rhs ) const
{
  return !( *this >= rhs );
}

bool uberzahl::operator > ( const uberzahl& rhs ) const
{
  return !( *this <= rhs );
}

bool uberzahl::operator == ( const uberzahl& rhs ) const
{
  return ( *this >= rhs ) && ( *this <= rhs );
}

/*bool uberzahl::operator == ( const mpz_class& rhs ) const
{
	mpz_class scale= mask+1;
	mpz_class rhstemp = rhs;
	for(int i=0;i<value_vector.size();i++) {
		if(rhstemp%scale != value_vector[i])
			return false;
		rhstemp/=scale;
	}
	return rhstemp == 0;
}*/


bool uberzahl::operator != ( const uberzahl& rhs ) const
{
  return !( *this == rhs );
}

// Bitwize operators done with true pass by reference
uberzahl uberzahl::operator | ( const uberzahl& rhs ) const
{
  if ( value_vector.size() > rhs.value_vector.size() )
    return rhs | *this;

  uberzahl retval = "0";
  retval.value_vector.pop_back();
  for ( size_t i=0; i < value_vector.size(); ++i ){
    smallType workbench = value_vector[i] | rhs.value_vector[i];
    retval.value_vector.push_back( workbench );
  }
  for ( size_t i=value_vector.size(); i < rhs.value_vector.size(); ++i )
    retval.value_vector.push_back( rhs.value_vector[i] );

  retval.clean_bits();
  return retval;
}

uberzahl uberzahl::operator & ( const uberzahl& rhs ) const
{
  if ( value_vector.size() > rhs.value_vector.size() )
    return rhs & *this;

  uberzahl retval = "0";
  retval.value_vector.pop_back();
  for ( size_t i=0; i < value_vector.size(); ++i ){
    smallType workbench = value_vector[i] & rhs.value_vector[i];
    retval.value_vector.push_back( workbench );
  }

  retval.clean_bits();
  return retval;
}

uberzahl uberzahl::operator ^ ( const uberzahl& rhs ) const
{
  if ( value_vector.size() > rhs.value_vector.size() )
    return rhs ^ *this;

  uberzahl retval = "0";
  retval.value_vector.pop_back();
  for ( size_t i=0; i < value_vector.size(); ++i ){
    smallType workbench = value_vector[i] ^ rhs.value_vector[i];
    retval.value_vector.push_back( workbench );
  }
  for ( size_t i=value_vector.size(); i < rhs.value_vector.size(); ++i )
    retval.value_vector.push_back( rhs.value_vector[i] );

  retval.clean_bits();
  return retval;
}

uberzahl uberzahl::random ( mediumType bits ){
  assert( bits > 0 );
  value_vector.clear();

  smallType shortbits = bits % maxBits;
  smallType longbits = bits / maxBits;

  for ( size_t i = 0; i <= longbits; ++i )
    value_vector.push_back( rand() );

  if ( shortbits ){
    value_vector[longbits] = value_vector[longbits] >> (maxBits - shortbits);
    value_vector[longbits] = value_vector[longbits] | (1 << shortbits - 1);
  }
  else{
    value_vector.pop_back();
    value_vector[longbits-1] = value_vector[longbits-1] | (1 << (maxBits - 1));
  }

  return *this;
}

uberzahl uberzahl::inverse ( const uberzahl& b) const
{
  // computes (*this)^(-1) mod b
  if(*this=="1")
    return "1";
  std::pair<std::pair<uberzahl,uberzahl>,bool> inv = inverse(*this,b);
  inv.first.first = inv.first.first % b;
  if(inv.second && inv.first.first!="0")
    return b-inv.first.first;
  else
    return inv.first.first;
}

uberzahl uberzahl::gcd( const uberzahl& that ) const
{
  // computes gcd(*this, that)
  if ( *this < that )
    return that.gcd( *this );
  else if( that == "0" )
    return *this;
  else 
    return that.gcd( *this % that );
}

std::pair<std::pair<uberzahl,uberzahl>,bool> uberzahl::inverse ( const uberzahl& a, const uberzahl& b) const
{
  // helper function for inverse
  uberzahl nexta = b%a, coeff = b/a;
  if(nexta == "0")
    return std::make_pair(std::make_pair("0","0"),true);
  else if(nexta == "1")
    return std::make_pair(std::make_pair(coeff,"1"),true);
  std::pair<std::pair<uberzahl,uberzahl>,bool> nextinv = inverse(nexta,a);
  return std::make_pair(std::make_pair(nextinv.first.second+nextinv.first.first*coeff,nextinv.first.first),!nextinv.second);
}

smallType uberzahl::bit ( mediumType n ) const
{
  // returns the nth bit (0 indexed)
  mediumType largeBit = n / maxBits;
  smallType smallBit = n % maxBits;

  if ( largeBit >= value_vector.size() )
    return 0;

  smallType bits = value_vector[largeBit];
  return (bits>>smallBit) & 1;
}

smallType uberzahl::bitLength ( void ) const
{
  // returns the bit length (size) of a uberzahl
  for(int i=value_vector.size()-1;i>=0;i--) {
    if(value_vector[i]!=0) {
      largeType k=1;
      int j=0;
      while(k<=value_vector[i]) {
        k = k<<1;
        j++;
      }
      return i*maxBits+j;
    }
  }
  return 0;
}

uberzahl uberzahl::exp ( const uberzahl& exponent ) const
{
  if ( exponent == "0" ) // exponent of 0
    return 1;
  else if ( exponent == "1" ) // exponent of 1
    return *this;

  if ( (exponent&1) == "1" ) // odd exponent
    return (this->exp(exponent^1)) * (*this);
  else {
    uberzahl tmp = this->exp(exponent >> 1);
    return tmp * tmp;
  }
}

// this is the recursive version, you need to implement the
// iterative version in your code!
uberzahl uberzahl::expm( const uberzahl& n, const uberzahl& mod ) const
{
  if ( n == "0" )
    return 1;
  else if ( n == "1" )
    return ( *this % mod );
  
  if ( (n&1) == "1" )
    return ( (this->expm(n^1, mod)) * (*this) ) % mod;
  else {
    uberzahl tmp = this->expm(n>>1, mod);
    return (tmp * tmp) % mod;
  }
}

uberzahl random ( const uberzahl& a, const uberzahl& b )
{
  if ( a > b ) return random( b, a );

  uberzahl retval;
  retval.random( b.bitLength() + 1 );
  retval = (retval % ( b-a )) + a;
  return retval;
}

bool rabinmiller ( const uberzahl& n, unsigned int k ){
  // take care of corner cases 1,2,3 and even
  if ( n < 2 ) return false;
  else if ( n < 4 ) return true;
  else if ( (n&1) == "0" ) return false;

  uberzahl d = (n - 1) >> 1;
  unsigned int s = 1;
  while ( (d&1) == "0" ){
    d = d >> 1;
    s = s + 1;
  }

  // witness loop
  uberzahl a;
  uberzahl x;
  for ( unsigned int i = 0; i < k; ++i ){
    a = random( 2, n-2 );
    x = a.expm(d,n);

    if ( x == "1" || x == n-1 ) continue;

    for ( unsigned int j = 0; j < s - 1; ++j ){
      x = x.expm(2,n);
      if ( x == "1" ) return false; // composite
      if ( x == n-1 ) break;
    }

    if ( x == n-1 ) continue;
    
    return false; // composite
  }

  return true; // probably prime
}

uberzahl nextprime ( const uberzahl& n, unsigned int accuracy = 50){
  // checks every odd number greater or equal to n for primality
  uberzahl retval = n;
  if ( (retval&1) == "0" ) retval = retval + 1;
  while ( !rabinmiller(retval, accuracy) )
    retval = retval + 2;
  return retval;
}
