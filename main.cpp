#include <iostream>
#include <cstdlib>
using namespace std;
#include "uberzahl.h"

int main( void ){
  srand(time(0));
  uberzahl x = 1;
  cout << x << " : " << x.bitLength() << endl;
  x = 3;
  cout << x << " : " << x.bitLength() << endl;
  x = 1023;
  cout << x << " : " << x.bitLength() << endl;
  x = 1025;
  cout << x << " : " << x.bitLength() << endl;
//  uberzahl y = "4";
  
  /*unsigned int x = -2;
  uberzahl a = x;
  uberzahl b = "5";
  cout << (a + 1);
  cout << (a + 1) * 2;

  uberzahl c = "1234";
  uberzahl d = "5678";
  uberzahl e = "123456789123456789";
  uberzahl f = "123456789123456789"; 

  if (c == d) cout << "c == d" << endl; //won't print
  if (e == f) cout << "e == f" << endl; //will print
  if (c < d) cout << "c < d" << endl; //will print
  if (c <= d) cout << "c <= d" << endl; //will print
  if (e > c) cout << "e > c" << endl; //will print
  if (f >= e) cout << "f >= e" << endl; //will print
  if (d >= f) cout << "d >= f" << endl; //won't print

  uberzahl g = "0000000000000000000000000000000000000000123456";
  uberzahl h = "123456"; 
  uberzahl i = "123457"; 

  if (g == h) cout << "g == h" << endl; //will print
  if (g >= h) cout << "g >= h" << endl; //will print
  if (i >= g) cout << "i >= g" << endl; //will print

  uberzahl j = "1239808213902318901089231283901238904392834982374"; 
  uberzahl k = "42949672964294967296"; 
  uberzahl l = "1239808213902318901089231283901238904392834982373";
  if (k <= j) cout << "k <= j" << endl;  
  if (j > k) cout << "j > k" << endl;
  if (k < j) cout << "k < j" << endl;
  if (l < j) cout << "l < j" << endl;
  if (l < l) cout << "l < l" << endl; 

  cout << j; 
  cout << k;*/
}

