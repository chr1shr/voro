#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gmpxx.h>
using namespace std;

int main() {
	mpz_class a,b,c;

	a="123456789101112";
	b="1234123412341234";
	c=a*b;
	c%=101;

	mpz_out_str(stdout,10,c.get_mpz_t());puts("");
}
