#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gmpxx.h>
using namespace std;

int main() {
	double q=3;
	mpq_class a,b;

	a="2/4";
	b=q;
	b*=b;
	q=b.get_d();
	printf("%g\n",q);

	mpz_out_str(stdout,10,b.get_num_mpz_t());puts("");
	mpz_out_str(stdout,10,b.get_den_mpz_t());puts("");
}
