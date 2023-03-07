#include <stdlib.h>
#include <stdio.h>

#define TRUE 1
#define FALSE 0

typedef unsigned long long int  ll;

ll input_vals[] = {
/*0x9bd7fb68, 0xd9023253, 0xe643d86c, 0xc4e09765,
0x6407abe0, 0x39df5ef7, 0x8814bcea, 0x5eec353c,
0xee03b859, 0x74cc8524, 0x2d255465, 0x682d3507,
0x4acfc99a, 0x6765868d, 0xd2060ecd, 0x1e5f2f94,
0x3aa98973, 0x65fd6759, 0x322ef3a6, 0x6941babc,
0x1843514d, 0x6649b2c3, 0xfd0bd93f, 0xc5c534bc,
0x8049987a, 0x867bccc3, 0x63865ad3, 0xe1b3fe56,
0x74480ed1, 0x75210657, 0x8c59a3d6, 0x84f9bbe5,
0xd43fb846, 0xd53fd2e0, 0x913c5d16, 0xdd621b45,
0xd49ea60f, 0x71def72f, 0x8a98b156, 0xba6964c2,
0xa8c6846f, 0xcf095848, 0x3890d905, 0xcb380eaf,
0x73c6683e, 0x43ef3bea, 0x5c541c30, 0xea7919eb,
0x5b2cca2d, 0x3831cdcf, 0x2d1ca408, 0xdf638b64,
0x4dc74c9b, 0xa9215c0a, 0x5530d35,  0x3bab720f
};
*/
0xad79d417, 0x87358ee5, 0x68c0ed8a, 0x6d73eb1d,
0x371d326d, 0xa4857029, 0xeeadc83c, 0x154c34d,
0x6cb328b4, 0xe0c7373d, 0xb4729020, 0x7d8e2d90,
0xcf19a0fb, 0xfbcefa9,  0x715070c1, 0xa9e467ad,
0xa31a11bf, 0x3c256122, 0xe2436417, 0x2ce03abc,
0x527a983f, 0xc77bdcbf, 0x2a62c559, 0xe9a5930a,
0x50fee11a, 0x54ff4b83, 0x44f1745b, 0x75886d16,
0xd1203b47, 0xd484195d, 0x31668f59, 0x13e6ef96,
0x12661ea,  0x19ef330e, 0x8e196b4e, 0x86cff959,
0x610d4537, 0x9926cb0c, 0xf42f64fd, 0x1714d2f3,
0xeaa625cf, 0x97f59d64, 0xc8bbce99, 0xa506eaf0,
0x2b3f2669, 0x9d961a35, 0x48183b35, 0x797cbe53,
0xb80ee164, 0xd3321493, 0xb4955195, 0xa0b4d41c,
0x901eaf81, 0xe77d7bdd, 0x2052f3a8, 0x7bb0d4f2,
0x6f5feda1, 0x6408c94e, 0x990f61b3, 0x13825d97}; 

ll mulmod(ll a, ll b, ll mod)
{
    ll x = 0;
    ll y = a % mod;

    while (b > 0) {
        if (b % 2 == 1) {
            x = (x + y) % mod;
        }
        y = (y * 2) % mod;
        b /= 2;
    }

    return x % mod;
}

ll modulo(ll base, ll exponent, ll mod)
{
    ll x = 1;
    ll y = base;

    while (exponent > 0) {
        if (exponent % 2 == 1)
            x = (x * y) % mod;
        y = (y * y) % mod;
        exponent = exponent / 2;
    }

    return x % mod;
}


int  miller(ll p, int iteration)
{
    int i;

    if (p < 2) {
        return FALSE;
    }

    if (p != 2 && p % 2==0) {
        return FALSE;
    }

    ll s = p - 1;
    while (s % 2 == 0) {
        s /= 2;
    }

    for (i = 0; i < iteration; i++) {
        ll a = rand() % (p - 1) + 1, temp = s;
        ll mod = modulo(a, temp, p);
        while (temp != p - 1 && mod != 1 && mod != p - 1) {
            mod = mulmod(mod, mod, p);
            temp *= 2;
        }

        if (mod != p - 1 && temp % 2 == 0) {
            return FALSE;
        }
    }

    return TRUE;
}

int main (int argc, char **argv)
{
    int i;
    int iteration = 5;
    //ll input_num = 85022837187228775553374836620586910990389778457989956277594338003496152540707;// 2910442519;
    
    for (i = 0; i < 56; i++) {
        if (miller(input_vals[i], iteration))
            printf(" %llu is a prime \n", input_vals[i]);
        else
            printf(" %llu is not a prime \n", input_vals[i]);
    }
/*
        if (miller(input_num[i], iteration))
            printf(" %llu is a prime \n", input_num);
        else
            printf(" %llu is not a prime \n", input_num);
*/    
    return 0;
}

