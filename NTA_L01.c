#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

//Realisation: tests on primeness, integer factorisation algorithms, "smart" function to determine best algorithm
//and find any prime

/*Part 1: primeness tests 
1.1 Solovey-Strassen test
1.2 Miller-Rabine test
*/

int EEA(int a, int n, int* u, int* v);
bool SS_test(uint64_t p, int k); //two args: number & test precision

static inline uint64_t ring_mul(uint64_t n1, uint64_t n2, uint64_t mod);
static inline uint64_t ring_pow(uint64_t base, uint64_t exponent, uint64_t mod);

int main(){
    printf("\nLehho Dorlw!");
    printf("\nTasks: \n1. Primeness tests \n2. Factorization algs");
    printf("\nRing mul & Ring pow tests");
    uint64_t a = ring_mul(11,11,7);
    uint64_t b = ring_pow(a,a,7);
    printf("\n11 * 11 (mod 7) = %llu \n %llu ^ %llu (mod 7) = %llu", a, a, a, b);
    
    return 0;
}

//Extended EA with Bezout`s coefficients (u,v saved in memory through pointer type var`s)
int EEA(int a, int n, int* u, int* v){
    if(n == 0){
        *u = 1;
        *v = 0;
        return a;
    }
    int u_i = 0, v_i = 0;
    int gcd = EEA(n, a%n, &u_i, &v_i);
    *u = v_i;
    *v = u_i - (a/n) * v_i;
    return gcd;
}

bool SS_test(uint64_t p, int k){
    if(p <= 2){
        return true;
    }
    return 0;
}

static inline uint64_t ring_mul(uint64_t n1, uint64_t n2, uint64_t mod){
    return (uint64_t)((__uint128_t) n1 * n2 % mod);
}

static inline uint64_t ring_pow(uint64_t base, uint64_t exponent, uint64_t mod){
    uint64_t rsd = 1;
    while(exponent){
        if(exponent & 1) rsd = ring_mul(rsd, base, mod);
        base = ring_mul(base, base, mod);
        exponent >>= 1;
    }
    return rsd;
}