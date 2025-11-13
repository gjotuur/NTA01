#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>

//Realization: tests on primeness, integer factorization algorithms, "smart" function to determine best algorithm
//and find any prime

/*Part 1: primeness tests 
1.1 Solovey-Strassen test
1.2 Miller-Rabine test
*/

uint64_t rnd_between(uint64_t x, uint64_t y);                           //simple RNG, to do it without srand() and rand
uint64_t EEA(uint64_t a, uint64_t n, uint64_t* u, uint64_t* v);         //Extended Euclidean alg 
bool SS_test(uint64_t p, int k);                                        //Solovey-Strassen test

static inline uint64_t ring_mul(uint64_t n1, uint64_t n2, uint64_t mod);            //actually, it`s mulmod, but i don`t like the name
static inline uint64_t ring_pow(uint64_t base, uint64_t exponent, uint64_t mod);    //actually, it`s powmod, but i don`t like the name as well :D

int main(){
    printf("\nLehho Dorlw!");
    printf("\nTasks: \n1. Primeness tests \n2. Factorization algs");
    printf("\nRing mul & Ring pow tests");
    uint64_t a = ring_mul(11,11,7);
    uint64_t b = ring_pow(a,a,7);
    printf("\n11 * 11 (mod 7) = %llu \n %llu ^ %llu (mod 7) = %llu", a, a, a, b);
    printf("\nTest for simple RNG");
    uint64_t* rnd_list = malloc(10 * sizeof(uint64_t));
    if (!rnd_list) {
        fprintf(stderr, "Ne piwlo, sprobuj shche!");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i<10; i++){
        rnd_list[i] = rnd_between(1, 65536);
    }
    printf("\n");

    for(int i = 0; i < 10; i++){
        printf("%2d. %6llu\n", i+1, rnd_list[i]);
    }

    return 0;
}

//Extended EA with Bezout`s coefficients (u,v saved in memory through pointer type var`s)
uint64_t EEA(uint64_t a, uint64_t n, uint64_t* u, uint64_t* v){
    if(n == 0){
        *u = 1;
        *v = 0;
        return a;
    }
    uint64_t u_i = 0, v_i = 0;
    uint64_t gcd = EEA(n, a%n, &u_i, &v_i);
    *u = v_i;
    *v = u_i - (a/n) * v_i;
    return gcd;
}

bool SS_test(uint64_t p, int k){
    uint64_t counter = 0;
    bool conclusion = false;
    if(p <= 2){
        return true;
    }
    uint64_t seed = rnd_between(1, p);
    uint64_t u,v;
    uint64_t n1 = EEA(seed, p, &u, &v);
    if(n1 != 1){
        return false;
    }
    //step 2: Eulers pseudoprime by base (sry, PC battery is almost empty...)


    return 0;
}

static inline uint64_t ring_mul(uint64_t n1, uint64_t n2, uint64_t mod){ //max 64 bit, less stress
    return (uint64_t)((__uint128_t) n1 * n2 % mod);
}

static inline uint64_t ring_pow(uint64_t base, uint64_t exponent, uint64_t mod){ //same, less stress
    uint64_t rsd = 1;
    while(exponent){
        if(exponent & 1) rsd = ring_mul(rsd, base, mod);
        base = ring_mul(base, base, mod);
        exponent >>= 1;
    }
    return rsd;
}

uint64_t rnd_between(uint64_t x, uint64_t y){
    uint64_t range = y - x + 1;
    uint64_t r, lim = UINT64_MAX - (UINT64_MAX % range);
    do {
        r = ((uint64_t)rand() << 32) ^ rand();
    } while (x >= lim);
    return x + r % range;
}