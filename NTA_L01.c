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
bool SS_test(uint64_t p, int k, uint64_t* counter);                     //Solovey-Strassen test
void test_number(uint64_t number, int precision);                       //Test with consol message
int8_t JCB(uint64_t a, uint64_t n);                                     //Jacobi symbol calculation

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
    int j1 = JCB(1719, 78936);
    int j2 = JCB(25799, 915133);
    printf("\n(1719 / 78936) = %d \n(25799 / 915133) = %d", j1, j2);
    uint64_t test1_num =  858135414131218189;
    test_number(test1_num, 100);
    uint64_t test2_num =  123456789012345671;
    test_number(test2_num, 200);
    uint64_t test3_num =  987654321098765429;
    test_number(test3_num, 200);
    test1_num = 123456789012345678;
    test_number(test1_num, 100);
    test2_num = 987654321098765433;
    test_number(test2_num, 100);

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

bool SS_test(uint64_t p, int k, uint64_t* counter){
    int c = 0;
    if (p == 2) {
        if(counter) *counter = k;  // всі k ітерацій "пройшли"
        return true;
    }
    if((p%2) == 0 || p <2) {
        if(counter) *counter = 0;  // не виконано жодної ітерації
        return false;
    }
    while(c < k){
        uint64_t seed = rnd_between(2, p - 1);
        uint64_t u,v;
        uint64_t n1 = EEA(seed, p, &u, &v);
        if(n1 != 1) {
            if(counter) *counter = c;
            return false;
        }
        uint64_t condition_a = ring_pow(seed, ((p - 1) / 2), p);
        int64_t condition_b = (int64_t)JCB(seed, p);

        //Normalize
        int64_t normalized_a;
        if (condition_a == p - 1) {
            normalized_a = -1;
        } else {
            normalized_a = (int64_t)condition_a;
        }
        if(normalized_a != condition_b) {
            if(counter) *counter = c;
            return false;
        }
        c++;
    }
    if(counter) *counter = (uint64_t)c;
    return true;
}

void test_number(uint64_t number, int precision){
    uint64_t k;
    bool conclusion = SS_test(number, precision, &k);
    double probability = 1 - pow(2, (double)(-k));
    if(conclusion){
        printf("\nNumber %20llu is prime with probability %2.2lf", number, probability*100);
    } else {
        printf("\nThe number %20llu is not prime", number);
    }
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

//simple RNG (built-in srand || rand sucks)
uint64_t rnd_between(uint64_t x, uint64_t y){
    uint64_t range = y - x + 1;
    uint64_t r, lim = UINT64_MAX - (UINT64_MAX % range);
    do {
        r = ((uint64_t)rand() << 32) ^ (uint64_t)rand();
    } while (r >= lim);
    return x + r % range;
}

//Jacobi sign (symbol?) calculation for SS test
int8_t JCB(uint64_t a, uint64_t n){
    int64_t jcb = 1;
    if(n <= 0 || (n % 2) == 0) return 0;    
    if(a == 1) return 1;                    
    a %= n;                                 
    while(a != 0){
        while((a % 2) == 0){
            a /= 2;
            if(n % 8 == 3 || n % 8 == 5){
                jcb = -jcb;
            }
        }
        uint64_t tmp = a;
        a = n;
        n = tmp;
        if(a % 4 == 3 && n % 4 == 3){
            jcb = -jcb;
        }
        a %= n;
    }
    return (n == 1) ? (int8_t)jcb : 0;
}