#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>
#include <time.h>

//Realization: tests on primeness, integer factorization algorithms, "smart" function to determine best algorithm
//and find any prime

/*Part 1: primeness tests 
1.1 Solovey-Strassen test
1.2 Miller-Rabine test
*/
//Structure to save factorized number (<=64 bit) and powers of primes, max size of array of those struct is ~50 (<2^64-1)
typedef struct{
    uint64_t prime;
    int power;
}fint64_t;

uint64_t rnd_between(uint64_t x, uint64_t y);                           //simple RNG, to do it without srand() and rand
uint64_t EEA(uint64_t a, uint64_t n, uint64_t* u, uint64_t* v);         //Extended Euclidean alg 
bool SS_test(uint64_t p, int k, uint64_t* counter);                     //Solovey-Strassen test
void test_number(uint64_t number, int precision);                       //Test with consol message
int8_t JCB(uint64_t a, uint64_t n);                                     //Jacobi symbol calculation

static inline uint64_t ring_mul(uint64_t n1, uint64_t n2, uint64_t mod);            //actually, it`s mulmod, but i don`t like the name
static inline uint64_t ring_pow(uint64_t base, uint64_t exponent, uint64_t mod);    //actually, it`s powmod, but i don`t like the name as well :D

//trial divisions
uint64_t* gen_primes(uint64_t n);                                                                      //prime numbers generation
uint64_t* td_method(uint64_t number, uint64_t* primes_array, size_t array_len, int* final_count);      //trial divisions method

//Rho-method
uint64_t rho_f(uint64_t number, uint64_t mod);                                                         //Internal function to find x_i, could be undeclared as prototype
uint64_t rho_g(uint64_t number, uint64_t mod);                                                         //internal function to find y_i, could be undeclared as prototype
uint64_t rho_factor(uint64_t number, bool method);                                                     //Pollard rho-method

//Brillhart - Morrison method
uint64_t* bm_factorbase(uint64_t number, uint64_t* count);               //factor-base generation
uint64_t* bm_cfrac(uint64_t number, uint64_t* b_size);                   //Chain fractions B-smooth array generation

int main(){
    srand(time(NULL));
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
    test2_num = 2500744714570633849;
    test_number(test2_num, 100);
    printf("\nMax value of int64_t is  %llu, 19 digits", INT64_MAX);
    //let`s generate some primes!
    uint64_t* some_primes = gen_primes(1000);
    for(int i = 0; i < 1000; i++){
        if(((i + 1) % 10) == 0 || i == 0) printf("\n");
        printf("%6llu",some_primes[i]);
    }
    int mistakes = 0;
    for(int i = 0; i < 1000; i++){
        uint64_t ctr;
        if(!SS_test(some_primes[i], 10, &ctr)) mistakes++;
    }
    printf("\nMistakes made: %d", mistakes);
    int max;
    uint64_t* lets_test = td_method(479001600, some_primes, 1000, &max);
    for(int i = 0; i < max; i++){
        printf("\nDivisor %d is %llu", i+1, lets_test[i]);
    }
    int max2;
    uint64_t* lets_test2 = td_method(6227020800, some_primes, 1000, &max2);
    for(int i = 0; i < max2; i++){
        printf("\nDivisor %d is %llu", i+1, lets_test2[i]);
    }

    printf("\nPollard rho-method, testing number: 1449863225586482579");
    uint64_t test4_num =  1449863225586482579;
    uint64_t factor = rho_factor(test4_num, true);
    if(factor) printf("\nFound factor: %llu, test: %llu", factor, test4_num % factor);
    printf("\nBrillhart-Morrison method: factor-base\n");
    uint64_t fb_count;
    uint64_t* fb = bm_factorbase(test4_num, &fb_count);
    for(uint64_t i = 0; i < fb_count; i++){
        printf("%10llu  ", fb[i]);
        (i == 9 || (i+1)%10 == 0) ? printf("\n") : 0;
    }
    printf("\nFactor base size is %llu", fb_count);
    /*uint64_t b_size;
    uint64_t* b_base = bm_cfrac(test4_num, &b_size);
    for(uint64_t i = 0; i < b_size; i++){
        printf("%10llu   ", b_base[i]);
        (i == 9 || (i+1)%10 == 0) ? printf("\n") : 0;    
    }*/

    return 0;
}

//Extended EA with Bezout`s coefficients (u,v saved in memory through pointer type var`s)
uint64_t EEA(uint64_t a, uint64_t n, uint64_t* u, uint64_t* v){
    if(n == 0){
        if(u) *u = 1;
        if(v) *v = 0;
        return a;
    }
    uint64_t u_i = 0, v_i = 0;
    uint64_t gcd = EEA(n, a%n, &u_i, &v_i);
    if(u) *u = v_i;
    if(v) *v = u_i - (a/n) * v_i;
    return gcd;
}

//Solovey-Strassen test body
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

//direct function for testing numbers with results
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

//For trial divisions method it would be good to generate some number of primes. We`ll limit a generation to first 10k primes
//using Erathosphene sieve
uint64_t* gen_primes(uint64_t n){
    if(n > 1000000 || n <= 0){
        printf("\nPrimes array not generated - limit exceeded");
        return NULL;
    }
    //estimate the potential size of sieve to create array
    uint64_t estimate = (uint64_t)((n * log(n)) * 1.2);
    bool* sieve = malloc(estimate*sizeof(bool));
    if(!sieve){
        printf("Memory allocation error");
        return NULL;
    }
    uint64_t* primes = malloc(n * sizeof(uint64_t));
    if(!primes){
        printf("Memory allocation error");
        return NULL;
    }
    //init sieve
    sieve[0] = false;
    sieve[1] = false;
    for(int i = 2; i < estimate; i++){
        sieve[i] = true;
    }
    //evaluate sieve
    uint64_t max_num = (uint64_t) sqrt(estimate);
    for(uint64_t p = 2; p <= max_num; p++){
        if(sieve[p]){
            // Викреслюємо всі кратні p, починаючи з p²
            for(uint64_t k = p * p; k < estimate; k += p){
                sieve[k] = false;
            }
        }
    }
    //insert values into an array
    uint64_t count = 0;
    for(uint64_t i = 2; i < estimate && count < n; i++){
        if(sieve[i]){
            primes[count] = i;
            count++;
        }
    }
    free(sieve);
    if(count < n) {
        if(count == 0){
            free(primes);
            printf("\nPrimes not found!");
        }
    }
    uint64_t* resized = realloc(primes, count * sizeof(uint64_t));

    if(resized){
        primes = resized;
    }
    return primes;
}

//trial divisions method
uint64_t* td_method(uint64_t number, uint64_t* primes_array, size_t array_len, int* final_count){
    uint64_t* factors = malloc(20 * sizeof(uint64_t));
    if(!factors){
        return NULL;
    }
    int k = 0;
    int count = 0;
    for(int i = 0; i < array_len; i++){
        if((number % primes_array[i]) == 0){
            factors[k] = primes_array[i];
            k++;
            count++;
            if(primes_array[i]*primes_array[i] > number) break;
        }
    }
    if(count == 0){
        free(factors);
        return NULL;
    }
    uint64_t* resize = realloc(factors, count*sizeof(uint64_t));
    if(!resize){
        return factors;
    }
    *final_count = count;
    factors = resize;
    return factors;
}

uint64_t rho_f(uint64_t number, uint64_t mod){
    return (ring_mul(number, number, mod) + 1) % mod;
}

uint64_t rho_g(uint64_t number, uint64_t mod){
    return rho_f(rho_f(number, mod), mod);
}

//Pollard rho-method single factor search
uint64_t rho_factor(uint64_t number, bool method){
    if(number % 2 == 0) return (uint64_t) 2;
    if(number <= 1) return 0;
    //calculate initial values (x, y, d)
    for(int shot = 0; shot < 5; shot++){
        uint64_t x_init = (method) ? rnd_between(2, number - 1) : 2;
        uint64_t y_init = x_init;
        uint64_t d = 1;
        uint64_t x_next, y_next;
        int cycler = 0;
        while(d == 1 && cycler < 10000000){
            cycler++;
            x_next = rho_f(x_init, number);
            y_next = rho_g(y_init, number);
            uint64_t check_val = (x_next > y_next) ? (x_next - y_next) : (y_next - x_next);
            d = EEA(check_val, number, NULL, NULL);
            x_init = x_next;
            y_init = y_next;
            if (d == number) break;
            if (d > 1) return d;
        }
    }
    return 0;
}

///Brillhart-Morrison method (CFRAC) and its friends

//factor-base generation
uint64_t* bm_factorbase(uint64_t number, uint64_t* count){
    long double L = expl(pow(log(number)*log(log(number)), 0.5));
    long double a = 1/sqrtl(2);
    long double L_a = powl(L, a);
    uint64_t prime_gen_lim = (uint64_t) ((L_a / (long double)log(L_a)) * 1.2);
    uint64_t* starting_base = gen_primes(prime_gen_lim);
    uint64_t* bm_base = malloc(prime_gen_lim * sizeof(uint64_t));
    if(!bm_base){
        //allocation error
        return NULL;
    }
    bm_base[0] = 2;
    uint64_t bm_count = 1;
    for(int i = 0; i < prime_gen_lim; i++){
        if(JCB(number, starting_base[i]) == 1 && starting_base[i] <= (uint64_t)L_a){
            bm_base[bm_count] = starting_base[i];
            bm_count++;
        }
    }
    free(starting_base);
    uint64_t* resized = realloc(bm_base, bm_count * sizeof(uint64_t));
    if(!resized){
        return bm_base;
    }
    bm_base = resized;
    if(count) *count = bm_count;
    return bm_base;
}

//BM B-smoothness
bool bm_smoothness(uint64_t x, const uint64_t* factorbase, const uint64_t fb_size){
    for(uint64_t i = 0; i < fb_size; i++){
        uint64_t p = factorbase[i];
        while(x % p == 0){
            x /= p;
        }
        if(x == 1) return true;
    }
    return (x == 1) ? true : false;
}

//BM method chain fractions. FUCK. THIS. SHIT.
uint64_t* bm_cfrac(uint64_t number, uint64_t* b_size){
    uint64_t fb_size;
    uint64_t* factorbase = bm_factorbase(number, &fb_size);
    uint64_t a_0 = (uint64_t)floorl(sqrtl(number));
    uint64_t u_init = a_0;
    uint64_t v_init = 1;
    uint64_t* b = malloc((fb_size + 1) * sizeof(uint64_t));
    if(!b) return NULL;
    uint64_t a_init = a_0;
    b[0] = 0;
    b[1] = 1;
    __int128_t v_next, u_next;
    uint64_t a_next;
    uint64_t found = 2;
    while(found < fb_size + 1){
        v_next = ((__int128_t)number - (__int128_t)u_init * (__int128_t)u_init) / v_init;
        a_next = (uint64_t)((uint64_t)sqrtl(number) + u_init) / v_next;
        u_next = a_next * v_next - u_init;
        u_init = u_next;
        v_init = v_next;
        a_init = a_next;
        uint64_t b_i = a_next * b[found - 1] + b[found - 2];
        uint64_t curr = ring_pow(b_i, 2, number);
        if(bm_smoothness(curr, factorbase, fb_size)){
            found++;
            b[found] = curr;
        }
    }
    *b_size = found;
    return b;
}