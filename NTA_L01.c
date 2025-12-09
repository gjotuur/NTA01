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
typedef struct {
    uint64_t prime;
    uint32_t power;
} PrimeFactor;

//Struct to save number in its canonical expansion
typedef struct {
    PrimeFactor* factors;
    uint32_t count;
    uint32_t capacity;
} Factorization;

//Struct to save Bismuth
typedef struct {
    uint64_t P_i;       
    uint64_t Q_i;
    Factorization fact;
} BSmooth;

//Struct for Gaussian elimination solutions
typedef struct {
    uint8_t* solution;  //usefull bits from equation
    uint64_t size;
} LinearDependency;

uint64_t rnd_between(uint64_t x, uint64_t y);                           //simple RNG, to do it without srand() and rand
uint64_t EEA(uint64_t a, uint64_t n, uint64_t* u, uint64_t* v);         //Extended Euclidean alg 
bool SS_test(uint64_t p, int k, uint64_t* counter);                     //Solovey-Strassen test
void test_number(uint64_t number, int precision);                       //Test with consol message
int8_t JCB(uint64_t a, uint64_t n);                                     //Jacobi symbol calculation

void factorization_init(Factorization* f);                                  //Init factorized numbah struct
void factorization_free(Factorization* f);                                  //Free factorized numbah struct
void factorization_add(Factorization* f, uint64_t prime, uint32_t power);   //Add new factor
bool factorize_with_base(uint64_t x, const uint64_t* factorbase, uint64_t fb_size, Factorization* result); //For cfrac
void print_factorization(const Factorization* f);                           //Let it be


static inline uint64_t ring_mul(uint64_t n1, uint64_t n2, uint64_t mod);            //actually, it`s mulmod, but i don`t like the name
static inline uint64_t ring_pow(uint64_t base, uint64_t exponent, uint64_t mod);    //actually, it`s powmod, but i don`t like the name as well :D

//trial divisions
uint64_t* gen_primes(uint64_t n);                                                                      //prime numbers generation
uint64_t* td_method(uint64_t number, uint64_t* primes_array, size_t array_len, int* final_count);      //trial divisions method

//Rho-method
uint64_t rho_f(uint64_t number, uint64_t mod);                                                         //Internal function to find x_i, could be undeclared as prototype
uint64_t rho_g(uint64_t number, uint64_t mod);                                                         //internal function to find y_i, could be undeclared as prototype
uint64_t rho_factor(uint64_t number, bool method);                                                     //Pollard rho-method (1 factor)

//Brillhart - Morrison method
uint64_t* bm_factorbase(uint64_t number, uint64_t* count);                                                             //factor-base generation
BSmooth* bm_cfrac(uint64_t N, uint64_t* smooth_count, const uint64_t* factorbase, uint64_t fb_size);                   //Chain fractions B-smooth array generation
uint8_t** bm_matrix(BSmooth* smooths, uint64_t smooths_count, const uint64_t* factorbase, uint64_t fb_size);           //Matrix generation
int bm_solvesys(uint8_t** M, int rows, int cols, uint8_t* solution);                                                    //Solve system with known matrix

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
    uint64_t test4_num =  1495056764861639599;
    //uint64_t factor = rho_factor(test4_num, true);
    //if(factor) printf("\nFound factor: %llu, test: %llu", factor, test4_num % factor);

    uint64_t res = test4_num;
    for(int i = 0; i < 3; i++){
        uint64_t factor = rho_factor(res, true);
        printf("\n%4llu", factor);
        res /= factor;
    }
    printf("\n%llu", res);

    /*printf("\nBrillhart-Morrison method: factor-base\n");
    uint64_t fb_count;
    uint64_t* fb = bm_factorbase(test4_num, &fb_count);
    for(uint64_t i = 0; i < fb_count; i++){
        printf("%10llu  ", fb[i]);
        (i == 9 || (i+1)%10 == 0) ? printf("\n") : 0;
    }
    printf("\nFactor base size is %llu", fb_count);

    uint64_t smooth_c;
    BSmooth* smooth = bm_cfrac(test4_num, &smooth_c, fb, fb_count);

    if(smooth && smooth_c > 0){
        printf("Found some");
        for(uint64_t i = 0; i < smooth_c; i++){
            printf("%2llu. Q_%llu = %llu =", i+1, i, smooth[i].Q_i);
            print_factorization(&smooth[i].fact);
            printf("\n");
        }
    }
    printf("\nfb_size = %llu,\nb_smooths = %llu", fb_count, smooth_c);
    
    return 0;*/
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

bool factorize_with_base(uint64_t x, const uint64_t* factorbase, uint64_t fb_size, Factorization* result) {
    factorization_init(result);
    uint64_t remaining = x;
    
    for(uint64_t i = 0; i < fb_size; i++) {
        uint64_t p = factorbase[i];
        uint32_t power = 0;
        
        while(remaining % p == 0) {
            remaining /= p;
            power++;
        }
        
        if(power > 0) {
            factorization_add(result, p, power);
        }
        
        if(remaining == 1) return true;
    }
    
    // Якщо remaining != 1, число не є B-smooth
    if(remaining != 1) {
        factorization_free(result);
        return false;
    }
    
    return true;
}

//factor-base generation
uint64_t* bm_factorbase(uint64_t number, uint64_t* count) {
    long double L = expl(sqrtl(logl(number) * logl(logl(number))));
    long double a = 1.0L / sqrtl(2.0L);
    long double L_a = powl(L, a);
    uint64_t prime_gen_lim = (uint64_t)((L_a / logl(L_a)) * 1.2);
    
    if(prime_gen_lim < 100) prime_gen_lim = 100;
    if(prime_gen_lim > 10000) prime_gen_lim = 10000;
    
    uint64_t* starting_base = gen_primes(prime_gen_lim);
    if(!starting_base) return NULL;
    
    uint64_t* bm_base = malloc(prime_gen_lim * sizeof(uint64_t));
    if(!bm_base) {
        free(starting_base);
        return NULL;
    }
    
    bm_base[0] = UINT64_MAX;
    uint64_t bm_count = 1;
    
    bm_base[bm_count++] = 2;
    
    for(uint64_t i = 1; i < prime_gen_lim; i++) {
        uint64_t p = starting_base[i];
        if(p > (uint64_t)L_a) break;
        
        if(JCB(number, p) == 1) {
            bm_base[bm_count++] = p;
        }
    }
    
    free(starting_base);
    
    uint64_t* resized = realloc(bm_base, bm_count * sizeof(uint64_t));
    if(resized) bm_base = resized;
    
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
BSmooth* bm_cfrac(uint64_t N, uint64_t* smooth_count, const uint64_t* factorbase, uint64_t fb_size) {
    printf("\n=== CFRAC Algorithm ===\n");
    printf("N = %llu\n", N);
    printf("Factor base size: %llu\n", fb_size);
    
    // Ініціалізація
    uint64_t sqrt_N = (uint64_t)sqrtl((long double)N);
    
    // P_0 = sqrt(N), Q_0 = 1, Q_1 = N - P_0²
    uint64_t P_prev = sqrt_N;
    uint64_t Q_prev = 1;
    uint64_t Q_curr = N - P_prev * P_prev;
    
    // Конвергенти: A_{-1} = 1, A_0 = P_0
    uint64_t A_prev2 = 1;
    uint64_t A_prev1 = sqrt_N;
    uint64_t A_curr;
    
    // Масив для B-smooth чисел
    uint64_t capacity = fb_size + 10;
    BSmooth* smooth_numbers = malloc(capacity * sizeof(BSmooth));
    if(!smooth_numbers) return NULL;
    
    uint64_t found = 0;
    uint64_t iterations = 0;
    uint64_t max_iterations = 100000; // Захист від зациклення
    
    printf("\nSearching for B-smooth numbers...\n");
    
    while(found < fb_size + 1 && iterations < max_iterations) {
        iterations++;
        
        uint64_t a_i = (sqrt_N + P_prev) / Q_curr;                  // a_i = floor((sqrt(N) + P_i) / Q_i)
        uint64_t P_curr = a_i * Q_curr - P_prev;                    // P_{i+1} = a_i * Q_i - P_i
        uint64_t Q_next = Q_prev + a_i * (P_prev - P_curr);         // Q_{i+1} = Q_{i-1} + a_i * (P_i - P_{i+1})

        A_curr = ring_mul(a_i, A_prev1, N) + A_prev2;               // A_i = a_i * A_{i-1} + A_{i-2}
        if(A_curr >= N) A_curr %= N;
        
        Factorization fact;
        bool is_smooth = factorize_with_base(Q_curr, factorbase, fb_size, &fact);
        
        if(is_smooth) {
            printf("  [%llu] Found: P=%llu, Q=%llu, A=%llu\n", found + 1, P_prev, Q_curr, A_prev1);
            
            smooth_numbers[found].P_i = P_prev;
            smooth_numbers[found].Q_i = Q_curr;
            smooth_numbers[found].fact = fact;
            found++;
            
            if(found >= capacity) {
                capacity *= 2;
                BSmooth* tmp = realloc(smooth_numbers, capacity * sizeof(BSmooth));
                if(!tmp) break;
                smooth_numbers = tmp;
            }
        }
        
        // Оновлюємо для наступної ітерації
        Q_prev = Q_curr;
        Q_curr = Q_next;
        P_prev = P_curr;
        A_prev2 = A_prev1;
        A_prev1 = A_curr;
        
        // Показуємо прогрес
        if(iterations % 1000 == 0) {
            printf("  Iterations: %llu, found: %llu/%llu\n", iterations, found, fb_size + 1);
        }
    }
    
    printf("\nCFRAC completed: %llu B-smooth numbers found in %llu iterations\n", found, iterations);
    
    if(smooth_count) *smooth_count = found;
    return smooth_numbers;
}

//Initialize struct for factorized numbah
void factorization_init(Factorization* f) {
    f->capacity = 20;
    f->count = 0;
    f->factors = malloc(f->capacity * sizeof(PrimeFactor));
}

//Free struct for factorized numbah
void factorization_free(Factorization* f) {
    free(f->factors);
    f->factors = NULL;
    f->count = 0;
}

//Add new factor to struct
void factorization_add(Factorization* f, uint64_t prime, uint32_t power) {
    if (f->count >= f->capacity) {
        f->capacity *= 2;
        f->factors = realloc(f->factors, f->capacity * sizeof(PrimeFactor));
    }
    f->factors[f->count].prime = prime;
    f->factors[f->count].power = power;
    f->count++;
}

void print_factorization(const Factorization* f) {
    for(uint32_t i = 0; i < f->count; i++) {
        if(i > 0) printf(" x ");
        printf("%llu", f->factors[i].prime);
        if(f->factors[i].power > 1) {
            printf("^%u", f->factors[i].power);
        }
    }
}

//build a matrix to solve system
uint8_t** bm_matrix(BSmooth* smooths, uint64_t smooths_count, const uint64_t* factorbase, uint64_t fb_size){
    uint8_t** M = malloc(smooths_count*sizeof(uint8_t*));
    if(!M){
        return NULL;
    }
    for(uint64_t i = 0; i < smooths_count; i++){
        M[i] = calloc(fb_size, sizeof(uint8_t));
        if(!M[i]) return NULL;

        Factorization* f = &smooths[i].fact;
        for(uint32_t k = 0; k < f->count; k++){
            uint64_t p = f->factors[k].prime;
            uint32_t pow = f->factors[k].power;

            for(uint64_t j = 0; j < fb_size; j++){
                if(factorbase[j] == p){
                    M[i][j] = (pow & 1);
                    break;
                }
            }
        }
    }
    return M;
}

int bm_solvesys(uint8_t** M, int rows, int cols, uint8_t* solution) {
    int r = 0;
    for (int c = 0; c < cols && r < rows; c++) {
        int pivot = -1;
        for (int i = r; i < rows; i++) {
            if (M[i][c]) {
                pivot = i;
                break;
            }
        }
        if (pivot == -1)
            continue;
        // swap
        if (pivot != r) {
            uint8_t* tmp = M[pivot];
            M[pivot] = M[r];
            M[r] = tmp;
        }
        // eliminate
        for (int i = 0; i < rows; i++) {
            if (i != r && M[i][c]) {
                for (int j = c; j < cols; j++) {
                    M[i][j] ^= M[r][j];
                }
            }
        }
        r++;
    }
    for (int i = 0; i < cols; i++) solution[i] = 0;

    for (int i = rows - 1; i >= 0; i--) {
        int leading = -1;
        for (int j = 0; j < cols; j++) {
            if (M[i][j]) {
                leading = j;
                break;
            }
        }
        if (leading == -1) continue;

        uint8_t val = 0;
        for (int j = leading + 1; j < cols; j++)
            val ^= (M[i][j] & solution[j]);

        solution[leading] = val;
    }

    return 1;
}