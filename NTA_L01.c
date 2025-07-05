#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int EEA(int a, int n, int* u, int* v);
int reduce(int a, int n);

int main(){
    //Tests for EEA
    int a = 0, n = 0;
    int x, y;
    printf("Enter the numbers to analyze");
    printf("\na: ");
    scanf("%d", &a);
    printf("\nb: ");
    scanf("%d", &n);
    int g = EEA(a, n, &x, &y);
    printf("\nGCD of %d and %d is %d", a, n, g);
    printf("\nBezout`s equation: %d + %d = %d", a*x, n*y, g);


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

int reduce(int a, int n){
    if(a < n){
        return a;
    }
    while(a > n){
        a -= n;
    }
    return a;
}