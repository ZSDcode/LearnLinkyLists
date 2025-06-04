#include <functional>
#include <random>
#include <iostream>
#include <vector>
using namespace std;
#include <bits/stdc++.h>

bool is_good_enuf(float guess, float x) {
    if (guess * guess - x < 0.001 && guess * guess - x > -0.01) {
        return true;
    }
    else {
        return false;
    }
};

float improve(float guess, float x) {
    return (guess + (x / guess))/2;
};

float sqrt_iter(float guess, float x) {
    if (is_good_enuf(guess, x)) {
        return guess;
    }
    else {
        return sqrt_iter(improve(guess, x), x);
    }
};

float sqrt_myfunc(float x) {
    return sqrt_iter(1, x);
};

int A(int x, int y) {
    if (y == 0) {
        return 0;
    }
    else if (x == 0) {
        return 2 * y;
    }
    else if (y == 1) {
        return 2;
    }
    else {
        return A(x-1, A(x, y-1));
    }
}

int fib_iter(int a, int b, int count) {
    if (count <= 3) {
        return a;
    }
    else {
        int total = a + b;
        b = a;
        a = total;
        return fib_iter(a, b, count - 1);
    }
    
}

int fib(int x) {
    return fib_iter(1, 1, x);
};

int f1_recurse(int n) {
    if (n < 3) {
        return n;
    }
    else {
        return f1_recurse(n-1) + 2 * f1_recurse(n-2) + 3 * f1_recurse(n-3);
    }
}

int f1_iter(int a, int b, int c, int count) {
    if (count < 4) {
        return 3 * a + 2 * b + c;
    }
    else {
        int new_c = 3 * a + 2 * b + c;
        return f1_iter(b, c, new_c, count-1);
    }
}

int f1_iterFull(int n) {
    if (n < 3) {
        return n;
    }
    else {
        return f1_iter(0, 1, 2, n); 
    }
}

vector<int> Pascal(int n) {
    if (n == 1) {
        return {1};
    }
    else if (n == 2) {
        return {1, 1};
    }
    else {
        vector<int> prev_row = Pascal(n-1);
        vector<int> new_row = {1};
        for (int i = prev_row.size()-1; i >= 1; --i) {
            new_row.push_back(prev_row[i]+prev_row[i-1]);
        }
        new_row.push_back(1);
        return new_row;
    }
}

int exp_recurse(int b, int n) {
    if (n == 1) {
        return b;
    }
    else {
        return b * exp_recurse(b, n-1);
    }
}

int exp_iter(int b, int count, int product) {
    if (count == 0) {
        return product;
    }
    else {
        return exp_iter(b, count-1, product * b);
    }
}

int exp(int b, int n) {
    return exp_iter(b, n, 1);
}

int fast_exp_recurse(int b, int n) {
    if (n == 1) {
        return b;
    }
    else if (n%2 == 0) {
        return fast_exp_recurse(b, n/2) * fast_exp_recurse(b, n/2);
    }
    else {
        return b * fast_exp_recurse(b, n-1);
    }
}

int fast_exp_iter(int b, int count, int product) {
    if (count == 0) {
        return product;
    }
    else if (count % 2 == 0) {
        return fast_exp_iter(b*b, count/2, product);
    }
    else {
        return fast_exp_iter(b, count - 1, product * b);
    }
}

int fast_exp(int b, int n) {
    return fast_exp_iter(b, n, 1);    
}

int double_myfunc(int a) {
    return a*2;
}

int halve(int a) {
    return a/2;
}

int multiply_recurse(int a, int b) {
    if (b == 0) {
        return 0;
    }
    else {
        if (b == 1) {
            return a;
        }
        else {
            return a + multiply_recurse(a, b-1);
        }
    }
}

int fast_multiply_recurse(int a, int b) {
    if (b == 0) {
        return 0;
    }
    else {
        if (b == 1) {
            return a;
        }
        else {
            if (b % 2 == 0) {
                return fast_multiply_recurse(double_myfunc(a), halve(b));
            }
            else {
                return a + fast_multiply_recurse(a, b-1);
            }
        }
    }
}

int fast_multiply_iter(int a, int count, int total) {
    if (count == 0) {
        return total;
    }
    else {
        if (count % 2 == 0) {
            return fast_multiply_iter(double_myfunc(a), halve(count), total);
        }
        else {
            return fast_multiply_iter(a, count-1, total + a);
        }
    }
}

int fast_multiply(int a, int b) {
    return fast_multiply_iter(a, b, 0);
}

int euclid_gcd(int a, int b) {
    if (b == 0) {
        return a;
    }
    else {
        return (euclid_gcd(b, a%b));
    }
}

bool check_divide(int a, int b) {
    if (a % b == 0) {
        return true;
    }
    else {
        return false;
    }
}

int check_divisibility(int n, int count, int no_of_divisors) {
    if (count * count > n) {
        return no_of_divisors;
    }
    else {
        if (check_divide(n, count)) {
            no_of_divisors += 1;
            return check_divisibility(n, count + 1, no_of_divisors);
        }
        else {
            return check_divisibility(n, count+1, no_of_divisors);
        }
    }
}

bool is_prime(int n) {
    if (check_divisibility(n, 2, 0) == 0) {
        return true;
    }
    else {
        return false;
    }
}

bool fermat_lil(int a, int n) {
    int modpow = fast_exp(a, n) % n;
    int mod = a % n;
    if (modpow == mod) {
        return true;
    }
    else {
        return false;
    }
}

random_device rd;
mt19937 gen(rd());

bool fermat_prime_test(int n, int repeat) {
    if (repeat <= 0) {
        return true;
    }
    else {
        uniform_int_distribution<int> distribution(1, n-1);
        int a = distribution(gen);
        if ( ! fermat_lil(a, n) ) {
            return false;
        }
        else {
            return fermat_prime_test(n, repeat - 1);
        }
    }
}

int summation (function<int(int)> term, int x, function<int(int)> next, int y) {
    if (x > y) {
        return 0;
    }
    else {
        return term(x) + summation(term, next(x), next, y);
    }
}

int cube(int x) {
    return x * x * x;
}

int inc(int x) {
    return x + 1;
}

float def_integral(function<float(float)> func, float x, float dx, float y, function<float(float, float)> change, float tot) {
    if (x > y) {
        return tot;
    }
    else {
        return def_integral(func, change(x, dx), dx, y, change, tot + func(x+dx/2) * dx);
    }
}

float change(float x, float dx) {
    return x + dx;
}

float random_math_func(float x) {
    return x * x + 2 * x;
}

int main() {
    cout << sqrt_myfunc(1000) << endl;
    cout << A(1, 10) << endl;
    cout << fib(10) << endl;
    cout << f1_recurse(4) << endl;
    cout << f1_iterFull(4) << endl;
    vector<int> row6 = Pascal(6);
    for (int i = 0; i < row6.size(); i++) {
        cout << row6[i] << ", ";
    }
    cout << endl;
    cout << exp_recurse(8, 3) << endl;
    cout << exp(8, 3) << endl;
    cout << fast_exp_recurse(2, 6) << endl;
    cout << fast_exp(2, 10) << endl;
    cout << multiply_recurse(2, 3) << endl;
    cout << fast_multiply_recurse(10, 20) << endl;
    cout << fast_multiply(2, 3) << endl;
    cout << euclid_gcd(100, 75) << endl;
    cout << is_prime(10) << endl;
    cout << is_prime(5) << endl;
    cout << fermat_lil(2, 11) << endl;
    cout << summation(cube, 1, inc, 3) << endl;
    cout << def_integral(random_math_func, 1, 0.01, 2, change, 0) << endl;
}