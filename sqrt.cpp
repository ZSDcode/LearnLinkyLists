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

float fast_exp_iter(float b, int count, float product) {
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

float fast_exp(float b, int n) {
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

bool is_prime(float n) {
    if (n <= 1) {
        return false;
    }
    if (check_divisibility(n, 2, 0) == 0) {
        return true;
    }
    else {
        return false;
    }
}

bool fermat_lil(int a, int n) {
    int modpow = static_cast<int>(fast_exp(a, n)) % n;
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
    return x * x * x + 2 * x;
}

float simpsons_est(function<float(float)> func, float initial, float final, float dx, float est) {
    if (initial + 2 * dx - 0.0000001 > final) {
        return est;
    }
    else {
        float term_sum_main = func(initial) + 4 * func(change(initial, dx)) + func(change(initial, 2 * dx));
        float term_sum_acc = term_sum_main * dx / 3;
        return simpsons_est(func, change(initial, 2 * dx), final, dx, est + term_sum_acc);
    }
}

float product_iter(function<float(float)> term, float initial, function<float(float)> next, float final, float result) {
    if (initial > final) {
        return result;
    }
    else {
        return product_iter(term, next(initial), next, final, result * term(initial));
    }
}

float product_recurse(function<float(float)> term, float initial, function<float(float)> next, float final) {
    if (initial >= final) {
        return term(initial);
    }
    else {
        return term(initial) * product_recurse(term, next(initial), next, final);
    }
}

float equiv(int n) {
    return (float)n;
}

float fact_abstraction(int n) {
    return product_iter(equiv, 1, inc, n, 1);
}

int inc_double(int x) {
    return x + 2;
}

float multiply_consec_eo_no(float n) {
    return product_iter(equiv, n, inc_double, n+2, 1);
}

float pi_approx(int n) {
    float top;
    float bottom;
    if (n % 2 == 0) {
        top = product_iter(multiply_consec_eo_no, 2, inc_double, (float)n, 1);
        bottom = 3 * (n+1) * product_iter(multiply_consec_eo_no, 3, inc_double, (float)n, 1);
    } else {
        top = product_iter(multiply_consec_eo_no, 2, inc_double, (float)n, 1) * (n+1);
        bottom = 3 * product_iter(multiply_consec_eo_no, 3, inc_double, (float)n, 1);
    }
    return top / bottom * 4;
}

float accumulate(function<float(float, float)> combiner, float null_value, function<float(float)> term, float a, function<float(float)> next, float b) {
    if (a > b) {
        return null_value;
    }
    else {
        return accumulate(combiner, combiner(null_value, term(a)), term, next(a), next, b);
    }
}

float accumulate_recursive(function<float(float, float)> combiner, float null_value, function<float(float)> term, float a, function<float(float)> next, float b) {
    if (a > b) {
        return null_value;
    }
    else {
        return combiner(term(a), accumulate_recursive(combiner, null_value, term, next(a), next, b));
    }
}

float filtered_accumulate(function<bool(float)> filter, function<float(float, float)> combiner, float null_value, function<float(float)> term, float a, function<float(float)> next, float b) {
    if (a > b) {
        return null_value;
    }
    else {
        if (filter(a)) {
            return filtered_accumulate(filter, combiner, combiner(null_value, term(a)), term, next(a), next, b);
        }
        else {
            return filtered_accumulate(filter, combiner, null_value, term, next(a), next, b);
        }
    }
};

float test_combiner(float a, float b) {
    return a + b;
}

float filtered_cubes(float x, float y) {
    return filtered_accumulate(
        [](float x) {return static_cast<int>(x) % 2 == 0;},
        [](float x, float y) {return x + y;}, 
        0,
        [](float x) {return x * x * x;},
        x,
        [](float x) {return x + 1;},
        y
    );
}

float search_for_root(function<float(float)> func, float neg_point, float pos_point) {
    auto good_enuf = [](float x, float y) {
        return abs(x - y) < 0.0001;
    };
    auto average = [](float x, float y) {return (x + y)/2;};
    if (good_enuf(neg_point, pos_point)) {
        return average(neg_point, pos_point);
    }
    else {
        if (func(neg_point) > 0 and func(pos_point) < 0) {
            return search_for_root(func, pos_point, neg_point);
        }
        else if (func(neg_point) < 0 and func(pos_point) > 0) {
            float avg = average(neg_point,pos_point);
            if (func(avg) > 0) {
                return search_for_root(func, neg_point, avg);
            }
            else if (func(avg) < 0) {
                return search_for_root(func, avg, pos_point);
            }
            else {
                return avg;
            }
        }
        else {
            return 92387091283128903;
        }
    }
}

float fixed_point(function<float(float)> func, float guess) {
    if (abs(guess - func(guess)) < 0.001) {
        return func(guess);
    }
    else if (abs(guess - func(guess)) > 1000 ) {
        return func(guess);
    }
    else if ( isnan(guess) || guess == INFINITY) {
        return guess;
    }
    else {
        cout << guess << endl;
        return fixed_point(func, func(guess));
    }
}

float fixed_damping(function<float(float)> func, float guess) {
    auto avg = [](float x, float y) -> float {
        return (x + y) / 2;
    };
    if (abs(guess - func(guess)) < 0.001) {
        return func(guess);
    }
    else if (abs(guess - func(guess)) > 1000 ) {
        return func(guess);
    }
    else if ( isnan(guess) || guess == INFINITY) {
        return guess;
    }
    else {
        cout << guess << endl;
        return fixed_point(func, func(0.5 * avg(guess, func(guess))));
    }
}

float log_approx(float base, float x) {
    auto log_mac_term = [](float n, float x) -> float {
        int power = static_cast<int>(n);
        if (power % 2 == 0) {
            return -1 * fast_exp(x, power) / power;
        } else {
            return fast_exp(x, power) / power;
        }
    };
    auto log_mac_use = [&log_mac_term, x](float n) -> float {
        if (x > 1) {
            return log_mac_term(n, (1 / x) - 1) * -1;
        } else{
            return log_mac_term(n, x-1);
        }
    };
    auto log_mac_base = [&log_mac_term, base](float n) -> float {
        return -1 * log_mac_term(n, (1 / base) - 1);
    };
    if (base <= 1 || x < -1) {
        return 21312323423212;
    }
    else {
        return accumulate(
            [](float a, float b){return a + b;}, 
            0,
            log_mac_use, 
            1, 
            [](float a){return a + 1;},
            15
        ) / accumulate(
            [](float a, float b){return a + b;}, 
            0,
            log_mac_base, 
            1, 
            [](float a){return a + 1;},
            15
        );
    }
}

float xx_solution(float equal) {
    return fixed_point([equal](float x){return log2(equal)/log2(x);}, 2);
}

float cont_frac_recursive(float i, function<float(float)> term_n, function<float(float)>term_d, float count) {
    if (i >= count) {
        return term_n(i)/term_d(i);
    }
    else {
        return term_n(i) / (term_d(i) + cont_frac_recursive(i+1, term_n, term_d, count));
    }
}

float cont_frac(function<float(float)> term_n, function<float(float)> term_d, float count, float result) {
    if (count <= 0) {
        return result;
    }
    else {
        return cont_frac(term_n, term_d, count-1, term_n(count)/(term_d(count) + result));
    }
}

float tan_cf(float x, float k, float result) {
    auto square = [](float x){
        return x * x;
    };
    if (k <= 0) {
        return result;
    }
    else if (k == 1) {
        return tan_cf(x, k-1, x/(k-result));
    }
    else {
        return tan_cf(x, k-1, square(x)/((2 * k - 1)-result));
    }
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
    cout << def_integral(random_math_func, 1, 0.001, 2, change, 0) << endl;
    cout << simpsons_est(random_math_func, 1, 2, 0.001, 0) << endl;
    cout << fact_abstraction(6) << endl;
    cout << pi_approx(7) << endl;
    cout << filtered_accumulate(is_prime, test_combiner, 0, sqrt_myfunc, 1, inc, 10) << endl;
    cout << filtered_cubes(1, 6) << endl;
    cout << search_for_root([](float x){return x * x * x + 2 * x - 3;}, 0, 100) << endl;
    cout << fixed_point([](float x){return x * x + 3 * x + 1;}, 0) << endl;
    cout << fixed_point([](float x) {return 1 + 1 / x;}, 1) << endl;
    cout << fixed_damping([](float x) {return 1 + 1 / x;}, 1) << endl;
    cout << log_approx(2, 4) << endl;
    cout << xx_solution(100) << endl;
    cout << cont_frac_recursive(1, 
        [](float i){
            return 1;
        }, 
        [](float i){
            if ((static_cast<int>(i) + 1) % 3 == 0) {
                return (i + 1)/3 * 2;
            } 
            else {
                return 1.0f;
            }
        }, 
        50) << endl;
    cout << cont_frac( 
        [](float i){
            return 1;
        }, 
        [](float i){
            if ((static_cast<int>(i) + 1) % 3 == 0) {
                return (i + 1)/3 * 2;
            } 
            else {
                return 1.0f;
            }
        }, 
        50, 0
    ) << endl;
    cout << tan_cf(1, 50, 1) << endl;
}