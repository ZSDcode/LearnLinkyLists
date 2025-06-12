#include <functional>
#include <random>
#include <iostream>
#include <vector>
using namespace std;
#include <bits/stdc++.h>

bool is_good_enuf(double guess, double x) {
    if (guess * guess - x < 0.001 && guess * guess - x > -0.01) {
        return true;
    }
    else {
        return false;
    }
};

double improve(double guess, double x) {
    return (guess + (x / guess))/2;
};

double sqrt_iter(double guess, double x) {
    if (is_good_enuf(guess, x)) {
        return guess;
    }
    else {
        return sqrt_iter(improve(guess, x), x);
    }
};

double sqrt_myfunc(double x) {
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

double fast_exp_iter(double b, int count, double product) {
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

double fast_exp(double b, int n) {
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

bool is_prime(double n) {
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

double def_integral(function<double(double)> func, double x, double dx, double y, function<double(double, double)> change, double tot) {
    if (x > y) {
        return tot;
    }
    else {
        return def_integral(func, change(x, dx), dx, y, change, tot + func(x+dx/2) * dx);
    }
}

double change(double x, double dx) {
    return x + dx;
}

double random_math_func(double x) {
    return x * x * x + 2 * x;
}

double simpsons_est(function<double(double)> func, double initial, double final, double dx, double est) {
    if (initial + 2 * dx - 0.0000001 > final) {
        return est;
    }
    else {
        double term_sum_main = func(initial) + 4 * func(change(initial, dx)) + func(change(initial, 2 * dx));
        double term_sum_acc = term_sum_main * dx / 3;
        return simpsons_est(func, change(initial, 2 * dx), final, dx, est + term_sum_acc);
    }
}

double product_iter(function<double(double)> term, double initial, function<double(double)> next, double final, double result) {
    if (initial > final) {
        return result;
    }
    else {
        return product_iter(term, next(initial), next, final, result * term(initial));
    }
}

double product_recurse(function<double(double)> term, double initial, function<double(double)> next, double final) {
    if (initial >= final) {
        return term(initial);
    }
    else {
        return term(initial) * product_recurse(term, next(initial), next, final);
    }
}

double equiv(int n) {
    return (double)n;
}

double fact_abstraction(int n) {
    return product_iter(equiv, 1, inc, n, 1);
}

int inc_double(int x) {
    return x + 2;
}

double multiply_consec_eo_no(double n) {
    return product_iter(equiv, n, inc_double, n+2, 1);
}

double pi_approx(int n) {
    double top;
    double bottom;
    if (n % 2 == 0) {
        top = product_iter(multiply_consec_eo_no, 2, inc_double, (double)n, 1);
        bottom = 3 * (n+1) * product_iter(multiply_consec_eo_no, 3, inc_double, (double)n, 1);
    } else {
        top = product_iter(multiply_consec_eo_no, 2, inc_double, (double)n, 1) * (n+1);
        bottom = 3 * product_iter(multiply_consec_eo_no, 3, inc_double, (double)n, 1);
    }
    return top / bottom * 4;
}

double accumulate(function<double(double, double)> combiner, double null_value, function<double(double)> term, double a, function<double(double)> next, double b) {
    if (a > b) {
        return null_value;
    }
    else {
        return accumulate(combiner, combiner(null_value, term(a)), term, next(a), next, b);
    }
}

double accumulate_recursive(function<double(double, double)> combiner, double null_value, function<double(double)> term, double a, function<double(double)> next, double b) {
    if (a > b) {
        return null_value;
    }
    else {
        return combiner(term(a), accumulate_recursive(combiner, null_value, term, next(a), next, b));
    }
}

double filtered_accumulate(function<bool(double)> filter, function<double(double, double)> combiner, double null_value, function<double(double)> term, double a, function<double(double)> next, double b) {
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

double test_combiner(double a, double b) {
    return a + b;
}

double filtered_cubes(double x, double y) {
    return filtered_accumulate(
        [](double x) {return static_cast<int>(x) % 2 == 0;},
        [](double x, double y) {return x + y;}, 
        0,
        [](double x) {return x * x * x;},
        x,
        [](double x) {return x + 1;},
        y
    );
}

double search_for_root(function<double(double)> func, double neg_point, double pos_point) {
    auto good_enuf = [](double x, double y) {
        return abs(x - y) < 0.0001;
    };
    auto average = [](double x, double y) {return (x + y)/2;};
    if (good_enuf(neg_point, pos_point)) {
        return average(neg_point, pos_point);
    }
    else {
        if (func(neg_point) > 0 and func(pos_point) < 0) {
            return search_for_root(func, pos_point, neg_point);
        }
        else if (func(neg_point) < 0 and func(pos_point) > 0) {
            double avg = average(neg_point,pos_point);
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

double fixed_point(function<double(double)> func, double guess) {
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
        return fixed_point(func, func(guess));
    }
}

double fixed_damping(function<double(double)> func, double guess) {
    auto avg = [](double x, double y) -> double {
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
        return fixed_point(func, func(0.5 * avg(guess, func(guess))));
    }
}

double log_approx(double base, double x) {
    auto log_mac_term = [](double n, double x) -> double {
        int power = static_cast<int>(n);
        if (power % 2 == 0) {
            return -1 * fast_exp(x, power) / power;
        } else {
            return fast_exp(x, power) / power;
        }
    };
    auto log_mac_use = [&log_mac_term, x](double n) -> double {
        if (x > 1) {
            return log_mac_term(n, (1 / x) - 1) * -1;
        } else{
            return log_mac_term(n, x-1);
        }
    };
    auto log_mac_base = [&log_mac_term, base](double n) -> double {
        return -1 * log_mac_term(n, (1 / base) - 1);
    };
    if (base <= 1 || x < -1) {
        return 21312323423212;
    }
    else {
        return accumulate(
            [](double a, double b){return a + b;}, 
            0,
            log_mac_use, 
            1, 
            [](double a){return a + 1;},
            15
        ) / accumulate(
            [](double a, double b){return a + b;}, 
            0,
            log_mac_base, 
            1, 
            [](double a){return a + 1;},
            15
        );
    }
}

double xx_solution(double equal) {
    return fixed_point([equal](double x){return log2(equal)/log2(x);}, 2);
}

double cont_frac_recursive(double i, function<double(double)> term_n, function<double(double)>term_d, double count) {
    if (i >= count) {
        return term_n(i)/term_d(i);
    }
    else {
        return term_n(i) / (term_d(i) + cont_frac_recursive(i+1, term_n, term_d, count));
    }
}

double cont_frac(function<double(double)> term_n, function<double(double)> term_d, double count, double result) {
    if (count <= 0) {
        return result;
    }
    else {
        return cont_frac(term_n, term_d, count-1, term_n(count)/(term_d(count) + result));
    }
}

double tan_cf(double x, double k, double result) {
    auto square = [](double x){
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

double deriv(function<double(double)> func, double x, double dx) {
    return (func(x + dx) - func(x))/dx;
}

double newton_method(function<double(double)> func, double guess) {
    if (abs(func(guess)) < 0.001) {
        return guess;
    }
    else {
        return newton_method(func, guess - func(guess)/deriv(func, guess, 0.0001));
    }
}

double cubic_solve(double a, double b, double c) {
    return newton_method([a, b, c](double x){return fast_exp(x, 3) + a * fast_exp(x, 2) + b * x + c;}, 1);
}

double composite2(function<double(double)> f, function<double(double)> g, double x) {
    return f(g(x));
}

double repeatedapp(function<double(double)> func, double x, int count) {
    if (count == 1) {
        return func(x);
    }
    else {
        return repeatedapp(func, func(x), count-1);
    }
}


double nfoldsmoothing(function<double(double)> func, double dx, double x, int count) {
    auto smooth = [func, dx](double x) -> double {
        double tot = func(x-dx) + func(x) + func(x+dx);
        return tot / 3;
    };
    return repeatedapp(smooth, x, count);
}

double damping_func_gen(double power, double x, double y) {
    if (power <= 2) {
        return 0.5 * (x + y);
    }
    else {

        return damping_func_gen(power - 1, x, 0.5 * (x + y));
    }
}

double solving_xnequals(double power, double equals) {
    auto damping_func = [power, equals](double x) -> double {
        return damping_func_gen(power, x, equals/fast_exp(x, power-1));
    };
    return fixed_point(damping_func, 1);
}

function<double(function<double(double)>, double)> iterative_improve(function<bool(function<double(double)>, double, double)> good_enough, function<double(double)> improvement) {
    function<double(function<double(double)>, double)> returning_func;
    returning_func = [&, good_enough, improvement](function<double(double)> solving, double x) -> double {
        if (good_enough(solving, x, 0.0001)) {
            return x;
        }
        else {
            return returning_func(solving, improvement(x));
        }
    };
    return returning_func;
}

double generic_equation_solver(function<double(double)> LHS, double RHS) {
    auto good_enough = [RHS](function<double(double)> solver, double x, double tolerance) -> bool {
        if (abs(solver(x) - RHS) < tolerance) {
            return true;
        }
        else {return false;}
    };
    auto newton_improve_func = [LHS, RHS](double x) -> double {
        return x - (LHS(x) - RHS)/deriv([LHS, RHS](double x){return LHS(x)-RHS;}, x, 0.0001);
    };
    return iterative_improve(good_enough, newton_improve_func)(LHS, 1);
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
    cout << search_for_root([](double x){return x * x * x + 2 * x - 3;}, 0, 100) << endl;
    cout << fixed_point([](double x){return x * x + 3 * x + 1;}, 0) << endl;
    cout << fixed_point([](double x) {return 1 + 1 / x;}, 1) << endl;
    cout << fixed_damping([](double x) {return 1 + 1 / x;}, 1) << endl;
    cout << log_approx(2, 4) << endl;
    cout << xx_solution(100) << endl;
    cout << cont_frac_recursive(1, 
        [](double i){
            return 1;
        }, 
        [](double i){
            if ((static_cast<int>(i) + 1) % 3 == 0) {
                return (i + 1)/3 * 2;
            } 
            else {
                return (double)1;
            }
        }, 
        50) << endl;
    cout << cont_frac( 
        [](double i){
            return 1;
        }, 
        [](double i){
            if ((static_cast<int>(i) + 1) % 3 == 0) {
                return (i + 1)/3 * 2;
            } 
            else {
                return (double)1;
            }
        }, 
        50, 0
    ) << endl;
    cout << tan_cf(1, 50, 1) << endl;
    cout << deriv([](double x){return x * x * x;}, 5, 0.00001) << endl;
    cout << newton_method([](double x){return x * x - 10;}, 1) << endl;
    cout << cubic_solve(0, 0, 1) << endl;
    cout << repeatedapp([](double x){return x + 1;}, 3, 2) << endl;
    cout << nfoldsmoothing([](double x){return x * x;}, 0.0001, 2, 3) << endl;
    cout << solving_xnequals(5, 32) << endl;
    cout << generic_equation_solver([](double x){return sin(x);}, 0.5);
}