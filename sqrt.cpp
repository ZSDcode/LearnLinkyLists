#include <limits>
#include <stdexcept>
#include <format>
#include <utility>
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

class Rational {
private:
    int numer;
    int denom;

    void simplify() {
        if (denom == 0) {
            cout << "Division by Zero Error";
            numer = 0;
            denom = 0;
        }
        else {
            int g = euclid_gcd(numer, denom);
            numer = numer / g;
            denom = denom / g;
        }
        if (denom < 0) {
            denom = -1 * denom;
            numer = -1 * numer;
        }
    }

public:
    Rational(int n, int d) : numer(n), denom(d) {
        simplify();
    }
    int get_numer() const {
        return numer;
    }
    int get_denom() const {
        return denom;
    }
    void print_rational() const {
        if (denom == 1) {
            cout << numer << endl;
        }
        else {
            cout << to_string(numer) + " / " + to_string(denom) << endl;
        }
    }

    Rational& operator+=(const Rational& other) {
        numer = numer * other.denom + other.numer * denom;
        denom = denom * other.denom;
        simplify();
        return *this;
    }

    Rational& operator*=(const Rational& other) {
        numer = numer * other.numer;
        denom = denom * other.denom;
        simplify();
        return *this;
    }

    Rational& operator-=(const Rational& other) {
        Rational neg(-1 * other.numer, other.denom);
        return *this += neg;
    }

    Rational& operator/=(const Rational& other) {
        Rational flip(other.denom, other.numer);
        return *this*=flip;
    }
};

Rational operator+(Rational r1, const Rational& r2) {
    r1 += r2;
    return r1;
}

Rational operator-(Rational r1, const Rational& r2) {
    r1 -= r2;
    return r1;
}

Rational operator*(Rational r1, const Rational& r2) {
    r1 *= r2;
    return r1;
}

Rational operator/(Rational r1, const Rational& r2) {
    r1 /= r2;
    return r1;
}

bool operator==(const Rational& r1, const Rational& r2) {
    return (r1.get_numer() == r2.get_numer() && r1.get_denom() == r2.get_denom());
}

bool operator!=(const Rational& r1, const Rational& r2) {
    return !(r1 == r2);
}

template<typename T>
Rational to_rationalise(T n) {
    function<Rational(T, int)> process;
    process = [&](T n, int denom) {
        if (abs(static_cast<int>(n * denom) - n * denom) < 0.000001) {
            return Rational(static_cast<int>(n*denom), denom);
        }
        else {
            return process(n, denom * 10);
        }
    };
    if constexpr (is_integral_v<T>) { 
        return Rational(n, 1);
    }
    else {
        return process(n, 1);
    }
}

class Point {
private:
    double x_coord;
    double y_coord;
public:
    Point(double x, double y) : x_coord(x), y_coord(y) {}
    string print_point() const {
        return ("(" + to_string(x_coord) + ", " + to_string(y_coord) + ")");
    }
    double get_x() const {
        return x_coord;
    }
    double get_y() const {
        return y_coord;
    }
};

class LineSegment {
private:
    Point start_pt;
    Point end_pt;
public:
    LineSegment (Point p1, Point p2) : start_pt(p1), end_pt(p2) {}
    string print_seg() const {
        return "Start pt: " + start_pt.print_point() + ", End pt: " + end_pt.print_point();
    }
    Point midpoint() const {
        return Point((start_pt.get_x() + end_pt.get_x())/2, (start_pt.get_y() + end_pt.get_y())/2);
    }
    double length() const {
        auto square = [](double x) -> double {
            return x * x;
        };
        double length = sqrt_myfunc(square(start_pt.get_x() - end_pt.get_x()) + square(start_pt.get_y() - end_pt.get_y()));
        return length;
    }
};

class Rectangle1 {
private:
    Point bLeft_corner;
    double length;
    double width;
    double rotation_angle;
    Point tLeft_corner;
    Point tRight_corner;
    Point bRight_corner;
    Rectangle1& normalize_angle() {
        if (rotation_angle > M_PI) {
            int k = ceil((rotation_angle - M_PI)/(2 * M_PI));
            rotation_angle -= 2 * k * M_PI;
            return *this;
        }
        else if (rotation_angle <= -1 * M_PI) {
            int k = floor((rotation_angle + M_PI)/(2 * M_PI));
            rotation_angle -= 2 * k * M_PI;
            return *this;
        }
        else {
            return *this;
        }
    }
    Rectangle1& calc_update_vertices() {
        auto pivot = [](Point pivot_pt, double l, double angle) -> Point {
            return Point(pivot_pt.get_x() + (l * cos(angle)), pivot_pt.get_y() + (l * sin(angle)));
        };
        tLeft_corner = pivot(bLeft_corner, width, rotation_angle + M_PI/2);
        bRight_corner = pivot(bLeft_corner, length, rotation_angle);
        tRight_corner = pivot(tLeft_corner, length, rotation_angle);
        return *this;
    }
public:
    Rectangle1(const Point& BL, double l, double w, double rot_angle, Point p_dummy) : bLeft_corner(BL), length(l), width(w), rotation_angle(rot_angle), tLeft_corner(p_dummy), tRight_corner(p_dummy), bRight_corner(p_dummy) {
        normalize_angle();
        calc_update_vertices();
    }
    string print_vertices() const {
        return ("P1: " + bLeft_corner.print_point() + ", P2: " + tLeft_corner.print_point()  + ", P3: " + bRight_corner.print_point() + ", P4: " + tRight_corner.print_point());
    }
    double area() const {
        return length * width;
    }
    double perimeter() const {
        return 2 * (length + width);
    }
};

function<double(function<double(double, double)>)> my_pair(double a, double b) {
    return [a, b](function<double(double, double)> m) -> double {
        return m(a, b);
    };
}

double head(function<double(function<double(double, double)>)> a) {
    auto z = [](double a, double b) -> double {
        return a;
    };
    return a(z);
}

double tail(function<double(function<double(double, double)>)> a) {
    auto z = [](double a, double b) -> double {
        return b;
    };
    return a(z);
}

double my_pair2(int a, int b) {
    return fast_exp(2, a) * fast_exp(3, b);
}

double head2(double x) {
    function<int(double, double)> counting2s;
    counting2s = [&](double a, double b) -> int {
        if (static_cast<int>(a) % 2 != 0) {
            return static_cast<int>(b);
        }
        else {
            return counting2s(a/2, b+1);
        }
    };
    return counting2s(x, 0);
}

function<function<int(int)>(function<int(int)>)> zero() {
    return [](function<int(int)> f) -> function<int(int)> {
        return [](int x) -> int {
            return x;
        };
    };
}

function<function<int(int)>(function<int(int)>)> one() {
    return [](function<int(int)> f) -> function<int(int)> {
        return [f](int x) -> int {
            return f(x);
        };
    };
}

function<function<int(int)>(function<int(int)>)> succ(function<function<int(int)>(function<int(int)>)> n) {
    return [n](function<int(int)> f) -> function<int(int)> {
        return [f, n](int x) -> int {
            return f((n(f))(x));
        };
    };
}

function<function<int(int)>(function<int(int)>)> two() {
    return [](function<int(int)> f) -> function<int(int)> {
        return [f](int x) -> int {
            return f(f(x));
        };
    };
}

class Interval {
    private:
        double upper;
        double lower;
    public:
        Interval(string method, double a, double b) {
            if (method == "low_high") {
                if (a > b) {
                    upper = a;
                    lower = b;
                } else {
                    upper = b;
                    lower = a;
                }
            } else if (method == "center_percent") {
                lower = a - (b / 100 * a);
                upper = a + (b / 100 * a);
            } else {
                throw invalid_argument("Error: String unrecognised by class Interval");
            }
        }

        Interval(double a, double b) : Interval("low_high", a, b) {}

        double get_upper() const {
            return upper;
        }

        double get_lower() const {
            return lower;
        }

        string print_interval() const {
            return "(" + to_string(lower) + ", " + to_string(upper) + ")";
        }

        Interval& operator += (const Interval& inter2) {
            upper += inter2.upper;
            lower += inter2.lower;
            return *this;
        }
        
        Interval& operator *= (const Interval& inter2) {
            double UU = upper * inter2.upper;
            double UL = upper * inter2.lower;
            double LU = lower * inter2.upper;
            double LL = lower * inter2.lower;
            vector<double> options = {UU, UL, LU, LL};
            double lowest = numeric_limits<double>::max();
            double highest = numeric_limits<double>::lowest();
            function<void(int)> z;
            z = [&](int n) -> void {
                if (n >= options.size()) {
                    return;
                } else {
                    if (options[n] > highest) {
                        highest = options[n];
                    } 
                    if (options[n] < lowest) {
                        lowest = options[n];
                    }
                }
                z(n+1);
            };
            z(0);
            upper = highest;
            lower = lowest;
            return *this;
        }
};

Interval operator + (Interval inter1, const Interval& inter2) {
    return inter1 += inter2;
}

Interval operator - (Interval inter1, const Interval& inter2) {
    Interval inter3("low_high", -1 * inter2.get_lower(), -1 * inter2.get_upper());
    return inter1 += inter3;
}

Interval operator / (const Interval& inter1, const Interval& inter2) {
    if (inter2.get_upper() == 0 or inter2.get_lower() == 0 or (inter2.get_upper() > 0 and inter2.get_lower() < 0)) {
        throw runtime_error("Error: Invalid 2nd interval, no 0 value");
    }
    Interval inter3("low_high", 1/inter2.get_upper(), 1/inter2.get_lower());
    Interval inter4(inter1.get_lower() * inter3.get_upper(), inter1.get_upper() * inter3.get_lower());
    return inter4;
}

Interval operator * (Interval inter1, const Interval& inter2) {
    return inter1 *= inter2;
}

Interval parallel_resist(const Interval& inter1, const Interval& inter2) {
    Interval numer = inter1 * inter2;
    Interval denom = inter1 + inter2;
    return numer / denom;
}

Interval parallel_resist2(const Interval& inter1, const Interval& inter2) {
    Interval base("low_high", 1, 1);
    Interval denom = (base / inter1) + (base / inter2);
    return base / denom;
}

struct Node {
    int x;
    shared_ptr<Node> y;
    Node(int data, shared_ptr<Node> next_node) : x(data), y(next_node) {}
};

void link_existing(shared_ptr<Node>& n1_ptr, shared_ptr<Node>& n2_ptr) {
    n2_ptr->y = n1_ptr->y;
    n1_ptr->y = n2_ptr;
}

void print_links(const Node& head) {
    if (head.y == nullptr) {
        cout << head.x << endl;
    } else {
        cout << head.x << " -> ";
        print_links(*head.y);
    }
}

int length_of_links(const Node& head) {
    if (head.y == nullptr) {
        return 1;
    } else {
        return 1 + length_of_links(*head.y);
    }
}

int get_after(const Node& head, int index) {
    if (index >= length_of_links(head)) {
        throw("Error: Index out of Range!");
    }
    if (index == 0) {
        return head.x;
    } else {
        return get_after(*head.y, index-1);
    }
}

shared_ptr<Node> reverse_list(shared_ptr<Node> curr_node) {
    if (curr_node == nullptr) {
        return nullptr;
    }
    if (curr_node->y == nullptr) {
        return curr_node;
    }
    shared_ptr<Node> new_head = reverse_list(curr_node -> y);
    curr_node->y->y = curr_node;
    curr_node->y = nullptr;
    return new_head;
}

int no_of_ways(int amount, shared_ptr<Node> head_ptr) {
    if (amount == 0) {
        return 1;
    } else if (amount < 0 or head_ptr == nullptr) {
        return 0;
    } else {
        return no_of_ways(amount - head_ptr->x, head_ptr) + no_of_ways(amount, head_ptr->y);
    }
}

shared_ptr<Node> reverse_iter(shared_ptr<Node>prev_ptr, shared_ptr<Node>curr_ptr) {
    if (curr_ptr == nullptr) {
        return prev_ptr;
    }
    shared_ptr<Node> temp_ptr = curr_ptr->y;
    if (temp_ptr == nullptr) {
        curr_ptr->y = prev_ptr;
        return curr_ptr;
    } else {
        curr_ptr->y = prev_ptr;
        return reverse_iter(curr_ptr, temp_ptr);
    }
}

function<int(int)> curried_sum(int x) {
    return [x](int y) -> int {
        return x + y;
    };
}

int brooks(function<function<int(int)>(int)> curried_f, shared_ptr<Node> head_ptr) {
    if (head_ptr->y == nullptr) {
        return head_ptr->x;
    } else {
        return curried_f(brooks(curried_f, head_ptr->y))(head_ptr->x);
    }
}

shared_ptr<Node> function_map(function<double(double)> f, shared_ptr<Node> head_ptr) {
    if (head_ptr == nullptr) {
        return nullptr;
    }
    function<shared_ptr<Node>(shared_ptr<Node>, shared_ptr<Node>)> z;
    z = [f, &z](shared_ptr<Node> start_ptr, shared_ptr<Node> curr_ptr) -> shared_ptr<Node> {
        if (curr_ptr->y == nullptr) {
            curr_ptr->x = f(curr_ptr->x);
            return start_ptr;
        } else {
            curr_ptr->x = f(curr_ptr->x);
            return z(start_ptr, curr_ptr->y);
        }
    };
    return z(head_ptr, head_ptr);
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
    cout << generic_equation_solver([](double x){return sin(x);}, 0.5) << endl;
    Rational r1(-3, 4);
    r1.print_rational();
    Rational r2(1, 4);
    Rational r3(2, 3);
    r1 += r2;
    r1 *= r3;
    r2 -= r3;
    Rational sum = r2 + r3;
    r1.print_rational();
    r2.print_rational();
    sum.print_rational();
    to_rationalise(10).print_rational();
    to_rationalise(0.5).print_rational();
    Point p1(1, 2);
    cout << p1.print_point() << endl;
    Point p2(4, 6);
    LineSegment L1(p1, p2);
    cout << L1.print_seg() << endl;
    cout << L1.midpoint().print_point() << endl;
    cout << L1.length() << endl;
    cout << head(my_pair(1,2)) << endl;
    cout << head2(my_pair2(1, 2)) << endl;
    cout << parallel_resist(Interval(1, 2), Interval(1, 3)).print_interval() << endl;
    cout << parallel_resist2(Interval(1, 2), Interval(1, 3)).print_interval() << endl;
    shared_ptr<Node> head = make_shared<Node>(1, nullptr);
    shared_ptr<Node> second = make_shared<Node>(2, nullptr);
    shared_ptr<Node> third = make_shared<Node>(3, nullptr);
    shared_ptr<Node> fourth = make_shared<Node>(4, nullptr);
    shared_ptr<Node> fifth = make_shared<Node>(5, nullptr);
    link_existing(head, second);
    link_existing(head, third);
    link_existing(second, fourth);
    link_existing(third, fifth);
    print_links(*head);
    print_links(*reverse_iter(nullptr, head));
    cout << length_of_links(*head) << endl;
    cout << get_after(*second, 1) << endl;
    print_links(*reverse_list(fourth));
    cout << no_of_ways(100, fourth) << endl;
    cout << brooks(curried_sum, head) << endl;
    print_links(*function_map([](double x){return x * 10;}, head));
    print_links(*function_map([](double x){ return x * x;}, head));
}