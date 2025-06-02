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

int main() {
    cout << sqrt_myfunc(1000) << endl;
    cout << A(1, 10) << endl ;
    cout << fib(10) << endl;
    cout << f1_recurse(4) << endl;
    cout << f1_iterFull(4) << endl;
    vector<int> row6 = Pascal(6);
    for (int i = 0; i < row6.size(); i++) {
        cout << row6[i] << ", ";
    }
    cout << endl;
}