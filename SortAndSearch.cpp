#include <iostream>
#include <bits/stdc++.h>
#include <vector>
using namespace std;

vector<int> AscInsertionSort(vector<int> data) {
    if (data.empty() || data.size() == 1) {
        return data;
    }
    else {
        vector<int> check = {data[0]};
        for (int i=1; i<data.size(); i++) {
            int key = data[i];
            int j=i-1; 
            while (j>=0 && key < check[j]) {
                j-=1;
            }
            if (j<0) {
                check.insert(check.begin(), key);
            }
            else if (j == i-1) {
                check.push_back(key);
            }
            else {
                check.insert(check.begin()+j+1, key);
            }
        }
        return check;
    }
}

vector<int> Merge(vector<int> sorted1, vector<int> sorted2) {
    vector<int> combined = {};
    int i = 0;
    int j = 0;
    while (i<sorted1.size() && j<sorted2.size()) {
        if (sorted1[i]<sorted2[j]) {
            combined.push_back(sorted1[i]);
            i += 1;
        }
        else if (sorted1[i]>=sorted2[j]) {
            combined.push_back(sorted2[j]);
            j += 1;
        }
    }
    while (i >= sorted1.size() && j < sorted2.size()) {
        combined.push_back(sorted2[j]);
        j += 1;
    }
    while (j >= sorted2.size() && i < sorted1.size()) {
        combined.push_back(sorted1[i]);
        i += 1;
    }
    return combined;
}

vector<int> MergeSort(vector<int> data) {
    int length = data.size();
    if (data.size() <= 1) {
        return data;
    }
    else {
        vector<int> subData1 = {};
        vector<int> subData2 = {};
        for (int i = 0; i <= floor(length/2)-1; i++) {
            subData1.push_back(data[i]);
        }
        for (int i = floor(length/2); i<length; i++) {
            subData2.push_back(data[i]);
        }
        subData1 = MergeSort(subData1);
        subData2 = MergeSort(subData2);
        return Merge(subData1, subData2);
    }
}

vector<int> MergeSortChanged(vector<int> data) {
    int length = data.size();
    if (data.size() <= 1) {
        return data;
    }
    else {
        vector<int> subData1 = {};
        vector<int> subData2 = {};
        for (int i = 0; i <= data.size()-2; i++) {
            subData1.push_back(data[i]);
        }
        subData2.push_back(data[data.size()-1]);
        subData1 = MergeSortChanged(subData1);
        subData2 = MergeSortChanged(subData2);
        return Merge(subData1, subData2);
    }
}

vector<int> BinarySearch(vector<int> data, int search) {
    if (data.size() == 1) {
        if (search != data[0]) {
            return {0};
        }
        else {
            return {1};
        }
    }
    else if (search == data[floor(data.size()/2)]) {
        return {1};
    }
    else {
        vector<int> newData = {};
        if (search > data[floor(data.size()/2)]) {
            for (int i = floor(data.size()/2); i < data.size(); i++) {
                newData.push_back(data[i]);
            }
        }
        else {
            for (int i = 0; i <= floor(data.size()/2); i++) {
                newData.push_back(data[i]);
            }
        }
        newData = BinarySearch(newData, search);
        return newData;
    }
}

int main() {
    vector<int>check = AscInsertionSort({13, 25, 62, 35, 14, 72});
    for (int i=0; i<check.size(); i++) {
        cout << check[i] << " ";
    }
    cout << endl;
    vector<int>checkMerge = Merge({1}, {1});
    for (int i=0; i<checkMerge.size(); i++) {
        cout << checkMerge[i] << " ";
    }
    cout << endl;

    vector<int>checkMergeSort = MergeSort({2, 1, 3, 4, 10, 5, 2});
    for (int i=0; i<checkMergeSort.size(); i++) {
        cout << checkMergeSort[i] << " ";
    }
    cout << endl;

    vector<int>checkMergeSort2 = MergeSortChanged({2, 1, 3, 4, 10, 5, 2});
    for (int i=0; i<checkMergeSort2.size(); i++) {
        cout << checkMergeSort2[i] << " ";
    }
    cout << endl;

    vector<int>check1 = BinarySearch({1, 2, 10, 102, 300}, 102);
    cout << check1[0];
}