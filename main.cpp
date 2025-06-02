#include <iostream>
#include <algorithm>
#include <chrono>
#include <iomanip>

int knapsackBruteForce(int W, int weights[], int values[], int n);
int knapsackDP(int w, int weights[], int values[], int n);
int knapsackBacktracking(int w, int weights[], int values[], int n, int i, int current_weight, int current_value);
int knapsackBranchAndBound(int w, int weights[], int values[], int n, int i, int current_weight, int current_value, int& best_value);

int main() {
    std::setlocale(0, "");

    const int n = 15;
    int values[n];
    int weights[n];
    int w = 2000;

    for (int i = 0; i < n; i++) {
        weights[i] = 10 + i % 20;
        values[i] = 50 + i % 100;
    }

    std::cout << std::fixed << std::setprecision(6);

    auto start = std::chrono::high_resolution_clock::now();
    int bruteForceResult = knapsackBruteForce(w, weights, values, n);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> bruteForceTime = end - start;
    std::cout << "����� ������ ����:                    " << bruteForceResult << " (�����: " << bruteForceTime.count() << " ���)\n";

    start = std::chrono::high_resolution_clock::now();
    int dpResult = knapsackDP(w, weights, values, n);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dpTime = end - start;
    std::cout << "����� ������������� ����������������: " << dpResult << " (�����: " << dpTime.count() << " ���)\n";

    start = std::chrono::high_resolution_clock::now();
    int backtrackingResult = knapsackBacktracking(w, weights, values, n, 0, 0, 0);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> backtrackingTime = end - start;
    std::cout << "����� ������ � ���������:             " << backtrackingResult << " (�����: " << backtrackingTime.count() << " ���)\n";

    int best_value = 0;
    start = std::chrono::high_resolution_clock::now();
    int branchAndBoundResult = knapsackBranchAndBound(w, weights, values, n, 0, 0, 0, best_value);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> branchAndBoundTime = end - start;
    std::cout << "����� ������ � ������:                " << branchAndBoundResult << " (�����: " << branchAndBoundTime.count() << " ���)\n";

    return 0;
}

int knapsackBruteForce(int W, int weights[], int values[], int n) {
    int max_value = 0;
    for (int mask = 0; mask < (1 << n); mask++) {
        int current_weight = 0, current_value = 0;
        for (int i = 0; i < n; i++) {
            if (mask & (1 << i)) {
                current_weight += weights[i];
                current_value += values[i];
            }
        }
        if (current_weight <= W && current_value > max_value) {
            max_value = current_value;
        }
    }
    return max_value;
}

int* knapsackDP(int w, int weights[], int values[], int n) {
    int** dp = new int* [w + 1];
    std::fill(dp, dp + w + 1, 0);

    for (int i = 0; i < n; i++) {
        for (int j = w; j >= weights[i]; j--) {
            dp[j] = std::max(dp[j], dp[j - weights[i]] + values[i]);
        }
    }

    return dp[w];
}

int knapsackBacktracking(int w, int weights[], int values[], int n, int i, int current_weight, int current_value) {
    if (i == n || current_weight >= w) {
        return current_value;
    }
    if (current_weight + weights[i] > w) {
        return knapsackBacktracking(w, weights, values, n, i + 1, current_weight, current_value);
    }
    return std::max(
        knapsackBacktracking(w, weights, values, n, i + 1, current_weight + weights[i], current_value + values[i]),
        knapsackBacktracking(w, weights, values, n, i + 1, current_weight, current_value)
    );
}

int knapsackBranchAndBound(int w, int weights[], int values[], int n, int i, int current_weight, int current_value, int& best_value) {
    if (current_weight <= w && current_value > best_value) {
        best_value = current_value;
    }
    if (i == n || current_weight >= w) {
        return best_value;
    }
    double bound = current_value;
    int remaining_weight = w - current_weight;
    int j = i;
    while (j < n && weights[j] <= remaining_weight) {
        remaining_weight -= weights[j];
        bound += values[j];
        j++;
    }
    if (j < n) {
        bound += (double)values[j] * remaining_weight / weights[j];
    }
    if (bound <= best_value) {
        return best_value;
    }
    knapsackBranchAndBound(w, weights, values, n, i + 1, current_weight + weights[i], current_value + values[i], best_value);
    knapsackBranchAndBound(w, weights, values, n, i + 1, current_weight, current_value, best_value);
    return best_value;
}
