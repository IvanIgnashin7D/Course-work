#include <iostream>
#include <algorithm>
#include <chrono>
#include <iomanip>

int knapsackBruteForce(int W, int weights[], int values[], int n);
int knapsackDP(int w, int weights[], int values[], int n);
int knapsackBacktracking(int w, int weights[], int values[], int n, int i, int current_weight, int current_value);
void sortItemsByValuePerWeight(double* valuePerWeight, int* indices, int n);
int knapsackBranchAndBound(int w, int* weights, int* values, int n);
int knapsackBranchAndBound(int w, int* weights, int* values, int n, int i, int current_weight, int current_value, int* best_value);
inline int max(int a, int b);

int main() {
    std::setlocale(0, "");

    const int n = 25;
    int values[n] = {
        45, 72, 18, 91, 33, 120, 25, 55, 80, 30,
        68, 95, 12, 50, 77, 60, 110, 40, 85, 65,
        42, 78, 22, 88, 37
    };

    int weights[n] = {
        8, 12, 5, 17, 9, 23, 6, 11, 15, 7,
        14, 19, 3, 10, 13, 11, 27, 8, 16, 12,
        7, 13, 4, 16, 8
    };
    int w = 50;


    std::cout << std::fixed << std::setprecision(6);

    auto start = std::chrono::high_resolution_clock::now();
    int bruteForceResult = knapsackBruteForce(w, weights, values, n);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> bruteForceTime = end - start;
    std::cout << "Метод грубой силы:                    " << bruteForceResult << " (Время: " << bruteForceTime.count() << " сек)\n";

    start = std::chrono::high_resolution_clock::now();
    int dpResult = knapsackDP(w, weights, values, n);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dpTime = end - start;
    std::cout << "Метод динамического программирования: " << dpResult << " (Время: " << dpTime.count() << " сек)\n";

    start = std::chrono::high_resolution_clock::now();
    int backtrackingResult = knapsackBacktracking(w, weights, values, n, 0, 0, 0);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> backtrackingTime = end - start;
    std::cout << "Метод поиска с возвратом:             " << backtrackingResult << " (Время: " << backtrackingTime.count() << " сек)\n";

    start = std::chrono::high_resolution_clock::now();
    int branchAndBoundResult = knapsackBranchAndBound(w, weights, values, n);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> branchAndBoundTime = end - start;
    std::cout << "Метод ветвей и границ:                " << branchAndBoundResult << " (Время: " << branchAndBoundTime.count() << " сек)\n";

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

int knapsackDP(int w, int weights[], int values[], int n) {
    int* prev = new int[w + 1]();
    int* curr = new int[w + 1];

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= w; j++) {
            if (weights[i - 1] <= j) {
                curr[j] = max(prev[j], values[i - 1] + prev[j - weights[i - 1]]);
            }
            else {
                curr[j] = prev[j];
            }
        }
        std::swap(prev, curr);
    }

    int result = prev[w];
    delete[] prev;
    delete[] curr;
    return result;
}

int knapsackBacktracking(int w, int weights[], int values[], int n, int i, int current_weight, int current_value) {
    if (i == n || current_weight >= w) {
        return current_value;
    }
    if (current_weight + weights[i] > w) {
        return knapsackBacktracking(w, weights, values, n, i + 1, current_weight, current_value);
    }
    return max(
        knapsackBacktracking(w, weights, values, n, i + 1, current_weight + weights[i], current_value + values[i]),
        knapsackBacktracking(w, weights, values, n, i + 1, current_weight, current_value)
    );
}

int knapsackBranchAndBound(int w, int* weights, int* values, int n) {
    double* valuePerWeight = new double[n];
    int* indices = new int[n];
    for (int j = 0; j < n; ++j) {
        valuePerWeight[j] = (double)values[j] / weights[j];
        indices[j] = j;
    }
    sortItemsByValuePerWeight(valuePerWeight, indices, n);

    int* sortedWeights = new int[n];
    int* sortedValues = new int[n];
    for (int j = 0; j < n; ++j) {
        sortedWeights[j] = weights[indices[j]];
        sortedValues[j] = values[indices[j]];
    }

    for (int j = 0; j < n; ++j) {
        weights[j] = sortedWeights[j];
        values[j] = sortedValues[j];
    }

    delete[] valuePerWeight;
    delete[] indices;
    delete[] sortedWeights;
    delete[] sortedValues;

    int best_value = 0;
    return knapsackBranchAndBound(w, weights, values, n, 0, 0, 0, &best_value);
}

int knapsackBranchAndBound(int w, int* weights, int* values, int n,
    int i, int current_weight, int current_value, int* best_value) {

    if (current_weight <= w && current_value > *best_value) {
        *best_value = current_value;
    }
    if (i == n || current_weight >= w) {
        return *best_value;
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

    if (bound <= *best_value) {
        return *best_value;
    }

    knapsackBranchAndBound(w, weights, values, n, i + 1,
        current_weight + weights[i], current_value + values[i], best_value);
    knapsackBranchAndBound(w, weights, values, n, i + 1,
        current_weight, current_value, best_value);

    return *best_value;
}

void sortItemsByValuePerWeight(double* valuePerWeight, int* indices, int n) {
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            if (valuePerWeight[j] < valuePerWeight[j + 1]) {
                double tempVal = valuePerWeight[j];
                valuePerWeight[j] = valuePerWeight[j + 1];
                valuePerWeight[j + 1] = tempVal;

                int tempIdx = indices[j];
                indices[j] = indices[j + 1];
                indices[j + 1] = tempIdx;
            }
        }
    }
}

inline int max(int a, int b) {
    return (a > b) ? a : b;
}