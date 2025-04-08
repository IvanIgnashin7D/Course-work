#include <iostream>
#include <vector>
#include <algorithm>


int knapsackBruteForce(int W, std::vector<int>& weights, std::vector<int>& values, int n); // grubaya sila
int knapsackDP(int W, std::vector<int>& weights, std::vector<int>& values, int n); // dinamicheskoe programmirovanie
int knapsackBacktracking(int W, std::vector<int>& weights, std::vector<int>& values, int i, int current_weight, int current_value); // poisk s vozvratom
int knapsackBranchAndBound(int W, std::vector<int>& weights, std::vector<int>& values, int i, int current_weight, int current_value, int& best_value); // vetvi i granitsi


int main() {
	std::cout << "Hello world" << '\n';
    return 0;
}


int knapsackBruteForce(int W, std::vector<int>& weights, std::vector<int>& values, int n) {
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


int knapsackDP(int W, std::vector<int>& weights, std::vector<int>& values, int n) {
    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(W + 1, 0));
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= W; j++) {
            if (weights[i - 1] <= j) {
                dp[i][j] = std::max(dp[i - 1][j], values[i - 1] + dp[i - 1][j - weights[i - 1]]);
            }
            else {
                dp[i][j] = dp[i - 1][j];
            }
        }
    }
    return dp[n][W];
}


int knapsackBacktracking(int W, std::vector<int>& weights, std::vector<int>& values, int i, int current_weight, int current_value) {
    if (i == weights.size() || current_weight >= W) {
        return current_value;
    }
    if (current_weight + weights[i] > W) {
        return knapsackBacktracking(W, weights, values, i + 1, current_weight, current_value);
    }
    return std::max(
        knapsackBacktracking(W, weights, values, i + 1, current_weight + weights[i], current_value + values[i]),
        knapsackBacktracking(W, weights, values, i + 1, current_weight, current_value)
    );
}


int knapsackBranchAndBound(int W, std::vector<int>& weights, std::vector<int>& values, int i, int current_weight, int current_value, int& best_value) {
    if (current_weight <= W && current_value > best_value) {
        best_value = current_value;
    }
    if (i == weights.size() || current_weight >= W) {
        return best_value;
    }
    double bound = current_value;
    int remaining_weight = W - current_weight;
    int j = i;
    while (j < weights.size() && weights[j] <= remaining_weight) {
        remaining_weight -= weights[j];
        bound += values[j];
        j++;
    }
    if (j < weights.size()) {
        bound += (double)values[j] * remaining_weight / weights[j];
    }
    if (bound <= best_value) {
        return best_value;
    }
    knapsackBranchAndBound(W, weights, values, i + 1, current_weight + weights[i], current_value + values[i], best_value);
    knapsackBranchAndBound(W, weights, values, i + 1, current_weight, current_value, best_value);
    return best_value;
}
