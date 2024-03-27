#include <utility>
#include <vector>
#include <iostream>
#include <cmath>
#include <random>
//#include <string>
#include <algorithm>
//#include <ctime>
#include <iomanip>

namespace mhe {
    struct Point {
        double x;
        double y;
    };

    class DistanceMatrix{
        std::vector<std::vector<double>> distanceMatrix;
        static double distance(const Point &p1, const Point &p2) {
            return std::sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y));
        }
    public:
        explicit DistanceMatrix(const std::vector<Point> &points){
            int n = int(points.size());
            std::vector<std::vector<double>> distances(n, std::vector<double>(n, 0.0));
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    distances[i][j] = distance(points[i], points[j]);
                }
            }
            this->distanceMatrix = distances;
        }

        explicit DistanceMatrix(std::vector<std::vector<double>> distances){
            this->distanceMatrix = std::move(distances);
        }

        DistanceMatrix(){
            this->distanceMatrix = {{0, 1, 1}, {1, 0, 1}, {1,1,0}};
        }

        void print(int precision = 5){
            for (const auto& row : distanceMatrix) {
                std::cout << std::fixed << std::setprecision(precision);
                for (const auto& dist : row) {
                    std::cout << dist << "|";
                }
                std::cout << "\n";
            }
        }

        double getDistance(int x, int y){
            return this->distanceMatrix[x-1][y-1];
        }

        //metoda celu
        [[nodiscard]]double cost(const std::vector<int>& solution){
            double cost = 0.0;
            for(int i = 0; i < solution.size()-1 ; i++){
                cost += getDistance(solution[i], solution[i+1]);
            }
            return cost;
        }

        int size(){
            return int(this->distanceMatrix.size());
        }
    };

    [[nodiscard]]std::vector<Point> generateRandomPoints(unsigned int num_points, int gen_seed = 42) {
        std::random_device rd;
        std::mt19937 gen(gen_seed);
        std::uniform_real_distribution<> dis(0.0, 1.0);
        std::vector<Point> points;
        for (int i = 0; i < num_points; ++i) {
            points.push_back({dis(gen), dis(gen)});
        }
        return points;
    }

    void drawSolution(const std::vector<int>& solution){
        std::string tekst = "{ ";
        for (auto point : solution){
            tekst += std::to_string(point);
            tekst += " ";
        }
        tekst += "}";
        std::cout << tekst << std::endl;
    }

    //funkcja zwracająca sąsiadów rozwiązania
    auto getNeighbors(const std::vector<int>& vector) -> std::vector<std::vector<int>>{
        std::vector<std::vector<int>> set_of_solutions;
        for(int i = 0; i < vector.size()-1; i++){
            for(int j = i+1; j<vector.size(); j++){
                auto v = vector;
                std::swap(v[i],v[j]);
                set_of_solutions.push_back(v);
            }
        }
        return set_of_solutions;
    }
    //generowanie losowego rozwiązania
    std::vector<int> generateRandomSolution(int size, int gen_seed=42){
        std::vector<int> result;
        std::random_device rd;
        std::mt19937 gen(gen_seed);
        std::uniform_int_distribution<> lop(0, size*size);
        std::uniform_int_distribution<> dis(0, size-1);
        for (int i = 1; i <= size; ++i) {
            result.push_back(i);
        }
        for(int i =0; i<lop(gen); i++){
            std::swap(result[dis(gen)], result[dis(gen)]);
        }
        return result;
    }

}

//algorytm pełnego przeglądu
namespace brute{
    std::vector<std::vector<int>> generatePermutations(std::vector<int>& nums, int start) {
        std::vector<std::vector<int>> result;
        if (start == nums.size() - 1) {
            result.push_back(nums);
            return result;
        }
        for (int i = start; i < nums.size(); ++i) {
            std::swap(nums[start], nums[i]);
            std::vector<std::vector<int>> permutations = generatePermutations(nums, start + 1);
            result.insert(result.end(), permutations.begin(), permutations.end());
            std::swap(nums[start], nums[i]);
        }
        return result;
    }

    [[nodiscard]]std::vector<std::vector<int>> generateAllPermutations(int x) {
        std::vector<int> result;
        for (int i = 1; i <= x; ++i) {
            result.push_back(i);
        }
        return generatePermutations(result, 0);
    }

    std::vector<int> solve(mhe::DistanceMatrix distances){
        auto permutations = generateAllPermutations(distances.size());
        double best_val = distances.cost(permutations[0]);
        int best_sol = 0;
        for(int i = 1; i < permutations.size(); i++){
            auto curr_val = distances.cost(permutations[i]);
            if(curr_val < best_val){
                best_val = curr_val;
                best_sol = i;
            }

        }
        return permutations[best_sol];
    }

}

//algorytm wspinaczkowy
namespace climb{
    std::vector<int> getBestNeighbor(const std::vector<int>& solution, mhe::DistanceMatrix distances){
        auto neighbors =  mhe::getNeighbors(solution);
        int best_sol = 0;
        double best_val = distances.cost(neighbors[0]);
        for(int i = 1; i < neighbors.size(); i++){
            auto curr_val = distances.cost(neighbors[i]);
            if(curr_val < best_val){
                best_sol = i;
                best_val = curr_val;
            }
        }
        if(best_val > distances.cost(solution)) return {0};
        return neighbors[best_sol];
    }
    std::vector<int> getRandomNeighbor(std::vector<int> solution){
        std::random_device rd;
        std::uniform_int_distribution<> dis(0, int(solution.size())-1);
        int i = dis(rd);
        int j = dis(rd);
        std::swap(solution[i],solution[j]);
        return solution;
    }

    //deterministyczna wersja
    std::vector<int> solve(mhe::DistanceMatrix distances, int loop_breaker = 1000){
        auto solution =  mhe::generateRandomSolution(distances.size());
        for(int i = 0; i < loop_breaker; i++){
            auto newSolution = getBestNeighbor(solution, distances);
            if(newSolution[0]==0) break;
            solution = newSolution;
        }
        return solution;
    }

    //wersja z losowym somsiadem
    auto solveRandom(mhe::DistanceMatrix distances, int  max_iterations = 10000) -> std::vector<int>{
        auto solution = mhe::generateRandomSolution(int(distances.size()));
        double best_value = distances.cost(solution);
        auto best_solution = solution;
        for(int i = 0; i < max_iterations; i++){
            solution = getRandomNeighbor(solution);
            double value = distances.cost(solution);
            if(value == best_value) break;
            else if (value < best_value){
                best_value = value;
                best_solution = solution;
            }
        }
        return best_solution;
    }
}

namespace tabu{
    std::vector<int> solve(mhe::DistanceMatrix distances, int max_iter=1000, int max_tab_size=5){
        //TODO cofanie sie do ostatniego punktu roboczego
        auto best_sol = mhe::generateRandomSolution(distances.size());
        auto best_candidate = best_sol;
        std::vector<std::vector<int>> tabu_list = {};
        for(int i = 0; i < max_iter; i++){
            auto neighbors = mhe::getNeighbors(best_candidate);
            double best_candidate_val = 0.0;
            for(auto candidate : neighbors){
                if(std::find(tabu_list.begin(), tabu_list.end(), candidate)==tabu_list.end()){
                    best_candidate = candidate;
                    best_candidate_val = distances.cost(best_candidate);
                }
                if(best_candidate_val == 0) break;
                if(best_candidate_val < distances.cost(best_sol)) best_sol = best_candidate;
                tabu_list.push_back(best_candidate);
                if(tabu_list.size() > max_tab_size) tabu_list.erase(tabu_list.begin());
            }
        }
        return best_sol;
    }
}


auto main() -> int{
    using namespace std;
    const unsigned int num_cities = 7;
    auto dis = mhe::DistanceMatrix(mhe::generateRandomPoints(num_cities));
    dis.print();

    auto brute_sol= brute::solve(dis);
    cout << "Bruteforce best solution: ";
    mhe::drawSolution(brute_sol);
    cout << "with cost: " << dis.cost(brute_sol) << endl;

    auto climb_sol = climb::solve(dis);
    cout << "Climb algorithm (deter. ver.) best solution: ";
    mhe::drawSolution(climb_sol);
    cout << "with cost: " << dis.cost(climb_sol) << endl;

    auto rand_climb_sol = climb::solveRandom(dis);
    cout << "Climb algorithm (random ver.) best solution: ";
    mhe::drawSolution(rand_climb_sol);
    cout << "with cost: " << dis.cost(rand_climb_sol) << endl;

    auto tabu_sol = tabu::solve(dis);
    cout << "Tabu algorithm best solution: ";
    mhe::drawSolution(tabu_sol);
    cout << "with cost: " << dis.cost(tabu_sol) << endl;
    return 0;
}
