#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <string>
#include <algorithm>
#include <ctime>
namespace mhe {
    struct Point {
        double x;
        double y;
    };


    double distance(const Point &p1, const Point &p2) {
        return std::sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y));
    }

    std::vector<Point> generateRandomPoints(int num_points) {
        std::random_device rd;
        std::mt19937 gen(42);
        std::uniform_real_distribution<> dis(0.0, 1.0);

        std::vector<Point> points;
        for (int i = 0; i < num_points; ++i) {
            points.push_back({dis(gen), dis(gen)});
        }
        return points;
    }

    std::vector<std::vector<double>> distanceMatrix(const std::vector<Point> &points) {
        int n = int(points.size());
        std::vector<std::vector<double>> distances(n, std::vector<double>(n, 0.0));

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                distances[i][j] = distance(points[i], points[j]);
            }
        }
        return distances;
    }

    std::string drawSolution(const std::vector<int>& solution){
        std::string tekst = "rozwiazanie: {";
        for (auto point : solution){
            tekst += std::to_string(point);
            tekst += ", ";
        }
        tekst += "}";
        return tekst;
    }


    auto target(const std::vector<int>& solution, const std::vector<std::vector<double>>& distances) -> double{
        double distance = 0.0;
        for (auto i = 0; i < solution.size()-1; i++){
            distance += distances[solution[i]][solution[i+1]];
        }
        return distance;
    }

    void drawDistances(const std::vector<std::vector<double>>& distances){
        for (const auto& row : distances) {
            for (const auto& dist : row) {
                std::cout << dist << " ";
            }
            std::cout << "\n";
        }
    }

    auto permute(const std::vector<int>& vector) -> std::vector<std::vector<int>>{
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
}
//rozwiązanie algorytmem wspinaczkowym
namespace wsp{
    auto change_random(std::vector<int> solution)->std::vector<int>{
        std::random_device rd;
        std::uniform_int_distribution<> dis(0, int(solution.size())-1);
        int i = dis(rd);
        int j = dis(rd);
        std::swap(solution[i],solution[j]);
        return solution;
    }

    auto getVector(int size) -> std::vector<int>{
        std::vector<int> vector;
        for (int i = 1; i < size; i++){
           vector.push_back(i);
        }
        std::mt19937 rng(std::time(nullptr));
        std::shuffle(vector.begin(), vector.end(), rng);
        return vector;
    }

    auto solve_random(const std::vector<std::vector<double>>& distances, int  max_iterations = 10000) -> std::vector<int>{
        auto solution = getVector(int(distances.size()));
        double best_value = mhe::target(solution,distances);
        auto best_solution = solution;
        for(int i = 0; i < max_iterations; i++){
            solution = change_random(solution);
            double value = mhe::target(solution, distances);
            if(value == best_value) break;
            else if (value < best_value){
                best_value = value;
                best_solution = solution;
            }
        }
        std::cout << "Best solution (random method): " << mhe::drawSolution(best_solution) << " with value: " << best_value << std::endl;
        return best_solution;
    }

    auto get_best_solution(const std::vector<int>& solution, const std::vector<std::vector<double>>& distances)->std::vector<int>{
        auto set = mhe::permute(solution);
        auto best_value = mhe::target(solution, distances);
        auto best_solution = solution;
        for(const auto & i : set){
            double tmp_val = mhe::target(i,distances);
            if(tmp_val < best_value){
                best_solution = i;
                best_value = tmp_val;
            }
        }
        if(solution==best_solution) return {-1};
        return best_solution;
    }

    auto solve(const std::vector<std::vector<double>>& distances) -> std::vector<int>{
        auto solution = getVector(int(distances.size()));
        double best_value = mhe::target(solution,distances);
        auto best_solution = solution;
        for(;;){
            solution = get_best_solution(solution, distances);
            if(solution[0]==-1) break;
            best_solution = solution;
            std::cout << ".";
        }
        std::cout << "Best solution (best neighbour): " << mhe::drawSolution(best_solution) << " with value: " << best_value << std::endl;
        return best_solution;
    }
}


auto main() -> int{
    using namespace mhe;
    auto dis_matrix = distanceMatrix(generateRandomPoints(200));
    wsp::solve_random(dis_matrix);
    wsp::solve(dis_matrix);
    return 0;
}
