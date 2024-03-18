#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <string>
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
        std::cout << "Wartosc funkcji celu dla " << drawSolution(solution) << ": " << distance << std::endl;
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
        return vector;
    }

    auto solve_random(const std::vector<std::vector<double>>& distances) -> std::vector<int>{
        auto solution = getVector(int(distances.size()));
        double best_value = mhe::target(solution,distances);
        auto best_solution = solution;
        while(true){
            solution = change_random(solution);
            double value = mhe::target(solution, distances);
            if(value == best_value) break;
            else if (value < best_value){
                best_value = value;
                best_solution = solution;
            }
        }
        std::cout << "Best solution: " << mhe::drawSolution(best_solution) << " with value: " << best_value << std::endl;
        return best_solution;
    }
}


auto main() -> int{
    using namespace mhe;
    auto dis_matrix = distanceMatrix(generateRandomPoints(5));
    drawDistances(dis_matrix);
    wsp::solve_random(dis_matrix);
    return 0;
}

