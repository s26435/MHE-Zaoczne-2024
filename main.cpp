#include <utility>
#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <iomanip>
#include <map>
#include <chrono>

#define SEED 42

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

    [[nodiscard]]std::vector<Point> generateRandomPoints(unsigned int num_points, int gen_seed = SEED) {
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

    std::vector<std::vector<int>> getNeighbors(const std::vector<int>& vector){
        std::vector<std::vector<int>> solutions;
        for(int i = 0; i < vector.size()-1; i++){
            for(int j = i+1; j<vector.size(); j++){
                auto v = vector;
                std::swap(v[i],v[j]);
                solutions.push_back(v);
            }
        }
        return solutions;
    }

    std::vector<int> generateRandomSolution(int size, int gen_seed=SEED){
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

    std::vector<int> solve(mhe::DistanceMatrix distances, int loop_breaker = 1000){
        auto solution =  mhe::generateRandomSolution(distances.size());
        for(int i = 0; i < loop_breaker; i++){
            auto newSolution = getBestNeighbor(solution, distances);
            if(newSolution[0]==0) break;
            solution = newSolution;
        }
        return solution;
    }

    std::vector<int> solveRandom(mhe::DistanceMatrix distances, int  max_iterations = 10000){
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
        auto best_sol = mhe::generateRandomSolution(distances.size());
        auto best_candidate = best_sol;
        std::vector<std::vector<int>> tabu_list = {};
        for(int i = 0; i < max_iter; i++){
            auto neighbors = mhe::getNeighbors(best_candidate);
            double best_candidate_val = 0.0;
            for(const auto& candidate : neighbors){
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

namespace ann{
    std::vector<int> solve(mhe::DistanceMatrix distances, double initial_temperature=25, int max_iterations = 1000){
        auto s = mhe::generateRandomSolution(distances.size());
        std::vector<int> best_s = s;
        auto best_val = distances.cost(s);
        std::random_device rd;
        std::uniform_int_distribution<> dis(0, 1);
        double temperature;
        for(int i = 0; i < max_iterations; i++){
            temperature = initial_temperature * (1 - static_cast<double>(i + 1) / max_iterations);
            auto news = climb::getRandomNeighbor(s);
            auto curr_val = distances.cost(news);
            if(curr_val<best_val||std::exp((best_val - curr_val) / temperature)>dis(rd)){
                s = news;
                if(distances.cost(s)<distances.cost(best_s)) best_s = s;
            }
        }
        return best_s;
    }
}

namespace gen {

    std::vector<int> uniformCrossover(std::vector<int> parent1, std::vector<int> parent2) {
        if (parent1.size() != parent2.size()) throw std::invalid_argument("Rodzice musza miec taka sama dlugosc!");
        double crossover_rate = 0.5;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        std::vector<int> child(parent1.size());
        for (size_t i = 0; i < parent1.size(); ++i) {
            if (dis(gen) < crossover_rate) child[i] = parent2[i];
            else child[i] = parent1[i];
        }
        return child;
    }

    std::vector<int> orderCrossover(std::vector<int> parent1, std::vector<int> parent2){
        if (parent1.size() != parent2.size()) throw std::invalid_argument("Rodzice musza miec taka sama dlugosc!");
        std::vector<int> child(parent1.size());
        for([[maybe_unused]] auto x : child) x = 0;
        std::map<int,int> missing;
        std::random_device rd;
        std::uniform_int_distribution<> dis(0, int(parent1.size()) - 1);
        auto a = dis(rd);
        auto b = dis(rd);
        while (a == b) b = dis(rd);
        if (a > b) std::swap(a, b);
        for(int i = 0; i < parent1.size(); i++){
            if(i>=a&&i<=b) child[i] = parent1[i];
            else missing[*std::find(parent1.begin(), parent1.end(), parent2[i])] = parent1[i];
        }
        for(int i = 0; i < parent1.size(); i++){
            if(child[i]==0) {
                int temp = missing.begin()->first;
                child[i] = missing[temp];
                missing.erase(temp);
            }
        }
        return child;
    }

    std::vector<int> rotation(std::vector<int> solution){
        std::random_device rd;
        std::uniform_int_distribution<> dis(0, int(solution.size()) - 1);
        std::uniform_int_distribution<> k(0, int(solution.size()));
        int a = dis(rd);
        int b = dis(rd);
        int c = k(rd);
        while (a == b) b = dis(rd);
        if (a > b) std::swap(a, b);
        std::vector<int> mutated(solution.size());
        std::vector<int> missing;
        for(int i = 0; i < solution.size(); i++){
            if(i>=a&&i<=b) mutated[i] = solution[i];
            else {
                missing.push_back(solution[i]);
                mutated[i] = 0;
            }
        }
        std::vector<int> rotated(missing.size());
        for (int i = 0; i < missing.size(); ++i) {
            rotated[(i + c) % missing.size()] = missing[i];
        }
        int ite = 0;
        for(int i = 0; i < solution.size(); i++){
            if(mutated[i]==0) {
                mutated[i] = rotated[ite];
                ite++;
            }
        }
        return mutated;
    }

    std::vector<int> inversion(std::vector<int> solution){
        std::random_device rd;
        std::uniform_int_distribution<> dis(0, int(solution.size()) - 1);
        int a = dis(rd);
        int b = dis(rd);
        while (a == b) b = dis(rd);
        if (a > b) std::swap(a, b);
        std::vector<int> mutated(solution.size());
        std::vector<int> missing;
        for(int i = 0; i < solution.size(); i++){
            if(i>=a&&i<=b) {
                missing.push_back(solution[i]);
                mutated[i] = 0;
            }
            else mutated[i] = solution[i];
        }
        for (int i = 0; i < missing.size()/2; ++i) {
            std::swap(missing[i], missing[missing.size()-i-1]);
        }
        int ite = 0;
        for(int i = 0; i < solution.size(); i++){
            if(mutated[i]==0) {
                mutated[i] = missing[ite];
                ite++;
            }
        }
        return mutated;
    }

    std::vector<int> tournamentSelection(const std::vector<std::vector<int>>& population, mhe::DistanceMatrix distances ){
        int tournament_size = 10;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, int(population.size()) - 1);
        std::vector<int> best_solution = population[dis(gen)];
        for(int i =0; i < tournament_size; i++){
            const std::vector<int>& random_solution = population[dis(gen)];
            if(distances.cost(random_solution)<distances.cost(best_solution)) best_solution = random_solution;
        }
        return best_solution;
    }

    std::vector<int> rouletteWheelSelection(const std::vector<std::vector<int>>& population, mhe::DistanceMatrix distances) {
        double totalFitness = 0;
        for (const auto& individual : population) {
            totalFitness += distances.cost(individual);
        }
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, totalFitness - 1);
        double selectedPoint = dis(gen);
        double cumulativeFitness = 0;
        for (const auto& individual : population) {
            cumulativeFitness += distances.cost(individual);
            if (cumulativeFitness > selectedPoint) {
                return individual;
            }
        }
    }

    int getBest(const std::vector<std::vector<int>>& population, mhe::DistanceMatrix distanceMatrix){
        int best = 0;
        double best_val = distanceMatrix.cost(population[0]);
        for(int i =0; i < population.size();i++){
            auto cur_val = distanceMatrix.cost(population[i]);
            if(best_val>cur_val){
                best = i;
                best_val = cur_val;
            }
        }
        return best;
    }


    std::vector<int> solve(int population_size, mhe::DistanceMatrix distances,int threshold,int max_iterations,const std::vector<std::string>& flags, bool debug=false){
        std::vector<int> (*selectionFunction)(const std::vector<std::vector<int>>&, mhe::DistanceMatrix);
        std::vector<int> (*mutationFunction)(std::vector<int>);
        std::vector<int> (*crossoverFunction)(std::vector<int>,std::vector<int>);

        if(max_iterations==-1) max_iterations = threshold*20;

        if(flags[0] == "tournament") selectionFunction = tournamentSelection;
        else if(flags[0]=="wheel") selectionFunction = rouletteWheelSelection;
        else throw std::invalid_argument("Blad we flagach algorytmu genetycznego. Podaj poprawyny rodzaj funkcji selekcji.");

        if(flags[1]=="inversion") mutationFunction = inversion;
        else if(flags[1]=="rotation") mutationFunction = rotation;
        else throw std::invalid_argument("Blad we flagach algorytmu genetycznego. Podaj poprawyny rodzaj funkcji mutacji.");

        if(flags[2]=="order") crossoverFunction = orderCrossover;
        else if(flags[2]=="uniform") crossoverFunction = uniformCrossover;
        else throw std::invalid_argument("Blad we flagach algorytmu genetycznego. Podaj poprawyny rodzaj funkcji krzyzowania.");

        std::vector<std::vector<int>> population;
        std::vector<std::vector<int>> new_generation;
        for(int i =0;i < population_size; i++) population.push_back(mhe::generateRandomSolution(distances.size()));
        int timer = 0;
        double temp_diff = 0;
        int counter=0;
        double diff;
        while(timer < max_iterations){
            new_generation.push_back(population[getBest(population,distances)]);
            for(int i = 1; i <population_size; i++){
                auto parent1 = selectionFunction(population, distances);
                auto parent2 = selectionFunction(population, distances);
                auto child = mutationFunction(crossoverFunction(parent1, parent2));
                if(distances.cost(parent2) < distances.cost(child)) new_generation.push_back(parent2);
                else if(distances.cost(parent1) < distances.cost(child))  new_generation.push_back(parent1);
                else new_generation.push_back(child);
            }
            population.clear();
            population = new_generation;
            new_generation.clear();
            timer++;
            diff = distances.cost(population[getBest(population,distances)]);
            if(temp_diff==diff) counter++;
            else counter = 0;

            if(timer%100==0&&debug) std::cout << timer << ": " << distances.cost(population[getBest(population,distances)]) << "\n";
            temp_diff = diff;
            if(counter >= threshold) break;
        }
        std::cout << "Iterations: " << timer <<  std::endl;
        return population[getBest(population,distances)];
    }

    void show_all(mhe::DistanceMatrix distanceMatrix){
        auto x = {"wheel","tournament"};
        auto y = {"inversion", "rotation"};
        auto z = {"uniform","order"};
        for(auto i : x){
            for(auto j : y){
                for(auto k : z){
                    std::cout << "Ga: " << i << " " << j << " "<< k << std::endl;
                    auto start = std::chrono::high_resolution_clock::now();
                    auto best = gen::solve(40, distanceMatrix, 10000, -1 ,{i,j,k}); //tournament, wheel; rotation, invertion; order, mapped
                    auto end = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> duration = end-start;
                    mhe::drawSolution(best);
                    std::cout << "with cost :" << distanceMatrix.cost(best) << " in " << duration.count() << std::endl;
                }
            }
        }
    }


}

auto main() -> int {
    using namespace std;
    const unsigned int num_cities = 50;
    auto dis = mhe::DistanceMatrix(mhe::generateRandomPoints(num_cities));
    gen::show_all(dis);
    return 0;
}