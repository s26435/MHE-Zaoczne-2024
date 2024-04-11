## Problem komiwojażera
### Opis problemu:
Dane jest n miast, które komiwojażer ma odwiedzić, oraz odległość pomiędzy każdą parą miast.
Celem jest znalezienie najkrótszej drogi łączącej wszystkie miasta. <br>
<a href="https://en.wikipedia.org/wiki/Travelling_salesman_problem">Wikipedia</a>

### Implementacja problemu
1. Implementacja problemu znajduje się w przestrzeni nazw ```mhe```
2. Do przechowywania odległości między miastami zastosowano klasę ```DistanceMatrix``` przyjmująca
w konstruktorze wektor punktów (```struct Point```). <br> 
Zawiera ona jeszcze dwie inne przydatne funkcje:
   *  ```double cost(const std::vector<int>& solution)``` - zwracającą odległości między miastami dla podanego rozwiązania
   * ```int size()``` - zwracającą ilość miast
3. Funkcje dodatkowe:
   * ```std::vector<Point> generateRandomPoints(unsigned int num_points, int gen_seed = SEED)``` - funkcja generująca losowe punkty. Używa stałej ```SEED``` 
    zdefiniowanej ```#define SEED 42```.
   * ```void drawSolution(const std::vector<int>& solution)``` - funkcja wypisująca w konsoli podane rozwiązanie.
   * ```std::vector<std::vector<int>> getNeighbors(const std::vector<int>& vector)``` - funkcja zwracająca __wszystkich__ sąsiadów podanego rozwiązania.
   * ```std::vector<int> generateRandomSolution(int size, int gen_seed=SEED)``` - funkcja zwracająca __losowego__ sąsiada podanego rozwiązania.

### Proste algorytmy:
1. Algorytm pełnego przeglądu
   * Służy do obliczania deterministycznego rozwiązania dla małej ilości rozwiązań.
   * przestrzeń nazw ```brute```
   * nagłówek funkcji rozwiązującej ```std::vector<std::vector<int>> generatePermutations(std::vector<int>& nums, int start)```
   * Reszta używanych funkcji:
      * ```std::vector<std::vector<int>> generatePermutations(std::vector<int>& nums, int start)``` - rekurencyjnie przeszukuje przestrzeń rozwiązań wyszukując wszystkie rozwiązania.
     Ze zbyt dużą ilością punktów może przekroczyć stos wywołań funkcji (dla 13 puntów powitał mnie niebieski uśmieszek <a href="https://filestore.community.support.microsoft.com/api/images/ca756db3-5559-4d35-8e50-6821a83060c0?upload=true">:) </a>)
      * ```std::vector<std::vector<int>> generateAllPermutations(int x)``` - funkcja wywołująca ```generatePermutations()```.
2. Algorytm wspinaczkowy:
   * przestrzeń nazw ```climb```
   * wersja algorytmu z wyborem najlepszego sąsiada:
     <br> funkcja rozwiązująca: ```std::vector<int> solve(mhe::DistanceMatrix distances, int loop_breaker = 1000)```
     przyjmuje argument ```loop_breaker``` określający maksymalną ilość iteracji głównej pętli algorytmu.
   * wersja algorytmu z wyborem losowego sąsiada: <br> funkcja rozwiązująca: ```std::vector<int> solveRandom(mhe::DistanceMatrix distances, int  max_iterations = 10000)```
     także przyjmuje argument ```loop_breaker``` określający maksymalną ilość iteracji głównej pętli algorytmu.
   * Funkcje użwyane przez algorytm:
       * ```std::vector<int> getBestNeighbor(const std::vector<int>& solution, mhe::DistanceMatrix distances)``` zwracająca najlepszego sąsiada rozwiązania lub wektor zerowy, jeśli rozwiązanie jest najlepsze w swoim sąsiedztwie.
       * ```std::vector<int> getRandomNeighbor(std::vector<int> solution)``` zwraca losowego sąsiada z sąsiedztwa podanego rozwiązania.
3. Algorytm Tabu:
    * przestrzeń nazw ```tabu```
    * funkcja rozwiązująca: ```std::vector<int> solve(mhe::DistanceMatrix distances, int max_iter=1000, int max_tab_size=5)```<br>
     argument ```max_iter``` przyjmuje maksymalną ilość iteracji głównej pętli algorytmu, a argument ```max_tab_size``` to wielkość tablicy pomocniczej przechowującej trasę algorytmu.
4. Algorytm symulowanego wyżarzania:
    * przestrzeń nazw ```ann```
    * funkcja rozwiązująca: ```std::vector<int> solve(mhe::DistanceMatrix distances, double initial_temperature=25, int max_iterations = 1000)``` <br>
    ```initial_temperature``` - określa początkową temperaturę układu <br>
    ```max_iterations``` - określa maksymalną ilość iteracji głównej pętli algorytmu
### Algorytm Genetyczny:
1. Przestrzeń nazw ```gen```
2. Główna funkcja algorytmu ```std::vector<int> solve(int population_size, mhe::DistanceMatrix distances,int threshold,int max_iterations,const std::vector<std::string>& flags, bool debug=false)```
   * ```population_size``` - określa wielkość populacji w algorytmie
   * ```distances``` - macierz dystansów problemu
   * ```threshold``` - określa, ile może nastąpić epok bez zmiany najlepszego rozwiązania, potem przerywa główną pętle algorytmu (pierwszy warunek stopu)
   * ```max_iterations``` - maksymalna ilość epok, zostaw "-1" aby przypisać domyślną ilość epok (threshold * 20) (drugi warunek stopu)
   * ```flags``` - wektor flag typu string ```{funkcja selekcji, funkcja mutacji, funkcja krzyżowania}```
     * funkcje selekcji:
       * "tournament" - selekcja turniejowa
       * "wheel" - selekcja kołem ruletki
     * funkcje mutacji:
       * "inversion" - mutacja przez inwersję 
       * "rotation" - mutacja przez rotację
     * funkcje krzyżowania:
       * "order" - order crossover OX1
       * "uniform" - uniform crossover
   * ```debug``` - funkcja w trybie debug - wypisuje wszystkich osobników populacji niektórych epok
3. Funkcje krzyżowania - przyjmują jako argumenty dwa rozwiązania```std::vector<int> parent1, std::vector<int> parent2``` i wykorzystują je jako rodziców
do zwracanego rozwiązania.
4. Funkcje mutacji przyjmują rozwiązanie ```std::vector<int> solution``` i zwracają jego zmutowaną wersję
5. Funkcje selekcji przyjmują populację ```const std::vector<std::vector<int>>& population``` i macierz dystansów``` mhe::DistanceMatrix distances ``` i zwracają jedno rozwiązanie
6. Funkcja ```int getBest(const std::vector<std::vector<int>>& population, mhe::DistanceMatrix distanceMatrix)``` zwraca indeks najlepszego osobnika w populacji.
7. Funkcja ```void show_all(mhe::DistanceMatrix distanceMatrix)``` porównuje działanie różnych algorytmów i przedstawia ich rozwiązania.




