#include<iostream>
#include<cmath>
#include<cstdlib>
#include<algorithm>
#include<sstream>
#include<cstring>
#include<map>
#include<list>
#include<queue>
#include<ctime>
#include<vector>
#include<fstream>
#include<chrono>
#define MAXGENS 100000000		//program restarts if no solution has been found after MAXGENS
#define MUTATION_RATE 0.8	//Rate of mutation (here 80%)
#define MATING_THRESHOLD 0.5	//Threshold for choosing Mating partners (here partners are chosen from top 50%)
#define POPULATIONSIZE 100	//Size of Population in each generation
#define SQUARESIZE 3		//Size of Default Square if not specified with --size
#define CROSSOVER_LIMIT_A 0.6	//Lower limit for Gene-Index in Crossover Function
#define CROSSOVER_LIMIT_B 0.9	//Upper limit for Gene-Index in Crossover Function
#define ELITE_MAX_AGE 5000	//Maximum Age of Elites
#define ELITE_PERCENTAGE 15	//Percantage of Individuals to be kept as elites until better ones are found (here 15%)
#define CHECK_DUPLICATES true	//If true forbid children to be duplicates of this or last Generations (Turn off for larger Population Sizes)

using namespace std;

typedef struct
{
    vector<vector<int>> board;
    int cost;
    int age;
} individual;

typedef struct
{
    individual* child1;
    individual* child2;
} children;


typedef vector<individual*> population_type;
population_type population;
int squareSize;

void print_square(vector<vector<int>> square) {
    std::cout << endl ;
    std::cout << "----------------------" << endl;
    for (int i = 0;i < squareSize;i++) {
        for (int j = 0;j < squareSize;j++) {
            std::cout << " " << square[i][j] << " |";
        }
        std::cout << endl;
        std::cout << "----------------------" << endl;
    }
}


int fitnessValue(vector<vector<int>> board,bool print=false) // ultra magic square
{
    int fitness = 0;
    int sum_rows = 0;
    int sum_columns = 0;
    int magicNumber = squareSize * (squareSize * squareSize + 1) / 2;
    int magicAssocNumber = (squareSize * squareSize) + 1;
    int sum_diagonals = 0;
    int sum_diagonals_reverse = 0;
    int sum_opposite_numbers = 0;
    for (int i = 0;i < squareSize;i++) {
        sum_diagonals = 0;
        sum_diagonals_reverse = 0;
        for (int j = 0;j < squareSize;j++) {
            sum_rows += board[i][j];
            sum_columns += board[j][i];
            //summing up all diagonals for pandiagonal MS
            sum_diagonals += board[(j + i + squareSize) % squareSize][j];
            sum_diagonals_reverse += board[(2 * squareSize - j - i - 1) % squareSize][j];
            //summing opposites for associative MS
            sum_opposite_numbers = board[i][j] + board[squareSize - 1 - i][squareSize - 1 - j];
            if (print) cout << "Opposite @ " << i << "," << j << " : " << abs(sum_opposite_numbers - magicAssocNumber) << endl;
            //Divide by 2 for duplicates
            if (sum_opposite_numbers != magicAssocNumber) fitness += ceil(float(abs(sum_opposite_numbers - magicAssocNumber))/1);
            sum_opposite_numbers = 0;
        }
        if (print) cout <<"Diagonale starting "<< (0 + i + squareSize) % squareSize <<" : " << abs(sum_diagonals - magicNumber) <<endl;
        if (print) cout << "Reverse Diagonale starting " << (2 * squareSize - 0 - i - 1) % squareSize << " : " << abs(sum_diagonals_reverse - magicNumber) <<endl;
        if (sum_diagonals != magicNumber) fitness += abs(sum_diagonals - magicNumber);
        if (sum_diagonals_reverse != magicNumber) fitness += abs(sum_diagonals_reverse - magicNumber);
        if (print)cout << "Zeile starting " << i <<" : "<< abs(sum_rows - magicNumber) << endl;
        if (print)cout << "Spalte starting " << i << " : " << abs(sum_columns - magicNumber) << endl;
        if (sum_rows != magicNumber) fitness +=abs(sum_rows-magicNumber);
        if (sum_columns != magicNumber) fitness +=abs(sum_columns-magicNumber);
        sum_rows = 0;
        sum_columns = 0;
    }
    return fitness;
}

individual* createNode()
{
    individual* newNode = new individual;
    vector<vector<int>> initBoard(squareSize, vector<int>(squareSize, 0));
    newNode->board = initBoard;
    newNode->age = 0;
    return newNode;
}

vector<vector<int>> create_shuffeled_board(vector<int> samples) {
    vector<vector<int>> sampleBoard(squareSize, vector<int>(squareSize, 0));
    for (int i = 0; i < squareSize * squareSize; i++)
    {
        sampleBoard[int(i / squareSize)][i % squareSize] = samples[i];
    }

    return sampleBoard;
}

void generatePopulation()
{
    population.clear();
    vector<vector<int>> sampleBoard(squareSize,vector<int>(squareSize,0));
    vector<int> samples(squareSize * squareSize, 0);
    individual* temp;
    for (int i = 0; i < squareSize*squareSize; i++)
    {
        samples[i] = (i + 1);
    }

    //adds entries to population list
    for (int i = 0; i < POPULATIONSIZE; i++)
    {
        random_shuffle(samples.begin(), samples.end());
        sampleBoard = create_shuffeled_board(samples);
        temp = createNode();
        temp->board = sampleBoard;
        temp->cost = fitnessValue(sampleBoard);
        population.push_back(temp);
    }
}

individual* reproduce_order1(individual* x, individual* y,int p1,int p2)
{
    individual* child = createNode();
    int n = squareSize*squareSize;
    //Order 1 Crossover
    //lower bound 0.6 upper bound 0.9 -> n*0.3
    int pos1 = min(p1, p2);
    int pos2 = max(p1, p2);
    vector<int>parent1_set(pos2-pos1,0);
    int counter = 0;
    for (int i = pos1;i < pos2; i++) {
        child->board[int(i / squareSize)][i % squareSize] = x->board[int(i / squareSize)][i % squareSize];
        parent1_set[counter] = x->board[int(i / squareSize)][i % squareSize];
        counter++;
    }
    int correct_i = 0;
    int child_pos_counter = pos2;
    int corrected_child_pos_counter = pos2;
    for (int i = pos2;i < (pos2 + n);i++) {
        correct_i = i % n;
        int insert_value = y->board[int(correct_i / squareSize)][correct_i % squareSize];
        bool used_flag = false;
        for (int j = 0;j < (pos2 - pos1);j++) {
            if (parent1_set[j] == insert_value) {
                used_flag = true;
                break;
            }
        }
        if (!used_flag) {
            corrected_child_pos_counter = child_pos_counter % (squareSize * squareSize);
            child->board[int(corrected_child_pos_counter / squareSize)][corrected_child_pos_counter % squareSize] = insert_value;
            child_pos_counter++;
        }
    }
    //cout << (pos2 - pos1) << " , " << (child_pos_counter - pos2) << endl;
    child->cost = fitnessValue(child->board);
    return child;
}

vector<int> board_inversion(individual* x) {
    int mat_size = squareSize * squareSize;
    vector<int> inv(mat_size, -1);
    for (int i = 0;i < mat_size;i++) {
        int bigger_counter = 0;
        for (int j = 0;j < mat_size;j++) {
            int val = x->board[int(j / squareSize)][j % squareSize];
            if (val == (i+1)) {
                inv[i] = bigger_counter;
                break;
            }
            else if (val > (i+1)) bigger_counter++;
        }
    }
    return inv;
}

vector<vector<int>> reverse_inversion(vector<int> inv) {
    int mat_size = squareSize * squareSize;
    vector<int> pos(mat_size, -1);
    vector<vector<int>> solution(squareSize,vector<int>(squareSize, -1));

    for (int i = (mat_size - 1);i >= 0;i--) {
        int bigger_counter = 0;
        int inv_i = inv[i];
        for (int j = (mat_size - 1);j >= 0;j--) {
            if (i == j) {
                pos[i] = inv_i;
                break;
            }
            else if (pos[j]>= inv_i){
                pos[j]++;
            }
        }
    }
    for (int i = 0;i < mat_size;i++) {
        int p = pos[i];
        solution[int(p/squareSize)][p%squareSize] = (i + 1);
    }
    return solution;
}

children reproduce_inversion(individual* x, individual* y) {
    vector<int> inv_x = board_inversion(x);
    vector<int> inv_y = board_inversion(y);
    int mat_size = squareSize * squareSize;
    int split_point = rand() % (mat_size);
    while (split_point<int(CROSSOVER_LIMIT_A * mat_size) || split_point>int(CROSSOVER_LIMIT_B * mat_size)) {
        split_point = rand() % (mat_size);
    }
    vector<int> inv_child1(mat_size,-1);
    vector<int> inv_child2(mat_size, -1);
    for (int i = 0;i < split_point;i++) {
        inv_child1[i] = inv_x[i];
        inv_child2[i] = inv_y[i];
    }
    for (int i = split_point;i < mat_size;i++) {
        inv_child1[i] = inv_y[i];
        inv_child2[i] = inv_x[i];
    }

    individual* child1 = createNode();
    individual* child2 = createNode();

    child1->board = reverse_inversion(inv_child1);
    child2->board = reverse_inversion(inv_child2);
    child1->cost = -1;
    child2->cost = -1;

    children ret_children;
    ret_children.child1 = child1;
    ret_children.child2 = child2;


    return ret_children;
}

void mutate(individual* child)
{

    int rand1 = rand() % (squareSize * squareSize);
    int rand2 = rand() % (squareSize * squareSize);
    int temp = 0;
    //swapping two random fields for mutation
    temp = child->board[int(rand1 / squareSize)][rand1 % squareSize];
    child->board[int(rand1 / squareSize)][rand1 % squareSize] = child->board[int(rand2 / squareSize)][rand2 % squareSize];
    child->board[int(rand2 / squareSize)][rand2 % squareSize] = temp;
    //return child;
}

int randomSelection()
{
    int randomPos = rand() % int(population.size()*MATING_THRESHOLD);
    return randomPos;
}

bool isFit(individual* test)
{
    if (test->cost == 0) 
        return true;
    return false;
}

bool in_pop(individual* ind1, individual* ind2, vector<individual*>pop) {
#pragma parallel for shared(found)
    for (int i = 0;i < pop.size();i++) {
        if (pop[i]->board == ind1->board) return true;
        else if (pop[i]->board == ind2->board) return true;
    }
    return false;;
}

bool comp(individual* a, individual* b)
{
    //sorting by lowest fitness value
    return(a->cost < b->cost);
}

void diffboard_printer(vector<vector<int>> old, vector<vector<int>> newer) {
    for (int i = 0;i < squareSize;i++) {
        for (int j = 0;j < squareSize;j++) {
            int entry_new = newer[i][j];
            int entry_old = old[i][j];
            if (entry_new == entry_old)printf(" %i |", entry_new);
            else printf(" \x1B[32m%i\033[0m |", entry_new);
        }
        cout << endl;
        cout << "--------------------------------" << endl;
    }
    //printf("\x1B[32mTexting\033[0m\t\t");
}

individual* GA()
{
    clock_t start_time, end_time;
    bool found = false;
    int gen_counter = 0;
    int prev_best_fitness = population[0]->cost;
    individual* prev_best_ind = population[0];
    int n = squareSize * squareSize;
    int best_fitness = 0;
    int best_fitness_counter = 0;
    start_time = clock();
    vector<vector<int>> last_board = population[0]->board;
    //while (!found&&(gen_counter<MAXGENS))
    while (!found)
    {
        gen_counter++;
        sort(population.begin(), population.end(), comp);
        if (population[0]->cost != best_fitness) {
            best_fitness = population[0]->cost;
            std::cout << "Best fitness: " << population[0]->cost << " Generation: " << gen_counter << endl;
        }
        if (gen_counter % 1000 == 0) {
            cout << "Generation: " << gen_counter << endl;
            diffboard_printer(last_board, population[0]->board);
            cout << "Fitness value of best individual: " << population[0]->cost << endl;
            end_time = clock();
            cout << "Time required for last 1000 Generations: \t" << 1000 * ((double)(end_time - start_time) / CLOCKS_PER_SEC) << " milliseconds." << "\n\n";
            start_time = clock();
            last_board = population[0]->board;
        }
        population_type new_population;
        
        //keep best ELITE_PERCENTAGE%
        for (int i = 0;i<int(population.size()*ELITE_PERCENTAGE/100);i++) {
            if (population[i]->age <= ELITE_MAX_AGE) {
                individual* temp = createNode();
                temp->board = population[i]->board;
                temp->cost = population[i]->cost;
                temp->age = population[i]->age + 1;
                new_population.push_back(temp);
            }
        }
        int missing_elements = population.size() - new_population.size();
        
        while(new_population.size()<population.size())
        {
            int randomNum1, randomNum2;
            individual* child1, * child2;
            individual* individualX, * individualY;
            randomNum1 = randomSelection();
            individualX = population[randomNum1];

            randomNum2 = randomSelection();
            individualY = population[randomNum2];
            children ch = reproduce_inversion(individualX, individualY);
            child1 = ch.child1;
            child2 = ch.child2;
            if (rand() % 100 < (MUTATION_RATE * 100))     
                mutate(child1);
            if (rand() % 100 < (MUTATION_RATE * 100))     
                mutate(child2);
            //prevent duplicates
            if (CHECK_DUPLICATES) {
                while (in_pop(child1, child2, new_population) || in_pop(child1, child2, population))
                {
                    delete child1;
                    delete child2;
                    randomNum1 = randomSelection();
                    individualX = population[randomNum1];
                    randomNum2 = randomSelection();
                    individualY = population[randomNum2];
                    children ch = reproduce_inversion(individualX, individualY);
                    child1 = ch.child1;
                    child2 = ch.child2;


                    if (rand() % 100 < (MUTATION_RATE * 100))     //random probability 50% chance too high
                        mutate(child1);
                    if (rand() % 100 < (MUTATION_RATE * 100))     //random probability 50% chance too high
                        mutate(child2);
                }
            }
            

            child1->cost = fitnessValue(child1->board);
            child2->cost = fitnessValue(child2->board);
            
            new_population.push_back(child1);
            if ((population.size() - new_population.size()) > 1)new_population.push_back(child2);
            else delete child2;
        }
        //freeing memory to avoid memory leak
        for (int i = 0; i < population.size();i++) {
            delete population[i];
        }
        //population.clear();
        population = new_population;
        for (int i = 0; i < population.size();i++) {
            if (isFit(population[i])) {
                return population[i];
            }
        }
        //temp.clear();
        new_population.clear();
    }
    individual* wrong = createNode();
    wrong->cost = -1;
    return wrong;
}

individual* GA_parallel()//tried parallelizing GA function
{
    clock_t start_time, end_time;
    bool found = false;
    int gen_counter = 0;
    int prev_best_fitness = population[0]->cost;
    individual* prev_best_ind = population[0];
    int n = squareSize * squareSize;
    int best_fitness = 0;
    int best_fitness_counter = 0;
    start_time = clock();
    vector<vector<int>> last_board = population[0]->board;
    while (!found && (gen_counter < MAXGENS))
    {
        gen_counter++;
        sort(population.begin(), population.end(), comp);
        if (population[0]->cost != best_fitness) {
            best_fitness = population[0]->cost;
            std::cout << "Best fitness: " << population[0]->cost << " Generation: " << gen_counter << endl;
        }
        if (gen_counter % 10 == 0) {
            cout << "Generation: " << gen_counter << endl;
            diffboard_printer(last_board,population[0]->board);
            cout << population[0]->cost << endl;
            cout << fitnessValue(population[0]->board) << endl;
            end_time = clock();
            cout << (end_time - start_time) << endl;
            start_time = clock();
            last_board = population[0]->board;
        }
        population_type new_population;

        //keep best ELITE_PERCENTAGE%
        //#pragma omp parallel for shared(new_population)
        for (int i = 0;i<int(population.size() * ELITE_PERCENTAGE / 100);i++) {
            if (population[i]->age <= ELITE_MAX_AGE) {
                individual* temp = createNode();
                temp->board = population[i]->board;
                temp->cost = population[i]->cost;
                temp->age = population[i]->age + 1;
                new_population.push_back(temp);
            }
        }
        int missing_elements = population.size() - new_population.size();
        population_type remaining_pop(missing_elements, NULL);
        for (int i =0; i<missing_elements;i++)
        {
            int randomNum1, randomNum2;
            individual* child1, * child2;
            individual* individualX, * individualY;
            randomNum1 = randomSelection();
            individualX = population[randomNum1];

            randomNum2 = randomSelection();
            individualY = population[randomNum2];
            children ch = reproduce_inversion(individualX, individualY);
            child1 = ch.child1;
            child2 = ch.child2;
            if (rand() % 100 < (MUTATION_RATE * 100))
                mutate(child1);
            if (rand() % 100 < (MUTATION_RATE * 100))
                mutate(child2);
            //prevent duplicates
            /*
            while (in_pop(child1, new_population) ||
                in_pop(child2, new_population) ||
                in_pop(child1, population) ||
                in_pop(child2, population))
            {
                delete child1;
                delete child2;
                randomNum1 = randomSelection();
                individualX = population[randomNum1];
                randomNum2 = randomSelection();
                individualY = population[randomNum2];
                children ch = reproduce_inversion(individualX, individualY);
                child1 = ch.child1;
                child2 = ch.child2;


                if (rand() % 100 < (MUTATION_RATE * 100))     //random probability 50% chance too high
                    mutate(child1);
                if (rand() % 100 < (MUTATION_RATE * 100))     //random probability 50% chance too high
                    mutate(child2);
            }
            */

            child1->cost = fitnessValue(child1->board);
            child2->cost = fitnessValue(child2->board);
            
            remaining_pop[i] = child1;

        }
        //freeing memory to avoid memory leak
        new_population.insert(new_population.end(), remaining_pop.begin(), remaining_pop.end());

        for (int i = 0; i < population.size();i++) {
            delete population[i];
        }
        population.clear();
        population = new_population;
        for (int i = 0; i < population.size();i++) {
            if (isFit(population[i])) {
                return population[i];
            }
        }
        //temp.clear();
        new_population.clear();
    }
    return prev_best_ind;
}

void write_square_to_file(string outfile, vector<vector<int>> board) {
    ofstream myfile;
    myfile.open(outfile);
    for (int i = 0;i < squareSize;i++) {
        for (int j = 0;j < squareSize;j++) {
            if (j != (squareSize - 1)) myfile << board[i][j] << ";";
            else myfile << board[i][j];
        }
        myfile << "\n";
    }
}

void initialize()
{
    srand(time(0));     //to ensure perfect randomness
    squareSize = SQUARESIZE;
}

void read_params_from_cmdline(int argc, char** argv) {
    vector<string> args(argv + 1, argv + argc);
    string infname, outfname, generations;
    for (auto i = args.begin(); i != args.end(); ++i) {

        if (*i == "--size") {
            if ((i + 1) != args.end()) {
                generations = *++i;
                squareSize = atoi(generations.c_str());
            };
        }
    }
};

int main(int argc, char** argv)
{
    initialize();
    read_params_from_cmdline(argc, argv);
    clock_t start_time, end_time;           //to keep a track of the time spent in computing
    map<vector<vector<int>>, int> solutionsFound;
    int maxSolutions = 1, numFound = 0;       
    start_time = clock();
    int sol_cost = -1;
    while (numFound != maxSolutions)
    {
        generatePopulation();
        sol_cost = -1;
        individual* solution = GA();
        sol_cost = solution->cost;
        while(sol_cost==-1) {
            generatePopulation();
            individual* solution = GA();
            sol_cost = solution->cost;
        }
        if (solution->board[0][0] == 0) continue;
        if (!solutionsFound[solution->board])
        {
            solutionsFound[solution->board] = 1;
            //cout << "Possible Solution #" << (++numFound) << ":\t" << solution->board << endl;
            long now = clock();
            string filename = "Result_square_";
            filename += to_string(squareSize);
            filename += "x";
            filename += to_string(squareSize);
            filename += "_";
            filename += to_string(now);
            filename += ".csv";
            numFound++;
            write_square_to_file(filename, solution->board);
            std::cout << "Found " << numFound << " unique Solutions" << endl;


        }
        cout << "Restarting Algorithm" << endl;
    }
 
    for (auto it = solutionsFound.begin(); it != solutionsFound.end(); it++)
    {
        print_square(it->first);    
    }
    end_time = clock();

    std::cout << "Time required for execution: \t" << 1000 * ((double)(end_time - start_time) / CLOCKS_PER_SEC) << " milliseconds." << "\n\n";
    return 0;
}
