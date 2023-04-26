/**
 * Projec : gtsp (voyageur de commerce)
 *
 * Date   : 07/04/2014
 * Author : Olivier Grunder
 */
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NBR_TOWNS 6

/* Distance matrix */
double dist[NBR_TOWNS][NBR_TOWNS];

/* Each edge has a starting and ending node */
int starting_town[NBR_TOWNS];
int ending_town[NBR_TOWNS];

/* no comment */
int best_solution[NBR_TOWNS];
double best_eval = -1.0;
int nb_nodes = 0;

/**
 * Berlin52 :
 *  6 towns : Best solution (2315.15): 0 1 2 3 5 4
 * 10 towns : Best solution (2826.50): 0 1 6 2 7 8 9 3 5 4
 */
float coord[NBR_TOWNS][2] =
    {
        {565.0, 575.0},
        {25.0, 185.0},
        {345.0, 750.0},
        {945.0, 685.0},
        {845.0, 655.0},
        {880.0, 660.0},
};

double distformula(float x1, float y1, float x2, float y2)
{
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

// compute the distance matrix from the coord array and returns a (n_coords, n_coords) matrix

void compute_dist_matrix()
{
    for (int i = 0; i < NBR_TOWNS; i++)
    {
        for (int j = 0; j < NBR_TOWNS; j++)
        {
            dist[i][j] = distformula(coord[i][0], coord[i][1], coord[j][0], coord[j][1]);
        }
    }
}
/*
 * print a matrix
 */
void print_matrix(double d[NBR_TOWNS][NBR_TOWNS])
{
    for (int i = 0; i < NBR_TOWNS; i++)
    {
        printf("%d:", i);
        for (int j = 0; j < NBR_TOWNS; j++)
        {
            printf("%6.1f ", d[i][j]);
        }
        printf("\n");
    }
}

/**
 * print a solution
 */
void print_solution(int *sol, double eval)
{
    int i;
    printf("(%.2f): ", eval);
    for (i = 0; i < NBR_TOWNS; i++)
        printf("%d ", sol[i]);
    printf("\n");
}

/**
 * evaluation of a solution
 */
double evaluation_solution(int *sol)
{
    double eval = 0;
    int i;
    for (i = 0; i < NBR_TOWNS - 1; i++)
    {
        eval += dist[sol[i]][sol[i + 1]];
    }
    eval += dist[sol[NBR_TOWNS - 1]][sol[0]];

    return eval;
}

/**
 * nearest neighbour solution
 */
double build_nearest_neighbour()
{
    /* solution of the nearest neighbour */
    int i, sol[NBR_TOWNS];

    /* evaluation of the solution */
    double eval = 0;

    // init solution
    for (i = 0; i < NBR_TOWNS; i++)
    {
        sol[i] = 0;
    }

    // find solution w/ nearest neighboor
    for (i = 1; i < NBR_TOWNS; i++)
    {
        double min = -1;
        int k;
        for (int j = 0; j < NBR_TOWNS; j++)
        {
            for (k = 0; k < i; k++)
            {
                if (sol[k] == j)
                {
                    break;
                }
            }
            if (k < i)
            {
                continue;
            }
            if (min < 0 || dist[sol[i - 1]][j] < min)
            {
                min = dist[sol[i - 1]][j];
                sol[i] = j;
            }
        }
    }

    eval = evaluation_solution(sol);
    printf("Nearest neighbour ");
    print_solution(sol, eval);

    for (i = 0; i < NBR_TOWNS; i++)
        best_solution[i] = sol[i];
    best_eval = eval;

    return eval;
}

/**
 *  Build final solution
 */
void build_solution()
{
    int i, solution[NBR_TOWNS];

    int indiceCour = 0;
    int villeCour = 0;

    while (indiceCour < NBR_TOWNS)
    {

        solution[indiceCour] = villeCour;

        // Test si le cycle est hamiltonien
        for (i = 0; i < indiceCour; i++)
        {
            if (solution[i] == villeCour)
            {
                /* printf ("cycle non hamiltonien\n") ; */
                return;
            }
        }
        // Recherche de la ville suivante
        int trouve = 0;
        int i = 0;
        while ((!trouve) && (i < NBR_TOWNS))
        {
            if (starting_town[i] == villeCour)
            {
                trouve = 1;
                villeCour = ending_town[i];
            }
            i++;
        }
        indiceCour++;
    }

    double eval = evaluation_solution(solution);

    if (best_eval < 0 || eval < best_eval)
    {
        best_eval = eval;
        for (i = 0; i < NBR_TOWNS; i++)
            best_solution[i] = solution[i];
        printf("New best solution: ");
        print_solution(solution, best_eval);
    }
    return;
}

/**
 *  Little Algorithm
 */
void little_algorithm(double d0[NBR_TOWNS][NBR_TOWNS], int iteration, double eval_node_parent)
{
    nb_nodes++;
    if (iteration == NBR_TOWNS)
    {
        build_solution();
        return;
    }

    /* Do the modification on a copy of the distance matrix */
    double d[NBR_TOWNS][NBR_TOWNS];
    memcpy(d, d0, NBR_TOWNS * NBR_TOWNS * sizeof(double));

    int i, j;

    double eval_node_child = eval_node_parent;

    /*
     * substract the min of the rows and the min of the columns
     * and update the evaluation of the current node
     */

    for (i = 0; i < NBR_TOWNS; i++)
    {
        double min = INFINITY;
        for (j = 0; j < NBR_TOWNS; j++)
        {
            if (d[i][j] < min)
            {
                if (d[i][j] != -1)
                    min = d[i][j];
            }
        }
        for (j = 0; j < NBR_TOWNS; j++)
        {
            if (d[i][j] != -1)
                d[i][j] -= min;
        }
        if (min != INFINITY)
            eval_node_child += min;
    }

    // Check if 0 in every column

    for (j = 0; j < NBR_TOWNS; j++)
    {
        double min = INFINITY;
        for (i = 0; i < NBR_TOWNS; i++)
        {
            if (d[i][j] < min)
            {
                if (d[i][j] != -1)
                    min = d[i][j];
            }
        }
        for (i = 0; i < NBR_TOWNS; i++)
        {
            if (d[i][j] != -1)
                d[i][j] -= min;
        }
        if (min != INFINITY)
            eval_node_child += min;
    }

    /* Cut : stop the exploration of this node */
    if (best_eval >= 0 && eval_node_child >= best_eval)
        return;

    /**
     *  Compute the penalities to identify the zero with max penalty
     *  If no zero in the matrix, then return, solution infeasible
     * row and column of the zero with the max penalty */

    int izero = -1, jzero = -1;
    double max_penalty = -1;
    for (i = 0; i < NBR_TOWNS; i++)
    {
        for (j = 0; j < NBR_TOWNS; j++)
        {
            if (d[i][j] == 0)
            {
                double penalty_x = INFINITY, penalty_y = INFINITY;
                double visited_x = 0, visited_y = 0;
                for (int k = 0; k < NBR_TOWNS; k++)
                {
                    if (k != i && d[k][j] < penalty_x && d[k][j] >= 0)
                    {
                        penalty_x = d[k][j];
                        visited_x = 1;
                    }
                    if (k != j && d[i][k] < penalty_y && d[i][k] >= 0)
                    {
                        penalty_y = d[i][k];
                        visited_y = 1;
                    }
                }
                double penalty = penalty_x + penalty_y;

                if (penalty > max_penalty && (visited_x || visited_y))
                {
                    max_penalty = penalty;
                    izero = i;
                    jzero = j;
                }
                else if (visited_x == 0 && visited_y == 0)
                {
                    max_penalty = INFINITY;
                    izero = i;
                    jzero = j;
                }
            }
        }
    }

    if (izero == -1 || jzero == -1)
    {
        return;
    }

    starting_town[iteration] = izero;
    ending_town[iteration] = jzero;

    /* Do the modification on a copy of the distance matrix */
    double d2[NBR_TOWNS][NBR_TOWNS];
    memcpy(d2, d, NBR_TOWNS * NBR_TOWNS * sizeof(double));

    /* Modify the dist matrix to explore the choice of the zero with the max penalty */
    for (i = 0; i < NBR_TOWNS; i++)
    {
        d2[izero][i] = -1;
        d2[i][jzero] = -1;
    }

    /* Explore left child node according to given choice */
    little_algorithm(d2, iteration + 1, eval_node_child);

    /* Do the modification on a copy of the distance matrix */
    memcpy(d2, d, NBR_TOWNS * NBR_TOWNS * sizeof(double));

    d2[izero][jzero] = -1;

    /* Explore right child node according to non-choice */
    little_algorithm(d2, iteration, eval_node_child);
}

/**
 *
 */
int main(int argc, char *argv[])
{

    clock_t t_start, t_end;
    t_start = clock();

    best_eval = -1;

    /* Print problem informations */
    printf("Points coordinates:\n");
    int i;
    for (i = 0; i < NBR_TOWNS; i++)
    {
        printf("Node %d: x=%f, y=%f\n", i, coord[i][0], coord[i][1]);
    }
    printf("\n");

    /* Compute the distance matrix */
    // use the compute_dist_matrix function, pass him the corrds and coord num
    for (int i = 0; i < NBR_TOWNS; i++)
    {
        int j = 0;
        while (j < i)
        {
            // Compute the distance between two towns
            dist[i][j] = sqrt(pow(coord[i][0] - coord[j][0], 2) + pow(coord[i][1] - coord[j][1], 2));
            // Symmetric matrix
            dist[j][i] = dist[i][j];
            j++;
        }
        dist[i][i] = -1;
    }

    // compute_dist_matrix(NBR_TOWNS, coord);

    printf("Distance Matrix:\n");
    print_matrix(dist);
    printf("\n");

    double nearest_neighbour = build_nearest_neighbour();

    int iteration = 0;
    double lowerbound = 0.0;
    little_algorithm(dist, iteration, lowerbound);
    printf("Best solution:");
    print_solution(best_solution, best_eval);

    t_end = clock();
    double time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
    printf("Time taken: %f\n", time);
    printf("Total number of nodes: %d\n", nb_nodes);

    printf("Hit RETURN!\n");
    getchar();

    return 0;
}
