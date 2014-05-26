#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<mpi.h>
#include<math.h>
#include<limits.h>

/**
 * CONSTANT'S
 */
#define EMPTY 0
#define WOLF 1
#define SQUIRREL 2
#define TREE 3
#define ICE 4
#define SQUIRREL_TREE 5
#define INVALIDO -1

#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3
#define WORLD 4

#define DOWN_TO_UP_TAG 0
#define UP_TO_DOWN_TAG 1

/**
 * Processor initial position, Position processor and Processor elements number
 */
#define PROC_INIT_LINE(procID, nproc, dim) ((dim / nproc) * procID)
#define PROC_LINES_NUM(procID, nproc, dim) ((procID != nproc - 1) ? (dim / nproc) : (dim / nproc) + (dim % nproc))
#define PROC_END_LINE(procID, nproc, dim) (((dim / nproc) * procID) + (PROC_LINES_NUM(procID, nproc, dim)))


/**
 * GLOBAL VARIABLES
 */
int NUM_GENERATIONS;
int WOLF_BREEDING_PERIOD;
int SQUIRREL_BREEDING_PERIOD;
int WOLF_STARVATION_PERIOD;
int DIM;
int current_generation;

int me;
int nprocs;

int INIT_LINE;
int N_LINES;
int END_LINE;

/**
 * DATA STRUCT DEFINITION
 */
typedef struct world {
    int last_generation;
    int type;
    int breeding_period;
    int starvation_period;
} world_t;

typedef struct pos {
    int row;
    int column;
} pos_t;

/**
 * Grids definition
 */

world_t ** world; // world grid

world_t *** generation_grid; // generation grid (current generation)


/**
 * MAIN FUNCTIONS DECLARATION
 */
void simulateGeneration();
void simulateReds();
void simulateBlacks();
void processCell(int row, int column);
void processWolf(pos_t coordinates);
void processSquirrel(pos_t coordinates);
void updateAnimal(world_t animal, pos_t coordinates, pos_t newCoordinates);
void removeAnimal(pos_t coordinates);
pos_t moveWolf(pos_t coordinates);
int updateGenerationGrid(world_t element, pos_t pos, int direction);

pos_t moveSquirrel(pos_t coordinates);

void removeAnimal(pos_t coordinates);
void breedWolf(pos_t coordinates);
void breedSquirrel(pos_t coordinates);

int findNearbySquirrels(pos_t coordinates, pos_t *result);
int findNearbyEmptyCells(pos_t coordinates, int considerTrees, pos_t *result);
int computeOptimalCellChoice(pos_t currentPosition, int numberOfChoices);
/*
 * AUX FUNCTIONS DECLARATION
 */
void processInput(int argc, char ** argv);
void processOutput();
void initGrids();
void updateWorldGrid();
void initBorderPositions();
void send_border_updates();
void send_result_of_update();

int isBorder(int me, int ln);

// Debug functions
void verifyAlloc(void * memoryPointer, char * errorMSG);
void printReceiveResult(int me, MPI_Status * status);
void printSendResult(int me, MPI_Status * status);
void printLine(world_t * line);
void printMyGrid(int me);
void printMyGridDetailed(int me);

world_t getWorldElement(pos_t pos) {
    return world[pos.row][pos.column];
}

void setWorldElement(pos_t pos, world_t new) {
    world[pos.row][pos.column] = new;
}

void setWorldElementType(pos_t pos, int type) {
    world[pos.row][pos.column].type = type;
}

/**
 * Receive a element and add it to the list of element in a pos (add to generation grid)
 */
int updateGenerationGrid(world_t element, pos_t pos, int dir) {
    generation_grid[pos.row][pos.column][dir] = element;
}

world_t getGenerationElement(pos_t pos, int dir) {
    return generation_grid[pos.row][pos.column][dir];
}

void send_border_updates() {

    world_t upBorderToSend[DIM];
    world_t downBorderToSend[DIM];

    int upBorderLine = INIT_LINE - 1;
    int downBorderLine = END_LINE;

    int i;

    MPI_Request requests[2];
    MPI_Status statuses[2];

    //send up border
    if (me != 0) {
        int destination = me - 1;

        for (i = 0; i < DIM; i++) {
            upBorderToSend[i] = generation_grid[upBorderLine][i][DOWN];
            generation_grid[upBorderLine][i][DOWN].type = INVALIDO;
        }

        MPI_Isend(upBorderToSend, DIM * sizeof (world_t), MPI_BYTE, destination, DOWN_TO_UP_TAG, MPI_COMM_WORLD, &requests[0]);
    }

    //send down border
    if (me != nprocs - 1) {
        int destination = me + 1;

        for (i = 0; i < DIM; i++) {
            downBorderToSend[i] = generation_grid[downBorderLine][i][UP];
            generation_grid[downBorderLine][i][UP].type = INVALIDO;
        }

        MPI_Isend(downBorderToSend, DIM * sizeof (world_t), MPI_BYTE, destination, UP_TO_DOWN_TAG, MPI_COMM_WORLD, &requests[1]);
    }

    if (me == 0) {
        MPI_Wait(&requests[1], &statuses[1]);

    } else if (me == nprocs - 1) {
        MPI_Wait(&requests[0], &statuses[0]);

    } else {
        MPI_Waitall(2, requests, statuses);
    }

}

void receive_border_updates() {

    world_t upBorderToReceive[DIM];
    world_t downBorderToReceive[DIM];

    int upBorderLine = INIT_LINE - 1;
    int downBorderLine = END_LINE;

    int i;

    MPI_Request requests[2];
    MPI_Status statuses[2];

    //receive down border
    if (me != nprocs - 1) {
        int source = me + 1;

        MPI_Irecv(downBorderToReceive, DIM * sizeof (world_t), MPI_BYTE, source, DOWN_TO_UP_TAG, MPI_COMM_WORLD, &requests[0]);
    }

    //receive up border
    if (me != 0) {
        int source = me - 1;

        MPI_Irecv(upBorderToReceive, DIM * sizeof (world_t), MPI_BYTE, source, UP_TO_DOWN_TAG, MPI_COMM_WORLD, &requests[1]);
    }


    if (me == 0) { // receive down border update
        MPI_Wait(&requests[0], &statuses[0]);
        for (i = 0; i < DIM; i++) {
            generation_grid[downBorderLine - 1][i][DOWN] = downBorderToReceive[i];
        }
    } else if (me == nprocs - 1) { // receive up border update
        MPI_Wait(&requests[1], &statuses[1]);
        for (i = 0; i < DIM; i++) {
            generation_grid[upBorderLine + 1][i][UP] = upBorderToReceive[i];
        }
    } else { // receive down and up border update
        MPI_Waitall(2, requests, statuses);
        for (i = 0; i < DIM; i++) {
            generation_grid[downBorderLine - 1][i][DOWN] = downBorderToReceive[i];
            generation_grid[upBorderLine + 1][i][UP] = upBorderToReceive[i];
        }
    }
}

void send_result_of_update() {
    int i;

    MPI_Request requests[2];
    MPI_Status statuses[2];

    //send up border
    if (me != 0) {
        int destination = me - 1;

        MPI_Isend(world[INIT_LINE], DIM * sizeof (world_t), MPI_BYTE, destination, DOWN_TO_UP_TAG, MPI_COMM_WORLD, &requests[0]);
    }

    //send down border
    if (me != nprocs - 1) {
        int destination = me + 1;

        MPI_Isend(world[END_LINE - 1], DIM * sizeof (world_t), MPI_BYTE, destination, UP_TO_DOWN_TAG, MPI_COMM_WORLD, &requests[1]);
    }

    if (me == 0) {
        MPI_Wait(&requests[1], &statuses[1]);
    } else if (me == nprocs - 1) {
        MPI_Wait(&requests[0], &statuses[0]);
    } else {
        MPI_Waitall(2, requests, statuses);
    }
}

void receive_result_of_update() {

    int upBorderLine = INIT_LINE - 1;
    int downBorderLine = END_LINE;

    int i;

    MPI_Request requests[2];
    MPI_Status statuses[2];

    //receive down border
    if (me != nprocs - 1) {
        int source = me + 1;

        MPI_Irecv(world[downBorderLine], DIM * sizeof (world_t), MPI_BYTE, source, DOWN_TO_UP_TAG, MPI_COMM_WORLD, &requests[0]);
    }

    //receive up border
    if (me != 0) {
        int source = me - 1;

        MPI_Irecv(world[upBorderLine], DIM * sizeof (world_t), MPI_BYTE, source, UP_TO_DOWN_TAG, MPI_COMM_WORLD, &requests[1]);
    }


    if (me == 0) {
        MPI_Wait(&requests[0], &statuses[0]);
    } else if (me == nprocs - 1) {
        MPI_Wait(&requests[1], &statuses[1]);

    } else {
        MPI_Waitall(2, requests, statuses);
    }
}

void printReceiveResult(int me, MPI_Status * status) {
    int msgLen;
    MPI_Get_count(status, MPI_BYTE, &msgLen);
    printf("Me: %d - Recebi do %d ", me, status->MPI_SOURCE);
    printf("%d bytes com a tag '%d'\n", msgLen, status->MPI_TAG);
}

void printSendResult(int me, MPI_Status * status) {
    int msgLen;
    MPI_Get_count(status, MPI_BYTE, &msgLen);
    printf("Me: %d - Enviei ao %d ", me, status->MPI_SOURCE);
    printf("%d bytes com a tag '%d'\n", msgLen, status->MPI_TAG);
}

void printLine(world_t * line) {
    int i;
    puts("Line elements type:");
    putchar('|');
    for (i = 0; i < DIM; i++) {
        printf("%d|", line[i].type);
    }
    putchar('\n');
}

/**
 * FUNTIONS DEFINITION
 */
int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    processInput(argc, argv);


    MPI_Barrier(MPI_COMM_WORLD);
    double start = omp_get_wtime();


    for (current_generation = 0; current_generation < NUM_GENERATIONS; current_generation++) {

        simulateGeneration();
    }


    int sendcounts;
    world_t * myBlockPointer = &(world[INIT_LINE][0]);

    sendcounts = N_LINES * DIM * sizeof (world_t);

    //printf("Me: %d sendCount: %d \n", me, sendcounts);
    int rcvcounts[nprocs];
    int displs[nprocs];

    int i;
    if (me == 0) {
        for (i = 0; i < nprocs; i++) {
            rcvcounts[i] = (PROC_LINES_NUM(i, nprocs, DIM)) * DIM * sizeof (world_t);
            displs[i] = (PROC_INIT_LINE(i, nprocs, DIM)) * DIM * sizeof (world_t);
            //printf("Proc: %d RcvCount: %d Displs: %d\n", i, rcvcounts[i], displs[i]);
        }
    }

    MPI_Gatherv(myBlockPointer, sendcounts, MPI_BYTE, &(world[0][0]), rcvcounts, displs, MPI_BYTE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    double end = omp_get_wtime();
    
    if (me == 0) {
        processOutput();
        fflush(stdout);
    }

    if (me == 0) {
        printf("\nExecution time: %f \n", end - start);
        fflush(stdout);
    }

    MPI_Finalize();
    return 0;
}

/**
 * Run the subgenerations (reds and blacks) and update world grid between and at end
 */
void simulateGeneration() {
    if (nprocs > 1) {
        simulateReds();

        send_border_updates();
        receive_border_updates();

        MPI_Barrier(MPI_COMM_WORLD);

        updateWorldGrid();

        send_result_of_update();
        receive_result_of_update();
        MPI_Barrier(MPI_COMM_WORLD);

        simulateBlacks();

        send_border_updates();
        receive_border_updates();

        MPI_Barrier(MPI_COMM_WORLD);
        updateWorldGrid();

        send_result_of_update();
        receive_result_of_update();
        MPI_Barrier(MPI_COMM_WORLD);
    } else {
        simulateReds();
        updateWorldGrid();
        simulateBlacks();
        updateWorldGrid();
    }
}

/*
 * In an even numbered row, red cells are the ones with an even column number, and, in an odd
 * numbered row, red cells are the ones with an odd column number
 */
void simulateReds() {
    int i, j;
    for (i = INIT_LINE; i < END_LINE; i++) {
        int oddRow = i % 2;
        for (j = 0; j < DIM; j++) {
            if ((!oddRow && (j % 2 == 0)) || /* even Row and even Column */
                    (oddRow && (j % 2 != 0))) /* odd Row and odd Column */
                processCell(i, j);
        }
    }
}

/*
 *
 */
void simulateBlacks() {
    int i, j;
    for (i = INIT_LINE; i < END_LINE; i++) {
        int oddRow = i % 2;
        for (j = 0; j < DIM; j++) {
            if (!((!oddRow && (j % 2 == 0)) || (oddRow && (j % 2 != 0))))
                processCell(i, j);
        }
    }
}

/**
 * Receive row and column and call specific method to process this element
 */
void processCell(int row, int column) {
    pos_t pos = {row, column};
    switch (getWorldElement(pos).type) {
        case WOLF:
            processWolf(pos);
            break;
        case SQUIRREL:
        case SQUIRREL_TREE:
            processSquirrel(pos);
    }
}

/**
 * Process the behavior of the wolf at the given coordinates
 * @param coordinates The coordinates for the wolf
 */
void processWolf(pos_t coordinates) {
    world_t currentWolf = getWorldElement(coordinates);
    pos_t newCoordinates;

    // Account for the starvation and breeding period
    if (currentWolf.last_generation != current_generation) {
        if (currentWolf.starvation_period <= 0) {
            removeAnimal(coordinates);
            return;
        }
        currentWolf.starvation_period--;
        currentWolf.breeding_period--;
        currentWolf.last_generation = current_generation;
    }

    newCoordinates = moveWolf(coordinates);

    // Only breed if the breeding period has ended and if the wolf moved
    if (newCoordinates.row != coordinates.row || newCoordinates.column != coordinates.column) {
        // The wolf moved, so remove it from the previous position
        removeAnimal(coordinates);
        // Breeding		
        if (currentWolf.breeding_period < 0) {
            // Spawn a new wolf
            breedWolf(coordinates);
            currentWolf.breeding_period = WOLF_BREEDING_PERIOD;
        }
    }

    // Update the wolf either at the new coordinates or on the same coordinates
    updateAnimal(currentWolf, coordinates, newCoordinates);
}

/**
 * Attempts to move a wolf at the given coordinates
 * @param coordinates The coordinates for the wolf
 * @return The coordinates in which the wolf ends up at
 */
pos_t moveWolf(pos_t coordinates) {
    int possibleCellCount = 0;
    pos_t targetCoordinates = coordinates;
    pos_t *possibleCells = (pos_t *) malloc(4 * sizeof (pos_t));

    // Find nearby squirrels
    possibleCellCount = findNearbySquirrels(coordinates, possibleCells);
    if (possibleCellCount > 0) {
        // hurray we can eat!
        targetCoordinates = possibleCells[computeOptimalCellChoice(coordinates,
                possibleCellCount)];
        // Add the squirrel to the generation update grid so that it can be accounted for during the conflict management
        updateGenerationGrid(getWorldElement(targetCoordinates),
                targetCoordinates, WORLD);

    } else {
        // Fine, we'll try to move then...
        possibleCellCount = findNearbyEmptyCells(coordinates, 0, possibleCells);
        if (possibleCellCount > 0) {
            targetCoordinates = possibleCells[computeOptimalCellChoice(
                    coordinates, possibleCellCount)];
        }
    }

    free(possibleCells);
    return targetCoordinates;
}

/**
 * Spawns a new wolf at the given coordinates
 * @param coordinates The coordinate in which the new wolf will be spawned
 */
void breedWolf(pos_t coordinates) {
    world_t newWolf;

    // Spawn a new wolf
    newWolf.type = WOLF;
    newWolf.breeding_period = WOLF_BREEDING_PERIOD;
    newWolf.starvation_period = WOLF_STARVATION_PERIOD;

    updateAnimal(newWolf, coordinates, coordinates);
}

/**
 * Process the behavior of the squirrel at the given coordinates
 * @param coordinates The coordinates for the squirrel
 */
void processSquirrel(pos_t coordinates) {
    world_t currentSquirrel = getWorldElement(coordinates);
    pos_t newCoordinates = moveSquirrel(coordinates);

    if (currentSquirrel.last_generation != current_generation) {
        currentSquirrel.breeding_period--;
        currentSquirrel.last_generation = current_generation;
    }

    // Only breed if the breeding period has ended and if the squirrel moved
    if (newCoordinates.row != coordinates.row || newCoordinates.column != coordinates.column) {
        // The squirrel moved, so remove it from the previous position
        removeAnimal(coordinates);
        // Breeding
        if (currentSquirrel.breeding_period < 0) {
            // Spawn a new squirrel
            breedSquirrel(coordinates);
            currentSquirrel.breeding_period = SQUIRREL_BREEDING_PERIOD;
        }
    }

    // Update the squirrel either at the new coordinates or on the same coordinates
    updateAnimal(currentSquirrel, coordinates, newCoordinates);
}

/**
 * Attempts to move the squirrel at the given coordinates
 * @param coordinates The coordinates for the squirrel
 * @return The coordinates in which the squirrel ends up at
 */
pos_t moveSquirrel(pos_t coordinates) {
    int possibleCellCount = 0;
    pos_t targetCoordinate;
    pos_t *possibleCells = (pos_t *) malloc(4 * sizeof (pos_t));

    // Find nearby empty cells or trees
    possibleCellCount = findNearbyEmptyCells(coordinates, 1, possibleCells);

    if (possibleCellCount > 0) {
        // hurray we can move
        targetCoordinate = possibleCells[computeOptimalCellChoice(coordinates,
                possibleCellCount)];

        free(possibleCells);
        return targetCoordinate;
    } else {
        // Ok, we'll stay put
        free(possibleCells);
        return coordinates;
    }

}

/**
 * Spawns a new squirrel at the given coordinates
 * @param coordinates The coordinate in which the new squirrel will be spawned
 */
void breedSquirrel(pos_t coordinates) {
    world_t newSquirrel;

    // Spawn a new squirrel
    newSquirrel.type = SQUIRREL;
    newSquirrel.breeding_period = SQUIRREL_BREEDING_PERIOD;

    updateAnimal(newSquirrel, coordinates, coordinates);
}

/**
 * Updates an animal from a given set of coordinates to a new set of coordinates
 * @param animal The animal to be updated
 * @param newCoordinates The coordinates where the animal will be moved to
 */
void updateAnimal(world_t animal, pos_t coordinates, pos_t newCoordinates) {

    // Replicate the animal in the new cell
    world_t actual_element = getWorldElement(newCoordinates);
    if (animal.type == SQUIRREL && actual_element.type == TREE || actual_element.type == SQUIRREL_TREE) {
        animal.type = SQUIRREL_TREE;
    } else if (animal.type == SQUIRREL_TREE && actual_element.type != TREE && actual_element.type != SQUIRREL_TREE) {
        animal.type = SQUIRREL;
    }

    int dir;

    if (newCoordinates.row > coordinates.row)
        dir = UP;
    else if (newCoordinates.row < coordinates.row)
        dir = DOWN;
    else if (newCoordinates.column > coordinates.column)
        dir = LEFT;
    else if (newCoordinates.row < coordinates.row)
        dir = RIGHT;
    else dir = WORLD;

    updateGenerationGrid(animal, newCoordinates, dir);
}

/**
 * Removes the given animal from the given position
 * @param coordinates The coordinates for the animal to be removed
 */
void removeAnimal(pos_t coordinates) {
    world_t animalFreeSlot = getWorldElement(coordinates);

    if (animalFreeSlot.type == SQUIRREL_TREE) {
        animalFreeSlot.type = TREE;
    } else {
        animalFreeSlot.type = EMPTY;
    }

    updateGenerationGrid(animalFreeSlot, coordinates, WORLD);
}

// Não foi copy paste 

/**
 * Finds nearby empty cells
 * @param coordinates The coordinates from which we're looking for empty cells
 * @param considerTrees If set to 1 trees will also be considered as empty cells
 * @param result A set of coordinates in which empty cells were found
 * @return the number of empty cells found nearby
 */
int findNearbyEmptyCells(pos_t coordinates, int considerTrees, pos_t *result) {
    int p = 0;
    world_t element;

    pos_t testPos;

    // Find nearby cells with squirrels in them
    // 12:00
    testPos.row = coordinates.row - 1;
    testPos.column = coordinates.column;


    if (coordinates.row - 1 >= 0) {
        element = getWorldElement(testPos);
        if ((element.type == EMPTY) || (considerTrees == 1)
                && (element.type == TREE)) {
            result[p].row = coordinates.row - 1;
            result[p].column = coordinates.column;
            p++;
        }
    }
    // 3:00
    testPos.row = coordinates.row;
    testPos.column = coordinates.column + 1;

    if (coordinates.column + 1 < DIM) {
        element = getWorldElement(testPos);
        if ((element.type == EMPTY) || (considerTrees == 1)
                && (element.type == TREE)) {
            result[p].row = coordinates.row;
            result[p].column = coordinates.column + 1;
            p++;
        }
    }

    // 6:00
    testPos.row = coordinates.row + 1;
    testPos.column = coordinates.column;

    if (coordinates.row + 1 < DIM) {
        element = getWorldElement(testPos);
        if ((element.type == EMPTY) || (considerTrees == 1)
                && (element.type == TREE)) {
            result[p].row = coordinates.row + 1;
            result[p].column = coordinates.column;
            p++;
        }
    }

    // 9:00
    testPos.row = coordinates.row;
    testPos.column = coordinates.column - 1;

    if (coordinates.column - 1 >= 0) {
        element = getWorldElement(testPos);
        if ((element.type == EMPTY) || (considerTrees == 1)
                && (element.type == TREE)) {
            result[p].row = coordinates.row;
            result[p].column = coordinates.column - 1;
            p++;
        }
    }
    return p;
}

/**
 * Finds nearby cells containing squirrels
 * @param coordinates The coordinates from which we're looking for squirrels
 * @param result A set of coordinates in which squirrels were found
 * @return the number of squirrels found nearby
 */
int findNearbySquirrels(pos_t coordinates, pos_t *result) {
    int p = 0;
    world_t element;

    pos_t testPos;

    // Find nearby cells with squirrels in them
    // 12:00
    testPos.row = coordinates.row - 1;
    testPos.column = coordinates.column;


    if (coordinates.row - 1 >= 0) {
        element = getWorldElement(testPos);
        if (element.type == SQUIRREL) {
            result[p].row = coordinates.row - 1;
            result[p].column = coordinates.column;
            p++;
        }
    }

    // 3:00
    testPos.row = coordinates.row;
    testPos.column = coordinates.column + 1;


    if (coordinates.column + 1 < DIM) {
        element = getWorldElement(testPos);
        if (element.type == SQUIRREL) {
            result[p].row = coordinates.row;
            result[p].column = coordinates.column + 1;
            p++;
        }
    }

    // 6:00
    testPos.row = coordinates.row + 1;
    testPos.column = coordinates.column;

    if (coordinates.row + 1 < DIM) {
        element = getWorldElement(testPos);
        if (element.type == SQUIRREL) {
            result[p].row = coordinates.row + 1;
            result[p].column = coordinates.column;
            p++;
        }
    }

    // 9:00
    testPos.row = coordinates.row;
    testPos.column = coordinates.column - 1;


    if (coordinates.column - 1 >= 0) {
        element = getWorldElement(testPos);
        if (element.type == SQUIRREL) {
            result[p].row = coordinates.row;
            result[p].column = coordinates.column - 1;
            p++;
        }
    }
    return p;
}

/**
 * Determines the option that should be taken for the animal movement according to the predefined heuristics
 * @param currentPosition The position in which the animal is at
 * @param numberOfChoices The number of choices presented to the animal
 * @return The choice that should be taken
 */
int computeOptimalCellChoice(pos_t currentPosition, int numberOfChoices) {
    return (currentPosition.row * DIM + currentPosition.column) % numberOfChoices;
}

/**
 * update the world grid with updates from current generation grid
 */
void updateWorldGrid() {
    int i, j, k;

    for (i = INIT_LINE; i < END_LINE; i++) {
        for (j = 0; j < DIM; j++) {
            pos_t pos = {i, j};
            int number_wolfs = 0, number_squirrels = 0, number_empty = 0, number_trees = 0;
            world_t final_squirrel = {-1, SQUIRREL, INT_MIN, INT_MIN};
            world_t final_wolf = {-1, WOLF, INT_MIN, INT_MIN};
            world_t element_grid;

            for (k = 0; k < 5; k++) {
                element_grid = getGenerationElement(pos, k);
                if (element_grid.type != INVALIDO) {
                    if (element_grid.type == SQUIRREL || element_grid.type == SQUIRREL_TREE) {
                        number_squirrels++;
                        if (element_grid.breeding_period > final_squirrel.breeding_period)
                            final_squirrel = element_grid;
                    } else if (element_grid.type == WOLF) {
                        number_wolfs++;
                        if (element_grid.starvation_period > final_wolf.starvation_period)
                            final_wolf = element_grid;
                        else if (element_grid.starvation_period == final_wolf.starvation_period) {
                            if (element_grid.breeding_period > final_wolf.breeding_period)
                                final_wolf = element_grid;
                        }
                    } else if (element_grid.type == EMPTY) {
                        number_empty++;
                    } else if (element_grid.type == TREE)
                        number_trees++;
                    element_grid.type = INVALIDO;
                    updateGenerationGrid(element_grid, pos, k);
                }
            }

            if (number_squirrels == 0 && number_wolfs == 0) {
                if (number_trees > 0)
                    setWorldElementType(pos, TREE);
                else if (number_empty > 0)
                    setWorldElementType(pos, EMPTY);
            } else if (number_squirrels > 0 && number_wolfs == 0) {// squirrels only
                if (getWorldElement(pos).type == TREE)
                    final_squirrel.type = SQUIRREL_TREE;
                setWorldElement(pos, final_squirrel);
            } else if (number_wolfs > 0 && number_squirrels == 0) { // wolves only
                setWorldElement(pos, final_wolf);
            } else if (number_wolfs > 0 && number_squirrels > 0) {// squirrels and wolves
                final_wolf.starvation_period = WOLF_STARVATION_PERIOD;
                setWorldElement(pos, final_wolf);
            }
        }
    }
}

/**
 * Receive the input of the program and initialize structs and variables
 */
void processInput(int argc, char ** argv) {
    if (argc < 6) {
        fprintf(stderr, "Número de argumentos insuficiente\n");
        exit(-1);
    }

    WOLF_BREEDING_PERIOD = atoi(argv[2]);
    SQUIRREL_BREEDING_PERIOD = atoi(argv[3]);
    WOLF_STARVATION_PERIOD = atoi(argv[4]);
    NUM_GENERATIONS = atoi(argv[5]);

    FILE * grid_file = fopen(argv[1], "r");

    if (!grid_file) {
        fprintf(stderr, "Erro na abertura do ficheiro: %s\n", argv[1]);
        exit(-1);
    }

    int dim;
    fscanf(grid_file, "%d\n", &dim);

    DIM = dim;


    INIT_LINE = PROC_INIT_LINE(me, nprocs, DIM);
    N_LINES = PROC_LINES_NUM(me, nprocs, DIM);
    END_LINE = PROC_END_LINE(me, nprocs, DIM);

    initGrids();
    initBorderPositions();

    while (!feof(grid_file)) {
        int ln, cn;
        char type;
        fscanf(grid_file, "%d %d %c\n", &ln, &cn, &type);

        int toMe = ln >= INIT_LINE && ln < END_LINE;

        if (toMe || isBorder(me, ln)) {

            switch (type) {
                case 's':
                    if (world[ln][cn].type == TREE) {
                        world[ln][cn].type = SQUIRREL_TREE;
                    } else
                        world[ln][cn].type = SQUIRREL;
                    world[ln][cn].breeding_period = SQUIRREL_BREEDING_PERIOD;
                    break;
                case 'w':
                    world[ln][cn].type = WOLF;
                    world[ln][cn].breeding_period = WOLF_BREEDING_PERIOD;
                    world[ln][cn].starvation_period = WOLF_STARVATION_PERIOD;
                    break;
                case 'i':
                    world[ln][cn].type = ICE;
                    break;
                case 't':
                    if (world[ln][cn].type == SQUIRREL)
                        world[ln][cn].type = SQUIRREL_TREE;
                    else
                        world[ln][cn].type = TREE;
                    break;
                default:
                    fprintf(stderr, "Ficheiro de input contêm dados inválidos\n");
            }
        }
    }

    if (fclose(grid_file) == EOF) {
        fprintf(stderr, "Erro no fecho do ficheiro: %s\n", argv[1]);
        exit(-1);
    }
}

void initGrids() {
    int i, j, k;

    //alloc pointers 
    world = (world_t **) malloc(DIM * sizeof (world_t *));
    generation_grid = (world_t ***) malloc(DIM * sizeof (world_t **));

    verifyAlloc(world, "world");
    verifyAlloc(generation_grid, "generation_grid");


    world_t * proc_world_block;

    if (me == 0) {
        proc_world_block = (world_t *) malloc(DIM * DIM * sizeof (world_t));
    } else {
        proc_world_block = (world_t *) malloc(N_LINES * DIM * sizeof (world_t));
    }

    verifyAlloc(proc_world_block, "proc_world_block");

    int procBlocksIndex = 0;

    for (i = 0; i < DIM; i++) {


        if (me == 0 || (i >= INIT_LINE && i < END_LINE)) {
            world[i] = &proc_world_block[procBlocksIndex * DIM];
            procBlocksIndex++;

            if (i >= INIT_LINE && i < END_LINE) {
                generation_grid[i] = (world_t **) malloc(DIM * sizeof (world_t *));


                for (j = 0; j < DIM; j++) {
                    world[i][j].type = EMPTY;
                    world[i][j].last_generation = -1;
                    generation_grid[i][j] = (world_t *) malloc(5 * sizeof (world_t));
                    for (k = 0; k < 5; k++) {
                        generation_grid[i][j][k].type = INVALIDO;
                    }
                }
            }
        } else {
            world[i] = NULL;
        }
    }
}

void initBorderPositions() {

    int upBorderLine = INIT_LINE - 1;
    int downBorderLine = END_LINE;

    int i, j, k;

    if (me != 0) {
        world[upBorderLine] = (world_t *) malloc(DIM * sizeof (world_t));
        generation_grid[upBorderLine] = (world_t **) malloc(DIM * sizeof (world_t *));

        verifyAlloc(world[upBorderLine], "world[upBorderLine]");
        verifyAlloc(generation_grid[upBorderLine], "generation_grid[upBorderLine]");

        for (i = 0; i < DIM; i++) {
            world[upBorderLine][i].type = EMPTY;
            world[upBorderLine][i].last_generation = -1;
            generation_grid[upBorderLine][i] = (world_t *) malloc(5 * sizeof (world_t));

            verifyAlloc(generation_grid[upBorderLine][i], "generation_grid[upBorderLine][i]");

            for (k = 0; k < 5; k++) {
                generation_grid[upBorderLine][i][k].type = INVALIDO;
            }
        }
    }

    if (me != nprocs - 1) {
        if (me != 0)
            world[downBorderLine] = (world_t *) malloc(DIM * sizeof (world_t));
        generation_grid[downBorderLine] = (world_t **) malloc(DIM * sizeof (world_t *));

        verifyAlloc(world[downBorderLine], "world[downBorderLine]");
        verifyAlloc(generation_grid[downBorderLine], "generation_grid[downBorderLine]");

        for (j = 0; j < DIM; j++) {
            world[downBorderLine][j].type = EMPTY;
            world[downBorderLine][j].last_generation = -1;
            generation_grid[downBorderLine][j] = (world_t *) malloc(5 * sizeof (world_t));

            verifyAlloc(generation_grid[downBorderLine][j], "generation_grid[downBorderLine][j]");

            for (k = 0; k < 5; k++) {
                generation_grid[downBorderLine][j][k].type = INVALIDO;
            }
        }
    }
}

int isBorder(int me, int ln) {
    int upBorder = -1;
    int downBorder = -1;
    if (me != 0)
        upBorder = INIT_LINE - 1;
//        upBorder = (PROC_END_LINE(me - 1, nprocs, DIM)) - 1;
    if (me != (nprocs - 1))
        downBorder = END_LINE;
//        downBorder = PROC_INIT_LINE(me + 1, nprocs, DIM);
    return ln == upBorder || ln == downBorder;
}

/**
 * Output the GRID (Standard)
 */
void processOutput() {
    if (me == 0) {
        fprintf(stdout, "%d\n", DIM);
    }

    //printf("proc: %d begin\n", me);
    int i, j;
    for (i = 0; i < DIM; i++) {
        //if (i >= INIT_LINE && i < END_LINE) {
        for (j = 0; j < DIM; j++) {
            switch (world[i][j].type) {
                case WOLF:
                    fprintf(stdout, "%d %d %c\n", i, j, 'w');
                    break;
                case SQUIRREL:
                    fprintf(stdout, "%d %d %c\n", i, j, 's');
                    break;
                case ICE:
                    fprintf(stdout, "%d %d %c\n", i, j, 'i');
                    break;
                case TREE:
                    fprintf(stdout, "%d %d %c\n", i, j, 't');
                    break;
                case SQUIRREL_TREE:
                    fprintf(stdout, "%d %d %c\n", i, j, '$');
                    break;
            }
        }
        //}
    }
    //printf("proc: %d end\n", me);
}

void verifyAlloc(void * memoryPointer, char * errorMSG) {
    if (!memoryPointer) {
        fprintf(stderr, "Erro na alocacao de memoria: %s\n", errorMSG);
        exit(-1);
    }
}

/**
 * Print the GRID graphically (helpful in debug)
 */
void printMyGrid(int me) {
    int i, j;
    printf("Proc: %d grid:\n", me);
    for (i = 0; i < DIM; i++) {
        if (i >= INIT_LINE && i < END_LINE) {
            for (j = 0; j < DIM; j++) {
                switch (world[i][j].type) {
                    case WOLF:
                        printf(" W ");
                        break;
                    case SQUIRREL:
                        printf(" S ");
                        break;
                    case ICE:
                        printf(" I ");
                        break;
                    case TREE:
                        printf(" T ");
                        break;
                    case SQUIRREL_TREE:
                        printf(" $ ");
                        break;
                    default:
                        printf(" - ");
                }
            }
            puts("\n");
        }
    }
    fflush(stdout);
}

/**
 * Print the the details of the GRID graphically (helpful in debug)
 */
void printMyGridDetailed(int me) {
    int i, j;
    printf("Proc: %d detailed grid:\n", me);
    for (i = 0; i < DIM; i++) {
        if (i >= INIT_LINE && i < END_LINE) {
            for (j = 0; j < DIM; j++) {
                switch (world[i][j].type) {
                    case WOLF:
                        printf("[W %2d %2d]", world[i][j].starvation_period,
                                world[i][j].breeding_period);
                        break;
                    case SQUIRREL:
                        printf("[S  - %2d]", world[i][j].breeding_period);
                        break;
                    case ICE:
                        printf("[   I   ]");
                        break;
                    case TREE:
                        printf("[   T   ]");
                        break;
                    case SQUIRREL_TREE:
                        printf("[$  - %2d]", world[i][j].breeding_period);
                        break;
                    default:
                        printf("[       ]");
                }
            }
            puts("\n");
        }
    }
}
