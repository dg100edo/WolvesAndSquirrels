#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <omp.h>

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
 * GLOBAL VARIABLES
 */
int NUM_GENERATIONS;
int WOLF_BREEDING_PERIOD;
int SQUIRREL_BREEDING_PERIOD;
int WOLF_STARVATION_PERIOD;
int DIM;
int current_generation;

world_t **world; // world grid

world_t ***generation_grid; // generation grid (current generation)

/**
 * MAIN FUNCTIONS DECLARATION
 */
void simulateGeneration();
void simulateReds();
void simulateBlacks();
void processCell(int row, int column);
void processWolf(pos_t coordinates);
void processSquirrel(pos_t coordinates);
void initGrids();
void resetGrids();
void updateGenerationGrid(world_t element, pos_t pos, int direction);
void updateWorldGrid();

/*
 * AUX FUNCTIONS DECLARATION
 */
void processInput(int argc, char ** argv);
void processOutput();
void printGrid();
void printGridDetailed();
pos_t moveSquirrel(pos_t coordinates);
pos_t moveWolf(pos_t coordinates);
void updateAnimal(world_t animal, pos_t coordinates, pos_t newCoordinates);
void removeAnimal(pos_t coordinates);
void breedWolf(pos_t coordinates);
void breedSquirrel(pos_t coordinates);
int findNearbySquirrels(pos_t coordinates, pos_t *result);
int findNearbyEmptyCells(pos_t coordinates, int considerTrees, pos_t *result);
int computeOptimalCellChoice(pos_t currentPosition, int numberOfChoices);

/**
 * FUNTIONS DEFINITION
 */
int main(int argc, char ** argv) {
	processInput(argc, argv);
	double start = omp_get_wtime();
	for (current_generation = 0; current_generation < NUM_GENERATIONS; current_generation++) {
		simulateGeneration();
	}	
	double end = omp_get_wtime();
	processOutput();
	resetGrids();
	printf("\nExecution time: %f \n",end-start);
	return 0;
}

/**
 * Run the subgenerations (reds and blacks) and update world grid between and at end
 */
void simulateGeneration() {
     	simulateReds();
        updateWorldGrid();  
	//printf("\nGereration %d Red\n", current_generation + 1); printGridDetailed();   
	simulateBlacks();
       	updateWorldGrid();
	//printf("\nGereration %d Black\n", current_generation + 1); printGridDetailed(); 
}

/*
 * In an even numbered row, red cells are the ones with an even column number, and, in an odd
 * numbered row, red cells are the ones with an odd column number
 */
void simulateReds() {
	int i, j;
	for (i = 0; i < DIM; i++) {
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
	for (i = 0; i < DIM; i++) {
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
	pos_t pos = { row, column };	
	switch (world[row][column].type) {
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
	world_t currentWolf = world[coordinates.row][coordinates.column];
	pos_t newCoordinates;

	// Account for the starvation and breeding period
	if(currentWolf.last_generation != current_generation) {	
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
	if (newCoordinates.row != coordinates.row || newCoordinates.column != coordinates.column){
		// The wolf moved, so remove it from the previous position
		removeAnimal(coordinates);
		// Breeding		
		if (currentWolf.breeding_period < 0){
			// Spawn a new wolf
			breedWolf(coordinates);
			currentWolf.breeding_period = WOLF_BREEDING_PERIOD;
		}
	}

	// Update the wolf either at the new coordinates or on the same coordinates
	updateAnimal(currentWolf, coordinates, newCoordinates);	
}

/**
 * Process the behavior of the squirrel at the given coordinates
 * @param coordinates The coordinates for the squirrel
 */
void processSquirrel(pos_t coordinates) {
	world_t currentSquirrel = world[coordinates.row][coordinates.column];
	pos_t newCoordinates = moveSquirrel(coordinates);

	if(currentSquirrel.last_generation != current_generation){
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
	pos_t *possibleCells = (pos_t *) malloc(4 * sizeof(pos_t));

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
 * Attempts to move a wolf at the given coordinates
 * @param coordinates The coordinates for the wolf
 * @return The coordinates in which the wolf ends up at
 */
pos_t moveWolf(pos_t coordinates) {
	int possibleCellCount = 0;
	pos_t targetCoordinates = coordinates;
	pos_t *possibleCells = (pos_t *) malloc(4 * sizeof(pos_t));

	// Find nearby squirrels
	possibleCellCount = findNearbySquirrels(coordinates, possibleCells);
	if (possibleCellCount > 0) {
		// hurray we can eat!
		targetCoordinates = possibleCells[computeOptimalCellChoice(coordinates,
				possibleCellCount)];
		// Add the squirrel to the generation update grid so that it can be accounted for during the conflict management
		updateGenerationGrid(
				world[targetCoordinates.row][targetCoordinates.column],
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
 * Updates an animal from a given set of coordinates to a new set of coordinates
 * @param animal The animal to be updated
 * @param newCoordinates The coordinates where the animal will be moved to
 */
void updateAnimal(world_t animal, pos_t coordinates, pos_t newCoordinates) {

	// Replicate the animal in the new cell
	world_t actual_element = world[newCoordinates.row][newCoordinates.column];
	if(animal.type == SQUIRREL && actual_element.type == TREE || actual_element.type == SQUIRREL_TREE){
		animal.type = SQUIRREL_TREE;
	}else if(animal.type == SQUIRREL_TREE && actual_element.type != TREE && actual_element.type != SQUIRREL_TREE){
        	animal.type = SQUIRREL;
        }

	int dir;

	if(newCoordinates.row > coordinates.row)
		dir = UP;
	else if(newCoordinates.row < coordinates.row)
		dir = DOWN;
	else if(newCoordinates.column > coordinates.column)
		dir = LEFT;
	else if(newCoordinates.row < coordinates.row)
		dir = RIGHT;
	else dir = WORLD;

	updateGenerationGrid(animal, newCoordinates, dir);
}

/**
 * Removes the given animal from the given position
 * @param coordinates The coordinates for the animal to be removed
 */
void removeAnimal(pos_t coordinates) {
	world_t animalFreeSlot = world[coordinates.row][coordinates.column];

	if (animalFreeSlot.type == SQUIRREL_TREE) {
		animalFreeSlot.type = TREE;
	} else {
		animalFreeSlot.type = EMPTY;
	}

	updateGenerationGrid(animalFreeSlot, coordinates, WORLD);
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
 * Finds nearby cells containing squirrels
 * @param coordinates The coordinates from which we're looking for squirrels
 * @param result A set of coordinates in which squirrels were found
 * @return the number of squirrels found nearby
 */
int findNearbySquirrels(pos_t coordinates, pos_t *result) {
	int p = 0;

	// Find nearby cells with squirrels in them
	// 12:00
	if ((coordinates.row - 1 >= 0)
			&& (world[coordinates.row - 1][coordinates.column].type == SQUIRREL)) {
		result[p].row = coordinates.row - 1;
		result[p].column = coordinates.column;
		p++;
	}
	// 3:00
	if ((coordinates.column + 1 < DIM)
			&& (world[coordinates.row][coordinates.column + 1].type == SQUIRREL)) {
		result[p].row = coordinates.row;
		result[p].column = coordinates.column + 1;
		p++;
	}
	// 6:00
	if ((coordinates.row + 1 < DIM)
			&& (world[coordinates.row + 1][coordinates.column].type == SQUIRREL)) {
		result[p].row = coordinates.row + 1;
		result[p].column = coordinates.column;
		p++;
	}
	// 9:00
	if ((coordinates.column - 1 >= 0)
			&& (world[coordinates.row][coordinates.column - 1].type == SQUIRREL)) {
		result[p].row = coordinates.row;
		result[p].column = coordinates.column - 1;
		p++;
	}

	return p;
}

/**
 * Finds nearby empty cells
 * @param coordinates The coordinates from which we're looking for empty cells
 * @param considerTrees If set to 1 trees will also be considered as empty cells
 * @param result A set of coordinates in which empty cells were found
 * @return the number of empty cells found nearby
 */
int findNearbyEmptyCells(pos_t coordinates, int considerTrees, pos_t *result) {
	int p = 0;

	// Find nearby cells with squirrels in them
	// 12:00
	if ((coordinates.row - 1 >= 0)
			&& ((world[coordinates.row - 1][coordinates.column].type == EMPTY)
					|| ((considerTrees == 1))
							&& (world[coordinates.row - 1][coordinates.column].type
									== TREE))) {
		result[p].row = coordinates.row - 1;
		result[p].column = coordinates.column;
		p++;
	}
	// 3:00
	if ((coordinates.column + 1 < DIM)
			&& ((world[coordinates.row][coordinates.column + 1].type == EMPTY)
					|| ((considerTrees == 1))
							&& (world[coordinates.row][coordinates.column + 1].type
									== TREE))) {
		result[p].row = coordinates.row;
		result[p].column = coordinates.column + 1;
		p++;
	}
	// 6:00
	if ((coordinates.row + 1 < DIM)
			&& ((world[coordinates.row + 1][coordinates.column].type == EMPTY)
					|| ((considerTrees == 1))
							&& (world[coordinates.row + 1][coordinates.column].type
									== TREE))) {
		result[p].row = coordinates.row + 1;
		result[p].column = coordinates.column;
		p++;
	}
	// 9:00
	if ((coordinates.column - 1 >= 0)
			&& ((world[coordinates.row][coordinates.column - 1].type == EMPTY)
					|| ((considerTrees == 1))
							&& (world[coordinates.row][coordinates.column - 1].type
									== TREE))) {
		result[p].row = coordinates.row;
		result[p].column = coordinates.column - 1;
		p++;
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
 * Receive a element and add it to the list of element in a pos (add to generation grid)
 */
void updateGenerationGrid(world_t element, pos_t pos, int dir) {	
	generation_grid[pos.row][pos.column][dir] = element; 
}

/**
 * update the world grid with updates from current generation grid
 */
void updateWorldGrid() {
	int i, j, k;
	for (i = 0; i < DIM; i++) {
		for (j = 0; j < DIM; j++) {
			int number_wolfs = 0, number_squirrels = 0, number_empty = 0, number_trees = 0;
			world_t final_squirrel = { -1, SQUIRREL, INT_MIN, INT_MIN };
			world_t final_wolf = { -1, WOLF, INT_MIN, INT_MIN };
			world_t element_grid;

			for(k= 0; k < 5; k++){
				element_grid = generation_grid[i][j][k];
				if(element_grid.type != INVALIDO){		
					if(element_grid.type == SQUIRREL || element_grid.type == SQUIRREL_TREE) {
						number_squirrels++;
						if(element_grid.breeding_period > final_squirrel.breeding_period)
							final_squirrel = element_grid;
					}else if (element_grid.type == WOLF) {
						number_wolfs++;
						if (element_grid.starvation_period > final_wolf.starvation_period)
							final_wolf = element_grid;
						else if(element_grid.starvation_period == final_wolf.starvation_period) {
							if(element_grid.breeding_period > final_wolf.breeding_period)
								final_wolf = element_grid;
						}
					}else if(element_grid.type == EMPTY){
		                            number_empty++;
		                        }else if(element_grid.type == TREE)
		                            number_trees++;		   
					generation_grid[i][j][k].type = INVALIDO;
				}
                        }
			                        
                        if( number_squirrels == 0 && number_wolfs == 0 ){
                            if(number_trees > 0)
                                world[i][j].type = TREE;
                            else if(number_empty > 0)
                                world[i][j].type = EMPTY;
                        } else if (number_squirrels > 0 && number_wolfs == 0) {// squirrels only
				if (world[i][j].type == TREE)
					final_squirrel.type = SQUIRREL_TREE;
				world[i][j] = final_squirrel;
			} else if (number_wolfs > 0 && number_squirrels == 0) {	// wolves only
				world[i][j] = final_wolf;
			} else if (number_wolfs > 0 && number_squirrels > 0) {// squirrels and wolves
				final_wolf.starvation_period = WOLF_STARVATION_PERIOD;
                                world[i][j] = final_wolf;
			}
		}
	}
}

/**world = (world_t **) malloc((DIM) * sizeof(world_t *));

 * Create a generation_grid in the start of a generation
 */
void initGrids() {
	int i, j, k;	
	
	//alloc pointers
	world = (world_t **) malloc(DIM * sizeof(world_t *));
	generation_grid = (world_t ***) malloc(DIM * sizeof(world_t **));
	
	//alloc blog of memory to world and generation grid
	world_t * block_world = (world_t *) malloc(DIM * DIM * sizeof(world_t));
	world_t ** block_generation_grid = (world_t **) malloc(DIM * DIM * sizeof(world_t *));

	//alloc blok of memory to vector of conflits
	world_t * block_conflito = (world_t *) malloc(5 * DIM * DIM * sizeof(world_t));
	
	for (i = 0; i < DIM; i++) {
		world[i] = &block_world[i * DIM];
		generation_grid[i] = &block_generation_grid[i * DIM];
		for (j = 0; j < DIM; j++){
			world[i][j].type = EMPTY;
                        world[i][j].last_generation = -1;

			generation_grid[i][j] = &block_conflito[i * 5 * DIM + j * 5];
			for( k = 0; k < 5 ; k++){
				generation_grid[i][j][k].type = INVALIDO;			
			}
		}
	}
}

/**
 * Free resources from generation_grid in the end of a generation
 */
void resetGrids() {
	free(world[0]);
	free(world);
	free(generation_grid[0][0]);
	free(generation_grid[0]);
	free(generation_grid);
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

	initGrids();
	
	while (!feof(grid_file)) {
		int ln, cn;
		char type;
		fscanf(grid_file, "%d %d %c\n", &ln, &cn, &type);
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

	if (fclose(grid_file) == EOF) {
		fprintf(stderr, "Erro no fecho do ficheiro: %s\n", argv[1]);
		exit(-1);
	}
}

/**
 * Output the GRID (Standard)
 */
void processOutput() {
	fprintf(stdout, "%d\n", DIM);
	int i, j;
	for (i = 0; i < DIM; i++) {
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
	}
}

/**
 * Print the GRID graphically (helpful in debug)
 */
void printGrid() {
	int i, j;
	for (i = 0; i < DIM; i++) {
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

/**
 * Print the the details of the GRID graphically (helpful in debug)
 */
void printGridDetailed() {
	int i, j;
	for (i = 0; i < DIM; i++) {
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
