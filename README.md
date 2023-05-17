# Simulation of Structure Change in Porous Media During Gas-Solid Reactions Using Cellular Automata Model

This repository contains a Python code for simulating the structure change in porous media during gas-solid reactions using a cellular automata model. The code is organized using object-oriented programming (OOP) principles, with various functionalities implemented in separate files to enhance readability and accessibility.

## Getting Started

To run the code, execute the "Main.py" file, which serves as the entry point for the simulation.

Upon running the code, you will be presented with two options:

1. Run a New Simulation: Choosing this option prompts you to input values for the required parameters. To assist you in getting started, sample input values are provided as comments after each input syntax in the "Initialize.py" file, giving you a rough idea of the expected input format.

2. Load Data from Previous Simulation: If you have previously saved simulation data, select this option to load the data and continue the simulation from where it left off.

## Obstacle Options

The code provides multiple options for defining obstacles within the simulation:

- Single Square Obstacle
- Single Round Obstacle
- Square Porous Obstacle (a combination of several square obstacles arranged in organized patterns)
- Single Round Porous Obstacle produced by QSGS method
- Loading from a previous simulation file

Feel free to choose the obstacle option that best suits your simulation requirements.

## Data Saving and Results

The code is structured to save the data during the simulation, enabling you to access it later or use it to continue the simulation from a specific point. You have the flexibility to choose the file name and location to save the simulation results.

## Collision Rules

The collision rules, which govern the interactions between different components in the simulation, are implemented using a dictionary named "CollisionTable1". To facilitate faster access, these collision rules are saved in a file named "CollisionRules.pkl". If you wish to modify these rules, update the "CollisionTable1" dictionary in the "Dictionaries.py" file. After making the changes, call the "MakeCollisionRules1" function from the "OtherFunctions.py" file to generate a new "CollisionRules.pkl" file with the updated rules.

## Binary and Decimal Conversion

The status of the cells is defined using decimal values, which need to be converted to their binary equivalents and vice versa during the simulation. To improve performance, we have precomputed these equivalents using the "MakeBin2DecDic" function in the "OtherFunctions.py" file. The results are saved in the "Bin2DecDic.pkl" file for faster access.

## File Structure

The repository contains the following files:

- Main.py: The main file for executing and starting the simulation.
- Initialize.py: Contains functions for initializing the simulation parameters and obtaining user input.
- Classes.py: Implements all the classes used to implement the logic using cellular automata model.
- OtherFunctions.py: Contains additional helper functions used in the simulation.
- Dictionaries.py: Defines dictionaries and collision rules used in the simulation.
- CollisionRules.pkl: A file containing the pre-computed collision rules for faster access.
- Bin2DecDic.pkl: A file containing the pre-computed binary string equivalents of all possible decimal statuses of the cells for faster access.
- README.md: This file, providing an overview of the code and its functionalities.

For more detailed comments and implementation details, please refer to the individual files.

Enjoy simulating the structure change in porous media during gas-solid reactions using this code! If you have any ideas to improve this code, please don't hesitate to reach out.