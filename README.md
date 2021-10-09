# The double traveling salesman problem with partial last‐in‐first‐out loading constraints

This repository contains the source code and data associated to the paper ["The double traveling salesman problem with partial last‐in‐first‐out loading constraints"](https://doi.org/10.1111/itor.12876) by Jonatas B. C. Chagas, Túlio A. M. Toffolo,  Marcone J. F. Souza, and Manuel Iori. The paper presents two Integer Linear Programming (ILP) formulations and a heuristic algorithm based on Biased Random-Key Genetic Algorithm (BRKGA) for solving the Double Traveling Salesman Problem with Partial Last-In-First-Out Loading Constraints (DTSPPL).

### Compiling the code

Before running our solution approaches, it is needed to compile their code. To this end, just run the following command:

```console
$ make
```

### Usage:

```console
$ ./dtsppl [parameters]

parameters:
                   --approach <approach_name> (options: ILP1, ILP2, or BRKGA)
                    --pickuparea <pickup_area_file_name> 
                    --deliveryarea <delivery_area_file_name> 
                    --n <number_of_items> 
                    --l <reloading_depth> 
                    --h <relocation_cost> 
                    --outputsolution <solution_file_name>
```
We provide a python script (see "src/run_all_experiments.py") for running each solution approach on instances described in the paper.
