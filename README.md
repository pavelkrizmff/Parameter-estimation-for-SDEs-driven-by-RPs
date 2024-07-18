This project contains source codes used for simulations and calculations, the results of which are presented in the paper "https://arxiv.org/abs/2403.12610" and should be considered as a supplement to the paper.
In particular, it consists of the following files:

1. simulateSDEwithRP.nb:

   Wolfram Notebook file for simulating the solutions to SDEs driven by an additive Rosenblatt process. Resulting trajectories are exported in *.csv format.

  
3. runEstimators.R:

    R file to be run. Reads the data from *.csv files, evaluates the estimates (using "estimatorFunctions.R") and generates files with vizualisation of the results. These figures are used in the paper above.

5. estimatorFunctions.R:

    R file with implemented estimators of  Hurst parameter H, noise intensity sigma and drift intensity lambda. These estimators are introduced and studied in the above mentioned paper.
