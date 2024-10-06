# 5GNRSimulation
Code used to generate the results presented in the article "Cross-Layer Latency Analysis for 5G NR in V2X Communications"

Prerequisites:
-To perform the simulation setup, access to the SHARPE (https://sharpe.pratt.duke.edu/) tool is necessary.
-Julia language. Packages: LinearAlgebra,MAT, Match, and Distributions.

How to use:
-The main code of the simulator is located in the file Simulation_5G_Main.jl. To execute this file, the following command should be used:
>julia Simulation_5G_Main.jl

The simulation settings, mainly the simulation duration, can be configured by modifying the event_no variable located in the Simulation_5G_Main.jl file. The behavior of the various protocols associated with 5G New Radio (NR) is defined within the Simulation_5G.jl file. It is permissible to modify the parameters of these protocols in order to evaluate diverse behaviors that may emerge in 5G networks. 
