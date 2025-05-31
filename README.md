# IntrinsicFrustration

Project: Higher-order Interactions Induce Chimera States in Globally Coupled Oscillators

This project reports the discovery of a novel Chimera state in a higher-order Kuramoto model.
We perform both numerical simulations and theoretical analysis of this state.
Details can be found in our paper:
"Higher-order Interactions Induce Chimera States in Globally Coupled Oscillators".

The project consists of three parts:

1. Anime
   Run the script to view an animation of the Chimera state we discovered.

2. Simulation
   This corresponds to the numerical simulations shown in Fig. 4 of the paper.

   a. Run the script `simulate_z.m` to obtain the partial order parameters `z_p` and `z_n`,
      representing the oscillators rotating in the positive and negative directions, respectively,
      under different values of K and omega_0.

   b. Simulation results are also included in this folder.
      Files are named like `z1ChimeraSim910D1.mat`, where:
         - `z1` denotes the partial order parameter for positive-frequency oscillators;
           `z2` would represent the negative-frequency ones;
         - `910` indicates the ratio of oscillators in the three synchronous clusters (9:1:0);
         - `D1` means each peak has a standard deviation D = 1.

3. Theory
   This corresponds to the theoretical analysis shown in Fig. 4 of the paper.

   a. `Theory_R31_vs_K_omega0.m` computes the boundary between the chimera state and traveling wave state.
      `Theory_R31_vs_K_omega0_Chimera.m` computes the boundary between the chimera and incoherent states.

   b. `R1_R3.m` is used to numerically determine the functional relationship between
      the partial order parameters R_{1,p} and R_{3,p} in each synchronous group.
      The outputs are saved in files like `R1_910.mat` and `R3_910.mat`.

   c. The results from `Theory_R31_vs_K_omega0.m` are saved in files named like
      `D1_TR910_150_200.mat`, where:
         - `D1` denotes the standard deviation;
         - `TR` stands for "theoretical result";
         - `910` is the same group ratio as above;
         - `150_200` indicates the range of omega_0 and K values.
      Note: Results from `Theory_R31_vs_K_omega0_Chimera.m` are not included, as we found
      the boundary between chimera and incoherent states is independent of omega_0.
      Therefore, we only calculated z_p vs. K to draw the theoretical boundary in Fig. 4.
      The time step used in generating results from `Theory_R31_vs_K_omega0.m` is dt = 0.0001 / 3.

Environment:
- MATLAB R2021a or later
- No additional toolboxes are required

Folder Structure:
- /Anime       : Chimera animation script
- /Simulation  : Simulation scripts and .mat result files
- /Theory      : Theoretical analysis scripts and .mat result files

How to Run:
1. View animation:
   Go to the `Anime` folder and run the animation script.

2. Run simulation:
   Go to the `Simulation` folder and run `simulate_z.m`.

3. Perform theoretical analysis:
   Go to the `Theory` folder and run `Theory_R31_vs_K_omega0.m` or `Theory_R31_vs_K_omega0_Chimera.m`.

File Naming Conventions:
- `z1ChimeraSim910D1.mat`: z1 = positive rotation, 910 = 9:1:0 ratio, D1 = standard deviation
- `D1_TR910_150_200.mat`: theoretical result for D = 1, ratio 910, omega0 range = 150, K range = 200

Citation:
If you use this code in your research, please cite our paper:
"Higher-order Interactions Induce Chimera States in Globally Coupled Oscillators".

License:
This project is distributed under the MIT License. See LICENSE file for details.
