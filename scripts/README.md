# Scripts:
Contains all scripts necessary to reproduce main text figures 2,3,4 and supplemental figures S2 and S3. Supplemental figure S2 is generated as a byproduct of the analysis used to generate figure 3, so is contained in figure_3.m. To run scripts with appropriate functions, please either run them in their current folder, or adjust the addpath() command at the top of each Matlab script to reflect where the functions are stored on your machine. 

## Scripts are as follows:

### figure_2.m 
Produces main text figure 2. Runs dynamics for no-feedback model under scarce iron and iron replete payoff scenarios. 

### figure_3_and_A2.m 
Produces main text figure 3, appendix figure 2. Runs dynamics for environmental feedback model and provides example simulations of all model behaviors. 

### figure_A3.m 
Produces supplemental figure 3. Demonstrates the difference between unstable interior fixed point (heteroclinic network) and no interior fixed point (2D orbits) behaviors using same parameters as heteroclinic network example in plotting_feedback_cases.m. For more details, see [supplemental information]

### figure_4.R (R code)
Produces main text figure 4. Identifies cycles in a heteroclinic network given model parameters. 