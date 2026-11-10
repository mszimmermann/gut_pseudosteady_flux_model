# gut\_pseudosteady\_flux\_model

Workflow for coarse-grain intestinal flux modelling of metabolite abundances along the gastrointestinal tract assuming pseudo-steady state



Requirements: scripts were developed with MATLAB 2023a. 



workflow\_beta.m script contains the example workflow applying the modelling framework to three different datasets. 



\## Data

Data folder contains three datasets for demonstration of the modelling pipeline:

\- modeled\_data\_annDCgut\_demo.csv



A subset of the original metabolomics dataset from mice colonized with a defined microbial community and germ-free mice. 

\- modeled\_data\_drugs.csv



Preprocessed data from the drug metabolism study (brivudine (BRV), sorivudine (SRV), clonazepam (CLZ)) from Zimmermann, Zimmermann-Kogadeeva, Wegman and Goodman, 2019. 

\- modeled\_data\_Meier.csv



Preprocessed data from Meier et al. study focusing on metabolite distribution along the gut of conventional SPF mice and germ-free mice. 



\## Code



Supporting functions ran from the main workflow:



* fit\_steadystate\_fluxes.m

Fit the model to the metabolites provided in the input table. Produces intermediary model files with solutions for the forward and reverse problems obtained with different optimization algorithms (Interior Point, Trust Region Reflective).



* plot\_and\_analyse\_modeling\_results.m

Reads combined solution from the model intermediary files and makes a clustergram. Output tables for model solutions are produced.



Combines modelling solutions based on the provided selection criteria (e.g. sel\_crit1="IP" means Interior Point algorithm solution, sel\_crit2="LI\_PCC\_within\_high\_total" means solution from the Trust Region Reflective Algorithm selected by the largest Pearson correlation (PCC) between original and restored data (by solving the reverse problem) in the large intestinal section, which is selected among the solutions with the total PCC above 0.7). Produces a clustergram, histograms of correlations between the original and restored metabolite profiles from model and random solutions, and combined model solution files as tables. 



\## Results

\- model\_\*\_intermediary folders contain model solutions from the fit\_steadystate\_fluxes.m script. 



Per each dataset run, three figures are generated: clustergram of model coefficients, and two histograms of correlations between original data and data restored in the reverse problem.

&#x20;

Per each dataset run, two tables are generated: model results for the forward problem, and model results for the reverse (\*reciprocal) problem.



