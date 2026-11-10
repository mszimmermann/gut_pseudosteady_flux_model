# gut\_pseudosteady\_flux\_model

Standalone modelling pipelines of the project

*Multiomics and quantitative modelling disentangle diet, host, and microbiota contributions to the host metabolome*

By Zimmermann-Kogadeeva, Maria; Henseler, Zachariah M.; Bencivenga-Barry, Natasha A; Gross, Matthias W.; Zimmermann, Michael; Bartmanski, Bartosz J.; Bork, Peer; Goodman, Andrew L.

🖥️**Software requirements:** 

The code is developed with MatLab2023a Optimization toolbox and requires MatLab version 2023a, although most scripts should be compatible with MatLab version 2019.

❕**License:**

The code is distributed under the terms of the GNU General Public License (please read copyright\_and\_license and LICENSE files for details).

☑️**Installation:**

There is no specific installation required to run the code.

📖 **Contents:**

Workflow for coarse-grain intestinal flux modelling of metabolite abundances along the gastrointestinal tract assuming pseudo-steady state

_workflow\_beta.m_ script contains the example workflow applying the modelling framework to three different datasets. 

📚 **Data:**

Data folder contains three datasets for demonstration of the modelling pipeline:

- modeled\_data\_annDCgut\_demo.csv

A subset of the original metabolomics dataset from mice colonized with a defined microbial community and germ-free mice. 

- modeled\_data\_drugs.csv

Preprocessed data from the drug metabolism study (brivudine (BRV), sorivudine (SRV), clonazepam (CLZ)) from Zimmermann, Zimmermann-Kogadeeva, Wegman and Goodman, 2019. 

- modeled\_data\_Meier.csv

Preprocessed data from _Meier et al. 2023_ study focusing on metabolite distribution along the gut of conventional SPF mice and germ-free mice. 

🖱️ **Code:**

Supporting functions ran from the main workflow:

- fit\_steadystate\_fluxes.m

Fit the model to the metabolites provided in the input table. Produces intermediary model files with solutions for the forward and reverse problems obtained with different optimization algorithms (Interior Point, Trust Region Reflective).

- plot\_and\_analyse\_modeling\_results.m

Reads combined solution from the model intermediary files and makes a clustergram. Output tables for model solutions are produced.

Combines modelling solutions based on the provided selection criteria (e.g. sel\_crit1="IP" means Interior Point algorithm solution, sel\_crit2="LI\_PCC\_within\_high\_total" means solution from the Trust Region Reflective Algorithm selected by the largest Pearson correlation (PCC) between original and restored data (by solving the reverse problem) in the large intestinal section, which is selected among the solutions with the total PCC above 0.7). Produces a clustergram, histograms of correlations between the original and restored metabolite profiles from model and random solutions, and combined model solution files as tables. 

🏁 **Results:**

- _model\_\*\_intermediary_ folders contain model solutions from the _fit\_steadystate\_fluxes.m_ script. 

- Per each dataset run, three figures are generated: clustergram of model coefficients, and two histograms of correlations between original data and data restored in the reverse problem.

- Per each dataset run, two tables are generated: model results for the forward problem, and model results for the reverse (\*reciprocal) problem.
