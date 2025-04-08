## Documentation is incomplete! (WIP)

# CINsim

## Installation
To install, use `devtools::install_github("Kaam-umcg/CINsim")`. Alternatively, download the repository into a folder, and install directly within the console.

## CINsim wrapper function
The main wrapper function in this package is `Cinsim()`, which can take many (optional) parameters as input. They are outlined below.

### Main parameters
* `karyotypes` a matrix with karyotypes (cells in rows, chromosomes in columns) actings as the founder population. If none is given, by default, then a founder population of 1 diploid mouse cell is used  (default = NULL).
* `euploid_ref` sets the euploid chromosome copy number, most importantly used to calculate aneuploidy score (default = 2).
* `g` the maximum number of generations (cycles) the simulation should run before it stops (default = 12)
* `max_num_cells` the maximum number of cells allowed (extrapolated from downsampled population) before the simulation stops (default = 2e+09).
* `pMisseg` sets the per chromatid probability of mis-segregation per mitosis (default = 0.0025).
* `pMissegG` a vector of generation/cycle numbers during which mis-segregation occurs. This parameter can be used to introduce periods of mis-segregation rather than continuous CIN (default = NULL).
* `pWGD` sets the probability of whole genome doubling (WGD) per cells (bypasses mitosis/mis-segregations) (default = 0).
* `pDivision` sets the baseline probability of division. Set to a value less than 1 for asynchronous cell division (default = 1).
* `copy_num_boundaries` sets the range (minimum and maximum) of viable copy numbers; karyotypes with copy numbers outside this range will die (default = 1-8).
* `qMods` a vector of length 3 to modify the internally calculate values of pDivision, pMisseg, and pSurvival respectively.

### Selection parameters
The following parameters relate to modes and strength of karyotype selection:
* `selection_mode` can be any of three possible modes: `"cn_based"`, `"rel_copy"`, or `"davoli"`. Default is `NULL` which yields no selection.

`"cn_based"` is the main mode of selection used in our study, relying on a copy number matrix that defines the relative fitness of chromosome copy number states (see below @ `selection_metric`).

`"rel_copy"` defines fitness on the degree of copy number deviation from the population median ploidy. Greater deviation from the modal ploidy will result in greater decreases in fitness.

`"davoli"` is based on the TSG-OG chromosome dosage scores to simulate the balance of tumour suppressor genes and oncogenes as a metric for oncogenicity of chromosomes (based on Davoli et al., 2013).

* `selection_metric` defines the selection metric used for the respective selection mode, such as a copy number matrix for the `"cn_based"` mode.

* `coef` a list of length three with karyotype fitness coefficients (default = NULL). These coeffecients can be calculated using the `make_cinsim_coefficients` function (see below for more details).

### Optional/additional parameters
The remaining parameters have not been explored in the main study, but could be of interest to users:

* `fit_misseg` determines whether fitness proportionally affects the rate of mis-segregation. More fit cells will mis-segregate less often (default = FALSE).
* `fit_division` determines whether fitness proportionally affects the rate of division. More fit cells will divide more often (default = FALSE).
* `chrom_weights` a vector defining the weights of individual chromosomes, if you would like specific chromosomes to contribute more the fitness metric than others (default = FALSE).
* `max_monosomy` the maximum number of monosomic chromosomes allowed before a cell dies (default = NULL).
* `min_euploid` the minimum number of chromosomes that must remain euploid before a cell dies (default = NULL).

### Systemic parameters
The following parameters are not related to biological concepts but technical aspects of the simulation.

* `down_sample` the maximum simulated population size (agents) before the population is down-sampled (default = 50.000).
* `down_sample_frac` the fraction of the population that is randomly selected for the next generation if down-sampling occurs (default = 0.25).
* `collect_fitness_scores` a boolean whether fitness scores are collected over time. This is a time-consuming process and can slow down the simulation rate significantly at larger scales (default = FALSE).

## Setting up a simulation (without karyotype-based selection).
To run a simple simulation, just call the `"Cinsim()"` function to create an object in which to store the results, and specify the parameters of interest. In this example we run a simulation with a fair amount of CIN for 25 cycles and no selection:

`sim_res <- Cinsim(pMisseg = 0.0025, g = 25, selection_mode = NULL)`

Depending on your computer's processing power, it may take a few seconds or minutes to complete the simulation. A message is displayed to communicate the current stage of the simulation run, including the amount of time it took to complete a step. To quickly inspect the resulting karyotype landscape call the `"cnvHeatmap()"` function on the results object:

`cnvHeatmap(sim_res)`

By default this function will sample 1,000 cells from the population to display. More cells can be displayed if desired at the expense of time and memory. We found that 1,000 cells are usually enough to get a good sense of the karytoypes. For a general summary of the copy number frequencies by chromosome at the end of the simulation, use the `"plot_cn()"` function on the results:

`plot_cn(sim_res)`

Alternatively, the copy number frequencies can be shown over time for each chromosome separately:

`plot_cn(sim_res, final_g = FALSE)`
