<b><p><h1 style="color:red;font-size:20px;">Project Description</b></p></h1>

Species-tailored sampling guidelines remain the most efficient method to conserve genetic diversity ex situ: a study on threatened oaks. Repository storing code, simulation, and parameter files that represent 14 IUCN Red List endangered oaks. This project is being completed in collaboration by Kaylee Rosenberger, Emily Schumacher, Dr. Sean Hoban, and Dr. Alissa Brown from the Morton Arboretum.
Repository storing code, simulation, and parameter files that represent IUCN Red List endangered oaks for my honor's capstone project at Northern Illinois University. This project is led by Kaylee Rosenberger and the code for most analyses is in her repo (found here: https://github.com/kayleejorose/Quercus_IUCN_samp_sims) but this folder contains contains code used to create maps for the final project, as well as first drafts of genetic diversity code implemented in the final work. 

<b><p><h1 style="color:red;font-size:20px;">Overview</b></p></h1>
The overall aim of this project is to contribute to practical seed sampling guidelines for creating and maintaining genetically diverse collections for botanic garden and arboreta. Informing these sampling guidelines is one way to ensure a genetically representative sample is obtained from wild populations. Prior work has found that it is important to consider species' traits like dispersal, mode of reproduction, population history, and more, when sampling from wild populations. For this project, we focused on the genus Quercus (oaks) and studied 14 species that the IUCN Red List describes as vulnerable. Oaks are keystone species for many environments and have a high ecological importance. In addition, oaks cannot be seed banked using traditional methods, so they must be conserved through living collections. Since maintaining living collections requires extensive space and energy for gardens, creating efficient collections that represent the diversity of wild oak populations is extremely important. Thus, creating and maintaining genetically diverse collections in botanic gardens and arboreta is essential for the future survival and restoration of these rare, endangered species, and this can be achieved through proper sampling techniques.

Here, we aim to determine if closely-related species (within the same genus) with similar biology, dispersal, and life history traits but different ranges and population sizes would have similar minimum sample sizes.

There were three main goals in this project:

Determine minimum sample sizes to capture 95% diversity for 14 threatened oaks
Determine if one minimum sample size can fit all 14 oaks studied
Determine the extent to which parameter values chosen for simulation impact the minimum sample size

<b><p><h1 style="color:red;font-size:20px;">Project Summary</b></p></h1>
We chose 14 species of oaks from the IUCN Red List of Endangered Species in the US. Each of these oaks the IUCN Red List describes as vulnerable. We created species-tailored parameter values to represent each species realistically in simulation, using the software fastsimcoal. In addition, we run alternative simulations for each species, varying the parameter value that we had the least amount of confidence in, so that we can account for estimations within our parameter values. We then created scripts with R that represent sampling from the simulated populations. Here, we tested a broad entire range of sampling for each species--from one individual, to 500 individuals. From this, we determine the minimum sample size required to capture 95% of the species’ genetic diversity (a common threshold for sufficient genetic diversity). With this data, we aim to recommend a minimum sample size to capture sufficient genetic diversity for each of these species, which would be directly useful to botanic gardens and arboreta. Furthermore, we aim to determine whether one minimum sample size can be recommended to sufficiently capture the diversity of all of these vulnerable oaks.

We also determine if changing the parameter values used for simulation significantly impacts the minimum sample size required, by creating 'alternative' simulations. In the alternative simulations for each species, we varied a parameter value that we had the least amount of confidence in, so that we can account for estimations within our parameter values and determine how much the parameters chosen for simulation impact the minimum sample size.

<b><p><h1 style="color:red;font-size:20px;">Analyses</b></p></h1>
We ran a generalized linear model (GLM) to determine the minimum sample sizes needed to adequately capture genetic diversity for each of the 14 oak species. We calculated pairwise contrasts using the emmeans package (Lenth, 2021) to determine whether the difference between the number of alleles captured between the two simulation types (original and alternative) was significant for each taxon.

<b><p><h1 style="color:red;font-size:20px;">Files Types</b></p></h1>

Parameter files:
<p>.par .txt</p>
<p>Edited in text editor Notepad++</p>
<p>These are input to the software fastsimcoal to create genetic datasets The .par signifies parameter files. They contain information to create the genetic datasets via a coalescent simulation, including population sizes and migration rates. Parameter files are written in the text editor Notepad++
You can run fastsimcoal through the command line with the prompt: fsc26 <\file_name> -g 1 -n 1000
Here we used 1 for the genetic data type of all simulations, representing diploid individuals and 1000 for 1000 simulation replicates.</p>

Genotype files:
<p>.arp .gen</p>
<p>Arlequin files are the outputs of the empirical results and simulations - .arp - which are then converted to genepop files (.gen) which are then read into R using adegenet as genind objects. </p>

Rscripts:
<p>.R .Rdata</p>
<p>For this project, R scripts were used to import .arp files into R for conversion to .gen files through adegenet package, convert .gen files to genind objects through adegenet package, analyze data through functions associated with the adegenet package. </p>


<b><p><h1 style="color:red;font-size:20px;">Directory Contents</b></p></h1>

<ul><li>Analyses</li></ul>
<ul><ul><li>Description: All input files utilized for generating butternut's species distribution model, habitat suitability maps, and projecting the habitat suitability maps into the past.</li></ul></ul>
<ul><ul><li>R_Scripts</li></ul></ul>
<ul><ul><ul><li>Reduce_Points_Mapping</li></ul></ul></ul>
<ul><ul><ul><li>IUCN_redlist_oaks_empirical_gendiv_summary_stats</li></ul></ul></ul>
<ul><ul><ul><li>IUCN_redlist_oaks_simulation_gendiv_summary_stats</li></ul></ul></ul>
<ul><ul><li>Results</li></ul></ul>
<ul><li>Data_Files</li></ul>
<ul><ul><li>Description: Genepop objects generated from empirical genetic diversity data collection on <i>Quercus acerifolia</i>, <i>Quercus hinckleyii</i>, <i>Quercus pacifica</i>, and <i>Quercus tomentella</i>. </li></ul></ul>
<ul><ul><li>QUAC</li></ul></ul>
<ul><ul><li>QUHI</li></ul></ul>
<ul><ul><li>QUPA</li></ul></ul>
<ul><ul><li>QUTO</li></ul></ul>

<ul><li>Figures</li></ul>

<ul><li>Mapping</li></ul>
