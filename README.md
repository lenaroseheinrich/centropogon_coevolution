# centropogon_coevolution

Datasets and code for the manuscript:

**Three sympatric species of hummingbird-pollinated *Centropogon*: A match for coevolution?**  
Asia Kaiser, Lena Heinrich, Sara Sofia Pedraza, Sebastian Forward, Raj Prasai, J. Edgardo Arévalo, Mauricio Fernández Otárola

---

## Short description of methods

This repository contains code and data for analyzing co-distributions of three hummingbird-pollinated *Centropogon* species and their pollinators across Costa Rica and Panama. We perform species distribution modeling (SDM) using occurrence data from GBIF and environmental variables, and quantify geographic overlap among species pairs. The workflow includes cleaning the occurrence points, performing species distribution modeling, visualizating the predicted distributions, and calculating range overlap percentages between species pairs.

---

## Files

**1. `data_cleaning.R`**  
Cleans and standardizes GBIF occurrence data for six focal species (*Centropogon* and hummingbirds). Removes duplicates, filters out ocean points, and saves a cleaned dataset (`clean.csv`) for use in modeling.

**2. `distribution_modeling.R`**  
Runs ensemble species distribution models (SDMs) using the cleaned occurrence data and climate predictors. Exports continuous habitat suitability maps (.tif), binarized presence/absence rasters, and polygon shapefiles for visual and quantitative overlap analyses.

**3. `plot_distributions.R`**  
Visualizes SDM outputs by plotting predicted distributions for each plant-pollinator pair. Calculates pairwise range overlap as a percentage of each species' total range and generates summary barplots and tables.

**4. `clean.csv`**  
Cleaned occurrence data for all six species, used as input to the SDM pipeline.

---

## Citations for GBIF datasets

- *Centropogon valerii*: GBIF.org (26 May 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.dumn22  
- *Centropogon talamancensis*: GBIF.org (26 May 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.rupmw9  
- *Centropogon costaricae*: GBIF.org (26 May 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.sytkyy  
- *Colibri cyanotus*: GBIF.org (26 May 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.gmxwp9  
- *Eugenes spectabilis*: GBIF.org (26 May 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.bjcvw6  
- *Panterpe insignis*: GBIF.org (26 May 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.454p9u  
