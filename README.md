# Productivity interacts with diversification rate in determining species richness and trait diversity of tetrapods in a global hotspot
Data from the article "Productivity interacts with diversification rate in determining species richness and trait diversity of tetrapods in a global hotspot"

## Folder Structure

The data used is available in the "data_from_figshare.zip" file and, once accepted, will be released for download on figshare.
Here is the folder structure of the project, along with a brief description of each folder:

- **01_processingdata**: This folder contains raw data that requires preprocessing for trait names, species composition, and phylogeny. The .qmd file details each step of the processing.

- **02_metrics**: In this folder, you will find the results of functional diversity metrics obtained with the data processed in the previous step.

- **03_sem**: Here are the structural equation analyses performed, generated plots, and examined spatial correlations.

- **Shapefiles**: This folder contains files used to obtain the presence/absence matrix (please note that it may be outdated due to size constraints).

- **assemblage_age**: The files in this folder are used to calculate diversification rate and assemblage age metrics using the Herodotools package. An additional script is required for ancestral range reconstruction using biogeobears.

**NOTE**: The "data_preparation.Rmd" file is out of date, to run the analysis flow go through the folders described above.