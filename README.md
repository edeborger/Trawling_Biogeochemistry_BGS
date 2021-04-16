# Trawling_Biogeochemistry_BGS

This repository contains the bottom trawling model used in the article:
"Impact of bottom trawling on sediment biogeochemistry: a modelling approach",
and instructions on how to run it.

The files in this repository are:

- 'fisheries_disturbance_model.R' contains the actual description of the trawling 
disturbance model, which is a wrapper around functions of the CNPDIA package.
This file also contains some example code to produce figures with the output of a 
simulation.
- 'run_simulation.Rmd' is a markdown file that shows how the model can be called
and parametrized, and how output may be visualized.
- /rda contains .rda files which contain output of the model. By loading these files
there is no need to run the models in the .rmd file yourself, as this takes quite 
some time.
- /nutrients contains output of the Copernicus Marine Environmental Monitoring Service
implementation of the ERSEM model, used to force the model.. These are annual bottom
water concentrations representative for locations used in the article.

Please read the materials and methods of the article and references therein to
understand the model parameters and process descriptions.

Please contact the authors (Emil De Borger - emil.de.borger@nioz.nl) for issues and
remarks, or to inquire about collaborations.
