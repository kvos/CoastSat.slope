
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3872442.svg)](https://doi.org/10.5281/zenodo.3872442)
[![Join the chat at https://gitter.im/CoastSat/community](https://badges.gitter.im/spyder-ide/spyder.svg)](https://gitter.im/CoastSat/community)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


# CoastSat.slope
#### Beach slope estimation from satellites
This toolkit is an extension of the main [CoastSat toolbox](https://github.com/kvos/CoastSat), it enables users to estimate the beach-face slope from satellite-derived shorelines and associated tide levels.

Visit the [CoastSat webGIS page](http://coastsat.wrl.unsw.edu.au/) to explore and download a regional-scale dataset of beach slopes.

![](./doc/intro_fig1.jpg)

The methodology is described in: Vos K., Harley M.D., Splinter K.D., Walker A., Turner I.L. (2019). Beach slopes from satellite-derived shorelines. Geophysical Research Letters. https://doi.org/10.1029/2020GL088365 (or preprint [here](https://www.essoar.org/doi/10.1002/essoar.10502903.1)).

Slides are also available [here](https://www.slideshare.net/KilianVos/beach-slopes-from-satellite-shorelines-coast2coast-presentation).

There are two Jupyter Notebooks in the repository, showing examples of beach slope estimation along transects at [Cable Beach](https://github.com/kvos/CoastSat.slope/blob/master/example_slope_Cable_beach.ipynb) and [Narrabeen-Collaroy](https://github.com/kvos/CoastSat.slope/blob/master/example_slope_Narrabeen.ipynb), Australia.

To run the examples you will need to install the `coastsat` environment (instructions in the main [CoastSat toolbox](https://github.com/kvos/CoastSat)).

If you want to use [FES2014](https://www.aviso.altimetry.fr/es/data/products/auxiliary-products/global-tide-fes/description-fes2014.html) global tide model to get the tide levels at the time of image acquisition, follow the [instructions](https://github.com/kvos/CoastSat.slope/blob/master/doc/FES2014_installation.md) provided to setup the model.

**If you like the repo put a star on it!**

Having a problem? Post an issue in the [Issues page](https://github.com/kvos/CoastSat.slope/issues) (please do not email).
