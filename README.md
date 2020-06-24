biosurvey: Tools for Biological Survey Planning
================
Claudia Nunez-Penichet, Marlon E. Cobos, A. Townsend Peterson, Narayani
Barve, Vijay Barve, Tomer Gueta

  - [Project description](#project-description)
  - [Package description](#package-description)
  - [Installing the package](#installing-the-package)

<br>

**This repository is for the project “Biological Survey Planning
Considering Hutchinson’s Duality” developed during the program GSoC
2020.**

<br>

## Project description

Student: *Claudia Nuñez-Penichet*

GSoC Mentors: *Narayani Barve, Vijay Barve, Tomer Gueta*

Motivation:

Given the increasing intensity of threats to biodiversity in the world,
one of the challenges in biodiversity conservation is to complete
inventories of existing species at distinct scales. Species
distributions depend on the relationships between accessible areas,
environmental conditions, and biotic interactions. As planning a survey
system only aims to register species in a region, biodiversity
interaction can be overlooked in this case. However, the relationship
between environmental conditions and the geographic configuration of an
area is of crucial importance when trying to identify key sites for
biodiversity surveys. Among the diverse packages in R for selecting
survey sites, such considerations are not implemented and are limited to
a random selection of sampling sites or analyses that allow detecting
potential sampling sites based on the environmental similarity between
sampled and unsampled areas. Given the need for more solutions, the
**biosurvey** package aimed for considering the relationship between
environmental and geographic conditions in a region when designing
survey systems that allow sampling of most of its biodiversity.

<br>

## Package description

The biosurvey R package implements multiple tools to allow users to
select sampling sites increasing efficiency of biodiversity survey
systems by considering the relationship of environmental and geographic
conditions in a region.

<br>

## Installing the package

biosurvey is in a GitHub repository and can be installed and/or loaded
using the code below (make sure to have Internet connection). If you
have any problem during installation, restart R session, close other
RStudio essions you may have open, and try again. If during the
installation you are asked to update packages, do so if you don’t need a
specific version of one or more of the packages to be installed. If any
of the packages gives an error when updating, please install it alone
using install.packages(), then try re-installing biosurvey again.

``` r
# Installing and loading packages
if(!require(remotes)){
  install.packages("remotes")
}
if(!require(biosurvey)){
  remotes::install_github("claununez/biosurvey")
  library(biosurvey)
}
```
