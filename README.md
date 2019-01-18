# R activities for the Evolution of Earth Systems course
These files to aid the R modeling component of the Earth System Evolution course -- Stanford Univ. 2019

## Lecture 0: Intro to R - Tutorial
This is the "Intro to R - Tutorial.Rmd" file. Before opening in RStudio, make sure you have installed the packages:

1. ggplot2
2. learnr
3. tidyverse

Then, open the file in RStudio and click "Run Document" (next to the green arrow right above the script). 

A window should open labeled RStudio at the top of the page. This window contains an interactive tutorial for you to follow. It covers the basics--meant for brand-new beginners to R. There is no quiz and you will not be graded on your performance.

## Lecture 1: Steady state C cycle
This tutorial walks through how to write a simple function in R for the steady state, long-term Carbon cycle. We will eventually derive the transient (time-dependent) solution to this model following the framework of Kump and Arthur, 1999. 

For now, we focus on an exercise inspired by **An Introduction to Isotopic Calculations** by *John Hayes (2004)*. We will derive the equations for the steady state isotopic balance of the carbon cycle, discuss assumptions that are often used, and evaluate the error associated with those assumptions. 

The same packages as those mentioned above (*Lecture 0*) are required to run this tutorial. 

## Lecture 3: C cycle box model
We're done with tutorials for a bit! Here, you will find the R script to run a one-box model of the geologic carbon cycle with formulations from and inspired by Kump and Arthur, 1999 (Chemical Geology) and the CLiBeSO model of Jeremy K Caves Rugenstein.

The R script is *KumpArthur_modified.R*

This should serve as an opportunity to test and troubleshoot changing box model parameters as well as a skeleton for the development of your own box model. The derivations required to arrive at the model equations are described in the PDF file: *C_cycle_BoxModel_deriv.PDF*.
