# Multispectral Image Fusion with PSO Algorithm

Description
----------
This code provides the fusion of PANchromatic (PAN) and MultiSpectral (MS) images using the Particle Swarm Optimization (PSO) algorithm. The steps for fusion is as follows: 

    1) Loading the dataset from its path
    
    2) Pre-processing steps including downsampling and normalization
    
    3) Initialization of PSO algorithm
    
    4) Obtaining the primitive detail map for each spectral band 
    
    5) Extracting the edge detectors for PAN and MS images
    
    6) Fine-tuning the gains of edge detectors using PSO algorithm

![PSO](https://user-images.githubusercontent.com/48659018/56169672-87da0b80-5fa4-11e9-9bad-0eec2ea1dfb8.gif)


Objective Function
--------------
For the purpose of optimizing the gains of edge detectors, the ERGAS metric is minimized. This metric is one of the widely used metrics for the objective evaluation of fusion results. 

Usage
------------
In order to run the code, first you need to specify the path of your dataset.
For example:

    Addpath QuickBird_Data
