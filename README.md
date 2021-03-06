# Multispectral Image Fusion with PSO Algorithm

[An adaptive multispectral image fusion using particle swarm optimization
](https://ieeexplore.ieee.org/abstract/document/7985325) is the paper for this MATLAB code.


Description
----------
This code provides the fusion of PANchromatic (PAN) and MultiSpectral (MS) images using the Particle Swarm Optimization (PSO) algorithm. The steps for fusion is as follows: 

    1) Loading the dataset from its path
    
    2) Pre-processing steps including downsampling and normalization
    
    3) Initialization of PSO algorithm
    
    4) Obtaining the primitive detail map for each spectral band 
    
    5) Extracting the edge detectors for PAN and MS images
    
    6) Fine-tuning the gains of edge detectors using PSO algorithm



Loss Function
--------------
For the purpose of optimizing the gains of edge detectors, the ERGAS metric is minimized. This metric is one of the widely used metrics for the objective evaluation of fusion results. 

Usage
------------
First you need to specify the path of your dataset.
For example:

    addpath QuickBird_Data
The Main_PSO.m is the main framework of the proposed fusion framework. The _**pre-processing**_ steps as well as the obtaining _**fusion outcome**_ is put into this M-file. The ERGAS_Index.m and ERGAS.m files are used for the purpose of optimization. The optimized gains of _**edge detectors**_ are computed as the output of of PSO algorithm. 

To _run_ the code, in the _command window_ use this: 

    Main_PSO.m

Objective Evaluation
----------
For objective assessment of the fusion result, first add the path of objective metric. For example: 

    addpath Objective_Evaluation

Sample Output
-----------
    
    The MS, PAN and pansharpened result of the QuickBird dataset from Sandarbans, Bangladesh. 

![LRMS](https://user-images.githubusercontent.com/48659018/56171542-5284ec00-5fab-11e9-93fb-a973ba1e8014.png)
![PAN](https://user-images.githubusercontent.com/48659018/56171559-603a7180-5fab-11e9-8626-c1103ca22e6d.png)
![Pansharpened](https://user-images.githubusercontent.com/48659018/56171570-6892ac80-5fab-11e9-86be-8f86e797e974.png)

Reference
--------
If you find this code helpful, please cite this paper: 

    A. Azarang and H. Ghassemian, "An adaptive multispectral image fusion using particle swarm optimization," in Proc. Iranian Conf. Elec. Eng. (ICEE), May 2017, pp. 1708-1712.

Contact
--------
If you have any question regarding the paper, codes and so on, do not hesitate to contact me. 

Arian Azarang  azarang@utdallas.edu
