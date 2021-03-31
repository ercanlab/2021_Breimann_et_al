<div align="center">
  
# Documentation for Breimann et al. 2021

</div>

_**Dynamic of condensin DC binding is regulated by a noncatalytic activity of the H4K20
demethylase DPY-21**_

 
Laura Breimann *, Ana Karina Morao *, Jun Kim, David Sebastian Jimenez, Nina Maryn, Krishna Bikkasani, Michael J Carrozza, Sarah E Albritton, Maxwell Kramer, Lena Annika Street, Kustrim Cerimi, Vic-Fabienne Schumann, Ella Bahry, Stephan Preibisch, Andrew Woehler, Sevinç Ercan 

bioRxiv:  ; doi: 

\* equal contribution  


### Content

* _**1.	Abstract**_
* _**2.	Requirements FRAP Matlab analysis**_
* _**3.	FRAP analysis tutorial**_
* _**4.	RNA-seq  analysis**_

<br />
<br />

<div style="text-align: justify">
 
### 1.	Abstract 
Abstract:
Condensin is a multi-subunit SMC complex that binds to and compacts chromosomes. Unlike cohesin, in vivo regulators of condensin binding dynamics remain unclear. Here we addressed this question using C. elegans condensin DC, which specifically binds to and represses transcription of both X chromosomes in hermaphrodites for dosage compensation.  Mutants of several chromatin modifiers that regulate H4K20me and H4K16ac cause varying degrees of X chromosome derepression. We used fluorescence recovery after photobleaching (FRAP) to analyze how these modifiers regulate condensin DC binding dynamics in vivo. We established FRAP using the SMC4 homolog DPY-27 and showed that a well-characterized ATPase mutation abolishes its binding. The greatest effect on condensin DC dynamics was in a null mutant of the H4K20me2 demethylase DPY-21, where the mobile fraction of the complex reduced from ~30% to 10%. Supporting a structural role for DPY-21, a catalytic mutant of dpy-21 did not regulate condensin DC mobility. Intriguingly, unlike FRAP, ChIP-seq and Hi-C detected little change in condensin DC binding and 3D DNA contacts in the dpy-21 null mutant. Together, our results suggest that DPY-21 has a non-catalytic role in regulating the dynamics of condensin DC binding, which is important for its function in transcription repression.
 

<br />
<br />

### 2.	Requirements for FRAP Matlab analysis

The analysis script was developed and tested in Matlab R2018a on Mac OS 10.15.7.

The following Matlab toolboxes are required to run “FRAP_analysis.m”: 
* curve_fitting_toolbox
* image_toolbox


The following scripts have to be in the same folder as “FRAP_analysis.m”

* tiffread2.m (by Francois Nedelec)
* struct2.csv.m (by James Slegers)
* matVIS.m (by S. Junek)
* dftregistration.m (by Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup)
* timesteps.m
* timestepsArray.m 


<br />
<br />


### 3.	FRAP analsysis tutorial

FRAP protocol: 

A FRAP dataset of a _C. elegans_ intestine nuclei can be found [here](https://github.com/ercanlab/2021_Breimann_et_al/tree/main/FRAP_analysis/Example_datasets )

_**Step 1: Select the data for analysis**_

The first step after running the script is to select the input folder with the raw FRAP images and an output folder for the analysis files. A window will pop up, and you can navigate to the respective folders (first input, then output). 

```Choose data_tif directory```
```Choose matlab_results directory``` 

Next, you need to select the FRAP dataset you want to analyze. The Leica SP8 creates two image stacks per FRAP experiment, one before the bleach point and one after the breakpoint. First, select the dataset before the bleach and then click on the post bleach image stack. 

```Load pre-stack```
```Load post-stack ```

_**Step 2: Manually outline the cell nucleus**_

To select only one nucleus for FRAP analysis, you can draw a **ROI**. For that, a window will open with a filtered post bleach image (mean of the first three images) and just start outlining the nucleus by clicking in the image. Once you are satisfied with the outline, double click in the middle of the selected ROI to accept it. 

<img src="https://github.com/ercanlab/2021_Breimann_et_al/blob/main/FRAP_analysis/Screenshots/nuc.png" alt="C. elegans intestine cell" width="400">

<img src="https://github.com/ercanlab/2021_Breimann_et_al/blob/main/FRAP_analysis/Screenshots/filtered_nuc.png" alt="C. elegans intestine cell with ROI" width="400">

_**Step 3: Automatically detect the bleach point**_

The next step automatically detects the bleach-point by subtracting the mean of the first couple post-bleach image from the mean of the last pre-bleach images. A window will appear with the pre-bleach and post-bleach images and the difference between the two images. The lower row depicts the selected mask from the previous step (thresholded) and the mask for the bleach point based on Fdiff above. 

At this point, there is the option to change the threshold for the bleach point selection. Simply press ```No``` in the second window and write a value between ```0-1``` in the ```Command Window``` (the starting point is 0.6), and press enter. If you are happy with the bleachpoint detection, press ```Yes```. The displayed overview image is saved as ```_mask.tif``` to the previously selected output folder. 


<img src="https://github.com/ercanlab/2021_Breimann_et_al/blob/main/FRAP_analysis/Screenshots/mask.png" alt="Masks for the nuclei and bleachpoint" width="400">

<img src="https://github.com/ercanlab/2021_Breimann_et_al/blob/main/FRAP_analysis/Screenshots/bleachpoint.png" alt="Window to select the bleachpoint" width="400">

_**Step 4: Inspect the results**_

The following steps are automatically executed and will save the results to the previously selected output folder. 

For correction acquisition bleaching, the rest of the nucleus is used. The file ```_bleaching_correction.``` shows the detected acquisition bleaching in the nuclear mask region for all post bleach images. The bleaching is fitted and used to correct the background fluorescence as well as the recovery fluorescence. The last timepoint intensity of the pre bleach detection in the bleach point is set to 1 to achieve comparable datasets between different replicates. 

<img src="https://github.com/ercanlab/2021_Breimann_et_al/blob/main/FRAP_analysis/Screenshots/DPY-27_bleaching_correction.tif" alt="Plots for the correction of the acquisition bleaching" width="400">


The fitted and normalized recovery curve is saved as ```_recovery.tif``` to the output folder. It displays the normalized fluorescence in the whole nucleus (red) and the bleach point recovery (black). The immobile and mobile fractions are displayed in the image, and the tau and t-half values from the fit of the curve. The fitting of the curve can only be changed directly in the code (Section 7). 

<img src="https://github.com/ercanlab/2021_Breimann_et_al/blob/main/FRAP_analysis/Screenshots/DPY-27_recovery.tif" alt="Plot for the FRAP recovery" width="400">


To check how well the fit describes the observation, a set of goodness of fit values is saved in the ```_gof.csv``` file in the results folder. It contains different statistics: 
* The sum of squares due to error (SSE) (values closer to 0 are good)
* R-square (values closer to 1 are good)
* Degrees of Freedom (DFE) 
* Adjusted R-square (values closer to 1 are good)
* Root mean squared error (RMSE) (values closer to 0 are good)
These values can be used to select the best fit or filter data.  


<img src="https://github.com/ercanlab/2021_Breimann_et_al/blob/main/FRAP_analysis/Screenshots/GOF.png" alt="Table for the GOF statistics" width="400">

For further analysis and averaging of different experiments, the normalized values for the FRAP curve, and the tau value and percent of the immobile and mobile fractions are saved to the file ```_pyan.txt```.  The data structure is as follows: the first value is the tau-value, then the mobile fraction and the immobile fraction. From the 4th value on, the normalized FRAP recovery values are listed. 

<img src="https://github.com/ercanlab/2021_Breimann_et_al/blob/main/FRAP_analysis/Screenshots/results.png" alt="Table of bleach valuess" width="400">


The t-half value is recorded by two different means. Firstly from the fit and the corresponding value is saved as ```_t_half_value_from_fit.txt```. The more direct way is to get the half-time of recovery from the recovery data. A visual representation of this can be found in the image ```_thalf_no_fit.tif``` and the estimated value in ```t_half_value_no_fit.txt```  


<img src="https://github.com/ercanlab/2021_Breimann_et_al/blob/main/FRAP_analysis/Screenshots/DPY-27_thalf_no_fit.tif" alt="Plot for the estimation of t-half" width="400">




<br />
<br />



### 4.	RNA-seq  analysis

The R script used for RNA-seq analysis can be found [here](https://github.com/ercanlab/2021_Breimann_et_al/tree/main/RNAseq_analysis )





<br />
<br />
