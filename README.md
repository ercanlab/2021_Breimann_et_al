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
* _**2.	Installation of FRAP Matlab analysis**_
* _**3.	FRAP analysis tutorial**_
* _**4.	RNA-seq  analysis tutorial**_

<br />
<br />

<div style="text-align: justify">
 
### 1.	Abstract 
Abstract:
Condensin is a multi-subunit SMC complex that binds to and compacts chromosomes. Unlike cohesin, regulators of condensin dynamics and function in vivo remain unclear. Here we addressed this question using C. elegans condensin DC, which specifically binds to and represses transcription of both X chromosomes in hermaphrodites for dosage compensation.  Mutants of several chromatin modifiers that increase and decrease H4K20me and H4K16ac on the X cause varying degrees of X chromosome derepression. We used fluorescence recovery after photobleaching (FRAP) to analyze how these mutations regulate condensin DC binding dynamics in vivo. We established the system using the SMC4 homolog DPY-27 and showed that a well-characterized ATPase mutation abolishes its binding as measured by FRAP and ChIP-seq. The greatest effect on condensin DC dynamics was in a null mutant of the H4K20me2 demethylase DPY-21, where the mobile fraction of the complex reduced from ~30% to 10%. Supporting a structural role for DPY-21, a catalytic mutant of dpy-21 did not regulate condensin DC mobility. Intriguingly, unlike FRAP, ChIP-seq and Hi-C detected little change in condensin DC binding and 3D DNA contacts in the dpy-21 null mutant. Together, our results suggest that DPY-21 has a non-catalytic role in regulating the dynamics of condensin DC binding, which is important for its function in transcription repression.
 

<br />
<br />

### 2.	Installation of FRAP Matlab analysis

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

Test dataset: 

Step 1: 

Select your data input directory and output directory

Step 2:



<br />
<br />



### 4.	RNA-seq  analysis tutorial







<br />
<br />
