-------------------------------------------------------------------------

REFERENCE

Repository containing the Matlab code for the scientific article :

Cardio-audio synchronization elicits prediction in auditory sequences during human wakefulness and sleep.

Andria Pelentritou, Christian Pfeiffer, Sophie Schwartz, Marzia De Lucia

https://www.biorxiv.org/content/10.1101/2022.03.03.482861v1

-------------------------------------------------------------------------

DEPENDENCIES

Matlab (version 2019b or later; The MathWorks, Natick, MA)

Fieldtrip (version used 20201205; https://www.fieldtriptoolbox.org/download/)

EEGLAB (version 13.4.4b; https://sccn.ucsd.edu/eeglab/download.php)

-------------------------------------------------------------------------

CONTENTS

Run_ControlAnalyses.m = performs quality control analyses related to experimental design from trigger information and ECG data

Run_ClusterAnalyses.m = performs condition contrast using cluster permutation statistical analyses on preprocessed ERP data

run_clustperm.m = function called in RunClusterAnalyses.m to run the cluster permutation statistical analysis step in fieldtrip

Extract_SlowOscillations.m = detects slow oscillations in EEG data

Extract_CardiacResponse.m = detects omission trial related RR intervals from ECG data

-------------------------------------------------------------------------

INPUT DATA

Input data paths have to be set within each script (will be modified to fit the data structure upon upon peer-reviewed publication).

Processed data will be made publicly available upon peer-reviewed publication.

-------------------------------------------------------------------------

AUTHORS

Author: Andria Pelentritou

PI: Marzia De Lucia

Laboratoire de Recherche en Neuroimagerie (LREN)

Department of Clinical Neurosciences (DNC)

Lausanne University Hospital (CHUV) and University of Lausanne (UNIL)

Mont-Paisible 16, CH-1011 Lausanne, Switzerland

Email: andria.pelentritou@gmail.com

Last updated: July 2022

-------------------------------------------------------------------------
