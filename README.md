# Mitometer
by Austin E.Y.T. Lefebvre

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4265280.svg)](https://doi.org/10.5281/zenodo.4265280)  
![GitHub](https://img.shields.io/github/license/aelefebv/Mitometer)

`Mitometer` is a MATLAB App containing tools for automated mitochondrial segmentation and tracking, and allows for visualization and preliminary data analysis of exportable quantitative data of individual mitochondrion motility and morphology features.

A version of this README is available with a step-by-step image walkthrough in the MitometerManual PDF after download.

# Installation
## DOWNLOAD
The Mitometer program is written in MATLAB (MathWorks). The MATLAB GUI Mitometer app and corresponding source code is freely available online at https://github.com/aelefebv/Mitometer.
## APP EXECUTION
1.	Open the MATLAB App Installer file titled “Mitometer”.  
2.	Once MATLAB opens, click “Install”.  
3.	Click the “Mitometer” app in the Apps menu.  
Note: Installation typically takes no more than 1 minute.  
## SOURCE CODE
All functions used in the Mitometer app are included in the Source Code folder. This includes the app’s App Designer GUI, titled “GUI”.

# UI Tutorial
## HOME TAB
An instance of Mitometer can take as an input either 2D or 3D images: 
 
### Start 2D: 
1.	Input pixel size and time between each temporal frame.
2.	Upload one or more Tif image stacks with temporal frames as the stack.  

Runs segmentation and tracking on each file in order of selection. 
Uploaded files will be shown in the right column under “Select file(s)”  
Note: You can analyze additional 2D files by pressing on Add 2D and following steps 1 and 2 again. 

### Start 3D: 
1.	Input pixel size, time between each temporal frame, number of spatial planes, and axial distance between the spatial planes.  
2.	Upload one or more Tif image stacks with the temporal and spatial frames as the stack. Order should follow z -> t, i.e. all sequential (bottom to top) z-planes for time 1, followed by all sequential z-planes for time 2 (Note: This is the default order for LSM files).  

Runs segmentation and tracking on each file in order of selection.  
Uploaded files will be shown in the right column under “Select file(s)”  

Note: You can analyze additional 3D files by pressing on Add 3D and following steps 1 and 2 again.   

### Reset: 
Reloads Mitometer.
 
## SEGMENTATION TAB
After the processing from the home page is finished, visualization of mitochondria segmentation for each individual file will be available.
 
### Visualize: 
1.	Select the file you wish to visualize under the “Select file(s)” menu.
2.	Press the “Select highlighted” button.  

The software allows visualization of the uploaded image, the image after diffuse background removal, randomly colored connected components (mitochondria), and the final segmented image.   

To scrub through frames of the images:  
1.	Select the image you wish you view under “Selected Image”   
2.	Use the frame spinner, slider, or play button to view temporal frames.  

In 3D, you can view different spatial frames by changing the spinner or slider next to “Z Plane”.  

### Save images:
To save the image stack to a tif file:  
1.	Select the image you wish to save under “Selected Image”.   
2.	Press the “Save .tif stack” button.  

To save a single spatiotemporal frame:  
1.	Hover over the image you wish to save.  
2.	Use the MATLAB menu that appears above the image.  

## TRACKING TAB
After the processing from the home tab is finished, visualization of mitochondria tracks for each individual file will be available.
 
### Visualize:
1.	Select the file you wish to visualize under the “Select file(s)” menu.  
2.	Press the “Select highlighted” button.  

The software allows visualization of all tracks, confident tracks, a single track (one track at a time), perinuclear tracks, and telenuclear tracks (as defined by below and above the Otsu’s threshold, respectively, of the histogram of mitochondrial distances from the mitochondrial aggregate center of mass).  
In 3D, you can change the aspect ratio of the Z plane in respect to the X-Y plane by changing the spinner next to “Z Aspect Ratio”.  
To remove tracks containing less than a specified number of points, change the spinner next to “Track length threshold”.  
When “Single track” is selected, you can change the selected track by changing the spinner value next to “Track to view”.  
Track starting points are shown in blue, track ending points are shown in red, fission events are shown in purple, and fusion events are shown in green.  
Outlines of mitochondria at every frame can be visualized by enabling the “Extrema” checkbox.  

### Save tracks:
To save the track data (refer to Section 3. Output, for information on the saved data):  
1.	Select the type of track you wish to save under “Track View” (default is all tracks).
2.	Select a track length threshold (default is 1).  

To save specific file(s):
1.	Select one or more files you wish to save under the “Select file(s)” menu.
2.	Press the “Select highlighted” button.
3.	Press the “Selected” button next to “Save tracks”.

To save all files:
1.	Press the “All” button next to “Save tracks”.  

Track data will be saved in a .mat file.

## ANALYSIS TAB
After the processing from the home tab is finished, preliminary analysis of basic mitochondrial morphological, motility, and dynamic features between files and conditions will be available.
 
### Visualize:
1.	Select the file(s) whose features you wish to visualize under the “Select file(s)” menu.  
a.	E.g. select all files of your control group.
2.	Press the “Select highlighted” button.
3.	Press the “Create group” button.
4.	Select whether you wish to analyze all tracks, confident tracks, perinuclear tracks, or telenuclear tracks (default is all tracks).
5.	Select a track length threshold (default is 1).
6.	Name your group.
7.	Repeat 1-6 for each condition/group you wish to compare.
8.	Select the first group to compare under the “Groups” dropdown.
9.	Select the feature you wish to analyze.   
a.	Note: Z Axis Length is only available in 3D.
10.	Press “Add to graph” to visualize the data.
11.	Repeat steps 8-10 for each condition/group you wish to compare.

The software allows visualization of either histogram distributions or bar graphs of median values of normal or logarithmically transformed data.  
Statistical comparison between groups is done via a One-way ANOVA.  
Clear off all groups from the graph by clicking the “Reset graph” button.

### Save feature data:
To save the feature data (refer to Section 3. Output, for information on the saved data):
1.	Select the file(s) whose features you wish to save under the “Select file(s)” menu.  
a.	E.g. select all files of your control group.
2.	Press the “Select highlighted” button.
3.	Press the “Create group” button.
4.	Select whether you wish to save all tracks, confident tracks, perinuclear tracks, or telenuclear tracks (default is all tracks).
5.	Select a track length threshold (default is 1).
6.	Name your group.
7.	Select the group to save under the “Groups” dropdown.
8.	Select the feature(s) you wish to save.  
a.	Note: Z Axis Length is only available in 3D.
9.	Press the “Save to txt” button  
a.	Your data will be saved to a comma delimited text file where each row is a track, and each column is a frame.
10.	Repeat 1-9 for each condition/group you wish to save.

## ADVANCED SETTINGS TAB
Before processing files from the Home tab, it is possible to adjust both segmentation and tracking parameters, though the default parameters are recommended.
 
### Segmentation settings
Minimum mitochondrial area is the minimum area threshold for a connected component to be analyzed as a mitochondrion.  
Maximum mitochondrial area is the maximum area threshold for a connected component to be analyzed as a mitochondrion.  
Custom gauss filt sigma is the standard deviation used in the 3x3 gaussian filter after the diffuse background subtraction, but before the mask creation.  
Custom threshold is the threshold level (1-256) used for creating a mask from the gaussian filtered image.  
  
### Tracking settings
Max velocity threshold is the maximum speed a mitochondrion can travel between frames.  
Max wait time threshold is the amount of time the algorithm will attempt to assign an unassigned track to a mitochondrion.  
Skip fission and fusion will skip fission and fusion assignments and speeds up processing if these features are not needed.  
Custom weights allows you to choose the cost matrix weighting for each mitochondrial morphological feature, and will speed up processing.  

## Output
### TRACKING OUTPUT
Saving a file’s data from the tracking tab will result in a .mat file of many essential and extra data from that file’s segmentation and tracking processes. Opening the .mat file in MATLAB will produce a variable named “trackList”. This variable is a 1xN structure where N is the number of tracks in the file. Each track has 1xM fields of each feature selected, where M is the number of frames in the track.

### ANALYSIS OUTPUT
Saving a file’s data from the analysis tab will result in a .txt file of the selected features. The file is comma delimited and is easily opened with programs such as excel for more complex analysis. Each row is a mitochondrion track, and each column is a frame.  
Tip: When pasting data from the text file into excel, you can split up the columns by going to Data -> Text to Columns -> Delimited -> Comma -> Finish.

## Requirements
### WINDOWS
Operating systems  
•	Windows 10 (version 1803 or higher)  
•	Windows 7 Service Pack 1  
•	Windows Server 2019  
•	Windows Server 2016  

Processors  
•	Minimum: Any Intel or AMD x86-64 processor  
•	Recommended: Any Intel or AMD x86-64 processor with four logical cores and AVX2 instruction set support  

RAM  
•	Minimum: 4 GB  
•	Recommended: 8 GB  

### MAC  
Operating systems  
•	macOS Catalina (10.15)  
•	macOS Mojave (10.14)  

Processors  
•	Minimum: Any Intel x86-64 processor  
•	Recommended: Any Intel x86-64 processor with four logical cores and AVX2 instruction set support

RAM  
•	Minimum: 4 GB  
•	Recommended: 8 GB  

### LINUX  
Operating systems    
•	Ubuntu 20.04 LTS  
•	Ubuntu 18.04 LTS  
•	Ubuntu 16.04 LTS  
•	Debian 10  
•	Debian 9  
•	Red Hat Enterprise Linux 8  
•	Red Hat Enterprise Linux 7 (minimum 7.5)  
•	SUSE Linux Enterprise Desktop 12 (minimum SP2)  
•	SUSE Linux Enterprise Desktop 15  
•	SUSE Linux Enterprise Server 12 (minimum SP2)   
•	SUSE Linux Enterprise Server 15  

Processors  
•	Minimum: Any Intel or AMD x86-64 processor  
•	Recommended: Any Intel or AMD x86-64 processor with four logical cores and AVX2 instruction set support  

RAM  
•	Minimum: 4 GB  
•	Recommended: 8 GB  

### ALL  
Mitometer has been tested only on MATLAB R2020a and R2020b.
