[This template readme file should be edited to be relevant to your dataset. The template proposes a basic set of information to be provided about a dataset. Sections 1-3 provide key information about the dataset and should be completed as fully as possible; Sections 4-5 provide information for interpretation and use of the dataset, and should be completed according to your judgement. Ask yourself in completing these sections: what information would the user of this dataset need in order to be able to understand it or replicate the results?

Use of the README plain text format for dataset documentation is not required, and may not be suitable for longer or more detailed documentation. In these cases, or if preferred, you can use PDF or MS Word.

Information provided here must correspond accurately with information provided in the dataset metadata record, e.g. the dataset title should match exactly, the same Creators should be listed, etc.

The readme file should be saved with the name README_[Creator surname]_[Publication year]. The file name should contain no spaces and not exceed 32 characters. Examples: README_Smith_2022.txt; README_Jones-etal_2022.txt.

Text within square brackets is instructional and should be deleted from the final version of the readme.]

1. ABOUT THE DATASET
------------

Title:	Back trajectories released from ocean catchment boundaries

Creator(s): Philip Craig (https://orcid.org/0000-0001-9213-4599)

Organisation(s): University of Reading

Rights-holder(s): Philip Craig

Publication Year: 2023

Description: This dataset contains a subset from a dataset of 14-day back trajectories released every 12 hours from points along the ocean drainage basin catchment boundaries on 17 vertical levels.

Cite as: Craig, P. (2023) Back trajectories released from ocean catchment boundaries. University of Reading. Dataset. [DOI]

Related publication: Craig (2018) The Atlantic/Pacific atmospheric moisture budget asymmetry: the role of atmospheric moisture transport. PhD Thesis, University of Reading (https://doi.org/10.48683/1926.00084855)
		     Craig, Ferreira & Methven (2023) <title to be decided> <to be submitted to JGR Atmospheres>

Acknowledgements: David Ferreira & John Methven also contributed to this dataset as PhD supervisors.


2. TERMS OF USE
------------

Copyright Philip Craig 2023. This dataset is licensed by the rights-holder(s) under a Creative Commons Attribution 4.0 International Licence: https://creativecommons.org/licenses/by/4.0/.


3. PROJECT AND FUNDING INFORMATION
------------

Title: The Atlantic/Pacific Atmospheric Moisture Budget Asymmetry: The Role of Atmospheric Moisture Transport

Dates: September 2014 to April 2018

Funding organisation: NERC

Grant no.: NE/L002566/1


4. CONTENTS
------------
File listing

utraj-df_afr2010073112
utraj-df_amr2010073112
utraj-df_ara2010073112
utraj-df_ari2010073112
utraj-df_arp2010073112
utraj-df_eaa2010073112
utraj-df_soa2010073112
utraj-df_soi2010073112
utraj-df_sop2010073112
utraj-tr_afr2010073112
utraj-tr_amr2010073112
utraj-tr_ara2010073112
utraj-tr_ari2010073112
utraj-tr_arp2010073112
utraj-tr_eaa2010073112
utraj-tr_soa2010073112
utraj-tr_soi2010073112
utraj-tr_sop2010073112
<Python script>

The file naming convention of the trajectory output files (utraj-df) is:

	utraj-df_<boundary>YYYYMMDDhh

The trajectory output files are human-readable text files & contain data at timesteps every 6 hours along each trajectory.
The metadata header contains the following information:

	TRAJECTORY BASE TIME (time of trajectories initialization)
	DATA BASE TIME (only relevant when using forecast data)
	DATA INTERVAL (size of timestep)
	TOTAL NUMBER OF TRAJECTORIES (number of trajectories contained in file)
	NUMBER OF ATTRIBUTES (number of variables interpolated along trajectory)
	ATTRIBUTE TYPES (variables interpolated along trajectory - represented by integers linked to ECMWF convention)
	NUMBER OF CLUSTERS (number of clusters, vertical levels in this case)
	CLUSTER POINTERS (trajectory numbers where each cluster begins)
	3D TRAJECTORY (T or F)
	FORECAST DATA (T or F)
	FORWARD TRAJECTORY (T or F)

The trajectory data are presented as follows:
	
	TRAJECTORY NUMBER     1 COMPRISES    56 INTERVALS

	STEP    HOURS  LAT          LON          P (MB)      ATTRIBUTES
    	0     0.00  0.40080E+02  0.41980E+02  0.76041E+03  0.29303E+03 -0.37576E-01  0.54442E-02  0.41508E+03  0.20706E+04
    	6    -6.00  0.40090E+02  0.43210E+02  0.76268E+03  0.28571E+03  0.30831E+00  0.57332E-02  0.48562E+03  0.95496E+03
   	12   -12.00  0.39938E+02  0.44454E+02  0.77076E+03  0.28728E+03  0.19342E+00  0.55529E-02  0.62563E+03  0.21712E+02

Header definitions:

	STEP (timestep along trajectory in hours)
	HOURS (time along trajectory, negative for back trajectory)
	LAT (latitude)
	LON (longitude)
	P (MB) (pressure in millibars)
	ATTRIBUTES (other variables interpolated along trajectory)

The file naming convention of the trajectory initialization files (utraj-tr) is:
	
	utraj-tr_<boundary>YYYYMMDDhh

The trajectory initialization files are human-readable text files & contain data at the trajectory initialization time.
The metadata headers are the same as the trajectory output files but without the DATA INTERVAL header.
The trajectory initialization data are presented in the same manner as the output files but only at the initialization time (STEP 0).

The <boundary> in the filenames is the catchment boundary of an ocean drainage basin along which trajectories are initialized.
These boundaries are:
	
	afr (Africa)
	amr (Americas)
	ara (Arctic Atlantic)
	ari (Arctic Indian)
	arp (Arctic Pacific)
	eaa (East Asia Australia*)
	soa (Southern Ocean Atlantic)
	soi (Southern Ocean Indian)
	sop (Southern Ocean Pacific)

*this boundary is referred to as "South-East Asia" in publications

The YYYYMMDDhh timestamp refers to the date and time at which the trajectories were initialized.
	
	YYYY (year)
	MM (month)
	DD (day)
	hh (hour)

The Python script is a small example of how to read the data.


5. METHODS
-----------

This dataset was generated using the ROTRAJ model (Methven, 1997; Methven et al., 2001; de Leeuw et al, 2017) and the ERA-Interim reanalysis (Dee et al, 2011) to calculate trajectories released from catchment boundaries surrounding ocean drainage basins (Craig, 2019; https://doi.org/10.17864/1947.195).

Back trajectories were released from points ~70 km apart along the ocean drainage basin catchment boundaries on 17 vertical levels. Each back trajectory is 14 days long and data are output every 6 hours.

Extensive information is included in chapter 2 of Philip Craig's PhD thesis (https://doi.org/10.48683/1926.00084855) and in Craig et al (2023, <to be submitted to JGR Atmospheres>).

References:
	Craig, P. (2018) The Atlantic/Pacific atmospheric moisture budget asymmetry: the role of atmospheric moisture transport. PhD Thesis, University of Reading (https://doi.org/10.48683/1926.00084855)
	Craig, P. (2019) Catchment boundaries of ocean drainage basins. University of Reading. Dataset. https://doi.org/10.17864/1947.195
	Craig, P. et al. (2023) <to be submitted to JGR Atmospheres>
	Dee, D. et al. (2011) The ERA-Interim Reanalysis: configuration and performance of the data assimilation system. Quarterly Journal of the Royal Meteorological Society, 137, 553-597
	de Leeuw, J. et al. (2017) Physical Factors Influencing Regional Precipitation Variability Attributed Using an Airmass Trajectory Method. Journal of Climate, 30, 7359-7378
	Methven, J. (1997) Offline trajectories: Calculation and accuracy (Tech. Rep. No. 44). Dept. of Meteorol., Univ. of Reading, U.K.: U.K. Univ. Global Atmos. Modelling Prog.
	Methven, J. et al. (2001) Estimating relationships between air mass origin and chemical composition. Journal of Geophysical Research: Atmospheres, 106, 5005-5019