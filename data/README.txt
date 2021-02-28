###############################################################################

File: README.txt for nyc_datasets directory
Author: Lawrence Chillrud
Email: lgc2139@columbia.edu
Date: 02/26/2021

Contents of nyc_datasets directory: 

    - 4 .csv files: [bronx, manhattan, queens, nyc]_daily_pm2.5_components.csv
    - 1 .png file: coverage-plots.png
    - 1 README.txt (this file you are reading right now)

###############################################################################

Overview of the 4 .csv files:

    - Each file has the following column names: Date, aluminum, ammonium ion, antimony, arsenic, barium, bromine, cadmium, calcium, chlorine, chromium, cobalt, copper, elemental carbon, iron, lead, magnesium, manganese, nickel, organic carbon, potassium ion, selenium, silicon, sodium, sulfur, titanium, total nitrate, vanadium, zinc

    - Each file spans: 2001 - 2020

    - For coverage of each pm 2.5 component for each file, see: coverage-plots.png

    - The data for each file was cleaned by Lawrence Chillrud (lgc2139@columbia.edu) using raw EPA AQS data files found at the following address:

    https://epa.maps.arcgis.com/apps/webappviewer/index.html?id=5f239fd3e72f424f98ef3d5def547eb5&extent=-146.2334,13.1913,-46.3896,56.5319

    ("PM2.5 Chemical Speciation Network - Active" layer was selected and the 3 subsequent NYC monitor sites were picked to compile the raw data, for the years 2001 - 2020.)

    (The "PM2.5 Chemical Speciation Network - Inactive" layer was also selected to supplement the manhattan dataset. See more below.)

    (Finally, the "IMPROVE (Interagency Monitoring of Protected Visual Environments) - Inactive" layer was selected in an attempt to supplement the bronx dataset's carbon gaps, but the IMPROVE monitor had overlapping data.) 

Discarded metadata of the 4 .csv files:

	- Units of measurement: Micrograms/cubic meter (LC)
	- Duration Description: 24 HOUR

###############################################################################

Notes: 
	
	- The manhattan dataset is actually made up of 2 different monitor sites. From 2003 - 2007, site number 62 is used, and from 2008 - 2020, site number 134 is used. Because these sites are so close to one another (they differ by .01º latitude) they were compiled together in the manhattan dataset.   

	- Cobalt and antimony were included on a whim. They could probably be removed for analyses or left in... 

	- Ionic potassium was preferred over elemental potassium. Elemental sodium was preferred over ionic sodium. Sulfur was preferred over sulfate. The coverage for each of these pairs was the same.

	- Elemental and organic carbon were gathered from the variables "EC CSN_Rev Unadjusted PM2.5 LC TOT" and "OC CSN_Rev Unadjusted PM2.5 LC TOT" respectively. They could have been compiled from other measurements on file, but these two were selected after conducting simple linear regressions and examining parameter coverages. 

	- Chlorine has a little blip in its coverage to keep an eye out for. See coverage-plots.png

###############################################################################

bronx_daily_pm2.5_components.csv discarded metadata:

	- State Name (Code): New York (36)
	- County Name (Code): Bronx (5)
	- Site Num: 110
	- POC(s) used: {1, 5}
	- Latitude, Longitude: 40.82, -73.9
	- Datum: WGS84
	- Address: IS 52 681 KELLY ST

manhattan_daily_pm2.5_components.csv discarded metadata:

	2003 - 2007:
		- State Name (Code): New York (36)
		- County Name (Code): New York (61)
		- Site Num: 62
		- POC(s) used: {1}
		- Latitude, Longitude: 40.72, -74
		- Datum: WGS84
		- Address: POST OFFICE, 350 CANAL STREET

	2008 - 2020:
		- State Name (Code): New York (36)
		- County Name (Code): New York (61)
		- Site Num: 134
		- POC(s) used: {5}
		- Latitude, Longitude: 40.71, -74
		- Datum: WGS84
		- Address: 40 DIVISION STREET, PS 124

queens_daily_pm2.5_components.csv discarded metadata:

	- State Name (Code): New York (36)
	- County Name (Code): Queens (81)
	- Site Num: 124
	- POC(s) used: {6}
	- Latitude, Longitude: 40.74, -73.82
	- Datum: WGS84
	- Address: Queens College 65-30 Kissena Blvd Parking Lot #6

nyc_daily_pm2.5_components.csv:

	- was compiled by averaging all available data between the above 3 files, to get the maximum amount of coverage. This means that for some observations in the nyc file, only one monitor contributed to the value, while for others, 2-3 monitors influenced the value. 

###############################################################################

END OF README.txt

Good job reading all that boring stuff now have fun with the data :)

###############################################################################