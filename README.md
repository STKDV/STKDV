In this Github repository, we conduct two case studies, which are (1) traffic accident hotspot detection in New York and (2) COVID-19 hotspot detection in Hong Kong.

**(1) Traffic accident hotspot detection in New York**

We show how we display STKDV as a time-evolving hotspot map in the New York traffic accident dataset [a], which contains nearly 1.3 million data points. Here, we adopt the 128x128x128 space-time cube for generating STKDV (i.e., the time-evolving hotspot map contains 128 timestamps).

[a] NYC open data. https://data.cityofnewyork.us/Public-Safety/Motor-Vehicle-Collisions-Crashes/h9gi-nx95

The first figure shows the time-evolving hotspot map in New York (from 2012 to 2019).
![](STKDV_New_York_Traffic_Accidents.gif)

We further zoom in to Upper Manhattan and generate STKDV again (with size 128x128x128). Here, we specify the time range in the space-time cube from 2015 to 2019. The second figure shows the time-evolving hotspot map in this region.
![](STKDV_New_York_Traffic_Accidents_zoom_in.gif)

In these two figures, we can observe that the traffic accident hotspots can change in different timestamps, which can be discovered by the STKDV tool.

**(2) COVID-19 hotspot detection in Hong Kong**

We also display a time-evolving hotspot map in the Hong Kong COVID-19 dataset [b]. Here, we also use 128x128x128 space-time cube for generating STKDV.

[b] Hong Kong GeoData Store https://geodata.gov.hk/gs/view-dataset?uuid=d4ccd9be-3bc0-449b-bd27-9eb9b615f2db&sidx=0

In the first figure, it shows the distribution of the COVID-19 cases, starting from 13th February 2020 to 5th February 2021. We can observe that STKDV tool clearly shows different waves in Hong Kong and indicates that Kowloon is the hotspot of COVID-19 cases in Hong Kong.
![](STKDV-Hong-Kong-COVID-19.gif)

In the second figure, it shows the distribution of the last wave of the COVID-19 cases (November 2020 to February 2021) in Kowloon. Observe that STKDV clearly shows that the hotspot changes from Mong Kok to Tsim Sha Tsui.
![](STKDV-Hong-Kong-Kowloon-COVID-19.gif)

**(3) Reproducing our results**

All the codes for generating these time-evolving hotspot map can be found in the file "Time_evolving_hotspot_map_code.zip". To run our code, you need to follow these three steps.

1: Install node.js

2: Use cmd (in Windows) to access the directory "./Time_evolving_hotspot_map_code/Time_evolving_hotspot_map_code"

3: Type these two commands (1) npm install and (2) npm start

Since the size of these result files are large, we only provide the time-evolving hotspot map for Upper Manhattan (i.e., "STKDV_New_York_Traffic_Accidents_zoom_in.gif").
