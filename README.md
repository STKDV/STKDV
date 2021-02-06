In this Github link, we conduct two case studies, which are (1) traffic accident hotspot detection in New York and (2) COVID-19 hotspot detection in Hong Kong.

Traffic accident hotspot detection in New York
We show how we display STKDV as a time-evolving hotspot map in the New York traffic accident dataset [a], which contains nearly 1.3 million data points. Here, we adopt the 128x128x128 space-time cube for generating STKDV (i.e., the time-evolving hotspot map contains 128 timestamps).

[a] NYC open data. https://data.cityofnewyork.us/Public-Safety/Motor-Vehicle-Collisions-Crashes/h9gi-nx95

The first figure shows the time-evolving hotspot map in New York (from 2012 to 2019).
![](STKDV_New_York_Traffic_Accidents.gif)

We further zoom in to Upper Manhattan and generate STKDV again (with size 128x128x128). Here, we specify the time range in the space-time cube from 2015 to 2019. The second figure shows the time-evolving hotspot map in this region.
![](STKDV_New_York_Traffic_Accidents_zoom_in.gif)

In these two figures, we can observe that the traffic accident hotspots can change in different timestamps, which can be discovered by the STKDV tool.

COVID-19 hotspot detection in Hong Kong
We also display a time-evolving hotspot map in the Hong Kong COVID-19 dataset [b]. Here, we also use 128x128x128 space-time cube for generating STKDV.

[b] Hong Kong GeoData Store https://geodata.gov.hk/gs/view-dataset?uuid=d4ccd9be-3bc0-449b-bd27-9eb9b615f2db&sidx=0

In the first figure, it shows the distribution of the COVID-19 cases, starting from 13th February 2020 to now. We can observe that STKDV tool clearly shows different waves in Hong Kong and indicates that Kowloon is the hotspot of COVID-19 cases in Hong Kong.

![](STKDV-Hong-Kong-COVID-19.gif)

In the second figure, it shows the distribution of the last wave of the COVID-19 cases (November 2020 to now) in Kowloon. Observe that STKDV clearly shows that the hotspot changes from Mong Kok to Tsim Sha Tsui.

![](STKDV-Hong-Kong-Kowloon-COVID-19.gif)
