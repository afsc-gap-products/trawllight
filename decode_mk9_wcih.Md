# Decoding and preparing MK9 files with Instrument Helper

This document contains instructions for decoding TDR-MK9 .wch files using Wildlife Computers Instrument Helper software, writing data to csv files, and preparing data for processing using the trawllight R package.

## 1. Download and install

Download and install Instrument Helper from the [Wildlife Computers website](https://wildlifecomputers.com/support/downloads/).

## 2. Import Tag Data (.wch) file

Open Instrument Helper and left click on 'Press Button to Begin.' Navigate to the wch file to be imported.

![Select data processing options](/assets/ih_1_import.png)

## 3. Select Export Variables

Click on the 'Archive' tab. In the lower left of the Archive Tab, click on the red 'C' button. Navigate to channels.

![Select data processing options](/assets/ih_2_archive.png)

In the channels window, check Depth, Temperature, and Light Level. Press 'Ok'

![Select data processing options](/assets/ih_3_channels.png)


## 4. Determine if tag is deck or trawl

Use plots to determine whether the tag was deployed on deck or on the trawl. A deck-mounted tag will show little variation in depth, while a trawl tag will show considerable variation in depth (blue lines below).


Deck:
![Select data processing options](/assets/ih_4_deck.png)

Trawl: 
![Select data processing options](/assets/ih_5_deck.png)

## 5. Save data to R Read.Table (.tab) format

Click on the red 'C' button in the lower left of the Archive Tab. Write data to R Read.Table Format (Save > Instrument Readings > R Read.Table Format). File names denote whether the tag was deployed on the trawl (trwl) or on deck (deck) and the tag serial number. An underscore is used between trwl/deck designation and the serial number. In the image below, there were two trawl files collected using different tags and the file open in Instrument Helper contains data from a trawl-mounted tag (serial: 1990030). Therefore, the output file will be named trwl2_1990030.tab

Note: If multiple tags were deployed or multiple data files exist from the same tag, add a number after trwl/deck noting that there are multiple tags. The order in which data were collected in the field does not matter (e.g., if one traw data file contains data from June 31-June 30 and another trawl file contains data from July 1-July 31, it doesn't matter which is trwl1 and which is trwl2).

