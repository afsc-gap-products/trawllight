# Decoding and preparing MK9 files with WC-DAP

This document contains instructions for decoding TDR-MK9 .wch files using Wildlife Computers DAP Processor (WC-DAP) software, writing data to csv files, and preparing data for processing using the trawllight R package.

## 1. Download and install

Download and install WC-DAP from the [Wildlife Computers website](https://wildlifecomputers.com/support/downloads/).

## 2. Import Tag Data (.wch) file

Open WC-DAP and use it to import data from a tag .wch using the menu to navigate to File > Import Tag Data > From File. 

![Import .wch using WC-DAP](/assets/dap_1_import.png)

Select the .wch file to decode.

![Select .wch file](/assets/dap_2_select.png)

## 3. Setup Data Processing Options

Change data processing options to export archived raw measurements as a .csv. First use the menu to navigate to Edit > Data Processing Options.

![Data processing options](/assets/dap_3_options.png)

In the CSV Files section of the Data Processing Options window, check 'Save Stored Archive Data in CSV files' and use the dropdown menu to set Default Date Format to 'Y-M-D h:m:s' (default character string format for POSIXct)

![Select data processing options](/assets/dap_4_sel_options.png)

## 4. Export Decoded Data to .csv

Export the decoded data as a .csv file. File > Export Decoded Data > Spreadsheet File (.csv).

![Select export](/assets/dap_5_export_menu.png)

## 5.  Use trawllight to convert Archive.csv to trawllight format

WC-DAP exports .csv archive files that are formatted and named differently than files that were previously exported by Wildlife Computers HexDecoder. To processs data in trawllight, the archive file must be converted to a format that is compatible with trawllight.

The file to be converted is the csv file ending in 'Archive.csv':
![Export directory](/assets/dap_7_export_ing_dir.png)

In R, use trawllight::mk9_convert_wcdap_archive() to convert the Archive csv file into a format that is compatible with trawllight:
```
library(trawllight)
trawllight::mk9_convert_wcdap_archive(dir_path = "C:/Users/sean.rohan/Work/wch/",
                                      file_name = "1190402_072821-Archive.csv",
                                      plot_check = TRUE,
                                      type = NULL)

```

If the argument type = NULL, The function detects whether the file is from a deck-mounted tag or a trawl-mounted tag. If the argument plot_check = TRUE, the function creates a plot of time versus depth in the specified filepath.

In the case above, mk9_convert_wcdap_archive writes two files, trwl_1190402.csv, which contains tag data formatted for processing with trawllight, and plotdat_1190402.png, which shows date/time versus depth (to provide a graphical check that the tag was collecting data).

![Export directory](/assets/dap_8_trwl_header.png)<br>
<p style="text-align:center"><b>Header of trwl_1190402.csv file</b>


![Export directory](/assets/dap_9_plotdat_1190402.png)<br>
<p style="text-align:center"><b>plotdat_1190402.png file</b>

## 6. First stage processing using trawllight

After repeating steps #1-5 for all .wch files (all deck and trawl tags), data can be processed using trawllight. Processing instructions are in [process_mk9.Rmd](/1_process_mk9.Rmd) ([html version](/1_process_mk9.html)).


