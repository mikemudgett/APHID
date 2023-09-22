# APHID
Arabidopsis Precise Homology-directed Insertion Designer

APHID is a tool that allows you to automatically design precise insertions in Arabidopsis thaliana using the system described in (TO BE PUBLISHED SOON).

The main purpose of APHID is to eliminate the tedious work of downloading genetic sequence and searching for guides, couting nucleotides, and manually picking out primers to make your construct. All you need to do is provide the transcript ID you are interested in, the desired insertion site, and the intended cargo sequence and APHID will output all the primers you need to order as well as reference genomic and plasmid sequences.

**USAGE**

First: If you are not comfortable with using Python scripts, the tool is available for in-browser use on Google Colab: APHID Online.
Follow the steps block-by-block to prepare the required files and run the tool.

If you would like to run APHID locally, download the three scripts APHID.py, aphid_utils.py, and aphid_prep.py

_File Descriptions_

APHID.py: The main script for running the APHID tool to design inseritons once setup is complete.
aphid_utils.py: File containing repeatedly-used functions in AHPID.
aphid_prep.py: Script for creating necessary files to run APHID. Allows for customization of PAM sequence and guide length.

_Dependencies_

&mdash; Biopython
Numpy
Pandas
Primer3-py
Regex

