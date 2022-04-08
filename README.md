# PDRMIP_summary_codes
Code to make summary of PDRMIP data

The code is commented and currently set up to work
on locally hosted data internal at cicero.

However, the various paths are marked, so the code should
be able to work on a different system if setup correctly

The code relies of the pdrmip data being organised in a folder
with subfolders per model (and nothing else), and subfolders for
each ocean setup (coupled or fsst) from there, with all the data
in those folders. If you want a different data organisation,
then more substantial rewrites are needed to make the code work,
but it should be doable

The code is written in python and should work with various python
version, but was last tested on python 3.7. It needs packages
sys, os, glob, netCD4, csv, math, pickle, numpy and operator from itemgetter
so make sure you have these installed.
