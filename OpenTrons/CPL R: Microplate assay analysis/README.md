# CPL R: Microplate assay analysis.

This sub repository contains all neccessary resources to analyse the results of "CPL7: Microplate assay" and produce a list of the wells containing phages. This list can be given to the OT-2 robot, so the wells can be separated for downstream processing or other assays. 

# Contents:

  -> ParseTecan.ipynb: A Jupyter Notebook which takes the Tecan output file from the microplate assay, and reformats it to be used in the R script.
  
  -> Dictionary_template.csv: A dictionary-style csv to define the contents of each well, to be used in the R script. By default the last two wells on the microplate will be for a negative control, and a blank (H11 and H12).
  
  -> CPLR_analysis.r: An R script which creates growth curve visualisations of the microplate assay and a list of the wells with a virulence index indicative of phage presence (default 0.5). This list can be inputted into "CPL 8: Recovery of phage positive wells".

# File usage in this repository
<img width="790" alt="image" src="https://github.com/citizenphage/protocols/assets/101196413/3262c30f-9575-42e5-a15a-95effe78aad1">

# Context within pipeline:
<img width="790" alt="image" src="https://github.com/citizenphage/protocols/assets/101196413/b0786863-7a7c-4229-b739-a0d5f862cabc">
