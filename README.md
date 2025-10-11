# ChEMBL_BacterialStrainClassifier

**Motivation:**<br>
A common issue when analysing antimicrobial activity data from ChEMBL is non-standardised strain names for entries in the database. 
for example when analysing data for Klebsiella pneumoniae, the strain ATCC 700603 may appear as
- '700603'
- 'ATCC 700603'
- 'ATCC700603'

While these diffreences can be resolved with a simple regular expression, the ATCC 700603 actually refers to a strain derived from a single isolate that is kept by multiple tissue repositories under different acession codes (e.g., LMG 20218, CCUG 45421, CCM 4985). 

All of these names refer to the same strain, with the same antimicobial susceptibility profile, however there is currently no standardised name for these bacteria in the ChEMBL database. 

**ChEMBL_BacterialStrainClassifier** <br>
This python script standardises the names from the 'assay_strain' column in chembl by creating a new column containing the BacDive IDs (1).

See example dataset and example usage files for further information on how to use it. 

<br><br><br><br>



**References** <br>
1. Schober I, Koblitz J, Sardà Carbasse J, Ebeling C, Schmidt ML, Podstawka A, Gupta R, Ilangovan V, Chamanara J, Overmann J, Reimer LC. Bac Dive in 2025:
 the core database for prokaryotic strain data. Nucleic acids research. 2025 Jan 6;53(D1):D748-56.
