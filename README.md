Locitracker is a shiny app to monitor mutations in PCR loci used in wastewater surveillance applications.

The app displays dynamic plots, with the ability to customize time periods including:

1) a heatmap of mutation frequency in clinical genomic sequences (see data sources below) by week at a given nucleotide position across PCR loci. For clarity, it is binary (position varies from reference sequence or not) and does show the nucelotide identity (i.e, A, C, G, T/U). The user can display two PCR loci at a time. **Currently there are 4 mapped PCR loci: CDC N1, CDC N2, E Sarbeco, N200.**

2) a stacked waterfall plot showing the average number of mutations across the aggregated weekly PCR locus over time. The user has the ability to add/remove PCR loci from the plot. Yellow highlighting cautions that <50 genomes were used in the analysis.

3) a line plot of the total sequence count by week (note log10 axis) used as the denominator in the mutation frequency calculations.

4) In a separate tab, Ottawa wastewater surveillance data can be explored with respect to differences between the two PCR loci targeted in testing. Namely N1 and N2. In the first plot, N1 and N2 signals normalized to PMMoV are plotted against each other over time.

5) The z-score which measures the number of standard deviations the N2 signal diverges away from (or towards) N1 serves as a method to assess accuracy of the signal. Deviation **may** imply mutation in one or both of the loci. A negative deviaton may indicate an underestimation in N1 signal.

6) the same mutational burden plot as above. Note that this is clinical data, and may not necessarily correlate with the local wastewater mutational/variant landscape.


DATA SOURCES
SARS-CoV-2 clinical genomic data: Worldwide mutation counts and total sequence counts are fetched on a weekly basis from Genbank via the API of Covspectrum. **Note that Canadian sequences are not represented in this database**

SARS-CoV-2 wastewater surveillance: This open data is from Ottawa, Canada. More information found here: https://github.com/Delatolla-lab/PHESD






