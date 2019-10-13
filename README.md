# TARP: Transposable Element Assembly Remapping Pipeline

Transposable elements are prolific throughout plant genomes and while originally regarded as nothing more than genetic parasites increasing number of studies have shown these elements to contribute significantly to the phenotypes of many significant agricultural species. Therefore, accurate and up to date information on the content of transposable elements in plant, genomes is an essential piece of information for many researchers.  

  
However, as plant genome assemblies continue to be improved and updates to old versions released the public libraries of transposable elements often remain the same. The repetitive nature of transposable elements means they can easily be left in unassembled scaffolds in early genome versions and then over time make their way into the chromosomal sequences as sequencing data becomes more robust. There is, therefore, a need for computational methods to reliably and accurately locate and remap transposable elements between outdated and current assemblies.  

  
Here, we present Transposable Element Assembly Remap Pipeline or TARP. TARP utilizes iterative global alignments of transposable element consensus sequences to locate known and novel elements in updated plant genome assemblies without searching for either one of these categories explicitly. This allows researchers interested in the transposable element content of updated plant assemblies to utilize previous transposable element libraries and data in their search instead of starting from scratch. 

  
Additionally, using flanking sequence comparisons TART can distinguish between novel and known elements in a new assembly even if there have been significant locational changes to known elements due to the incorporation of scaffolds into updated assemblies. This allows researchers to compare where elements of interest are located between two or more assembly versions.   
