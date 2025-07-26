A while ago, we decided to apply the silicon drug network algorithm written by Guney and colleagues to our list of candidate genes for PSC disease. Upon reviewing the paper and running the code, I identified several issues related to parameter settings and the code itself. I reviewed two subsequent papers that referenced Guney's work and found that they had acknowledged one of the problems, attempting to resolve it through the cherry-picking method. However, this solution didn't meet my standards, so I decided to update the code, correcting the errors and explaining the causes and how they occurred.
Furthermore, I incorporated a method to score drug target genes, enhancing the algorithm's performance. This method was published in DGIdb, four years after Guney's paper was published. To complete the update, I utilized a method to score the drugs and disease target genes based on protein-protein interactions (PPIs). I've named this updated method TheraNet. 
The TharaNet.pdf document provides a scientific documentation of the method, explaining what Guney's method entails, identifying the problems encountered, and outlining how they were addressed. I hope the updated version helps other researchers generate more accurate results. 
Enjoy!

* Flowcharts explain the steps of the method.

* The data folder must be downloaded.

>To run this code, you need a list of disease target genes with Ensembl IDs in numerical format. 
To obtain such data, remove the first four letters ('ENSG') from the Ensembl gene IDs, name the column 'ENSEMBL ID', and store the data as a .csv file. 
There is a file named disease_target_genes.csv. You can use it as a sample.

* The file TheraNet.pfd explains the details of the method and includes all discussions, related figures, and tables.
