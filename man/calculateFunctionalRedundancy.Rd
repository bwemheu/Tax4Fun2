\name{calculateFunctionalRedundancy}
\alias{calculateFunctionalRedundancy}
\title{Predictiong functional redundancies}
\description{Predicts (multi-)functional redundancies based on the phyologentic distribution and releative abundance of KEGG functions}
\arguments{
  \item{path_to_otu_table}{The path to the otu table. The table must be tab-delimited with sample names in the first row and OTU ids in the first column. No additional information should be included (e.g. taxonomy). The OTU ids must correspond to the sequence ids in the fasta file}
  \item{path_to_reference_data}{The path to folder containg the reference data ("Tax4Fun2_ReferenceData_v1.1")}
	\item{path_to_temp_folder}{The same folder as provided in the runReferenceBlast function}
	\item{database_mode}{Reference sequences in the default dataset were preclustered at 99\% and 100\% identity. The 'Ref99NR' database is smaller and requires less memory than the 'Ref100NR. (Valid choices are 'Ref99NR' and 'Ref100NR'). Use the same mode as provided in the runReferenceBlast function}
	\item{min_identity_to_reference}{The minimum identity to consider hits from the BLAST (default is 97\%). A log file will document the amount of OTUs and seqeunces not to be considerd in the prediction}
  \item{use_uproc}{UProC was used to functionally anotate the genomes in the reference data. In addition, diamond was used and predictions can be made based on the diamond results. (default is TRUE)}
}