\name{makeFunctionalPrediction}
\alias{makeFunctionalPrediction}
\title{Predictiong functional profiles}
\description{This function predicts functional profiles based on 16S rRNA gene data using the results of the runReferenceBlast function.}
\arguments{
  \item{path_to_otu_table}{The path to the otu table. The table must be tab-delimited with sample names in the first row and OTU ids in the first column. No additional information should be included (e.g. taxonomy). The OTU ids must correspond to the sequence ids in the fasta file}
  \item{path_to_reference_data}{The path to folder containg the reference data ("Tax4Fun2_ReferenceData_v1.1")}
	\item{path_to_temp_folder}{The same folder as provided in the runReferenceBlast function}
	\item{database_mode}{Reference sequences in the default dataset were preclustered at 99\% and 100\% identity. The 'Ref99NR' database is smaller and requires less memory than the 'Ref100NR. (Valid choices are 'Ref99NR' and 'Ref100NR'). Use the same mode as provided in the runReferenceBlast function}
	\item{normalize_by_copy_number}{Normalize the OTU abundance by 16S rRNA gene copy numbers (default is TRUE)}
	\item{min_identity_to_reference}{The minimum identity to consider hits from the BLAST (default is 97\%). A log file will document the amount of OTUs and seqeunces not to be considerd in the prediction}
	\item{include_user_data}{Include user data in the prediction (default is FALSE).}
	\item{path_to_user_data}{The path to the user data. Identical to the parameters provided by the generateUserData and generateUserDataByClustering.}
	\item{name_of_user_data}{The name of the user data. Identical to the parameters provided by the generateUserData and generateUserDataByClustering.}
	\item{use_uproc}{UProC was used to functionally anotate the genomes in the reference data. In addition, diamond was used and predictions can be made based on the diamond results. (default is TRUE)}
	\item{normalize_pathways}{Different to Tax4Fun, when converting from KEGG functions to KEGG pathways, Tax4Fun2 does not equally split KO gene abundances between pathways a functions is affiliated to. The full predicted abundance is affiliated to each pathway. Use TRUE to split the abundances (default is FALSE)}
}
\author{Written by Bernd and Franziska Wemheuer}