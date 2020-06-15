\name{runRefBlast}
\alias{runRefBlast}
\title{Running BLAST against reference data sets}
\description{Aligns 16S rRNA gene sequences against a reference database. Requires BLAST.}
\arguments{
  \item{path_to_otus}{The path to the fasta file containing your OTU/ASV sequences}
  \item{path_to_reference_data}{The path to folder containg the reference data ("Tax4Fun2_ReferenceData_v1.1")}
	\item{path_to_temp_folder}{The folder where the output is stored}
	\item{database_mode}{Reference sequences in the default dataset were preclustered at 99\% and 100\% identity. The 'Ref99NR' database is smaller and requires less memory than the 'Ref100NR. (Valid choices are 'Ref99NR' and 'Ref100NR')}
	\item{use_force}{Overwrite existing results (default is FALSE)}
	\item{num_threads}{Number of threads to use during BLAST (default is 1)}
	\item{include_user_data}{Include user data in the prediction (default is FALSE).}
	\item{path_to_user_data}{The path to the user data. Identical to the parameters provided by the generateUserData and generateUserDataByClustering.}
	\item{name_of_user_data}{The name of the user data. Identical to the parameters provided by the generateUserData and generateUserDataByClustering.}
}
\author{Written by Bernd and Franziska Wemheuer}