\name{generateUserData}
\alias{generateUserData}
\title{Generating user-defined reference profiles for Ta4Fun2}
\description{A function to generate a simple user-defined reference data set for Tax4Fun2. Recommended for smaller sets of genomes.}
\arguments{
  \item{path_to_reference_data}{The path to folder containg the reference data ("Tax4Fun2_ReferenceData_v1.1")}
	\item{path_to_user_data}{The path to the user data. Identical to the parameters provided by the generateUserData and generateUserDataByClustering.}
	\item{name_of_user_data}{The name of the user data. Identical to the parameters provided by the generateUserData and generateUserDataByClustering.}
	\item{SSU_file_extension}{The file extension of the file/files containing the extracted SSU sequences}
	\item{KEGG_file_extension}{The file extension of the file/files containing the functional profiles}
	\item{use_force}{Overwrite existing results (default is TRUE)}
}
\author{Written by Bernd and Franziska Wemheuer}