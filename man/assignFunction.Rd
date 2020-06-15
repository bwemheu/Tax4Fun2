\name{assignFunction}
\alias{assignFunction}
\title{Functionally annotaing prokaryotic genomes}
\description{This function uses diamond to funtionally annoate prokaryotic genomes. Protein sequences are predicted using prodigal. A database for eukaryotes is provided as well but is currently unused because of the diffuclties in gene prediction for Eukaryotes. Mac users need to compile diamond in order to use this function. See Tax4Fun2 Wiki for help}
\arguments{
  \item{genome_file}{The path to a single prokaryotic genomes}
  \item{genome_folder}{The path to a folder with prokaryotic genomes. One file per genome.}
  \item{path_to_reference_data}{The path to folder containg the reference data ("Tax4Fun2_ReferenceData_v1.1")}
	\item{num_threads}{Number of threads to use by diamond (default is 1)}
	\item{file_extension}{The file extension of the genome file/files SSU sequences should be extracted from}
	\item{fast}{Diamond runs in default mode but if set to FALSE is called in '--senstive' mode.}
	\item{path_to_diamond_binary_mac}{For Mac users only: the path to the complied diamond binary.}
}
\author{Written by Bernd and Franziska Wemheuer}