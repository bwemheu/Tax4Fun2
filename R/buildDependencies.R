buildDependencies = function(path_to_reference_data = "Tax4Fun2_ReferenceData_v2", install_suggested_packages = TRUE, use_force=FALSE)
{
  message("Install Tax4Fun2 dependencies")
  for (d in c("KEGG","Ref100NR","Ref99NR","SILVA","TOOLS"))
  {
    if(!dir.exists(file.path(path_to_reference_data, d))) stop("The folder you specified might be incorrect.")
  }
  if(dir.exists(file.path(path_to_reference_data, "blast_bin")))
  {
    if(use_force) unlink(file.path(path_to_reference_data, "blast_bin"), recursive = T)
    if(!use_force) stop("Blast binary folder detected (use_force=TRUE to overwrite)")
  }
  
  # Test for packages ape and seqinr, warning if missing
  is.installed = function(mypkg) is.element(mypkg, installed.packages()[,1])
  if(install_suggested_packages)
  {
    if(!is.installed('ape')) install.packages("ape")
    if(!is.installed('seqinr')) install.packages("seqinr")
  } else {
    if(!is.installed('ape')) warning("R package ape is not installed!")
    if(!is.installed('seqinr')) warning("R package seqinr is not installed!")
  }
  
  if(!dir.exists(path_to_reference_data)) stop("Refernce folder not found! Please run the buildReferenceData() command first.")

  # BLAST FTP
  if (tolower(Sys.info()[["sysname"]]) == "linux") blast_ftp = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz'
  if (tolower(Sys.info()[["sysname"]]) == "darwin") blast_ftp = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-macosx.tar.gz'
  if (tolower(Sys.info()[["sysname"]]) == "windows") blast_ftp = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-win64.tar.gz'
  
  dwl_file = file.path(path_to_reference_data, "ncbi-blast-2.9.0.tar.gz")
  download.file(url = blast_ftp, destfile = dwl_file, mode = 'wb')
  blast_exdir = untar(tarfile = dwl_file, list = T)[1]
  untar(tarfile = dwl_file)
  file.rename(from = blast_exdir, to = file.path(path_to_reference_data, "blast_bin"))
  unlink(x = dwl_file)
  
  blast_bin = file.path(path_to_reference_data, "blast_bin/bin/blastn -help")
  if (tolower(Sys.info()[["sysname"]]) == "windows") blast_bin = file.path(path_to_reference_data, "blast_bin/bin/blastn.exe -help")
  res = system(command = blast_bin, intern = T)
  if(length(res) != 0) message("blastn installation successfull!")
}
