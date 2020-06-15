buildReferenceData = function(path_to_working_directory = ".", use_force = F, install_suggested_packages = T)
{
  require(tools)
  
  #path_to_reference_data = "./Tax4Fun2_ReferenceData_v2"
  path_to_reference_data = file.path(path_to_working_directory, "Tax4Fun2_ReferenceData_v2")
  if(dir.exists(path_to_reference_data))
  {
    if(use_force) unlink(x = path_to_reference_data, recursive = T)
    if(!use_force) stop("Folder exists. Use use_force=T to overwrite or delete manually!")
  }
  dir.create(path_to_reference_data)
  
  list_file = file.path(path_to_reference_data, "fileList.txt")
  
  out = tryCatch(
    {
    message("Testing the downloadURL")
    expr = download.file(url = "https://cloudstor.aarnet.edu.au/plus/s/EpS8jMFeLlbueyM/download", destfile = list_file, method = "auto", mode = "wb", quiet = T)
    
    },
    error=function(cond) {
      message(paste("URL does not seem to exist:", url))
      message("Here's the original error message:")
      message(cond)
      return(NA)
    },
    warning=function(cond) {
      message(paste("URL caused a warning:", url))
      message("Here's the original warning message:")
      message(cond)
      return(NULL)
    },
    finally={
      message("Test download successful!")
      message("Downloading the Tax4Fun2 reference data!")
      
      #KEGG
      dwl_file = file.path(path_to_reference_data, "KEGG.zip")
      n = 1
      repeat {
        download.file(url = "https://cloudstor.aarnet.edu.au/plus/s/MbE1nNcBwka0Y71/download", destfile = dwl_file, method = "auto", mode = "wb")
        if(as.character(md5sum(dwl_file)) == "6f72b5e02ba9d9498a1d7f02ec6e6f63") break
        if(n == 3) stop("Download error!")
        n = n+1
      }
      unzip(zipfile = dwl_file, exdir = path_to_reference_data, overwrite = T)
      file.remove(dwl_file)
      
      #Ref99NR
      dwl_file = file.path(path_to_reference_data, "Ref99NR.zip")
      n = 1
      repeat {
        download.file(url = "https://cloudstor.aarnet.edu.au/plus/s/DkoZIyZpMNbrzSw/download", destfile = dwl_file, method = "auto", mode = "wb")
        if(as.character(md5sum(dwl_file)) == "4cb3b479176d1a451588d6cc5f758a0a") break
        if(n == 3) stop("Download error!")
        n = n+1
      }
      unzip(zipfile = dwl_file, exdir = path_to_reference_data, overwrite = T)
      file.remove(dwl_file)
      
      #Ref100NR
      dwl_file = file.path(path_to_reference_data, "Ref100NR.zip")
      n = 1
      repeat {
        download.file(url = "https://cloudstor.aarnet.edu.au/plus/s/jIByczak9ZAFUB4/download", destfile = dwl_file, method = "auto", mode = "wb")
        if(as.character(md5sum(dwl_file)) == "dc81bf68e5cae609ec62dfd4f1790d2a") break
        if(n == 3) stop("Download error!")
        n = n+1
      }
      unzip(zipfile = dwl_file, exdir = path_to_reference_data, overwrite = T)
      file.remove(dwl_file)
      
      #SILVA
      dwl_file = file.path(path_to_reference_data, "SILVA.zip")
      n = 1
      repeat {
        download.file(url = "https://cloudstor.aarnet.edu.au/plus/s/EieEdRT8gEm49Ed/download", destfile = dwl_file, method = "auto", mode = "wb")
        if(as.character(md5sum(dwl_file)) == "8ed4e4d5c7c230c7330fdcc60fd50ab1") break
        if(n == 3) stop("Download error!")
        n = n+1
      }
      unzip(zipfile = dwl_file, exdir = path_to_reference_data, overwrite = T)
      file.remove(dwl_file)
      
      #TOOLS
      dwl_file = file.path(path_to_reference_data, "TOOLS.zip")
      n = 1
      repeat {
        download.file(url = "https://cloudstor.aarnet.edu.au/plus/s/5QRas9TTIogiSvo/download", destfile = dwl_file, method = "auto", mode = "wb")
        if(as.character(md5sum(dwl_file)) == "60ee0de351a4a016b58d3fe3b7e60dd0") break
        if(n == 3) stop("Download error!")
        n = n+1
      }
      unzip(zipfile = dwl_file, exdir = path_to_reference_data, overwrite = T)
      file.remove(dwl_file)
      
      message("Download finished! Checking the database for missing files...")
    }
  )
  #return(out)
  
  unlink(x = list_file)
  testReferenceData(path_to_reference_data)
  
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
  
  message("All done! Have fun with Tax4Fun2 version 1.2")
  message("If you have any issues here, try to download the entire reference data using the link provided on the Tax4Fun2 homepage.")
}

testReferenceData = function(path_to_reference_data = "Tax4Fun2_ReferenceData_v2")
{
  require(tools)
  list_file = file.path(path_to_reference_data, "fileList.txt")
  download.file(url = "https://cloudstor.aarnet.edu.au/plus/s/EpS8jMFeLlbueyM/download", destfile = list_file, method = "auto", mode = "wb", quiet = T)
  
  faulty_install = F
  test_reference_file_list = read.delim(list_file)
  for(n in 1:nrow(test_reference_file_list))
  {
    cat(paste0("Testing ", test_reference_file_list$file_name[n]))
    f = file.path(path_to_reference_data, test_reference_file_list$file_name[n])
    if(md5sum(f) != test_reference_file_list$md5sum[n]) 
    {
      warning(paste0("File corrupted: ", f))
      faulty_install = T
    } else {
      cat(" OK\n")
    }
  }
  unlink(x = list_file)
  if(!faulty_install) message("Looks good. All files are present.")
}
