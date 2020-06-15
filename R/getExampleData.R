getExampleData = function(path_to_working_directory = ".", use_force = F)
{
  require(tools)
  
  out = tryCatch(
    {
    message("Testing the downloadURL")
    dwl_file = file.path(path_to_working_directory, "test")
    expr = download.file(url = "https://cloudstor.aarnet.edu.au/plus/s/EpS8jMFeLlbueyM/download", destfile = dwl_file, method = "auto", mode = "wb", quiet = T)
    file.remove(dwl_file)
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
      n = 1
      dwl_file = file.path(path_to_working_directory, "ExampleData.zip")
      repeat {
        download.file(url = "https://cloudstor.aarnet.edu.au/plus/s/WC2t7AkLpqzwblE/download", destfile = dwl_file, method = "auto", mode = "wb")
        if(as.character(md5sum(dwl_file)) == "baaec067ccafdf0621e5bf1b1fcc704f") break
        if(n == 3) stop("Download error!")
        n = n+1
      }
      unzip(zipfile = dwl_file, overwrite = T, exdir = path_to_working_directory)
      Sys.sleep(1)
      file.remove(dwl_file)
    }
  )
  #return(out)
}
