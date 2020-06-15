assignFunction = function(genome_folder = "", genome_file = "", file_extension = "fasta", path_to_reference_data = "Tax4Fun2_ReferenceData_v2", num_of_threads = 1, fast = F, path_to_diamond_binary_mac = "")
{
  if(file_extension == "") stop("No file extension specified. Must be provided, e.g. [default = fasta]!") 
  if(nchar(genome_folder) > 0 && nchar(genome_file) > 0) stop("You can only specify a single genome file or a folder with several genome files!")
  if(nchar(genome_folder) == 0 && nchar(genome_file) == 0) stop("Neither genome folder nor file given!")
  if(!file.exists(genome_file) && nchar(genome_file) > 0) stop("Genome file not found!")
  if(!dir.exists(genome_folder) && nchar(genome_folder) > 0) stop("Genome folder not found!")
  if(length(dir(genome_folder, pattern = file_extension)) == 0 && nchar(genome_folder) > 0) stop("Genome folder is either empty or wrong extension provided!")
  if(substr(x = file_extension, 1, 1) != '.') file_extension = paste0('.', file_extension)
  
  if (tolower(Sys.info()[["sysname"]]) == "linux") path_to_prodigal = file.path(path_to_reference_data, "TOOLS", "prodigal_v2.6.3", "prodigal.linux")
  if (tolower(Sys.info()[["sysname"]]) == "darwin") path_to_prodigal = file.path(path_to_reference_data, "TOOLS","prodigal_v2.6.3", "prodigal.macos")
  if (tolower(Sys.info()[["sysname"]]) == "windows") path_to_prodigal = file.path(path_to_reference_data, "TOOLS", "prodigal_v2.6.3", "prodigal.exe")
  
  if (tolower(Sys.info()[["sysname"]]) == "linux") path_to_diamond = file.path(path_to_reference_data, "TOOLS", "diamond_v0.9.24", "diamond.linux")
  if (tolower(Sys.info()[["sysname"]]) == "darwin") path_to_diamond = file.path(path_to_reference_data, "TOOLS", "diamond_v0.9.24", "diamond.macos")
  if (tolower(Sys.info()[["sysname"]]) == "windows") path_to_diamond = file.path(path_to_reference_data, "TOOLS", "diamond_v0.9.24", "diamond.exe")

  if (tolower(Sys.info()[["sysname"]]) != "windows")
  {
    system(paste("chmod +x", path_to_prodigal))
    system(paste("chmod +x", path_to_diamond))
  }
  
  if (tolower(Sys.info()[["sysname"]]) == "darwin")
  {
    res = system(command = paste(path_to_diamond, "help"), intern = T, ignore.stderr = T)
    if(length(res) == 0) stop("Mac users might need to re-compile diamond. Please check the Tax4Fun2 wiki for help.")
  }

  res = system(command = paste(path_to_diamond, "help"), intern = T, ignore.stderr = T)
  if(length(res) == 0) warning("Please check your diamond binary file.")
  
  path_to_kegg = file.path(path_to_reference_data, "KEGG", "prokaryotes.dmnd")

  genome_list = genome_file
  if(nchar(genome_folder) > 0) genome_list = dir(genome_folder, pattern = file_extension, full.names = T)
  
  for(g in genome_list)
  {
    message("Currently annotating ", basename(g))
    if(tolower(Sys.info()[["sysname"]]) == "windows")
    {
      system(paste(path_to_prodigal, "-q -i", g, "-a", gsub(file_extension, ".faa", g)), show.output.on.console = F)
      if(fast) system(paste(path_to_diamond, "blastp -p", num_of_threads, "-e 1e-10 -d", path_to_kegg, "-q", gsub(file_extension, ".faa", g), "-o", gsub(file_extension, ".blast", g)), show.output.on.console = F)
      if(!fast) system(paste(path_to_diamond, "blastp --sensitive -p", num_of_threads, "-e 1e-10 -d", path_to_kegg, "-q", gsub(file_extension, ".faa", g), "-o", gsub(file_extension, ".blast", g)), show.output.on.console = F)
    }
    if(tolower(Sys.info()[["sysname"]]) != "windows")
    {
      system(paste(path_to_prodigal, "-q -i", g, "-a", gsub(file_extension, ".faa", g)), ignore.stdout = T, ignore.stderr = T)
      if(fast) system(paste(path_to_diamond, "blastp -p", num_of_threads, "-e 1e-10 -d", path_to_kegg, "-q", gsub(file_extension, ".faa", g), "-o", gsub(file_extension, ".blast", g)), ignore.stdout = T, ignore.stderr = T)
      if(!fast) system(paste(path_to_diamond, "blastp --sensitive -p", num_of_threads, "-e 1e-10 -d", path_to_kegg, "-q", gsub(file_extension, ".faa", g), "-o", gsub(file_extension, ".blast", g)), ignore.stdout = T, ignore.stderr = T)
    }
    if(tolower(Sys.info()[["sysname"]]) != "windows") system(paste0("rm ", gsub(file_extension, ".faa", g)))
    if(tolower(Sys.info()[["sysname"]]) == "windows") shell(paste0("del ", normalizePath(gsub(file_extension, ".faa", g), "\\")))
    blast_tab = read.table(file = gsub(file_extension, ".blast", g), sep = "\t", header = F)
    ko_list1 = unlist(lapply(strsplit(x = unique(paste(blast_tab$V1, blast_tab$V2)), split = " "), `[[`, 2))
    ko_list2 = unlist(strsplit(ko_list1, split = ":"))
    write.table(x = data.frame(table(ko_list2)), file = gsub(file_extension, "_funPro.txt", g), append = F, quote = F, sep = ",", row.names = F, col.names = F)
    unlink(gsub(file_extension, ".blast", g))
  }
}
