extractSSU = function(genome_folder = "", genome_file = "", file_extension = "fasta", path_to_refernce_data = "Tax4Fun2_ReferenceData_v2", use_vsearch = F)
{
  require(seqinr)

  blast_bin = "blastn"
  res = system(command = paste(blast_bin, "-help"), intern = T)
  if(length(res) == 0)
  {
    blast_bin = file.path(path_to_reference_data, "blast_bin/bin/blastn")
    if (tolower(Sys.info()[["sysname"]]) == "windows") blast_bin = file.path(path_to_reference_data, "blast_bin/bin/blastn.exe")
    res = system(command = paste(blast_bin, "-help"), intern = T)
    if(length(res) == 0) stop("blastn not found! Consider to use the buildDependencies() command.")
    system(paste("chmod +x", blast_bin))
  }
  
  makeblastdb_bin = "makeblastdb"
  res = system(command = paste(blast_bin, "-help"), intern = T)
  if(length(res) == 0)
  {
    blast_bin = file.path(path_to_reference_data, "blast_bin/bin/makeblastdb")
    if (tolower(Sys.info()[["sysname"]]) == "windows") blast_bin = file.path(path_to_reference_data, "blast_bin/bin/makeblastdb.exe")
    res = system(command = paste(blast_bin, "-help"), intern = T)
    if(length(res) == 0) stop("blastn not found! Consider to use the buildDependencies() command.")
    system(paste("chmod +x", makeblastdb_bin))
  }
  
  path_to_silva_database = file.path(path_to_refernce_data, "SILVA", "SILVA_132_SSURef_NR90.fasta")
  if(!file.exists(paste0(path_to_silva_database, ".nsq")))
  {
    message("Building blast database!")
    system(paste0(makeblastdb_bin, " -dbtype nucl -in ", path_to_silva_database))
  }
  
  if(nchar(file_extension) == 0) stop("No file extension specified. Must be provided, e.g. [fasta]!") 
  if(nchar(genome_folder) > 0 && nchar(genome_file) > 0) stop("You can only specify a single genome file or a folder with several genome files!")
  if(nchar(genome_folder) == 0 && nchar(genome_file) == 0) stop("Neither genome folder nor file given!")
  if(!file.exists(genome_file) && nchar(genome_file) > 0) stop("Genome file not found!")
  if(!dir.exists(genome_folder) && nchar(genome_folder) > 0) stop("Genome folder not found!")
  if(length(dir(genome_folder, pattern = file_extension)) == 0 && nchar(genome_folder) > 0) stop("Genome folder is either empty or wrong extension provided!")
  if(substr(x = file_extension, 1, 1) != '.') file_extension = paste0('.', file_extension)
  
  genome_list = genome_file
  if(nchar(genome_folder) > 0) genome_list = dir(genome_folder, pattern = file_extension, full.names = T)
 
  for(g in genome_list)
  {
    message("Currently extracting from ", basename(g))

    system(paste(blast_bin, "-db", path_to_silva_database, "-query", g, "-evalue 1e-10 -outfmt 6 -out", gsub(file_extension, ".blast", g)))
    blast_tab = read.delim(file = gsub(file_extension, ".blast", g), header = F)
    blast_tab = blast_tab[which(blast_tab$V4 > 300),]

    j = NULL
    pos = NULL
    for(i in which(blast_tab$V2 == blast_tab$V2[1]))
    {
      ctg = as.character(blast_tab$V1[i])
      if(blast_tab$V7[1] < blast_tab$V8[1])
      {
        pos1 = blast_tab$V7[i]
        pos2 = blast_tab$V8[i]
      } else {
        pos1 = blast_tab$V8[i]
        pos2 = blast_tab$V7[i]
      }
      if(is.null(j))
      {
        j = pos1:pos2
        pos = rbind(pos, c(ctg, as.character(pos1), as.character(pos2)))
      } else {
        if(length(intersect(j, pos1:pos2)) == 0)
        {
          j = sort(c(j, pos1:pos2))
          pos = rbind(pos, c(ctg, as.character(pos1), as.character(pos2)))
        }
      }
    }
    n = 1
    
    seq_name = paste0(">", gsub(file_extension, "_", g))
    genome_fasta = read.fasta(g)
    fasta_out = file(description = gsub(file_extension, "_16SrRNA.ffn", g), open = "w")
    for(i in 1:nrow(pos))
    {
      pos1 = as.numeric(pos[i,2])
      pos2 = as.numeric(pos[i,3])
      write(x = paste0(seq_name, n), file = fasta_out)
      write(x = paste(genome_fasta[[which(names(genome_fasta) == pos[i,1])]][pos1:pos2], collapse = ""), file = fasta_out)
      n = n + 1
    }
    close(fasta_out)
    
    unlink(gsub(file_extension, ".blast", g))
  }
}