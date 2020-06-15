generateUserDataByClustering = function(path_to_reference_data = "Tax4Fun2_ReferenceData_v2", path_to_user_data, name_of_user_data = '', SSU_file_extension = '_16SrRNA.ffn', KEGG_file_extension = '_funPro.txt', use_force = T, similarity_threshold = 0.99)
{
  require(seqinr)

  read.profile = function(filename, kegg_data) 
  {
    ko_profile = read.table(file = filename, sep = ',')
    ko_profile2 = NULL
    for(i in 1:nrow(ko_profile))
    {
      for(j in strsplit(as.character(ko_profile$V1[i]), split = ":")[[1]])
      {
        ko_profile2 = c(ko_profile2, rep(j, ko_profile$V2[i]))
      }
    }
    ko_profile2 = data.frame(table(ko_profile2))
    names(ko_profile2) = c("ko","Freq")
    ko_profile = merge(kegg_data, ko_profile2, by = "ko", all.x = T, all.y = F)
    ko_profile$Freq[is.na(ko_profile$Freq)] = 0
    return(ko_profile)
  }
  
  # Check for file presence
  if(SSU_file_extension == '') stop('No file extension for 16S rRNA genes provided !')
  if(!file.exists(path_to_user_data)) stop('User folder not found!')
  if(!dir.exists(path_to_reference_data)) stop('Reference folder not found!')
  if(SSU_file_extension == KEGG_file_extension) stop('File extensions of 16S rRNA genes and profiles identical!')
  
  if (tolower(Sys.info()[["sysname"]]) == "linux") path_to_vsearch = file.path(path_to_reference_data, "TOOLS/vsearch_v2.14.1/vsearch.linux")
  if (tolower(Sys.info()[["sysname"]]) == "darwin") path_to_vsearch = file.path(path_to_reference_data, "TOOLS/vsearch_v2.14.1/vsearch.macos")
  if (tolower(Sys.info()[["sysname"]]) == "windows") path_to_vsearch = file.path(path_to_reference_data, "TOOLS/vsearch_v2.14.1/vsearch.exe")
  if (tolower(Sys.info()[["sysname"]]) != "windows")
  {
    system(paste("chmod +x", path_to_vsearch))
  }
  
  if(KEGG_file_extension != "")
  {
    if(length(dir(path = path_to_user_data, pattern = SSU_file_extension)) != length(dir(path = path_to_user_data, pattern = KEGG_file_extension))) stop('Unequal number of rRNA files and uproc profile files detected!')
  }

  path_to_user_profiles = file.path(path_to_user_data, name_of_user_data)
  # message(path_to_user_profiles)
  # Generate new files
  if(dir.exists(path_to_user_profiles))
  {
    if(use_force) unlink(x = path_to_user_profiles, recursive = T)
    if(!use_force) stop('Output folder exist! Use \'use_force = T\' to overwrite')
  }
  dir.create(path_to_user_profiles)
  
  copy_numbers = NULL
  fasta_file_out = file.path(path_to_user_profiles, "joined.fasta")
  vsearch_file_out = file.path(path_to_user_profiles, "joined.uc")
  
  if(file.exists(fasta_file_out)) unlink(x = fasta_file_out)
  file.create(fasta_file_out)
  fasta_files1 = dir(path = path_to_user_data, pattern = SSU_file_extension, full.names = F)
  fasta_files2 = dir(path = path_to_user_data, pattern = SSU_file_extension, full.names = T)
  for(i in 1:length(fasta_files1))
  {
    genome_name = gsub(SSU_file_extension, "", fasta_files1[i])
    fasta_file = read.fasta(fasta_files2[i], as.string = T)
    copy_numbers = rbind(copy_numbers, c(genome_name, length(fasta_file)))
    names(fasta_file) = paste(genome_name, "_", 1:length(fasta_file), sep = "")
    write.fasta(sequences = fasta_file, names = names(fasta_file), file.out = fasta_file_out, open = "a")
  }
  copy_numbers = data.frame(gid = copy_numbers[,1], cn = as.numeric(copy_numbers[,2]))
  
  if(similarity_threshold > 1) similarity_threshold = similarity_threshold / 100
  cmd = paste(path_to_vsearch, "--cluster_fast", fasta_file_out, "--uc", vsearch_file_out, "--id", similarity_threshold)
  
  if (tolower(Sys.info()[["sysname"]]) != "windows") system(cmd, ignore.stdout = T, ignore.stderr = T)
  if (tolower(Sys.info()[["sysname"]]) == "windows") system(cmd, intern = T)
  
  path_to_user_tab = file.path(path_to_user_profiles, 'id2seq.txt')
  vsearch_file = read.delim(file = vsearch_file_out, header = F)
  vsearch_file = vsearch_file[which(vsearch_file$V1 != "C"),]
  vsearch_file_tab = NULL
  for(i in unique(vsearch_file$V2))
  {
    rep_seq = as.character(vsearch_file$V9[which(vsearch_file$V2 == i)])[1]
    genomes_in_cluster = gsub("_[0-9]*$", "", as.character(vsearch_file$V9[which(vsearch_file$V2 == i)]))
    vsearch_file_tab = rbind(vsearch_file_tab, c(i, rep_seq, paste(genomes_in_cluster, collapse = ",")))
  }
  vsearch_file_tab = data.frame(vsearch_file_tab)
  names(vsearch_file_tab) = c("cluster_id", "rep_seq", "assigned_genome_ids")
  write.table(x = vsearch_file_tab, file = path_to_user_tab, append = F, quote = F, sep = "\t", row.names = F, col.names = T)
  
  path_to_new_rRNA = paste(path_to_user_profiles, '/USER.fasta', sep = '')
  fasta_file = read.fasta(file = fasta_file_out, as.string = T)
  user_tab = read.delim(file = path_to_user_tab, header = T)
  reference_fasta_file = file(description = path_to_new_rRNA, open = "w")
  for(i in 1:nrow(user_tab))
  {
    w1 = paste(">user", as.character(user_tab$cluster_id)[i], "_1", sep = "")
    w2 = as.character(fasta_file[[which(names(fasta_file) == as.character(user_tab$rep_seq[i]))]])
    write(x = w1, file = reference_fasta_file)
    write(x = w2, file = reference_fasta_file)
  }
  close(reference_fasta_file)
  
  unlink(fasta_file_out)
  unlink(vsearch_file_out)

  kegg_data = data.frame(ko = read.delim(paste(path_to_reference_data, '/KEGG/ko.txt', sep = ''))[,1])
  for(i in 1:nrow(user_tab))
  {
    du = NULL
    dn = NULL
    uu = NULL
    un = NULL
    
    genome_ids = strsplit(as.character(user_tab$assigned_genome_ids[i]), split = ",")[[1]]
    for(gid in genome_ids)
    {
      if(KEGG_file_extension != "")
      {
        cn = copy_numbers$cn[which(copy_numbers$gid == gid)]
        file.path(path_to_user_data, paste0(gid, KEGG_file_extension))
        genome_profile = read.profile(file.path(path_to_user_data, paste0(gid, KEGG_file_extension)), kegg_data)
        uu = rbind(uu, genome_profile$Freq)
        un = rbind(un, genome_profile$Freq / cn)
      }
    }
    if(is.null(du)) du = rbind(rep(0, nrow(kegg_data)),rep(0, nrow(kegg_data)))
    if(is.null(dn)) dn = rbind(rep(0, nrow(kegg_data)),rep(0, nrow(kegg_data)))
    if(is.null(uu)) uu = rbind(rep(0, nrow(kegg_data)),rep(0, nrow(kegg_data)))
    if(is.null(un)) un = rbind(rep(0, nrow(kegg_data)),rep(0, nrow(kegg_data)))
    
    functional_profile = data.frame(du = colMeans(du), dn = colMeans(dn), uu = colMeans(uu), un = colMeans(un))
    uid = paste("user", as.character(user_tab$cluster_id)[i], sep = "")
    path_to_functional_profile = file.path(path_to_user_profiles, paste0(uid, ".tbl"))
    write.table(x = functional_profile, file = path_to_functional_profile, append = F, quote = F, sep = "\t", row.names = F, col.names = T)
  }
}
