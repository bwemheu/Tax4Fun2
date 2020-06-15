generateUserData = function(path_to_reference_data = "Tax4Fun2_ReferenceData_v2", path_to_user_data, name_of_user_data = '', SSU_file_extension = '_16SrRNA.ffn', KEGG_file_extension = '_funPro.txt', use_force = T)
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
  if(KEGG_file_extension == '') stop('No file extension for functional profiles provided !')
  if(!file.exists(path_to_user_data)) stop('User folder not found!')
  if(!dir.exists(path_to_reference_data)) stop('Reference folder not found!')
  if(SSU_file_extension == KEGG_file_extension) stop('File extensions of 16S rRNA genes and profiles identical!')
  
  if(KEGG_file_extension != "")
  {
    if(length(dir(path = path_to_user_data, pattern = SSU_file_extension)) != length(dir(path = path_to_user_data, pattern = KEGG_file_extension))) stop('Unequal number of rRNA files and uproc profile files detected!')
  }

  path_to_user_profiles = file.path(path_to_user_data, name_of_user_data)
  path_to_user_tab = file.path(path_to_user_profiles, 'id2seq.txt')
  path_to_new_rRNA = file.path(path_to_user_profiles, 'USER.fasta')
  
  # Generate new files
  if(dir.exists(path_to_user_profiles))
  {
    if(use_force) unlink(x = path_to_user_profiles, recursive = T)
    if(!use_force) stop('Output folder exist! Use \'use_force = T\' to overwrite')
  }
  dir.create(path_to_user_profiles)
  
  copy_numbers = NULL
  fasta_file_out = file(description = path_to_new_rRNA, open = "w")
  fasta_files1 = dir(path = path_to_user_data, pattern = SSU_file_extension, full.names = F)
  fasta_files2 = dir(path = path_to_user_data, pattern = SSU_file_extension, full.names = T)
  for(i in 1:length(fasta_files1))
  {
    genome_id = paste("user", i, sep = "")
    genome_name = gsub(SSU_file_extension, "", fasta_files1[i])
    fasta_file = read.fasta(file = fasta_files2[i], as.string = T)
    copy_numbers = rbind(copy_numbers, c(genome_id, genome_name, length(fasta_file)))
    
    seq_num = 1
    for(fasta_seq in fasta_file)
    {
      w1 = paste(">", genome_id, "_", seq_num, sep = "")
      w2 = as.character(fasta_seq)
      write(x = w1, file = fasta_file_out)
      write(x = w2, file = fasta_file_out)
      seq_num = seq_num + 1
    }
  }
  close(fasta_file_out)
  user_tab = data.frame(gid = copy_numbers[,1], genome_name = copy_numbers[,2], cn = as.numeric(copy_numbers[,3]))
  names(user_tab) = c("cluster_id", "assigned_genome_ids", "cn")
  write.table(x = user_tab, file = path_to_user_tab, append = F, quote = F, sep = "\t", row.names = F, col.names = T)
  
  kegg_data = data.frame(ko = read.delim(paste(path_to_reference_data, '/KEGG/ko.txt', sep = ''))[,1])
  for(i in 1:nrow(user_tab))
  {
    du = rep(0, nrow(kegg_data))
    dn = rep(0, nrow(kegg_data))
    uu = rep(0, nrow(kegg_data))
    un = rep(0, nrow(kegg_data))
    user_tab
    gid = strsplit(as.character(user_tab$assigned_genome_ids[i]), split = ",")[[1]]
    cn = user_tab$cn[i]
    if(KEGG_file_extension != "")
    {
      file.path(path_to_user_data, paste0(gid, ".", KEGG_file_extension))
      genome_profile = read.profile(file.path(path_to_user_data, paste0(gid, KEGG_file_extension)), kegg_data)
      uu = genome_profile$Freq
      un = genome_profile$Freq / cn
    }

    functional_profile = data.frame(du, dn, uu, un)
    uid = as.character(user_tab$cluster_id)[i]
    path_to_functional_profile = paste(path_to_user_profiles, "/", uid, ".tbl", sep = "")
    write.table(x = functional_profile, file = path_to_functional_profile, append = F, quote = F, sep = "\t", row.names = F, col.names = T)
  }
}
