calculateFunctionalRedundancy = function(path_to_otu_table, path_to_reference_data = "Tax4Fun2_ReferenceData_v2", path_to_temp_folder, database_mode = 'Ref100NR', min_identity_to_reference = 97, use_uproc = T)
{
  require(package = 'ape')
  # Read old log file to control the database mode
  path_to_log_file = file.path(path_to_temp_folder, 'logfile1.txt')
  log_file = read.delim(path_to_log_file, header = F)
  if(database_mode != as.character(log_file$V1[2])) warning('Logfile indicates other database mode')
  
  # Set path to reference files
  if(database_mode == 'Ref99NR') {
    path_to_ref_profiles = file.path(path_to_reference_data, 'Ref99NR')
    path_to_ref_tree = file.path(path_to_ref_profiles, 'Ref99NR.tre')
  } else if (database_mode == 'Ref100NR') {
    path_to_ref_profiles = file.path(path_to_reference_data, 'Ref100NR')
    path_to_ref_tree = file.path(path_to_ref_profiles, 'Ref100NR.tre')
  } else {
    stop('Database mode unknown! valid choces are Ref99NR and Ref100NR')
  }

  # Write to new log file
  path_to_log_file = file.path(path_to_temp_folder, 'logfile3.txt')
  write(x = "Tax4fun2 beta", file = path_to_log_file, append = T)
  write(x = date(), file = path_to_log_file, append = T)

  # Reading and reducing the blast file
  ref_blast_result = read.delim(file.path(path_to_temp_folder, 'ref_blast.txt'), h = F)
  ref_blast_result_reduced = ref_blast_result[which(ref_blast_result$V3 >= min_identity_to_reference), 1:2]

  # Reading and filtering the otu table
  otu_table = read.delim(path_to_otu_table)
  otu_table_reduced = merge(x = ref_blast_result_reduced, y = otu_table, by.x = 'V1', by.y = names(otu_table)[1])[,-1]
  otu_table_reduced_aggregated = aggregate(x = ifelse(otu_table_reduced[,-1]>0,1,0), by = list(otu_table_reduced[,1]), sum)
  
  # Write unknown fraction to log file
  if((ncol(otu_table) - 1) == 1)
  {
    unknown_fraction1 = as.data.frame(round(1 - sum(ifelse(otu_table_reduced[,-1]>0,1,0)) / sum(ifelse(otu_table[,-1]>0,1,0)), digits = 5))
    write(x = 'Unknown fraction (amount of otus unused in the prediction) for each sample:', file = path_to_log_file, append = T)
    write.table(x = unknown_fraction1, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
    unknown_fraction2 = as.data.frame(round(1 - sum(otu_table_reduced_aggregated[,-1]) / sum(otu_table[,-1]), digits = 5))
    write(x = 'Unknown fraction (amount of sequences unused in the prediction) for each sample:', file = path_to_log_file, append = T)
    write.table(x = unknown_fraction2, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
  } else {
    unknown_fraction1 = as.data.frame(round(1 - colSums(ifelse(otu_table_reduced[,-1]>0,1,0)) / colSums(ifelse(otu_table[,-1]>0,1,0)), digits = 5))
    write(x = 'Unknown fraction (amount of otus unused in the prediction) for each sample:', file = path_to_log_file, append = T)
    write.table(x = unknown_fraction1, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
    unknown_fraction2 = as.data.frame(round(1 - colSums(otu_table_reduced_aggregated[,-1]) / colSums(otu_table[,-1]), digits = 5))
    write(x = 'Unknown fraction (amount of sequences unused in the prediction) for each sample:', file = path_to_log_file, append = T)
    write.table(x = unknown_fraction2, file = path_to_log_file, append = T, quote = F, sep = ': ', row.names = T, col.names = F)
  }
  
  # Reading reference tree
  reference_tree = read.tree(path_to_ref_tree)
  distance_matrix = cophenetic(reference_tree)
  
  rownumber = NULL
  for(reference_id in otu_table_reduced_aggregated$Group.1)
  {
    rownumber = c(rownumber, which(row.names(distance_matrix) == reference_id))
  }
  distance_matrix_reduced = distance_matrix[rownumber,rownumber]
  distance_matrix_mean = mean(as.dist(distance_matrix))
  distance_matrix_reduced_mean = mean(as.dist(distance_matrix_reduced))
  rm(distance_matrix)
  
  # normalize or not normalize is the question
  n = 1
  if(use_uproc) n = 3

  # Generate reference profile
  print('Generating reference profile')
  reference_profile = NULL
  for(reference_id in otu_table_reduced_aggregated$Group.1)
  {
    print(reference_id)
    path_to_profile = path_to_ref_profiles
    reference_file_path = file.path(path_to_profile, paste0(reference_id, '.tbl.gz'))
    reference_file_path
    reference_file = read.delim(file = reference_file_path)
    reference_profile = rbind(reference_profile, as.numeric(reference_file[,n]))
  }
  dim(reference_profile)

  ko_list = read.delim(file.path(path_to_reference_data, 'KEGG/ko.txt'))
  # Calculate functional redundancy sample-wise
  abs_functional_redundancy_tab = NULL
  rel_functional_redundancy_tab = NULL
  for(sample in 2:ncol(otu_table_reduced_aggregated))
  {
    print(paste('Calculate functional redundancy for sample', names(otu_table_reduced_aggregated[sample])))
    functional_prediction_sample = reference_profile * as.numeric(otu_table_reduced_aggregated[,sample])
    sum(functional_prediction_sample)
    functional_prediction_sample_mod = ifelse(functional_prediction_sample>=1,1,0)
    sum(functional_prediction_sample_mod)
    
    abs_functional_redundancy_sample = NULL
    rel_functional_redundancy_sample = NULL
    for(i in 1:nrow(ko_list))
    {
      ko_count = functional_prediction_sample_mod[,i]
      ko_count
      aFRI = (mean(as.dist(distance_matrix_reduced * ko_count)) * (sum(ko_count) / length(ko_count))) / distance_matrix_mean
      rFRI = (mean(as.dist(distance_matrix_reduced * ko_count)) * (sum(ko_count) / length(ko_count))) / distance_matrix_reduced_mean
      abs_functional_redundancy_sample = c(abs_functional_redundancy_sample, aFRI)
      rel_functional_redundancy_sample = c(rel_functional_redundancy_sample, rFRI)
    }
    abs_functional_redundancy_tab = cbind(abs_functional_redundancy_tab, abs_functional_redundancy_sample)
    rel_functional_redundancy_tab = cbind(rel_functional_redundancy_tab, rel_functional_redundancy_sample)
  }
  
  abs_functional_redundancy_tab = data.frame(abs_functional_redundancy_tab)
  rel_functional_redundancy_tab = data.frame(rel_functional_redundancy_tab)
  
  colnames(abs_functional_redundancy_tab) = names(otu_table)[2:ncol(otu_table_reduced_aggregated)]
  colnames(rel_functional_redundancy_tab) = names(otu_table)[2:ncol(otu_table_reduced_aggregated)]

  abs_functional_redundancy_final = data.frame(KO = ko_list$ko, abs_functional_redundancy_tab, description = ko_list$description)
  rel_functional_redundancy_final = data.frame(KO = ko_list$ko, rel_functional_redundancy_tab, description = ko_list$description)
  write.table(x = abs_functional_redundancy_final, file = file.path(path_to_temp_folder, 'absolute_functional_redundancy.txt'), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(x = rel_functional_redundancy_final, file = file.path(path_to_temp_folder, 'relative_functional_redundancy.txt'), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
}
