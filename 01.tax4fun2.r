rm(list = ls())
####获得当前路径#### 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(11111)
library(Tax4Fun2)

runRefBlast(#blast_tool_path = "C:/Program Files/NCBI/blast-2.15.0+/bin",
            path_to_otus = "02.ASV.all.no.NC.fasta", 
            path_to_reference_data = "Tax4Fun2_ReferenceData_v2", 
            path_to_temp_folder = "Kelp_Ref99NR", 
            database_mode = "Ref99NR", 
            use_force = T, num_threads = 8)

makeFunctionalPrediction(path_to_otu_table = '02.ASV.all.even.txt',
                         path_to_reference_data = 'Tax4Fun2_ReferenceData_v2', 
                         path_to_temp_folder = 'Kelp_Ref99NR', 
                         database_mode = 'Ref99NR', 
                         normalize_by_copy_number = TRUE, 
                         min_identity_to_reference = 0.97,
                         normalize_pathways = TRUE)

calculateFunctionalRedundancy(path_to_otu_table = '02.ASV.all.even.txt', 
                              path_to_reference_data = 'Tax4Fun2_ReferenceData_v2', 
                              path_to_temp_folder = 'Kelp_Ref99NR', 
                              database_mode = 'Ref99NR', 
                              min_identity_to_reference = 0.97)
