#Interaction Tests

#Convert response to binary number
only_data$response_binary <- ifelse(only_data$response == "response", 1,
                                    ifelse(only_data$response == "no_response", 0,NA))
table(only_data$response_binary)

signatures <- cols(only_data)[[13:20, 22, 35, 43:46, 50]]
clusters <- cols(only_data)[[1:10]]


#Run test in a loop
clusters <- c("cluster_0.x" ,"cluster_1.x","cluster_2.x","cluster_3.x","cluster_4.x","cluster_5.x","cluster_6.x","cluster_7.x","cluster_8.x" , "cluster_9.x")
signatures <- c("TMB", "Clonal_TMB", "Subclonal_TMB", "Indel_TMB", "NMD_escape_TMB" , "Signature.4", "Signature.7","signature_apobec","MUC16", "SDI" , "CD8A", "T_inflam_GEP" ,"CXCL9", "CD274", 'master_purity_R1')

final_dat <- matrix(nrow = 10, ncol = 15) %>% as.data.frame()
rownames(final_dat) <- clusters
colnames(final_dat) <- signatures

for(i in clusters) {
  cluster <- i %>% as.character()
  
  for (j in signatures) {
    signature <- j %>% as.character()
    tmp <- only_data %>% select(response_binary, cluster, signature, cancer_type.x)
    colnames(tmp) <- c("response", "cluster", "signature", "type")
    res <- glm(response~ signature * cluster + type, family=binomial(link="logit"),data=tmp)
    sum <- summary(res)
    sum <- sum$coefficients
    sum <- sum[6, 4]
    
    final_dat[i, j] <- sum
    
  }
}
