cat("\n")
suppressMessages(library("taxize"))
suppressMessages(library("usethis"))
suppressMessages(library("myTAI"))
suppressMessages(library("tidyverse"))
suppressMessages(library("seqinr"))
suppressMessages(library("stringr"))
suppressMessages(library("data.table"))

PosArgs <- as.character(commandArgs(trailingOnly = TRUE))

#PosArgs <- c("-q","../data/Proteins_fasta/All_queries.faa","-s","../data/Proteins_fasta/All_Results_fastacmd3.faa","-P","/home/abelardo/get_homologues/db/Pfam-A.hmm")

stopifnot(length(which(PosArgs == "-q")) != 0)
stopifnot(length(which(PosArgs == "-s")) != 0)
stopifnot(length(which(PosArgs == "-P")) != 0)

Query.pwd <- PosArgs[(which(PosArgs == "-q")+1)]
Subject.pwd <- PosArgs[(which(PosArgs == "-s")+1)]
Path_to_Pfam <- PosArgs[(which(PosArgs == "-P")+1)]

Query.fasta <- read.fasta(Query.pwd, seqtype = "AA",as.string = TRUE)
Subject.fasta <- read.fasta(Subject.pwd, seqtype = "AA",as.string = TRUE)

system("mkdir -p ../results/hmmscan")

Query.pfam.file <- tail(stringr::str_split(Query.pwd,"/")[[1]],1)
Subject.pfam.file <- tail(stringr::str_split(Subject.pwd,"/")[[1]],1)

Query.hmmscan.out <- paste(Query.pfam.file,".pfam", sep="")
Subject.hmmscan.out <- paste(Subject.pfam.file,".pfam", sep="")

#Perform hmmscan
Query.pfam.command <- paste("hmmscan -o",paste("../results/hmmscan/",Query.hmmscan.out,sep = ""),"--noali --acc --cut_ga --cpu 1",Path_to_Pfam,Query.pwd,sep = " ")
Subject.pfam.command <- paste("hmmscan -o",paste("../results/hmmscan/",Subject.hmmscan.out,sep = ""),"--noali --acc --cut_ga --cpu 1",Path_to_Pfam,Subject.pwd,sep = " ")

if(length(which(list.files("../results/hmmscan/") == Query.hmmscan.out)) == 0){
  cat("No previous hmmscan found for query... calling:\n\n")
  cat(Query.pfam.command)
  #cat("\n")
  system(Query.pfam.command)
  system(paste("sed -i \'1 i\\# raw hmmscan command: ",Query.pfam.command,"\'"," ",paste("../results/hmmscan/",Query.hmmscan.out,sep = ""),sep = ""))
  system(paste("sed -i \'1 i\\# parameters:\'",paste("../results/hmmscan/",Query.hmmscan.out,sep = "")))
  Query.pfam <- system(paste("cat",paste("../results/hmmscan/",Query.hmmscan.out,sep = ""),sep = " "),intern = TRUE)
  }else{
  cat("Previous hmmscan found for query!\n")
  Query.pfam <- system(paste("cat",paste("../results/hmmscan/",Query.hmmscan.out,sep = ""),sep = " "),intern = TRUE)
  }

if(length(which(list.files("../results/hmmscan/") == Subject.hmmscan.out)) == 0){
  print("No previous hmmscan found for subject... calling:\n")
  cat(Subject.pfam.command)
  #cat("\n")
  system(Subject.pfam.command)
  system(paste("sed -i \'1 i\\# raw hmmscan command: ",Subject.pfam.command,"\'"," ",paste("../results/hmmscan/",Subject.hmmscan.out,sep = ""),sep = ""))
  system(paste("sed -i \'1 i\\# parameters:\'",paste("../results/hmmscan/",Subject.hmmscan.out,sep = "")))
  Subject.pfam <- system(paste("cat",paste("../results/hmmscan/",Subject.hmmscan.out,sep = ""),sep = " "),intern = TRUE)
}else{
  cat("Previous hmmscan found for subject!\n\n")
  Subject.pfam <- system(paste("cat",paste("../results/hmmscan/",Subject.hmmscan.out,sep = ""),sep = " "),intern = TRUE)
}

total_hits = length(which(grepl("^>> PF",Query.pfam)))

Query.domain.annot.df <- as.data.frame(matrix(nrow = total_hits,ncol = 16))
colnames(Query.domain.annot.df) <- c("Query_ID","Num_q_local","aa_length","Pfam_id","Pfam_name","HMMER_score","biasq","c_value","i_value","hmmfrom","hmmto","alifrom","alito","envfrom","envto","acc")

locount = 0
hitcount = 0
for (l in Query.pfam){
  if (l == Query.pfam[2]){
    locount = 0
  }
  locount = locount + 1
  if (any(grepl("# raw hmmscan command:",l))){
    qcount = 0
  }
  if (any(grepl("^Query:",l))){
    qcount = qcount + 1
    num_query_total = gsub("ID","",strsplit(gsub("[^[:alnum:] ]", "", l), " +")[[1]][2])
    aa_length = as.numeric(gsub("L","",strsplit(gsub("[^[:alnum:] ]", "", l), " +")[[1]][3]))
  }
  if (any(grepl("^>> PF",l))){
    PFid = strsplit(l, " +")[[1]][2]
    PF = strsplit(l, " +")[[1]][2]
    endline = length(strsplit(l, " +")[[1]])
    PFindex = grep(PF,strsplit(l, " +")[[1]])
    PFname = paste(strsplit(l, " +")[[1]][(PFindex + 1):endline], collapse = ' ')
    PFname = gsub(" ","_",PFname)
    #print(paste(t,num_query_total,qcount,aa_length,PFid,PFname))
    parametros <- strsplit(gsub("\\]","",gsub("\\[","",Query.pfam[locount+4]))," +")
    parametros2 <- parametros[[1]][-which(gsub("\\.","",parametros[[1]])=="")]
    scoreq = parametros2[3]
    biasq = parametros2[4]
    c_value = parametros2[5]
    i_value = parametros2[6]
    hmmfrom = parametros2[7]
    hmmto = parametros2[8]
    alifrom = parametros2[9]
    alito = parametros2[10]
    envfrom = parametros2[11]
    envto = parametros2[12]
    acc = parametros2[13]
    hitcount = hitcount + 1
    #print(paste(num_query_total,qcount,aa_length,PFid,PFname,scoreq,biasq,c_value,i_value,hmmfrom,hmmto,alifrom,alito,envfrom,envto,acc))
    Query.domain.annot.df$Query_ID[hitcount] <- num_query_total
    Query.domain.annot.df$Num_q_local[hitcount] <- qcount
    Query.domain.annot.df$aa_length[hitcount] <- aa_length
    Query.domain.annot.df$Pfam_id[hitcount] <- PFid
    Query.domain.annot.df$Pfam_name[hitcount] <- PFname
    Query.domain.annot.df$HMMER_score[hitcount] <- scoreq
    Query.domain.annot.df$biasq[hitcount] <- biasq
    Query.domain.annot.df$c_value[hitcount] <- c_value
    Query.domain.annot.df$i_value[hitcount] <- i_value
    Query.domain.annot.df$hmmfrom[hitcount] <- hmmfrom
    Query.domain.annot.df$hmmto[hitcount] <- hmmto
    Query.domain.annot.df$alifrom[hitcount] <- alifrom
    Query.domain.annot.df$alito[hitcount] <- alito
    Query.domain.annot.df$envfrom[hitcount] <- envfrom
    Query.domain.annot.df$envto[hitcount] <- envto
    Query.domain.annot.df$acc[hitcount] <- acc
  }
}

total_hits = length(which(grepl("^>> PF",Subject.pfam)))

Subject.domain.annot.df <- as.data.frame(matrix(nrow = total_hits,ncol = 16))
colnames(Subject.domain.annot.df) <- c("Subject_ID","Num_q_local","aa_length","Pfam_id","Pfam_name","HMMER_score","biasq","c_value","i_value","hmmfrom","hmmto","alifrom","alito","envfrom","envto","acc")

locount = 0
hitcount = 0
for (l in Subject.pfam){
  if (l == Subject.pfam[2]){
    locount = 0
  }
  locount = locount + 1
  if (any(grepl("# raw hmmscan command:",l))){
    qcount = 0
  }
  if (any(grepl("^Query:",l))){
    qcount = qcount + 1
    num_query_total = gsub("ID","",strsplit(gsub("[^[:alnum:] ]", "", l), " +")[[1]][2])
    aa_length = as.numeric(gsub("L","",strsplit(gsub("[^[:alnum:] ]", "", l), " +")[[1]][3]))
  }
  if (any(grepl("^>> PF",l))){
    PFid = strsplit(l, " +")[[1]][2]
    PF = strsplit(l, " +")[[1]][2]
    endline = length(strsplit(l, " +")[[1]])
    PFindex = grep(PF,strsplit(l, " +")[[1]])
    PFname = paste(strsplit(l, " +")[[1]][(PFindex + 1):endline], collapse = ' ')
    PFname = gsub(" ","_",PFname)
    #print(paste(t,num_query_total,qcount,aa_length,PFid,PFname))
    parametros <- strsplit(gsub("\\]","",gsub("\\[","",Subject.pfam[locount+4]))," +")
    parametros2 <- parametros[[1]][-which(gsub("\\.","",parametros[[1]])=="")]
    scoreq = parametros2[3]
    biasq = parametros2[4]
    c_value = parametros2[5]
    i_value = parametros2[6]
    hmmfrom = parametros2[7]
    hmmto = parametros2[8]
    alifrom = parametros2[9]
    alito = parametros2[10]
    envfrom = parametros2[11]
    envto = parametros2[12]
    acc = parametros2[13]
    hitcount = hitcount + 1
    #print(paste(num_query_total,qcount,aa_length,PFid,PFname,scoreq,biasq,c_value,i_value,hmmfrom,hmmto,alifrom,alito,envfrom,envto,acc))
    Subject.domain.annot.df$Subject_ID[hitcount] <- num_query_total
    Subject.domain.annot.df$Num_q_local[hitcount] <- qcount
    Subject.domain.annot.df$aa_length[hitcount] <- aa_length
    Subject.domain.annot.df$Pfam_id[hitcount] <- PFid
    Subject.domain.annot.df$Pfam_name[hitcount] <- PFname
    Subject.domain.annot.df$HMMER_score[hitcount] <- scoreq
    Subject.domain.annot.df$biasq[hitcount] <- biasq
    Subject.domain.annot.df$c_value[hitcount] <- c_value
    Subject.domain.annot.df$i_value[hitcount] <- i_value
    Subject.domain.annot.df$hmmfrom[hitcount] <- hmmfrom
    Subject.domain.annot.df$hmmto[hitcount] <- hmmto
    Subject.domain.annot.df$alifrom[hitcount] <- alifrom
    Subject.domain.annot.df$alito[hitcount] <- alito
    Subject.domain.annot.df$envfrom[hitcount] <- envfrom
    Subject.domain.annot.df$envto[hitcount] <- envto
    Subject.domain.annot.df$acc[hitcount] <- acc
  }
}


#Nested
Nested_q_dannot <- Query.domain.annot.df %>% group_nest(Pfam_id)
Nested_s_dannot <- Subject.domain.annot.df %>% group_nest(Pfam_id)

#Subset
Subject_subset_by_query <- Nested_s_dannot[which(Nested_s_dannot$Pfam_id %in% Nested_q_dannot$Pfam_id),] %>% unnest(data)

#Build fasta out
# indexes of Subject.fasta that are also in query
Subject_hits <- unique(Subject_subset_by_query$Subject_ID)
Subject_hits_indexes <- which(gsub("[^[:alnum:] ]", "", sort(names(Subject.fasta))) %in% Subject_hits)
Subject_hits_names <- names(Subject.fasta)[Subject_hits_indexes]
Subject_hits_annots <- seqinr::getAnnot(Subject.fasta[Subject_hits_indexes])
Subject_hits_seq <- seqinr::getSequence(Subject.fasta[Subject_hits_names], as.string = TRUE) %>% unlist()

#Output names

Query_name <- head(strsplit(tail(strsplit(Query.pwd,"/")[[1]],1),"\\.")[[1]],1)
Subject_name <- head(strsplit(tail(strsplit(Subject.pwd,"/")[[1]],1),"\\.")[[1]],1)

fasta_name <- paste(Subject_name,"_subset_by_",Query_name,"_pfams.fasta",sep = "")
csv_name <- paste(Subject_name,"_subset_by_",Query_name,"_pfams.csv",sep = "")
out_pwd <- "../results/subset_by_pfam/"

#Conditional output prepair

if(!("-o" %in% PosArgs)){system("mkdir -p ../results/subset_by_pfam")}else{out_pwd <- paste(head(strsplit(PosArgs[(which(PosArgs == "-o")+1)],"/")[[1]],(length(strsplit(PosArgs[(which(PosArgs == "-o")+1)],"/")[[1]])-1)),collapse = "/")}
if("-o" %in% PosArgs){fasta_name <- paste(PosArgs[(which(PosArgs == "-o")+1)],".fasta",sep = "")}
if("-o" %in% PosArgs){csv_name <- paste(PosArgs[(which(PosArgs == "-o")+1)],".csv",sep = "")}
if(str_sub(out_pwd,start = -1) != "/"){out_pwd <- paste(out_pwd,"/",sep = "")}
if("-o" %in% PosArgs){export_fasta <- fasta_name}else{export_fasta <- paste(out_pwd,fasta_name,sep = "")}
if("-o" %in% PosArgs){export_csv <- csv_name}else{export_csv <- paste(out_pwd,csv_name,sep = "")}


#Write out
# fasta
sink(export_fasta)
for (i in 1:length(Subject_hits_seq)){
  cat(paste(Subject_hits_annots[i],Subject_hits_seq[i],"\n",sep = "\n"))
}
sink()
# table
write.csv(Subject_subset_by_query, file = export_csv)
#write dannots
write.csv(Query.domain.annot.df,file = paste("../results/hmmscan/",Query_name,"_pfam_annotation.csv",sep = ""))
write.csv(Query.domain.annot.df,file = paste("../results/hmmscan/",Subject_name,"_pfam_annotation.csv",sep = ""))

#stdout
N_prot_q <- length(Query.fasta)
N_prot_s <- length(Subject.fasta)

N_pfams_q <- length(unique(Query.domain.annot.df$Pfam_id))
N_pfams_s <- length(unique(Subject.domain.annot.df$Pfam_id))

#cat("\n\n\n")
cat(paste(paste("Number of proteins in query \"",Query_name,"\" = ",N_prot_q, sep = "")),paste("Number of Pfams =",N_pfams_q),sep = "\t|\t")
cat("\n")
cat(paste(paste("Number of proteins in subject \"",Subject_name,"\" = ",N_prot_s, sep = "")),paste("Number of Pfams =",N_pfams_s),sep = "\t|\t")
cat("\n\n")
cat("Subject counts of sequences with Pfams matching query:\n")
system(paste("cat",export_csv," | cut -d, -f6 | sort | uniq -c | grep -v \"Pfam_name\""))
cat("\n")


if("--stdout" %in% PosArgs){
  system("clear")
  for (i in 1:length(Subject_hits_seq)){
    cat(paste(Subject_hits_annots[i],Subject_hits_seq[i],"\n",sep = "\n"))
  }
}

