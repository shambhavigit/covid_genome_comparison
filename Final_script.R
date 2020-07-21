trimers <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "ATA", "ATC", "ATG", "CAA", 
           "CAC", "CAG", "CCA","CCC","CCG","CGA","CGC","CTA","CTC","GAA","GAC","GCA","GCC","GGA","GTA","TAA","TCA")

tetramers <- c('AAAA','AAAC','AAAG','AAAT','AACA','AACC','AACG','AACT','AAGA','AAGC','AAGG','AAGT','AATA','AATC',
             'AATG','AATT','ACAA','ACAC','ACAG','ACAT','ACCA','ACCC','ACCG','ACCT','ACGA','ACGC','ACGG','ACGT',
             'ACTA','ACTC','ACTG','AGAA','AGAC','AGAG','AGAT','AGCA','AGCC','AGCG','AGCT','AGGA','AGGC','AGGG',
             'AGTA','AGTC','AGTG','ATAA','ATAC','ATAG','ATAT','ATCA','ATCC','ATCG','ATGA','ATGC','ATGG','ATTA',
             'ATTC','ATTG','CAAA','CAAC','CAAG','CACA','CACC','CACG','CAGA','CAGC','CAGG','CATA','CATC','CATG',
             'CCAA','CCAC','CCAG','CCCA','CCCC','CCCG','CCGA','CCGC','CCGG','CCTA','CCTC','CGAA','CGAC','CGAG',
             'CGCA','CGCC','CGCG','CGGA','CGGC','CGTA','CGTC','CTAA','CTAC','CTAG','CTCA','CTCC','CTGA','CTGC',
             'CTTA','CTTC','GAAA','GAAC','GACA','GACC','GAGA','GAGC','GATA','GATC','GCAA','GCAC','GCCA','GCCC',
             'GCGA','GCGC','GCTA','GGAA','GGAC','GGCA','GGCC','GGGA','GGTA','GTAA','GTAC','GTCA','GTGA','GTTA',
             'TAAA','TACA','TAGA','TATA','TCAA','TCCA','TCGA','TGAA','TGCA','TTAA')


# Read the file. File is read as multiple lines
gs <- scan("D:/work/projects/R/genomic analysis/LC5282832/LC528232.txt", what="character", sep=NULL)


# Merge multiple lines into 1 single line
gs_single_line <- paste(gs, collapse="")

# Split the single line into characters
gs_list <- strsplit(gs_single_line, "")

OR library(seqinr)
 mysequence <- s2c(gs)


oligo_counter <- function(gs_vector,oligo_collection,oligo_len){
	i <- oligo_len
	j <- 0
      counter <- c()
	len_gs_vector <- length(gs_vector)
	for (i in i:len_gs_vector)
		{     j <- 0
			temp_extract <- c()
			for ( j in 1:oligo_len-1) 
				{     
					temp_extract <- c(temp_extract,gs_vector[c(i-j)])
					#print(temp_extract)
				}               

 			curr_seq <- paste(temp_extract ,collapse='')
			#print(curr_seq)
 			k <- 1   
 			# Merge cs into single string   
 			for (k in k:length(oligo_collection))
  			  {
     				if (curr_seq == oligo_collection[k])
     				{
         			 counter <- c(counter, oligo_collection[k])
				}
     			  } 
		}

	a <- table(counter) 
	Seq <- as.data.frame(a)
	
}

ebola_tetra <- oligo_counter(ebola_seq ,tetramers,4)

# Read and vectorize sequences one by one
covid <- read.fasta("D:/work/projects/R/genomic analysis/datasets_fasta/cov2.fasta",forceDNAtolower = FALSE)
SARS <-  read.fasta("D:/work/projects/R/genomic analysis/datasets_fasta/sars.fasta",forceDNAtolower = FALSE)
Mers <-  read.fasta("D:/work/projects/R/genomic analysis/datasets_fasta/mers.fasta",forceDNAtolower = FALSE)
Bat_sars <-  read.fasta("D:/work/projects/R/genomic analysis/datasets_fasta/BAT_SARS.fasta",forceDNAtolower = FALSE)
Bat_RaTG13 <- read.fasta("D:/work/projects/R/genomic analysis/datasets_fasta/bat_RaTG13.txt",forceDNAtolower = FALSE)
camel <- read.fasta("D:/work/projects/R/genomic analysis/datasets_fasta/camelus.fasta",forceDNAtolower = FALSE)
civet <- read.fasta("D:/work/projects/R/genomic analysis/datasets_fasta/Civet-SARS.fasta",forceDNAtolower = FALSE)
ebola <- read.fasta("D:/work/projects/R/genomic analysis/datasets_fasta/EBOLAV.fasta",forceDNAtolower = FALSE)
hedgehog <- read.fasta("D:/work/projects/R/genomic analysis/datasets_fasta/hedgehog.fasta",forceDNAtolower = FALSE)
hiv <- read.fasta("D:/work/projects/R/genomic analysis/datasets_fasta/hiv2.fasta",forceDNAtolower = FALSE)
malaria <- read.fasta("D:/work/projects/R/genomic analysis/datasets_fasta/malariae.fasta",forceDNAtolower = FALSE)
Pangolin <- read.fasta("D:/work/projects/R/genomic analysis/datasets_fasta/pangolin.txt",forceDNAtolower = FALSE)
Bat_SL_CoVZC45 <- read.fasta("D:/work/projects/R/genomic analysis/datasets_fasta/bat-SL-CoVZC45.txt",forceDNAtolower = FALSE)





# Convert first sequence into vector
covid_seq <- unlist(covid[[1]][1:length(covid[[1]])])
sars_seq <- unlist(SARS[[1]][1:length(SARS[[1]])])
mers_seq <- unlist(Mers[[1]][1:length(Mers[[1]])])
bat_seq <- unlist(Bat_sars[[1]][1:length(Bat_sars[[1]])])
bat_ratg13_seq<- unlist(Bat_RaTG13[[1]][1:length(Bat_RaTG13[[1]])])
camel_seq <- unlist(camel[[1]][1:length(camel[[1]])])
civet_seq <- unlist(civet[[1]][1:length(civet[[1]])])
ebola_seq <- unlist(ebola[[1]][1:length(ebola[[1]])])
hedgehog_seq <- unlist(hedgehog[[1]][1:length(hedgehog[[1]])])
hiv_seq <- unlist(hiv[[1]][1:length(hiv[[1]])])
malaria_seq <- unlist(malaria[[1]][1:length(malaria[[1]])])
Pangolin_seq <- unlist(Pangolin[[1]][1:length(Pangolin[[1]])])
Bat_SL_CoVZC45_seq <- unlist(Bat_SL_CoVZC45[[1]][1:length(Bat_SL_CoVZC45[[1]])])

# Calculate trimer count of all
covid_tri <- oligo_counter(covid_seq ,trimers,3)
sars_tri <- oligo_counter(sars_seq ,trimers,3)
mers_tri <- oligo_counter(mers_seq ,trimers,3)
bat_tri <- oligo_counter(bat_seq ,trimers,3)
bat_rat_tri <- oligo_counter(bat_ratg13,trimers,3)
camel_tri <- oligo_counter(camel_seq ,trimers,3)
civet_tri <- oligo_counter(civet_seq ,trimers,3)
ebola_tri <- oligo_counter(ebola_seq ,trimers,3)
hedgehog_tri <- oligo_counter(hedgehog_seq ,trimers,3)
hiv_tri <- oligo_counter(hiv_seq ,trimers,3)
malaria_tri <- oligo_counter(malaria_seq ,trimers,3)
pangolin_tri<- oligo_counter(Pangolin_seq,trimers,3)
Bat_CovZC45_tri <- oligo_counter(Bat_SL_CoVZC45_seq,trimers,3)



#Calculate tetramer count of all
covid_tetra <- oligo_counter(covid_seq ,tetramers,4)
sars_tetra <- oligo_counter(sars_seq ,tetramers,4)
mers_tetra <- oligo_counter(mers_seq ,tetramers,4)
bat_tetra <- oligo_counter(bat_seq ,tetramers,4)
bat_rat_tetra <- oligo_counter(bat_ratg13,tetramers,4)
ebola_tetra <- oligo_counter(ebola_seq ,tetramers,4)
hiv_tetra <- oligo_counter(hiv_seq ,tetramers,4)
malaria_tetra <- oligo_counter(malaria_seq ,tetramers,4)
pangolin_tetra <- oligo_counter(Bat_SL_CoVZC45_seq,tetramers,4)
Bat_CovZC45_tetra <- oligo_counter(Bat_SL_CoVZC45_seq,tetramers,4)


# create a dataframe for plotting graph
install.packages(gdata)
library(gdata)

## merge tetramers

tetramer_merge<- cbind(covid_tetra,sars_tetra$Freq,mers_tetra$Freq,bat_tetra$Freq,bat_rat_tetra$Freq,hiv_tetra$Freq,malaria_tetra$Freq,pangolin_tetra$Freq,Bat_CovZC45_tetra$Freq)

tetramer_plot <- merge(tetramer_merge,ebola_tetra,by="tetramers",all.x=TRUE)

#convert NA to zero
tetramer_plot[is.na(tetramer_plot)] <- 0


## only 5 viruses
tetramer_merge_5<- cbind(covid_tetra,sars_tetra$Freq,bat_rat_tetra$Freq,pangolin_tetra$Freq,Bat_CovZC45_tetra$Freq)



names(tetramer_merge)[1]<- "tetramers"
names(tetramer_merge)[2]<- "cov"
names(tetramer_merge)[3]<- "sars"
names(tetramer_merge)[4]<- "mers"
names(tetramer_merge)[5]<- "bat"
names(tetramer_merge)[6]<- "bat_ratg"
names(tetramer_merge)[7]<- "hiv"
names(tetramer_merge)[8]<- "mal"
names(tetramer_merge)[9]<- "pangolin"
names(tetramer_merge)[10]<- "Bat_CovZC45"
names(tetramer_plot)[11]<- "ebola"




#Merge Trimers

trimers_merge <- cbind(covid_tri,sars_tri$Freq,mers_tri$Freq,bat_tri$Freq,bat_rat_tri$Freq,ebola_tri$Freq,hiv_tri$Freq,malaria_tri$Freq,pangolin_tri$Freq,Bat_CovZC45_tri$Freq)

# Normalised freqeuncy of k? =
# Number of occurrences of k? / total number of k-mers
trimers_merge[,-1] = apply(trimers_merge[,-1],2,function(x){x/sum(x)})
names(trimers_merge)[1]<- "trimers"
names(trimers_merge)[2]<- "cov"
names(trimers_merge)[3]<- "sars"
names(trimers_merge)[4]<- "mers"
names(trimers_merge)[5]<- "bat"
names(trimers_merge)[6]<- "bat_ratg"
names(trimers_merge)[7]<- "ebola"
names(trimers_merge)[8]<- "hiv"
names(trimers_merge)[9]<- "mal"
names(trimers_merge)[10]<- "pangolin"
names(trimers_merge)[11]<- "Bat_CovZC45"




## only 5 viruses

trimer_merge_5 <- data.frame(trimers=trimers,cov=trimers_merge$cov,sars=trimers_merge$sars,bat_ratg=trimers_merge$bat_rat,pangolin=trimers_merge$pangolin,Bat_CovZC45=trimers_merge$Bat_CovZC45)
##Draw the stacked barchart
#trimers
library(reshape2)
plot<- melt(trimers_merge,id.var = "trimers")
# melt for only 3 genomes

plot5 <- melt(trimer_merge_5 ,id.var = "trimers")

library(ggplot2)

#Bar plot
trimer_bar_plot <- ggplot(plot, aes(x = trimers, y = value, fill = variable)) + 
  geom_tile(stat = "identity")+ theme(axis.text.x = element_text(angle = 90))
trimer_bar_plot

#Line Plot
trimer_line_plot <- ggplot(plot, aes(x = trimers, y = value,group = variable, color = variable)) + 
  geom_line()+ theme(axis.text.x = element_text(angle = 90))
trimer_line_plot


# line plot for 5 genomes
trimer_line_plot5 <- ggplot(plot5, aes(x = trimers, y = value,group = variable, color = variable)) + 
  geom_line()+ theme(axis.text.x = element_text(angle = 90))
trimer_line_plot5

#tetramers

# Normalised freqeuncy of k? =
# Number of occurrences of k? / total number of k-mers
tetramer_plot[,-1] = apply(tetramer_plot[,-1],2,function(x){x/sum(x)})

plot_tetra <- melt(tetramer_plot,id.var = "tetramers")

## for 5 genomes

tetramer_merge_5[,-1] = apply(tetramer_merge_3[,-1],2,function(x){x/sum(x)})
plot_tetra_5 <- melt(tetramer_merge_5,id.var = "tetramers")
names(tetramer_merge_5)[1]<- "tetramers"
names(tetramer_merge_5)[2]<- "cov"
names(tetramer_merge_5)[3]<- "sars"
names(tetramer_merge_5)[4]<- "bat_ratg"
names(tetramer_merge_5)[5]<- "pangolin"
names(tetramer_merge_5)[6]<- "Bat_CovZC45"


## line plot
tetra_line_plot5 <- ggplot(plot_tetra_5, aes(x = tetramers, y = value,group = variable, color = variable)) + 
  geom_line()+ theme(axis.text.x = element_text(angle = 90))
tetra_line_plot5


tetra_cov_mal_1 <- data.frame(tetramer_count$tetramers,tetramer_count$cov,tetramer_count$mal)
plot1 <- melt(tetra_cov_mal_1,id.var = "tetramer_count.tetramers")


#Line Plot
tetra_line_plot <- ggplot(plot_tetra, aes(x = tetramers, y = value,group = variable, color = variable)) + 
  geom_line()+ theme(axis.text.x = element_text(angle = 90))
tetra_line_plot

# for GC content
#Option 1
covid_GC <- GC(covid_seq)
sars_GC <- GC(sars_seq)
mers_GC <- GC(mers_seq)
bat_GC <- GC(bat_seq)
bat_ratg13 <- GC(bat_ratg13)
ebola_GC <- GC(ebola_seq)
hiv_GC <- GC(hiv_seq)
malaria_GC <- GC(malaria_seq)
pangolin_GC <- GC(Pangolin_seq)
Bat_CovZC45_GC <- GC(Bat_SL_CoVZC45_seq )


#plot for GC content
GC_content <- c(covid_GC,sars_GC,mers_GC,bat_GC,bat_ratg13,ebola_GC,hiv_GC,malaria_GC,pangolin_GC,Bat_CovZC45_GC)
GC_df <- data.frame("Names"= character(10),"GC"=GC_content)

GC_df<- edit(GC_df)
GC_df

plot_GC <- melt(GC_df,id.var = "Names")
plot_GC <- ggplot(GC_df, aes(x = Names, y = GC)) + 
  geom_bar(stat = "identity")+ theme(axis.text.x = element_text(angle = 90))




## Pairwise Sequence alignment
# First use a fixed substitution matrix

covid_str <- paste(covid_seq , collapse="")
hiv_str <- paste(hiv_seq , collapse="")
ebola_str <- paste(ebola_seq , collapse="")
malaria_str <- paste(malaria_seq , collapse="")
mers_str <- paste(mers_seq , collapse="")
sars_str <- paste(sars_seq , collapse="")
ratg_str <- paste(bat_ratg13_seq , collapse="")
pangolin_str <- paste(Pangolin_seq, collapse="")
Bat_CoVZC45_str <- paste(Bat_SL_CoVZC45_seq ,collapse="")

corona_dnastring <- DNAString(covid_str)
hiv_dnastring <- DNAString(hiv_str)
ebola_dnastring <- DNAString(ebola_str)
malaria_dnastring <- DNAString(malaria_str )
mers_dnastring <- DNAString(mers_str)
sars_dnastring <- DNAString(sars_str)
bat_ratg_dnastring <- DNAString(ratg_str)
pangolin_dnastring <- DNAString(pangolin_str)
Bat_CoVZC45_dnastring <- DNAString(Bat_CoVZC45_str)


mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
corona_hiv <-
    pairwiseAlignment(corona_dnastring,hiv_dnastring, substitutionMatrix = mat,
                      gapOpening = 5, gapExtension = 2)
corona_malaria <-
    pairwiseAlignment(corona_dnastring,malaria_dnastring, substitutionMatrix = mat,
                      gapOpening = 5, gapExtension = 2)
corona_sars <-
    pairwiseAlignment(corona_dnastring,sars_dnastring, substitutionMatrix = mat,
                      gapOpening = 5, gapExtension = 2)
corona_ebola <-
    pairwiseAlignment(corona_dnastring,ebola_dnastring, substitutionMatrix = mat,
                      gapOpening = 5, gapExtension = 2)
corona_mers <-
    pairwiseAlignment(corona_dnastring,mers_dnastring, substitutionMatrix = mat,
                      gapOpening = 5, gapExtension = 2)
 corona_bat_ratg <-
    pairwiseAlignment(corona_dnastring,bat_ratg_dnastring, substitutionMatrix = mat,
                      gapOpening = 5, gapExtension = 2)
corona_pangolin <-
    pairwiseAlignment(corona_dnastring,pangolin_dnastring, substitutionMatrix = mat,
                      gapOpening = 5, gapExtension = 2)
corona_bat_CoVZC45 <-
    pairwiseAlignment(corona_dnastring,Bat_CoVZC45_dnastring , substitutionMatrix = mat,
                      gapOpening = 5, gapExtension = 2)







### Perform a multiple sequence alignment
### read in Rel sequences from the file
# from package Biostrings
library(Biostrings)
dna <- readDNAStringSet("D:/work/projects/R/genomic analysis/datasets_fasta/cov2.fasta", format="fasta",
               use.names=TRUE)


## Perform a multiple sequence alignment
library(msa)
msa <- msaClustalW(dna,cluster="upgma", gapOpening=1, 
                gapExtension= 1, maxiters= 16, 
                substitutionMatrix="clustalw",type="dna", 
                order=c("aligned", "input"), verbose=FALSE,
                help=FALSE)
### Turn your alignment into a tree
# convert the alignment for the seqinr package
myseq <- msaConvert(msa, type="seqinr::alignment")


##From the manual for the seqinr package
# generate a distance matrix using seqinr package
library(seqinr)
d <- dist.alignment(myseq, "identity")

# have a look at the output
as.matrix(d)


# generate the tree with the ape package
# the nj() function allows neighbor-joining tree estimation
library(ape)
myTree <- nj(d)

# plot the tree
plot(myTree, main="Phylogenetic Tree")


# for upgma 
install.packages("phangorn")
library(phangorn)

tree <- upgma(d)
plot(tree)




