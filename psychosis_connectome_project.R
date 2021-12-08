# ==============================================================================
# title           : psychosis connectome
# description     : connectivity of lesions associated with psychosis
# author          : team
# date            : 2021/10/31
# version         : 1
# usage           : RStudio version 1.4.1106
# notes           : # But: répliquer les études de "lesion network mapping" avec un dataset largement nouveau, et avec un connectome structurel plutôt que fonctionnel
# Hypothèse générale: les différentes lésions cérébrales associées à un même symptômes peuvent être retracées à une région commune par leur réseau de connectivité
# ==============================================================================

rm(list = ls())
dev.off()
setwd("/Users/raphaelbourque/desktop/projet psychosis connectome")

# Load
library("tm")
library("SnowballC")
library("wordcloud")
library("RColorBrewer")
library("car")
library("ggplot2")

# Load patient data: n = 59 patients with stroke associated with psychotic symptoms
# Source: Stangeland, H., Orgeta, V., & Bell, V. (2018). Poststroke psychosis: a systematic review. Journal of Neurology, Neurosurgery & Psychiatry, 89(8), 879-885.
patients <- read.delim2("patient_list.csv",header=TRUE, sep=";", dec = ",",na.strings = c("NA", ""), comment.char = "", check.names = FALSE)

# Load connectome: average structural connectome from n = 842 general population subjects
# The connectome is parcellated in 65 regions based on AAL parcellation
# Source: Yeh, F. C., Panesar, S., Fernandes, D., Meola, A., Yoshino, M., Fernandez-Miranda, J. C., ... & Verstynen, T. (2018). Population-averaged atlas of the macroscale human structural connectome and its network topology. Neuroimage, 178, 57-68.
connectome <- read.delim2("connectogram.csv",header=TRUE, sep=";", dec = ",",na.strings = c("NA", ""), comment.char = "", check.names = FALSE)
rownames(connectome) <- connectome[,1] 
connectome <- connectome[,-1] 

# For sensitivity analysis: varying the threshold of binarization, i.e. at which regions are considered as "connected"
# T value 5,7,9,11 (cf. )
# "When deciding which threshold to use for the ‘primary’ analysis, we choose the highest threshold that still demon- strates the relevant finding, as this will improve specificity."

# ==============================================================================
# Textual analysis of psychotic symptoms
# ==============================================================================
# Reference: http://www.sthda.com/english/wiki/text-mining-and-word-cloud-fundamentals-in-r-5-simple-steps-you-should-know
docs <- Corpus(VectorSource(patients$symptoms[-c(1,2)]))
docs <- tm_map(docs, removeWords, stopwords("english"))
docs <- tm_map(docs, removeWords, c("and","with","like","syndrome","symptoms")) 
docs <- tm_map(docs, stripWhitespace)
dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
head(d, 10)
set.seed(1234)
wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words=200, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"),scale = c(3, 0.2))

barplot(d[1:10,]$freq/59*100, las = 2, names.arg = d[1:10,]$word,
        col ="lightblue", main ="Most prevalent symptoms",
        ylab = "Symptom prevalence (%)",cex.names=0.7)

# ==============================================================================
# Connectivity of lesions
# ==============================================================================
# For each patient, project lesion to connectome
patients_brain <- patients[3:61,4:68]
patients_brain[] <- lapply(patients_brain,as.numeric)
connectome <- as.matrix(connectome)
projections <- c() # to store the results
for (i in 1:59){
  print(paste("working on patient",i))
  this_patient_lesion <- as.numeric(patients_brain[i,])
  this_patient_projections <- connectome %*% this_patient_lesion # multiply matrix by vector
  C1 <- c(this_patient_projections)
  projections <- rbind(projections,C1)
  colnames(projections) <- c(rownames(this_patient_projections))
}
projections <- data.frame(projections)

# ==============================================================================
# Descriptive statistics
# ==============================================================================

# Patients age
age <- as.numeric(patients$age[-c(1,2)])
hist(age,col ="lightblue",ylab="Number of patients")
min(patients$age[-c(1,2)])
max(patients$age[-c(1,2)])
median(patients$age[-c(1,2)])

# Lesion localization
lesions <- apply(patients_brain,2,sum)
index <- order(lesions,decreasing = TRUE)
lesions <- lesions[index]
region_names <- colnames(patients_brain)
region_names <- region_names[index]

barplot(lesions/59*100, las = 2, names.arg = region_names,
        col ="lightblue", main ="Most prevalent lesion localizations",
        ylab = "Localization prevalence (%)",cex.names=0.3)

# Lesion connectivity
binarized_projections <- projections
binarized_projections[projections>0] <- 1
lesions <- apply(binarized_projections,2,sum)
index <- order(lesions,decreasing = TRUE)
lesions <- lesions[index]
region_names <- colnames(patients_brain)
region_names <- region_names[index]

barplot(lesions[1:65]/59*100, las = 2, names.arg = region_names[1:65],
        col ="lightblue", main ="Most prevalent lesion connectivity",
        ylab = "Connectivity prevalence (%)",cex.names=0.3)
# ==============================================================================
# Test association of lesion connectivity with specific symptoms
# ==============================================================================
# Add data concerning symptoms and control variables

projections$delusion <- rep(0,length(nrow(projections)))
projections$delusion[grep("delusion",patients$symptoms[-c(1,2)])] <- 1
projections$delusion <- as.factor(projections$delusion)

projections$hallucination <- rep(0,length(nrow(projections)))
projections$hallucination[grep("hallucination",patients$symptoms[-c(1,2)])] <- 1
projections$hallucination <- as.factor(projections$hallucination)

projections$othello <- rep(0,length(nrow(projections)))
projections$othello[grep("Othello",patients$symptoms[-c(1,2)])] <- 1
projections$othello <- as.factor(projections$othello)

projections$visual <- rep(0,length(nrow(projections)))
projections$visual[grep("visual",patients$symptoms[-c(1,2)])] <- 1
projections$visual <- as.factor(projections$visual)
# 
# projections$parasit <- rep(0,length(nrow(projections)))
# projections$parasit[grep("parasit",patients$symptoms[-c(1,2)])] <- 1
# projections$parasit <- as.factor(projections$parasit)



projections$age <- patients$age[-c(1,2)]

# Modèle linéaire prédisant l'atteinte pour chaque région cérébrale en fonction des symptômes les plus communs et corrigeant pour les autres
results <- c()
for (i in c(1:65)){ # for each of 65 brain regions
  projections$outcome <- projections[,i] # select region
  mod1<- lm(outcome ~
              delusion+
              hallucination+
               othello+
               visual+
               as.numeric(age),
             data = projections)
  
 # VIF <- vif(mod1) # Multicolinearity check and print if applies
 # if (any(VIF>=10) == TRUE) {print (paste("Check for multicolinearity, VIF =",VIF))}
  estimate <- coef(mod1)
  intervals <- confint(mod1)
  p_values <- summary(mod1)$coefficients[,4]
  labels <- names(estimate)
  C1 <- cbind(estimate,intervals,labels,colnames(projections)[i],p_values)
  results <- rbind(results,C1)
  colnames(results) <- c("effect_size","lower","upper","symptoms","brain_region","p_values")
}
results <- data.frame(results)
results[,c(1:3,6)] <- lapply(results[,c(1:3,6)],as.numeric)
# results$p_values <- p.adjust(results$p_values, method = "bonf", n = length(results$p_values))
results$significance[results$p_values>0.05] <- NA
results$significance[results$p_values<0.05] <- "p<0.05"
results$significance[results$p_values<0.01] <- "p<0.01"
results$significance[results$p_values<0.001] <- "p<0.001"
results$significance <- as.factor(results$significance)

results <- results[(results$symptoms == "(Intercept)")==FALSE,]
results <- results[(results$symptoms == "as.numeric(age)")==FALSE,]
# Plot
ggplot() +
  geom_bar(data=results,aes(x=brain_region, y=effect_size), color="grey", position=position_dodge(width=0.9),stat="identity") +
  geom_point(data=results,aes(x=brain_region, y=3, color=significance), position=position_dodge(width=0.9),stat="identity") +
  geom_errorbar(data=results,aes(x=brain_region, ymin=lower, ymax=upper), position=position_dodge(width=0.9),stat="identity")+
  facet_wrap(~ symptoms, nrow = 2, scales="free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Interprétation: parmi les patients ayant des symptômes psychotiques après une lésion causée par un AVC (stroke), 
# la connectivité de lésion avec la région Y est significativement plus élevée chez les patients ayant le symptôme X que chez ceux n'ayant pas ce symptôme,
# en maintenant constant l'âge, et les autres symptômes plus communs.

