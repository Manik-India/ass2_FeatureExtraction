#To stop warnings that come during library() runs use following 2 lines
defaultW <- getOption("warn")
options(warn = -1)

#import required libraries
library(Peptides)
library(peptider)
library(log4r)

#Creating Log File
file.create("101903077.log")
cat("File created")
my_logfile = "101903077.log"
my_console_appender = console_appender(layout = default_log_layout())
my_file_appender = file_appender(my_logfile, append = TRUE,layout = default_log_layout())

my_logger <- log4r::logger(threshold = "INFO",appenders = list(my_console_appender,my_file_appender))

log4r_info <- function() {
  log4r::info(my_logger, "Number of Arguments will be 2. Example: Rscript <filename> input.csv")
}

log4r_error <- function() {
  log4r::error(my_logger, "File not found!!!!")
}

#to take file name which we pass in command line as an argument
args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 1){
  log4r_info()
}

tryCatch({mydata <- read.csv(args[1])
data("AAdata")  #A list with a collection of properties, scales and indices for the 20 naturally occurring amino acids from various sources
i = 0
for(Sequence in mydata$Peptide.Sequence)
{
  i = i+1
  cat("Completed: ",i,"/", nrow(mydata), "\n")
  
  Length <- lengthpep(Sequence) #for length of peptide sequence
  
  No_of_neighbour<-getNofNeighbors(Sequence, blosum = 1, method = "peptide", libscheme = NULL)

  Aliphatic_index <- aIndex(Sequence) #For Aliphatic_Index
  
  Boman_index <- boman(Sequence)  #For Boman_Index
  
  Instability_index <- instaIndex(Sequence)   #For Instability_Index
  
  Hmoment_index <- hmoment(Sequence, angle = 160, window = 11) #for Hydorphobic_moment
  
  Prob_peptide <- ppeptide(Sequence, libscheme = "NNK", N = 10^8) #probability of detection of peptide Sequence
  
  Mass_shift <- massShift(Sequence, label = "silac_13c")   #For Mass_difference
  
  Auto_correlation_index <- autoCorrelation(Sequence, lag = 1, property = AAdata$Hydrophobicity$KyteDoolittle, center = TRUE)
  
  Aacomp <- as.numeric(unlist(aaComp(seq = Sequence)))
  
  Tiny_class = Aacomp[1]
  
  Small_class = Aacomp[2]
  
  Aliphatic_class = Aacomp[3]
  
  Aromatic_class = Aacomp[4]
  
  Nonpolar_class = Aacomp[5]
  
  Polar_class = Aacomp[6]
  
  Charged_class = Aacomp[7] #names(Aacomp)<- paste("Aacomp_",c(1:length(Aacomp)),sep="")
  
  Basic_class = Aacomp[8]
  
  Acidic_class = Aacomp[9]
  
  Peptide_Charge <- charge(Sequence, pH = 7, pKscale = "Dawson") #for peptide_charge (we can change value of pKscale)
  
  Hydrophobicity <- hydrophobicity(Sequence, scale = "BlackMould")  #For Hydrophobicity (we can change the value of scale)
  
  Molecular_weight <- mw(Sequence, monoisotopic = FALSE, avgScale = "expasy", label = "none", aaShift = NULL) #For molecular_weight
  
  Isoelecric_Point <- pI(Sequence, pKscale = "EMBOSS")  #For iso-electric point (We can change pKscale)
  
  kidera <- as.numeric(unlist(kideraFactors(seq = Sequence)))
  
  Kidera_Factor = kidera[1]
  
  Target <- mydata$Target[i]
  
  result <- data.frame(Sequence,Length, No_of_neighbour, Aliphatic_index, Boman_index, Instability_index, Hmoment_index, Prob_peptide,Mass_shift, Auto_correlation_index,
                       Tiny_class, Small_class, Aliphatic_class, Aromatic_class, Nonpolar_class, Polar_class, Charged_class, Basic_class, Acidic_class, Peptide_Charge,
                       Hydrophobicity, Molecular_weight, Isoelecric_Point, Kidera_Factor,Target)
  finall <- c(result)
  if(i==1){
    write.table(finall, "output-101903077.csv", sep = ",", row.names = F, col.names = T, append = T)}
  else
  {
    write.table(finall, "output-101903077.csv", sep = ",", row.names = F,col.names = F, append = T)
  }
  
}
cat("Done")},error = function(e) {
  log4r_error()
})