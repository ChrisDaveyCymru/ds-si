#Selection index R-script

#IMPORTANT
#Note that the code to generate the G and P matrices used for the index calculations
#are in the block of code that calculates rg from the genetic covariance. This
#applies to the situation where traits were measured in the same year.
#To use this code for the index calculations ensure all traits are from the same year.




#Remove all the objects currently in R to ensure there are no clashing variable names.
#The print should display "character(0)" if this has worked.
#print("Check that character(0) is displayed below to ensure R has been reset correctly.")
#rm(list=ls(all=TRUE))
#print(ls())


#Load packages used.
library(lme4)

 
 
#---------------------------------------------------------------------------------------------
#The data file has 2 types of columns. There are columns that contain support info such as genotype, block,
#column, row, altitude, lat. etc and the columns containing actual trait data. There is a separate row for each plant monitored, 
#the top row contains R-friendly column names.
#
#The support info columns can be in any order but must be a contiguous block of columns. 
#However, if the following info is provided use the standardised column names given in brackets: 
#genotype (geno), block/rep (block), block row (row), block column (col), species (species).
#This makes specifying the mixed models more standardised, particularly important for geno where this name is explicitly
#used in down-stream calcs.
#
#The actual trait data columns must also be contiguous.
#This block of trait data can occur anywhere in the file's column set: i.e. at the start, in the middle or at the end of
#the set of columns. The trait data can be of any type with any name but each year's data must have a separate column. 
#Thus, the same trait measured on the same plant in different years will have separate columns. This is done because 
#rg uses a different formula for cross-year data. 
#
#Missing data must be represented by NA in both the trait and support info columns.
#
#Note that the code to generate the G and P matrices used for the index calculations
#are in the block of code that calculates rg from the genetic covariance. This
#applies to the situation where traits were measured in the same year.
#To use this code for the index calculations ensure all traits are from the same year.

#---------------------------------------------------------------------------------------------
#To use the code.
#There are demarcated blocks of code to do all these jobs below.
#(1) Input file name of the trait data file.
#(2) Set the variables holding the column numbers of the start and end of the block of support info data in the input file.
#(3) Set the variables holding the column numbers of the start and end of the block of trait data in the input file.
#(4) Set the file names for the output data files (will need to uncomment the code lines that 
#    write the data as well, these are at the very end of the code).
#(5) There are 3 places below which use linear mixed models (5 models in all, marked by a row of "*"s).
#    It is important to ensure the  model specification is the same for all the models.
#(6) Edit the input data set to include only the traits needed for the index and to put the trait of
#    interest for selection as the first trait column (then adjust (2) and (3) appropriately).
#(7) Update the vector which holds the years of each trait.
#--------------------------------------------------------------------------------------------- 
#Input trait data file name.
trait.support.input <- "si8_pheno_data_08.csv"

#----------------------------------------------------------------------------------------------
#Input the column number of the start and end of the block of support info data in the input file.
#If you remove unneeded columns of trait data below then the values here must be modified to
#reflect this so that the column numbers reflect the data actually analysed.
info.start <- 1
info.end <- 13

#Input the column number of the start and end of the block of trait data in the input file.
#If you remove unneeded columns of trait data below then the values here must be modified to
#reflect this so that the column numbers reflect the data actually analysed.
trait.start <- 14
trait.end <- 27

#----------------------------------------------------------------------------------------------
#Set names of output data files.
#You will need to un-comment the actual write statements at the end of the code.

#Matrix of rg values. 
#Diagonal are sqrt(H2). Off diagonal are rg. If years are not the same
#rg calculated from rp. If years are the same rg is from the 
#genetic covar.
#r.g.output <- "rg.csv"

#Matrix of rg and rp values.
#Diagonal are sqrt(H2). rp above diagonal. Below diagonal are rg. 
#If years are not the same rg calculated from rp. 
#If years are the same rg is from the genetic covar.
#r.g.p.output <- "rgrp.csv"

#Matrix of rg blup values.
#rg from correlations of BLUPs.
#rg.blup.output <- "rgblup.csv"


#Matix of all blup values.
#Has index BLUPs also.
blup.all.output <- "si8_xxx_blup_all.csv"

#Matrix of H2 values.
#Has H2 for index also.
h2.output <- "si8_xxx_h2.csv"


#Matrix of the info columns from the trait data
#along with the equivalent index values.
#trait.info.i.output <- "info.csv"

#--------------------------------------------------------------------------------------------- 
#Load trait data into program. 
#Do not touch this code.
data.trait <- read.table(file = trait.support.input, header = TRUE, sep = ",", stringsAsFactor=FALSE)

#---------------------------------------------------------------------------------------------
#At this point remove trait columns not needed for the index you are creating but
#leave in info columns that are needed.
#Ensure all traits were measured in the same year.
#Also put the trait to be selected (eg dry matter) as the first trait column.
#Do this by using the c() of column indexes and their order in c(). 
#Then adjust the values for the start and end of the info and trait columns above
#so that the column numbers reflect the data actually analysed.
data.trait <- data.trait[, c(1:13, 20, 14, 15, 16, 17)]


#Ensure the year associated with each of the remaining traits is entered in the data frame below.
#The dates are in the same order as the trait columns. Thus the first trait corresponds to
#the first year in c().
data.year <- data.frame(year = as.integer( c(2008, 2008, 2008, 2008, 2008) ))

#---------------------------------------------------------------------------------------------
#Check input data and that year data is int or other easily comparable type.
print("Check input data, ensure year data is int.")
print("Input year data")
str(data.year)
View(data.year)

print("Trait data to use. Ensure trait of interest for selection is the 1st trait column.")
str(data.trait)
View(data.trait)

#---------------------------------------------------------------------------------------------
#This is the end of what the user needs to change other than un-commenting the file save code
#and checking the linear mixed models.
#---------------------------------------------------------------------------------------------

#Calc trait dimensions.
num.traits <- as.integer( trait.end - trait.start + 1 )
#Includes all the reps of each genotype
num.rows <- nrow(data.trait)

#Create vector of unique genotypes
unique.geno <- unique(data.trait$geno)
#Find number of unique genotypes
unique.geno.num <- length(unique.geno)

#Generate a vector of trait names.
traits <- colnames(data.trait)[trait.start:trait.end]

#Create a blank matrix for H2 value of each trait.
H2 <- matrix(NA, num.traits, 1)
rownames(H2) <- traits
colnames(H2) <- "H2"

#Make a copy of the original input data and add a column to store
#the column of data to currently analyse.
data.trait.1 <- cbind(data.trait, curr.trait = data.trait[,trait.start])

#Create matrices for the blup values.
#Note these matrices have a single row for each unique genotype.
blup <- matrix(NA, unique.geno.num, num.traits)
blup.all <- data.frame(geno = unique.geno)



#Calculate the H2 value for each trait.
for (i in 1:num.traits){
   #Put current trait data in special column of data frame.
   #Dynamic data typing deals with changes in data (eg int to num).
   data.trait.1$curr.trait <- data.trait.1[, (trait.start - 1) + i ]
   
   #******************************
   #Generate linear mixed effects model. If you change the model here then change the
   #equivalent models used later in the code.
   lmem.1 <- lmer(curr.trait ~ (1|geno) + (1|block) + (1|col) + (1|row), data = data.trait.1)   
   
   #Important to use geno as col name.
   var.g.1 <- VarCorr(lmem.1)$geno[1,1]
   var.e.1 <- ( summary(lmem.1)$sigma )^2
   
   H2[i] <- var.g.1/(var.g.1 + var.e.1)
 
   #Get current traits blup values. 
   blup.i <- cbind( rownames( ranef(lmem.1)$geno ), ranef(lmem.1)$geno )
   colnames(blup.i) <- c("geno", traits[i])
   
   #Combine data with blup.all, including all the (unique) genotypes in blup.all in the result.
   blup.all <- merge(blup.all, blup.i, by="geno", all.x=T)
   
}#end for


#View(blup.i)
#View(blup.all)
#View(H2)


#Get genetic correlations as correlations of blups
r.g.blup <- cor( blup.all[,2:(num.traits+1)], use="pairwise.complete.obs" )
rownames(r.g.blup) <- traits
colnames(r.g.blup) <- traits

#View(r.g.blup)

#Initialise the G and P matrices specified in the index selection paper:
#Mackay 2015.
G <- matrix(NA, num.traits, num.traits)
rownames(G) <- traits
colnames(G) <- traits
P <- matrix(NA, num.traits, num.traits)
rownames(P) <- traits
colnames(P) <- traits



#Looped code below calculates rg from the genetic covariance or rp.
#Create the symmetrical matrices to hold the trait x trait results.
#Contains both rg (below diagonal) and rp (above diagonal) data.
r.g.p <- matrix(NA, num.traits, num.traits)
rownames(r.g.p) <- traits
colnames(r.g.p) <- traits

#Contains only rg both above and below diagonal.
r.g <- matrix(NA, num.traits, num.traits)
rownames(r.g) <- traits
colnames(r.g) <- traits

#To the copy of the original data add three extra columns to hold the trait
#data indexed by the loop indexes i and j.
data.trait.1 <- cbind(data.trait.1, trait.i = data.trait.1[,trait.start], trait.j = data.trait.1[,trait.start], trait.i.plus.j = data.trait.1[,trait.start] )



#Loop through the traits creating pairs of traits to test. Two loop counters
#are used (i and j) to select the traits along the matrix rows and columns to
#give a calculated value to store at matrix location i,j. The loops select matrix
#locations only in one of the triangles without including the diagonal. The equivalent
#location for the same pair of traits in the opposite triangle is referenced by swooping
#the indexes i and j when referencing the matrix. Thus [i,j] and [j,i] are matrix locations
#for the same combination of two traits but in opposite triangles.
#Once i and j have been used to store traits in $trait.i and $trait.j ($trait.i.plus.j is also calculated)
#then rp (phenotypic correlation) for the pair of traits is calculated and stored above the diagonal in matrix r.g.p.
#The calculation of rg of the two traits depends on if the year of measurement is the same or not.
#Finally the rg is stored below the diagonal in matrix r.g.p and above and below the diagonal
#in matrix r.g.
#For the loops: i counts through the rows and j the columns for that row, below the diagonal.
#Thus, the first loop location misses the top left corner as this is a diagonal.
for (i in 2:num.traits){
   #Store trait i and year it was measured.
   data.trait.1$trait.i <- data.trait.1[, (trait.start - 1) + i ]
   year.trait.i <- data.year[i, 1]
     
   for  (j in 1:(i-1)){
      #Store trait j and year it was measured.
      data.trait.1$trait.j <- data.trait.1[, (trait.start - 1) + j ]
      year.trait.j <- data.year[j, 1]
   
      #Calculate rp and store above diagonal.
	  rp.value <- cor(data.trait.1$trait.i, data.trait.1$trait.j, use = "pairwise.complete.obs")
      r.g.p[j, i] <- rp.value
   
   
      #Calculate rg depending on whether the two traits were measured in the same year or not.   
      if (year.trait.i != year.trait.j){
         #Two years different: rg = rp/[sqrt(H2i) sqrt(H2j)]: Howe et al. (2000) TAG 101.
		 #Note you will not get the matrices needed for the index if the calculations are
		 #done using this block of code.
         rg.value <- rp.value /( sqrt(H2[i]) * sqrt(H2[j]) )
		 
		 #Store below the diagonal in r.g.p and on both sides in r.g.
		 r.g.p[i, j] <- rg.value
		 r.g[i, j] <- rg.value
		 r.g[j, i] <- rg.value		 
      } else {
         #Two years the same. rg is calculated using genetic covariance.	  
		 #This code block is used to calculate the matrices for the index.
         #******************************
		 #Three linear mixed effects models are used here. The basic models must be the same as for the H2 calcs above.
         lmem.i <- lmer(trait.i ~ (1|geno) + (1|block) + (1|col) + (1|row), data = data.trait.1)   	  
         lmem.j <- lmer(trait.j ~ (1|geno) + (1|block) + (1|col) + (1|row), data = data.trait.1)   	  	  
         data.trait.1$trait.i.plus.j <- data.trait.1$trait.i + data.trait.1$trait.j 	  
         lmem.i.plus.j <- lmer(trait.i.plus.j ~ (1|geno) + (1|block) + (1|col) + (1|row), data = data.trait.1)   	  	  	  
		 
         #Extract the genotypic variances from the models.
         var.g.i <- VarCorr(lmem.i)$geno[1,1]		 
         var.g.j <- VarCorr(lmem.j)$geno[1,1]	  
         var.g.i.plus.j <- VarCorr(lmem.i.plus.j)$geno[1,1]	  
		 
		 #Extract the environmental variances from the models.
		 var.e.i <- (summary(lmem.i)$sigma)^2
		 var.e.j <- (summary(lmem.j)$sigma)^2
		 var.e.i.plus.j <- (summary(lmem.i.plus.j)$sigma)^2
		 
         #Calculate covariances.
         covar.g.ij <- (var.g.i.plus.j - var.g.i - var.g.j)/2		
		 covar.e.ij <- (var.e.i.plus.j - var.e.i - var.e.j)/2

         #Calculate rg value (Howe et al. (2000) TAG 101).
	     rg.value <- covar.g.ij / ( sqrt(var.g.i) * sqrt(var.g.j) )

		 #Store below the diagonal in r.g.p and on both sides in r.g.
		 r.g.p[i, j] <- rg.value
		 r.g[i, j] <- rg.value
		 r.g[j, i] <- rg.value		 	

         #Update the matrices needed for the index.
		 G[i,j] <- covar.g.ij
		 G[j,i] <- covar.g.ij
		 G[i,i] <- var.g.i
		 G[j,j] <- var.g.j
		 P[i,j] <- covar.g.ij + covar.e.ij
		 P[j,i] <- covar.g.ij + covar.e.ij
		 P[i,i] <- var.g.i + var.e.i
		 P[j,j] <- var.g.j + var.e.j 	 
	  }#End if.
   
   }#End j for.
}#End i for.


#Set the diagonals in the two matrices.
diag(r.g.p) <- sqrt(H2) 
diag(r.g) <- sqrt(H2)

#View(r.g.p)
#View(r.g)


#Calculate the selection indices as in Mackay 2015.
#Note this uses the formulae in Falconer p327-328 (eq 19.15).
#It is assumed the trait of interest for selection is the 1st column
#of trait data in data.trait after it has been modified to remove 
#traits not needed by the index.
#Create matrix e of economic values. Only the 1st trait (i.e. the one for
#selection) has a value of 1, the others are 0.
e <- matrix(0, num.traits, 1)
e[1,1] <- 1
print("Matrix e of economic values")
print(e)

#RHS intermediate matrix
R <- G %*% e
print("Matrix R: RHS intermediate matrix")
print(R)

#Find b matrix of weighting factors for the actual index.
b <- solve(P,R)
print("Matrix b of index weighting factors")
print(b)


#Calculate relative response of index selection vs single trait selection.
#Calc SD of index (Falconer Eq 19.18, p329).
rhs.sum <- 0
for  (ti in 1:num.traits){ 
   rhs.sum <- rhs.sum + (b[ti] * R[ti])
}#End for.
sigma.i <- sqrt(rhs.sum)
print("sigma i: SD of index")
print(sigma.i)

#Genotypic SD of selected trait.
sigma.a <- sqrt(G[1,1])
print("sigma a: genotypic SD of selected trait")
print(sigma.a)

#Relative response to same selection intensity. Combines Falconer
#eq 19.17 (p326) and eq 11.4 (p189).
R.rel <- sigma.i / (sqrt(H2[1]) * sigma.a)
print("-----------")
print("Relative response of index selection vs single trait selection.")
print(R.rel)



#Calc index for each plant in the trait data.
#This will give a separate index for each genotype rep.
Index <- c(rep(NA, times=num.rows))

#Make a data frame containing only the trait data.
only.traits <- data.trait[, c(trait.start:trait.end)]

#Loop around each plant and calc its index.
for  (pl.cnt in 1:num.rows){ 
   #This plants index.
   pl.I <- 0
   
   #Loop across all the traits to get b by trait terms.
   for  (tr.cnt in 1:num.traits){ 
      pl.I <- pl.I + ( b[tr.cnt] * only.traits[pl.cnt, tr.cnt] )
   }#End for tr.cnt

   Index[pl.cnt] <- pl.I   
}#End for pl.cnt.


#Augment data frames with indexes.
trait.info.I <- cbind(data.trait[,c(info.start:info.end)], index=Index) 
data.trait <- cbind(data.trait, index=Index)


#******************************
#Get H2 and BLUPs for index.
lmem.bi <- lmer(index ~ (1|geno) + (1|block) + (1|col) + (1|row), data.trait)
var.g.bi <- VarCorr(lmem.bi)$geno[1,1]
var.e.bi <- (summary(lmem.bi)$sigma)^2

H2.bi <- var.g.bi/(var.g.bi + var.e.bi)
H2 <- rbind(H2, index=H2.bi)

blup.bi <- cbind(rownames(ranef(lmem.bi)$geno), ranef(lmem.bi)$geno)
colnames(blup.bi) <- c("geno", "index")
blup.all <- merge(blup.all, blup.bi, all.x=T)



#Write data files.
#Un-comment the write statements to save the files.
#Matrix of rg
#write.csv(r.g, r.g.output)

#Matrix of rg and rp values.
#write.csv(r.g.p, r.g.p.output)

#Matrix of rg blup values.
#write.csv(r.g.blup ,rg.blup.output) 


#These 2 files are the ones to uncomment to save the data.
#Matix of all blup values.
#write.csv(blup.all, blup.all.output) 

#Matrix of H2 values.
#write.csv(H2, h2.output) 




#Matrix of trait info data and equivalent index values.
#write.csv(trait.info.I, trait.info.i.output)

print("-----Run complete-----")


