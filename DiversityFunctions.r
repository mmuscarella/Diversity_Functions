### Tracy and Mario's Codes of Subsampling Stuff
### Created 22 Aug 2012
### Last Update 0 Oct 2012 - Muscarella

# system( ) #mothur command version

# Example files used for testing
# input_file <- "KBS.an.shared" 

# Needed packages
require(vegan)

# Read OTU File
read.otu <- function(shared = " ", cutoff = "0.03"){
  matrix <- read.delim(shared, head=T)
  matrix.cutoff <- subset(matrix, matrix$label == cutoff)
  matrix.out <- as.matrix(matrix.cutoff[1:dim(matrix.cutoff)[1],4:(3+mean(matrix.cutoff$numOtus))])
  row.names(matrix.out) <- matrix.cutoff$Group
  return(matrix.out)
  } 

# Count All Groups 
count.groups <- function(otu.matrix = " "){
  counts <- rowSums(otu.matrix)
  return(counts)
  }  

# Subsampling wrapper
sub.sample <- function(otu.matrix = " ", sample.size = "min(count.groups(test))"){
  counts <- count.groups(otu.matrix) 
  statement <- counts > sample.size  # Add warning message
  otu.matrix <- subset(otu.matrix, rowSums(otu.matrix)>sample.size)
  x <- rrarefy(otu.matrix, sample.size)
  return(x)
  }

# Example function for calculating species richness  
richness.iter <- function(shared = " ", cutoff = " ", size = " ", iters = " "){
  otu.matrix <- read.otu(shared, cutoff)
  counts <- count.groups(otu.matrix) 
  statement <- counts > size  # Add warning message
  otu.matrix <- subset(otu.matrix, rowSums(otu.matrix)>size)
  rich.matrix <- matrix(NA, dim(otu.matrix)[1], iters)
  rownames(rich.matrix) <- rownames(otu.matrix)
  for (i in 1:iters){
    temp.matrix <- sub.sample(otu.matrix, size)
    rich.matrix[,i] <- rowSums((temp.matrix>0)*1)
    }
  return(rich.matrix)
  } 
  
# Example function for calculating diversity
diversity.iter <- function(shared = " ", index = "shannon", cutoff = " ", size = " ", iters = " "){
  otu.matrix <- read.otu(shared, cutoff)
  counts <- count.groups(otu.matrix) 
  statement <- counts > size  # Add warning message
  otu.matrix <- subset(otu.matrix, rowSums(otu.matrix)>size)
  div.matrix <- matrix(NA, dim(otu.matrix)[1], iters)
  rownames(div.matrix) <- rownames(otu.matrix)
  for (i in 1:iters){
    temp.matrix <- sub.sample(otu.matrix, size)
    div.matrix[,i] <- diversity(temp.matrix, index)  
    }
  return(div.matrix)
  }
  
# Example function for calculating species evenness (Pielou's Evenness = J' = H'/H'max = H'/lnS, H' = Shannon)
evenness.iter <- function(shared = " ", cutoff = " ", size = " ", iters = " "){
  otu.matrix <- read.otu(shared, cutoff)
  counts <- count.groups(otu.matrix)
  statement <- counts > size  # Add warning message
  otu.matrix <- subset(otu.matrix, rowSums(otu.matrix)>size)
  even.matrix <- matrix(NA, dim(otu.matrix)[1], iters)
  rownames(even.matrix) <- rownames(otu.matrix)
  for (i in 1:iters){
    temp.matrix <- sub.sample(otu.matrix, size)
    even.matrix[,i] <- diversity(temp.matrix, "shannon")/log(rowSums((temp.matrix>0)*1))
    }
  return(even.matrix)
  } 
  
