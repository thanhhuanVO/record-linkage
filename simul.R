
# Packages
library(RecordLinkage)
library(stringr)
library(ggplot2)

rm(list = ls())

# Initial data
df1 <-  read.csv("https://query.data.world/s/fophzgeexhxmh66cyjscdsusgagote", header=TRUE, stringsAsFactors=FALSE);
fname_list <- df1$name[1:10000];
pcode_list <- c(35000:35050)

# This function using to create error for string (name, address,...)
str_error <- function(string, s){
  # string
  # s : maximum proportion of characters having errors for each string
  ns <- nchar(string) # length of string
  nse <- round(s*ns) # maximum number of error-character
  if (nse <1){
    nse <- 1
  }
  nse <- sample(1:nse,1) # random number of error-character
  
  # We account for 4 types of string error
  # 1: insertion (huan to huaan)
  # 2: deletion (huan to hun)
  # 3: subsitution (huan to huam)
  # 4: transposition (huan to haun)
  # err_type <- sample(c(1,2,3,4),1,prob = c(0.2,0.2,0.3,0.3))
  
  pos <- sample(1:ns, nse) # random index of error-character
  
  for (i in 1:nse){
    
    err_type <- sample(c(1,2,3,4),1,prob = c(0.2,0.2,0.3,0.3)) 
    
    if (err_type == 1){
      temp <- str_sub(string, pos[i],ns)
      str_sub(string, pos[i],pos[i]) <- sample(LETTERS,1)
      str_sub(string, pos[i]+1,ns+1) <- temp
    }else if (err_type == 2){
      str_sub(string,pos[i], pos[i]) <- ""
      
    }else if (err_type == 3){
      str_sub(string,pos[i],pos[i]) <- sample(LETTERS,1)
      
    }else if (err_type == 4){
      temp <- str_sub(string,pos[i],pos[i])
      
      if (pos[i] == ns){
        str_sub(string,pos[i],pos[i]) <- str_sub(string,pos[i]-1, pos[i]-1)
        str_sub(string,pos[i]-1,pos[i]-1) <- temp
      }else{
        str_sub(string,pos[i],pos[i]) <- str_sub(string,pos[i]+1, pos[i]+1)
        str_sub(string,pos[i]+1,pos[i]+1) <- temp
      }
    }
  }
  return(string)
}

make_error <- function(field, type, e, s=0.2){
  # Fields: fname, lname, sex, age,...
  # Type of records: string(=1), binary (=2), category (=3), continuous (=4)
  # e: maximum proportion of records having errors
  # s: maximum proportion of characters having errors for each string !
  
  n <- length(field);
  ne <- round(e*n); # maximum number of records containing errors
  
  ne <- sample(1:ne,1) # random number of records containing errors
  
  pos <- sample(1:n, ne) # random error index
  if (type==1){
    for (i in 1:ne){
      temp <- field[pos[i]];
      field[pos[i]] <- str_error(temp, s)
    }
  }else if (type==2){
    
    field[pos] <- 1-field[pos]
  }else if (type== 3){
    
    field[pos] <- sample(unique(field),ne, replace = TRUE)
  }else if (type ==4){
    # This rule for age!
    field[pos] <- field[pos] + sample(c(-3,-2,-1,1,2,3), ne, replace = TRUE)
  }
  
  Z <- rep(0,n)
  Z[pos] <- 1
    
  return(list("field"= field,"err_id" = Z))
  #return(field)
}

# Generate data function
gen_data <- function(nA, nB, e, s=0.2 ){
  # nA, nB: sample size of A and B
  # e: maximum proportion of records containing error
  # s: (for string only)  maximum proportion of characters having errors for each string 
  
  # Records of A
  id <- 1:nA
  fname <- sample(fname_list, size = nA, replace = TRUE)
  sex <- sample(c(0,1), size = nA, replace = TRUE)
  age <- round(rnorm(nA,50,6))
  pcode <- sample(pcode_list, size = nA, replace = TRUE)  
  
  A <- data.frame(id, fname,sex,age,pcode, stringsAsFactors = FALSE)
  
  # Records of B

  B_id <- c(A$id,(nA+1):nB)
  B_fname <- c(A$fname, sample(fname_list, size = nB-nA, replace = TRUE))
  
  B_sex <- c(A$sex, sample(c(0,1), size = nB-nA, replace = TRUE))
  B_age <- c(A$age, round(rnorm(nB-nA,50,6)))
  B_pcode <- c(A$pcode, sample(pcode_list, size = nB-nA, replace = TRUE))
  
  B <- data.frame(B_id, B_fname, B_sex, B_age, B_pcode, stringsAsFactors = FALSE)
  colnames(B) <- c("id", "fname", "sex", "age", "pcode")
  #Shuffle
  B <- B[sample(nrow(B)),]
  
  # Create errors for A
  
  temp <- make_error(A$fname,1,e,s)
  A$fname <- temp$field
  A$err_fname <- temp$err_id
  
  temp <- make_error(A$sex,2,e,s)
  A$sex <- temp$field
  A$err_sex <- temp$err_id
  
  temp <- make_error(A$age,4,e,s)
  A$age <- temp$field
  A$err_age <- temp$err_id
  
  temp <- make_error(A$pcode,3,e,s)
  A$pcode <- temp$field
  A$err_pcode <- temp$err_id
  
  # Create errors for B
  temp <- make_error(B$fname,1,e,s)
  B$fname <- temp$field
  B$err_fname <- temp$err_id
  
  temp <- make_error(B$sex,2,e,s)
  B$sex <- temp$field
  B$err_sex <- temp$err_id
  
  temp <- make_error(B$age,4,e,s)
  B$age <- temp$field
  B$err_age <- temp$err_id
  
  temp <- make_error(B$pcode,3,e,s)
  B$pcode <- temp$field
  B$err_pcode <- temp$err_id
  
  return(list("A" = A,"B" = B))
}


# Record linkage using Fellegi-Sunter model function
FS <- function(nA,nB,e,s=0.2,lambda, mu){
 
  data <- gen_data(nA,nB,e,s)
  
  A <- data$A
  B <- data$B
  
  # Using package record linake
  rPairs <- compare.linkage(A,B,exclude = c(1,6,7,8,9), strcmp = T, strcmpfun = levenshteinSim, identity1 = A$id, identity2 = B$id )
  #summary(rPairs)
  
  w1 <- emWeights(rPairs, cutoff = 0.7)
  #summary(w1)
  
  result <- emClassify(w1, my = 0.05, ny = 0.1)
  #summary(result)
  
  temp1 <- result$prediction == "L"
  temp2 <- result$pairs
  
  accuracy <- temp1 == temp2$is_match
  
  return(sum(accuracy)/length(accuracy))
}


# Evaluation
# Vary sample size, fixed e
nA <- 5
nB <- 6
e <- 0.2
s <- 0.2

L <- 10 # number of simulated databases
p <- matrix(0,30,2)

for (j in 1:3){
  nA <- 2*nA
  nB <- 2*nB
  p[((j-1)*L):(j*L),1] <- nA
  for (i in 1:L){
    p[L*(j-1)+i,2] <- FS(nA, nB,e,s)
  }
}

p <- data.frame(p)
colnames(p) <- c("nA","accuracy")

ggplot(p, mapping = aes(x = nA, y = accuracy, fill = nA, group = "x")) + 
  geom_boxplot(outlier.color = "black", outlier.shape = 8, outlier.size = 2)



### test