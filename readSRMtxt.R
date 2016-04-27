read_txt_SMR<-function(file){
  spectra_file <- file(file, open="r")
  TI_matrix <- matrix(data=NA, nrow=1, ncol=1)
  while (length(oneLine <- readLines(spectra_file, n=1)) > 0){
    if (str_trim(unlist(strsplit(oneLine, ":"))[1]) == "index"
        & suppressWarnings(as.numeric(unlist(strsplit(oneLine, ":"))[2])!=0)){
      TI_matrix <- cbind(TI_matrix, NA);
      numCol <- dim(TI_matrix)[2];
      id <- unlist(strsplit(str_trim(oneLine<-readLines(spectra_file, n=1)), " "))[2:5];
      dimnames(TI_matrix)[[2]][numCol]<-paste(id[1], id[3], id[4], sep=" ");
      times <- suppressWarnings(
        as.numeric(unlist(strsplit(unlist(strsplit(oneLine <- (readLines(spectra_file, n=6))[6], "] "))[2], " "))));#time
      Intensity <- suppressWarnings(
        as.numeric(unlist(strsplit(unlist(strsplit(oneLine <- (readLines(spectra_file, n=3))[3], "] "))[2], " "))));#intensity
      if (max(Intensity) > 0){
        tempMi <- matrix(c(times,  Intensity), nrow=length(times), ncol=2)
        tempM <- tempMi[tempMi[,2] != 0,]
        TI_matrix <- Build_Int_M(TI_matrix, tempM);
      }
    }
  }
  colnames(TI_matrix)[1] <- "Time"
  TI_matrix <- TI_matrix[-1,]
  close(spectra_file);
  
  return(list(sample_name = scan_features(file)[1],
              timestamp = scan_features(file)[2],
              transitions = dim(TI_matrix)[2]-1,
              TI_matrix = TI_matrix[order(TI_matrix[,1]),]))
}

Build_Int_M <- function(poolM, tempM){
  numCol <- dim(poolM)[2]
  if (length(tempM) == 2){
    if (tempM[1] %in% poolM[,1]){
      poolM[which(tempM[1] == poolM[,1]), numCol] <- tempM[2];
      return (poolM)}
    else {
      poolM <- rbind(poolM, c(tempM[1], rep(NA, numCol-2), tempM[2]));
      return(poolM)
    }
  }
  else {
    if (tempM[1,1] %in% poolM[,1]){
      poolM[which(tempM[1,1] == poolM[,1]), numCol] <- tempM[1,2];
      return(Build_Int_M(poolM, tempM[-1,]))
    }
    else {
      poolM <- rbind(poolM, c(tempM[1,1], rep(NA, numCol-2), tempM[1,2]))
      return(Build_Int_M(poolM, tempM[-1,]))} 
  }
}

scan_features <- function(files){
  spectra_file <- file(files, open="r")
  sample_name <- str_trim(unlist(strsplit(readLines(spectra_file,n=2)[2], ":")))[2]
  timestamp <- str_trim(unlist(strsplit(readLines(spectra_file,n=48)[48], ":")))[2]
  close(spectra_file);
  return (c(sample_name,timestamp))
}

is.SMRM <- function(file){
  spectra_file <- file(file, open="r")
  a <- FALSE
  i <- 1
  while (i < 20){
    oneLine <- readLines(spectra_file, n = 1)
    i <- 1 + i
    if (str_trim(unlist(strsplit(oneLine, ":"))[1])=="cvParam"){
      if( str_trim(unlist(strsplit(oneLine, ":"))[2])=="selected reaction monitoring chromatogram")
      {a <- TRUE}
    }
  }
  close(spectra_file)
  return(a)
}


getSRM <- function(x) UseMethod("getSRM")
getSRM.default <- function(x){
  if (is.SMRM(x) == F){
    print ("This is not a SRM text file")
  }
  
  else {
    read_chr <- read_txt_SMR(x)
    read_chr$call <- match.call()
    class(read_chr) <- "getSRM"
    read_chr
  }
  
}

print.getSRM <- function(x){
  cat("Call:\n")
  print(x$call)
  cat(paste("\nSample name: ", x$sample_name))
  cat(paste("\nNumber of transitions: ", x$transitions))
  cat(paste("\nTime range: ", as.character(round(x$TI_matrix[1,1],2)), "-", as.character(round(x$TI_matrix[ dim(x$TI_matrix)[1],1],2)))) 
  cat("\n")
}


