#Load drm package which is used for the 5-parameter logistic curve optimization
library(nCal)
library(vioplot)
library(gridExtra)
library(grid)

######################################################################################
####The functions used to build the bigger algorithm##################################

#Calculate coefficient of variation for each concentration/dosage
cov <- function(data, conc, mfi){
  #Get the target variable to calculate %CV
  target.no <- grep(colnames(data), pattern = conc)
  mfi.index <- grep(colnames(data), pattern = mfi)
  
  #Check to see if the other variable is a numeric or factor
  check <- is.numeric(data[,target.no])
  if (check == TRUE){
    #If it was numeric find unique values
    dose <- unique(data[,target.no])
  } else {
    #If it was factor find the unique levels
    dose <- levels(data[,target.no])
  }
  
  #calculate each %CV for a concentration/dosage
  nlvl <- length(x = dose)
  
  #Declare the return matrix that contains the means, stdev, and %CV
  ret <- data.frame(array(dim=c(nlvl,4)))
  name <- colnames(data)[mfi.index]
  colnames(ret) <- c("Conc", "mean", "stdev", "CV")
  
  #Calculate the %CV
  i = 1
  for (i in 1:nlvl){
    #put the dosage/concentration to the first column of return matrix
    ret[i,1] <- dose[i]
    
    #Grab the replicates for each dosage/concentration
    replicates <- data[which(data[, target.no] == dose[i]), mfi.index]
    replicates.val <- as.numeric(unlist(replicates))
    #Calculate the mean and stdev and %CV
    mean <- mean(replicates.val)
    std <- sd(replicates.val)
    cv <- std/mean*100
    ret[i,c(2:4)] <- c(mean, std, cv)
  }
  
  #Return the matrix as function
  return(ret)
}

#Function to prune the standard had too much variance, user-defined or algorithm picked
#data = data that used in the coev function, and coev = result from the coev function
#Limit = CV% allowed for each standards, for CAP 20% precision or CoV is standard limit
prune.stds.cv <- function(data, coev, limit=20){
  
  #Get the average of the %CV of the standard, limit is not chosen
  
  if (is.null(limit) == TRUE) {
    mean.cv <- mean(coev$CV) + 3 * sd(coev$CV) #Find the outlier values for CoV based on the CV of MFI measured
  } else {
    mean.cv <- limit
  }
  
  #Get rid of the standard had %CV above the limit
  outliers <- coev[which(coev$CV > mean.cv),1]
  
  #Get the name of dosage/concentration from coev and find it in data
  dosage <- grep(pattern = colnames(coev)[1], colnames(data))
  
  #Eliminate all standards that had CV above the limit
  if (length(outliers) == 0) {
    data <- data
  } else {
    for (j in 1:length(outliers)){
      data <- data[-which(data[,dosage] == outliers[j]),]
    }
  }
  
  data <- data[which(data[, dosage] != 0), ]
  data <- data[order(-data[,dosage]),]
  return(data)
}

#Detect single replicates meaning only one standard sample was found
# Rufei's code
#singles <- function(list){
#  uniques = unique(list)
#  count = unlist(lapply(X = list, FUN = function(x) sum(list == x)))
#  sin = cbind(uniques, count)
#  singles.out = sin[sin[, "count"] == 1, "uniques"]
#  return (singles.out)
#}

# Hua's code 
singles <- function(list){
  duplicates<-list[duplicated(list)]
  singles.out<-list[!list %in% duplicates]
  return (singles.out)
}

#Prune standards that are lower than LOD which is defined as LOD function
#If one of the duplicate of standards fall out of range, then all duplicates were eliminated
prune.stds.lod <- function(data, expected = "Conc", measured, lod, plate.id = "plate.id"){
  exp.id = grep(pattern = expected, x = colnames(data))
  mea.id = grep(pattern = measured, x = colnames(data))
  
  #eliminate the samples below LOD
  if (length(which(data[, mea.id] < lod)) == 0) {
    data = data
  } else {
    data = data[-which(data[, mea.id] < lod), ]
  }
 
  
  #eliminate single standard replicate after prune
  plates = unique(data[, plate.id])
  output = data.frame()
  for (i in 1:length(plates)){
    temp = data[data[, plate.id] == plates[i], ]
    list = temp[temp[, plate.id] == plates[i], exp.id]
    single = singles(list = list)
    if (length(single) == 0) {
      temp = temp
    } else {
      single.index = unlist(lapply(X = single, FUN = function(x) grep(pattern = x, x = temp[, exp.id])))
      temp = temp[-c(single.index), ]
    }
    
    output = rbind(output, temp)
  }

  return(output)
}

#This Function gets rid of the uneven standard ranges. So the most common range is preserved whereas
#the high ranges were eliminated
std.range <- function(data, x, y, group){
  #Get the x, y, and group column number
  if (is.numeric(x) == TRUE) {
    n.x <- x
  } else { n.x <- grep(pattern = x, x = colnames(data))}
  
  if (is.numeric(y) == TRUE) {
    n.y <- y
  } else { n.y <- grep(pattern = y, x = colnames(data))}
  
  if (is.numeric(group) == TRUE) {
    group <- group
  } else { group <- grep(pattern = group, x = colnames(data))}
  
  plate.id <- unique(data[, group])
  
  hi <- array(dim = c(length(plate.id), 2))
  for (i in 1:length(plate.id)){
    hi[i, 1] <- plate.id[i]
    hi[i, 2] <- max(data[which(data[, group] == plate.id[i]), n.x])
  }
  
  ret <- data.frame()
  cutoff <- min(hi[,2]) * 1.5
  for (i in 1:length(plate.id)){
    std.temp <- data[which(data[, group] == plate.id[i] & data[, n.x] < cutoff),]
    ret <- rbind(ret, std.temp)
  }
  
  return(ret)
}

#Calculating LOB, according to CAP LOB = mean(blanks) + 1.645 * sd(blanks)
LOB.single <- function(blanks, n = 1.645){
  mean.blanks = mean(blanks) #Average of blanks, the blank sample across plates can be used here
  sd.blanks = sd(blanks) #Standard deviation of blanks
  LOB = mean.blanks + (n * sd.blanks) #LOB based on CAP definition, it's 1.645*SD
  return(c(blanks.ave = mean.blanks, blanks.sd = sd.blanks, LOB = LOB)) 
}

#Calculate LOB across plates
# Rufei's Code
#LOB.acrossplates <- function(data, id, assay, blank.pattern = "Blank.*$"){
  #Hua's Code
 LOB.acrossplates <- function(data, id, assay, blank.pattern){
  blanks.rowindex <- grep(pattern = blank.pattern, x = data[,id])
  assay.colindex <- grep(pattern = assay, x = colnames(data))
  blanks = na.omit(data[blanks.rowindex, assay])
  LOB = LOB.single(blanks = blanks)
  return(LOB)
}

#Extract standards for each assay for each plate or all, can be used to extract measure and expected 
#from its perspective file
stds.extract <- function(data, id, assay, std = "Std", plate.id = NULL, plate.name = NULL) {
  pattern.std = paste0(std, ".*$")
  id.index = grep(pattern = id, x = colnames(data))
  plateid.colindex= grep(pattern = "plate.id", x = colnames(data))
  stds.index = grep(pattern = pattern.std, x = data[, id])
  assay.colindex = grep(pattern = assay, x = colnames(data))
  stds = data[stds.index, c(id.index, assay.colindex, plateid.colindex)]
  if (is.null(plate.id) == FALSE){
    plateid.index = grep(pattern = plate.id, x = stds[, "plate.id"])
    stds = stds[plateid.index, ]
  }
  
  if (is.null(plate.name) == FALSE){
  stds = stds[stds[, "plate.id"] == plate.name, ]
  }
  
  return(stds)
}

#Combine expected and known values in one file
StandardCurve.file <- function(data.master, id, assay, std = "Std", platenames){
  combined.master = data.frame()
  master = data.frame()
  for (i in 1:length(platenames)){
    #extract individual plate measured and expected values
    measured = stds.extract(data = data.frame(data.master[[1]]), id = id, std = std, assay = assay, plate.name = platenames[i])
    expected = stds.extract(data = data.frame(data.master[[2]]), id = id, std = std, assay = assay, plate.name = platenames[i])
    
    combined.master = measured
    Conc = c(1:nrow(combined.master))
    combined.master = cbind(combined.master, Conc)
    for (j in 1:nrow(combined.master)){
      combined.master[j, "Conc"] = expected[c(expected[, id] == as.character(combined.master[j, id])), assay]
    }
    
    master = rbind(master, combined.master) 
  }

  return(master)
}

#Determine LOD based on CAP approved method which is LOB + 1.645 * sd(low concentration samples)
LOD.cal <- function(std.files, expected = "Conc", measured, LOB, sd.blank = "yes"){
  exp.id = grep(pattern = expected, x = colnames(std.files))
  mea.id = grep(pattern = measured, x = colnames(std.files))
  
  expected.min = min(std.files[, exp.id])
  lowConc = std.files[c(std.files[, exp.id] == expected.min), mea.id]
  if (sd.blank == "yes") {
    LOD = LOB[[3]] + 1.645 * LOB[[2]]
  } else {
    LOD = LOB[[3]] + 1.645 * sd(lowConc)
  }

  return(c(LOD = LOD, LowConc.sd = sd(lowConc)))
}

#Combine expected known standards with measured MFI from those known values
stdCurve <- function(stdCurve.mod){
  colnames(stdCurve.mod)[2] = "MFI"
  colnames(stdCurve.mod)[4] = "Conc"
  colnames(stdCurve.mod)[3] = "Plateid"
  stdCurve.mod = data.frame(stdCurve.mod)
  stdCurve.mod[, "MFI"] = log10(stdCurve.mod[, "MFI"])
  stdCurve.mod[, "Conc"] = log10(stdCurve.mod[, "Conc"])

  model<-drm(MFI~Conc, data=stdCurve.mod, curveid = stdCurve.mod[, "Plateid"], fct=LL.5(), type="continuous", robust="mean", logDose = 10)
  return(model)
}

#Method 1 of LOQ, whihc is Signal to Noise ratio of 10, which means LOQ = LOB + 10 * sd(blanks)
LOQ.1 <- function(LOB, LOD){
  LOQ = LOB[3] + 10 * LOD[2]
  return(LOQ)
}

#Method 2 of LOQ, which is LOQ = LOD + 1.645 * sd(low conc)
LOQ.2 <- function(LOD){
  LOQ = LOD[1] + 1.645 * LOD[2]
  return(LOQ)
}
  
#A specific function to extract the IDs from the columan names
extract <- function (pattern, x){
  ret <- sub(pattern, "\\1", x)
  return(ret)
}

#Calculate the concentration of the assay based on the MFI
calculate.conc <- function(model, MFI, platename){
  #Extract calculated parameter from the model
  coefficients = model$coefficients
  coeff.names = names(coefficients)
  platename.index =c(grep(pattern = paste0(".*?",platename), coeff.names))
  plate.coeff = coefficients[platename.index]  
  b = plate.coeff[1]
  c = plate.coeff[2]
  d = plate.coeff[3]
  e = plate.coeff[4]
  f = plate.coeff[5]
  
  #Calculate the concentration based on the equation provided by the model 
  #f(x) = c + \frac{d-c}{(1+\exp(b(\log(x)-\log(e))))^f}
  MFI.log = log10(MFI)
  conc = exp((log(((d - c)/(MFI.log - c))^(1/f) - 1))/b)*e
  return(conc)
}

#Calculate the MFI based on the 5-parameter logisitc curve
calculate.mfi <- function(b, c, d, e, f, Conc){
  MFI = c + (d - c)/((1 + exp(b*log(Conc/e)))^f)
  return(MFI)
}

#Calculate coefficient of variance
coev <- function(values){
  coev = sd(values) / mean(values)
  return(coev)
}

#Method 3 of LOQ, where CAP accept the <20% CV%
LOQ.3 <- function(stdsCurve.files, platenames, cv.limit = 0.2){
  stdsCurve.files[, "calculated.conc"] = log10(stdsCurve.files[, "calculated.conc"])

  #Calculate CV% for each unique expected concentration
  exp.conc = unique(stdsCurve.files[, "Conc"])
  CV = c(1:nrow(stdsCurve.files))
  stdsCurve.files = cbind(stdsCurve.files, CV)
  for (i in 1:length(exp.conc)){
    cal.conc.temp = stdsCurve.files[stdsCurve.files[, "Conc"] == exp.conc[i], "calculated.conc"]
    if (length(cal.conc.temp) > 2) {
      cv = abs(sd(cal.conc.temp)/mean(cal.conc.temp)) #Calculate the CV%
      stdsCurve.files[stdsCurve.files[, "Conc"] == exp.conc[i], "CV"] = round(cv, 4)
    } else {
      stdsCurve.files[stdsCurve.files[, "Conc"] == exp.conc[i], "CV"] = NA
    }
  }
  
  #Determine the lowest expected concentration as LOQ 
  temp <- na.omit(stdsCurve.files)
  temp = temp[temp[, "CV"] < cv.limit, ]
  lloq = min(temp[, "Conc"])
  uloq = max(temp[, "Conc"])
  return(list(lloq, uloq, stdsCurve.files)) 
}

#Method 4 of LOQ, where LOQ is defined by an acceptable difference between groups, since the threshold is
#flexible, this can be adapted into any recursive analysis for downstream
LOQ.4 <- function(stdsCurve.files, model.l5, platenames, LOQ, a.2 = 1.645, b = 1.645, delta, n = 1){
  Calculated.Conc = c(1:nrow(stdsCurve.files))
  stdsCurve.files = cbind(stdsCurve.files, Calculated.Conc)
  
  #Back calculate all concentration based on the 5-parameter logisitc curve
  for (i in 1:length(platenames)){
    stdCurve.temp = stdsCurve.files[stdsCurve.files$plate.id == platenames[i], 2]
    stdCurve.conc = lapply(X = stdCurve.temp, FUN = function(x) calculate.conc(model = model.l5, platename = platenames[i], MFI = x))
    stdCurve.conc = log10(unlist(stdCurve.conc))
    names(stdCurve.conc) = NULL
    stdsCurve.files[stdsCurve.files$plate.id == platenames[i], ncol(stdsCurve.files)] = stdCurve.conc
  }
  stdsCurve.files[, 5] = log10(stdsCurve.files[, 5])
  stdsCurve.files[, 2] = log10(stdsCurve.files[, 2])
  
  #Calculate CV% for each unique expected concentration
  exp.conc = unique(stdsCurve.files[, 5])
  CV = c(1:nrow(stdsCurve.files))
  Accuracy = c(1:nrow(stdsCurve.files))
  sd.conc = c(1:nrow(stdsCurve.files))
  stdsCurve.files = cbind(stdsCurve.files, CV, Accuracy, sd.conc)
  for (i in 1:length(exp.conc)){
    cal.conc.temp = stdsCurve.files[stdsCurve.files[, 5] == exp.conc[i], "Calculated.Conc"]
    cv = abs(coev(values = cal.conc.temp)) #Calculate the CV%
    stdsCurve.files[stdsCurve.files[, 5] == exp.conc[i], "CV"] = cv
    accu = abs(sd(cal.conc.temp)/exp.conc[i]) #Accuracy based on sd(calculated concentrations)/expected values
    stdsCurve.files[stdsCurve.files[, 5] == exp.conc[i], "Accuracy"] = accu
    sd.conc = sd(cal.conc.temp) #Accuracy based on sd(calculated concentrations)/expected values
    stdsCurve.files[stdsCurve.files[, 5] == exp.conc[i], "sd.conc"] = sd.conc
  }
  
  #Determine the lowest expected concentration as LOQ  
  sigma = sqrt(n)*delta/(sqrt(2)*((a.2 + b)^2))
  min.temp = na.omit(stdsCurve.files[stdsCurve.files[, "CV"] < sigma, ])
  min.temp = min.temp[min.temp[, 2] > log10(LOQ), ]
  LOQ.exp.conc = min(min.temp[, 5])
  return(c(LOQ.exp.conc, stdsCurve.files)) 
}

#This function is used to read in all the data in to one single file with sample MFI and expected concentrations based on assay
#Directory = directory all the raw files are stored in
#filepattern.samples = a particular pattern that sample files are stored in, for examples: if all files of samples are stored as "Study 1 FI.csv", "Study 2 FI.csv" ... then pattern is "xxx FI.csv", where xxx stands for wild cards
#filepattern.knownValues = a particular pattern that known values or expected values files are stored in, 
#for examples: if all files of samples are stored as "Study 1 Expected Conc.csv", "Study 2 Expected Conc.csv" ... then pattern is "xxx FI.csv", where xxx stands for wild cards

read.to.masterfile <- function(dir, filepattern.samples, filepattern.knownValues){
  setwd(dir) #Set the working directory
  filenames <- dir(dir) #Grab the file names of ALL files in the folder or directory
  pattern.samples <- paste0(sub("(.*?)xxx(.*?)$", "\\1", filepattern.samples),"(.*?)",
                            sub("(.*?)xxx(.*?)$", "\\2", filepattern.samples),"$") #Change the pattern in to R code command, e.g. "(.*?) FI.csv$" 
  pattern.knownValues <- paste0(sub("(.*?)xxx(.*?)$", "\\1", filepattern.knownValues),"(.*?)",
                                sub("(.*?)xxx(.*?)$", "\\2", filepattern.knownValues),"$") #Change the pattern in to R code command, e.g. "(.*?) Expected Values.csv$"
  samples <- grep(pattern = pattern.samples, filenames) #get sample files indexes can be converted to filenames later
  knownValues <- grep(pattern = pattern.knownValues, filenames) #get known value files indexes and can be converted to filenames later
  platenames <- sapply(X = filenames[samples], FUN = extract, pattern = pattern.samples) #Extract plate names as it appears in the directory, e.g. "Study 1 FI.csv", would be extracted as "Study 1" as plate 1 name
  
  masterfile.samples <- data.frame() #declare master file to store all sample data in
  masterfile.knownValues <- data.frame() #declare master file to store all known values in
  
  #Loop to combine all sample files into a masterfile to be used later for calculation
  for (i in 1:length(samples)){
    tempfile.samples <- read.table(file = filenames[samples[i]], header = T, sep = ",") #Extract in to a temp file to be combined into master file
    plate.id = rep(platenames[i], nrow(tempfile.samples))
    tempfile.samples <- cbind(tempfile.samples, plate.id)
    masterfile.samples <- data.frame(rbind(masterfile.samples, tempfile.samples)) #Loop to combine all files
  }
  rm(tempfile.samples) #Clear the memory
  
  #Loop to combine all known values files into a masterfile to be used later for calculation
  for (i in 1:length(knownValues)){
    tempfile.knownValues <- read.table(file = filenames[knownValues[i]], header = T, sep = ",") #Extract in to a temp file to be combined into master file
    plate.id = rep(platenames[i], nrow(tempfile.knownValues))
    tempfile.knownValues <- cbind(tempfile.knownValues, plate.id)
    masterfile.knownValues <- data.frame(rbind(masterfile.knownValues, tempfile.knownValues)) #Loop to combine all files
  }
  rm(tempfile.knownValues) #Clear the memory
  
  return(list(masterfile.samples, masterfile.knownValues, platenames))
}

#Main function to do that pre-stat QC
#filepattern.samples = "xxx FI.csv"
#filepattern.knownValues = "xxx Expected Conc.csv"
#sampleid =  "Barcode"
#dir = "/Users/rufeilu/Desktop/OMRF Analysis/BOLD cytokine data/Combined/Plate A"
#pattern.s = "Std"
#pattern.blank = "Blank.*$"
# Rufei's code
# pre.analysis.qc <- function(filepattern.samples, filepattern.knownValues, dir, sampleid, pattern.s, pattern.blank, no.lob.limit = "yes", lob = "yes"){
# Hua's code  
pre.analysis.qc <- function(filepattern.samples, filepattern.knownValues, dir, sampleid, pattern.s, pattern.blank, no.lob.limit = "yes", lob = "yes", blank.pattern, control.pattern){
  files <- read.to.masterfile(dir = dir, filepattern.samples = filepattern.samples, filepattern.knownValues = filepattern.knownValues)
  #turn all variable into numeric
  temp <- files[[1]][, c(2:(ncol(files[[1]])-1))]
  files[[1]][, c(2:(ncol(files[[1]])-1))] <- apply(X = temp, MARGIN = 2, function(x) as.numeric(as.character(x)))
  temp <- data.frame(files[[2]][, (2:(ncol(files[[1]])-1))])
  files[[2]][, c(2:(ncol(files[[2]])-1))] <- apply(X = temp, MARGIN = 2, function(x) as.numeric(as.character(x)))
  
  
  #####Use method 3 of LOQ, which is lowest concentration < 20% CV%
  platenames = files[[3]]
  assay.names = colnames(files[[1]]) #Getting all assay names
  assay.names = assay.names[-c(grep(pattern = sampleid, x = assay.names), grep(pattern = "plate.id", x = assay.names))] #Keep only assay names, but eliminate sample id and plate id columns
  n.assay = length(assay.names)
  
  #Keep just samples
  # Rufei's code
  # samples = files[[1]]
  # samples.legacy = files[[1]] #keep original files for all MFI values.
  # stds.n = grep(pattern = paste0(pattern.s, ".*$"), samples[, sampleid])
  # samples = samples[-c(stds.n), ]
  # conc.s = samples #Pass to be calculated
  # conc.all = samples.legacy #pass to legacy
  
  # Hua's code
  #Keep just samples
  samples = files[[1]]
  samples.legacy = files[[1]]
  stds.n = grep(pattern = paste0(pattern.s, ".*$"), samples[, sampleid])
  blank.n = grep(pattern = paste0(blank.pattern, ".*$"), samples[, sampleid])
  control.n = grep(pattern = paste0(control.pattern, ".*$"), samples[, sampleid])
  samples = samples[-c(stds.n, blank.n, control.n), ]
  conc.s = samples #Pass to be calculated
  conc.all = samples.legacy #pass to legacy
  
  #Generate blank QC report
  dir.create(file.path(dir, "QC Output"), showWarnings = FALSE) #Create a directory for output
  setwd(file.path(dir, "QC Output"))
  QC.Steps = c("Step 1", "No. Stds > 20% CV", "No. Stds < LOD", "Step 2", "No. Stds outside of 80% to 120% accuracy",
               "Problem Plates", "Step 3", "No. Samples outside of LOQ", "No. Samples < LOD")
  table.report = data.frame(matrix(ncol = n.assay + 1, nrow = 9))
  colnames(table.report) = c("QC Steps", assay.names)
  table.report[, 1] = QC.Steps
  
  pdfname = "QC plots.pdf"
  pdf(file = pdfname, onefile = TRUE, paper = "letter", height = 12)
  for (i in 1: n.assay){
    assayid = assay.names[i]
    
    #Determine the LOB
    # Rufei's code
#    LOB <- LOB.acrossplates(data = data.frame(files[1]), id = sampleid, assay = assayid, blank.pattern = pattern.blank) #Extract the LOB for each assay
    
    #Hua's code
    LOB <- LOB.acrossplates(data = data.frame(files[1]), id = sampleid, assay = assayid, blank.pattern = blank.pattern) #Extract the LOB for each assay
    
    LOB.val = LOB[3]
    
    #Get master standards file
    master.stds <- StandardCurve.file(data.master = files, std = pattern.s, id = sampleid, assay = assayid, platenames = files[[3]])
    no.start.stds = length(unique(master.stds[, sampleid]))
    
    #Determine the LOD
    LOD <- LOD.cal(std.files = master.stds, LOB = LOB, measured = assayid, expected = "Conc", sd.blank = "no")
    
    if (lob == "yes"){
      LOD.val = LOB.val
    } else if (lob == "no"){
      LOD.val = LOD[1]
    }
    
    
    master.stds.prune = data.frame()
    if (no.lob.limit == "yes"){
      master.stds.prune = master.stds #Default option to leave the below dectection low concentration in for curve estimate purposes, only use for bad curves
    } else if (no.lob.limit == "no"){
      master.stds.prune <- prune.stds.lod(data = master.stds, expected = "Conc", measured = assayid, lod = LOD.val) #Cleaner way of eliminate stds that are below the detection, but does eliminate a lot in some cytokines
    }
    
    #Estimate how many samples are below limit of blank or detection
    n.sample = nrow(samples) #Estimate total sample size
    n_lod = round(length(which(samples[, assayid] < LOD.val)) / n.sample, 3) #Estimate fraction of samples are below limit of blank or detection
    
    samples[which(samples[, assayid] < LOD.val), assayid] = NA #Change one below LOB or LOD to NA so won't be calculated
    conc.s[which(conc.s[, assayid] < LOD.val), assayid] = NA #Change one below LOB or LOD to NA so won't be calculated
    table.report[9, i + 1] = n_lod
    
    #Elminate stds with > 20% CV
    no.cv.drop = NULL
    no.lod.drop = NULL
    master.stds.prune.2 = data.frame()
    for (j in 1:length(platenames)){
      plate.stds = master.stds.prune[master.stds.prune[, "plate.id"] == platenames[j], ]
      cv <- cov(data = plate.stds, conc = "Conc", mfi = assayid) #Output matrix contains mean and %CV for each concentration
      plate.stds.prune <- prune.stds.cv(data = plate.stds, coev = cv, limit = 20) #Eliminates %CV > 20%
      
      lod.drop = round(((no.start.stds * 2) - nrow(plate.stds)) / (no.start.stds * 2), 4)
      no.lod.drop = c(no.lod.drop, lod.drop)
      
      drop = round((nrow(plate.stds) - nrow(plate.stds.prune)) / (no.start.stds * 2), digits = 2)
      no.cv.drop = c(no.cv.drop, drop)
      
      master.stds.prune.2 = rbind(master.stds.prune.2, plate.stds.prune)
    }
    table.report[2, i + 1] = paste0(no.cv.drop, collapse = ", ")
    table.report[3, i + 1] = paste0(no.lod.drop, collapse = ", ")
    
    
    #Eliminate stds fall outside of 80% to 120% range of accuracy
    suppressWarnings(model <- stdCurve(master.stds.prune.2))
    recovery = master.stds.prune.2 #Transfer over
    calculated.conc = c(1:nrow(recovery))
    accuracy = c(1:nrow(recovery))
    recovery = cbind(recovery, calculated.conc, accuracy)
    for (j in 1:length(platenames)){
      plate = recovery[which(recovery[, "plate.id"] == platenames[j]), ]
      mfi.rec = recovery[which(recovery[, "plate.id"] == platenames[j]), assayid]
      suppressWarnings(
        conc.rec <- unlist(lapply(X = mfi.rec, FUN = function(x) calculate.conc(model = model, MFI = x, platename = platenames[j])))
      )
      recovery[c(recovery[, "plate.id"] == platenames[j]), "calculated.conc"] = conc.rec
      accuracy = abs((conc.rec - plate[, "Conc"]) / plate[, "Conc"])
      plate[, "accuracy"] = accuracy
      
      #stds.unique = unique(plate[, "Conc"])
      #for (k in 1:length(stds.unique)){
      #  rec.ave = mean(plate[plate[, "Conc"] == stds.unique[k], "accuracy"])
      #  plate[plate[, "Conc"] == stds.unique[k], "accuracy"] = rec.ave
      #}
      
      recovery[c(recovery[, "plate.id"] == platenames[j]), "accuracy"] = plate[, "accuracy"]
    }
    
    #take the stds with average accuracy outside of the 80% to 120% range
    master.stds.prune.3 = na.omit(recovery)
    master.stds.prune.3 = master.stds.prune.3[master.stds.prune.3[, "accuracy"] < 0.2, ]
    
    #Make sure the plate is not dropping out too many stds
    no.drop.final = NULL
    for (j in 1:length(platenames)){
      plate = master.stds.prune.3[c(master.stds.prune.3[, "plate.id"] == platenames[j]), ]
      no.drop.temp = 1 - round(length(unique(plate[, "Conc"])) / (no.start.stds), 2)
      no.drop.final = c(no.drop.final, no.drop.temp)
    }
    table.report[5, i + 1] = paste0(no.drop.final, collapse = ", ")
    
    #All the plate with > 60% drop rate
    p = which(no.drop.final > 0.6)
    problem = platenames[p]
    table.report[6, i + 1] = paste0(problem, collapse = ", ")
    
    #After pruning make the final model
    
    model.final <- stdCurve(stdCurve.mod = master.stds.prune.3)
    
    loq <- LOQ.3(stdsCurve.files = master.stds.prune.3, platenames = platenames, cv.limit = 0.2)
    lloq = loq[[1]]
    uloq = loq[[2]]
    
    #Calculate concentration based on each plates standard curves
    for (j in 1:length(platenames)){
      samples.assay = conc.s[c(conc.s[, "plate.id"] == platenames[j]), assayid]
      samples.legacy.assy = conc.all[c(conc.all[, "plate.id"] == platenames[j]), assayid]
      suppressWarnings(
        sample.concs <- unlist(lapply(X = samples.assay, FUN = function(x) calculate.conc(model = model.final, MFI = x, platename = platenames[j])))  
      )
      suppressWarnings(
        sample.concs.all <- unlist(lapply(X = samples.legacy.assy, FUN = function(x) calculate.conc(model = model.final, MFI = x, platename = platenames[j])))  
      )
      conc.s[which(conc.s[, "plate.id"] == platenames[j]), assayid] = sample.concs
      conc.all[which(conc.all[, "plate.id"] == platenames[j]), assayid] = sample.concs.all
    }
    
    #Calculate how many sample fall outside of 80% to 120%
    n.low = length(which(conc.s[, assayid] < lloq))
    n.hi = length(which(conc.s[, assayid] > uloq))
    n_loq = (n.low + n.hi) / n.sample
    table.report[8, i + 1] = n_loq
    
    conc.s[which(conc.s[, assayid] < lloq), assayid] = lloq
    conc.s[which(conc.s[, assayid] > uloq), assayid] = uloq
    
    #Set up plot panels
    par(mfrow = c(2, 1))
    #Dynamic range set by max and min of calculated concentration and expected concentration
    y.max = log10(max(master.stds[, assayid]) * 1.4)
    y.min = log10(min(master.stds[, assayid]) * 0.8)
    x.max = log10(max(master.stds[, "Conc"]) * 1.4)
    x.min = log10(min(master.stds[, "Conc"]) * 0.8)
    col.pal <- rainbow(length(unique(platenames)))
    plot(model.final, col=col.pal, xlab=paste("log10(", assayid, "Concentration)"), 
         ylab=paste("log10(", assayid, ")MFI"), type="all", xlim = c(x.min, x.max), ylim = c(y.min, y.max))
    #graph points below the detection
    graph.below.lod <- na.omit(data.frame(Conc = log10(conc.all[, assayid]), MFI = log10(samples.legacy[, assayid]), plate.id = conc.all[, "plate.id"]))
    graph.below.lod <- graph.below.lod[which(graph.below.lod$MFI < log10(LOD.val)), ]
    points(x = graph.below.lod[, 1], y = graph.below.lod[, 2], cex = 1, pch = 16, col="black")
    #Plot points calucated based on 5-parameter logitis curve
    graph = na.omit(data.frame(Conc = log10(conc.s[, assayid]), MFI = log10(samples[, assayid]), plate.id = conc.s[, "plate.id"]))
    col.pt <- col.pal[c(as.numeric(graph[, "plate.id"]))] #Set same colors as the standard curves
    points(x = graph[, 1], y = graph[, 2], cex = 1, pch = 16, col=col.pt)
    
    #Calculate LLOQ for LOD for each plate
    lloq.at.lod = NULL
    for (k in platenames){
      suppressWarnings(
        lloq.cal <- calculate.conc(model = model.final, MFI = LOD.val, platename = k)  
      )
      lloq.at.lod <- c(lloq.at.lod, lloq.cal)
    }
    lloq.graph <- min(lloq.at.lod)
    
    abline(v = log10(lloq), col="green", lwd = 2, lty = 2)
    text(log10(lloq), 2, "LLOQ", adj = c(0,0))
    abline(v = log10(uloq), col="green", lwd = 2, lty = 2)
    text(log10(uloq), 2, "ULOQ", adj = c(0,0))
    abline(h = log10(LOD.val), col = "red", lwd = 1, lty = 2)
    text(log10(lloq.graph), log10(LOD.val), "LOD", adj = c(0,0))
    abline(h = log10(LOB[1]), col = "red", lwd = 1.2, lty = 1)
    #text(log10(lloq.graph), log10(LOB[1]), "Blank Average", adj = c(0,0))
    mtext(text = paste0(n_lod*100, "% below LOD"), side = 3)

    #Generate violin plot for each plate
    if (length(na.omit(log10(conc.s[, assayid]))) == 0){
      range = c(0, 1)
    } else {
      range = range(na.omit(log10(conc.s[, assayid])))
    }
    plot(1, 1, xlim = c(0, length(platenames) + 1), ylim = range, xlab = "", ylab = paste0("log10(", assayid, ")"), xaxt = 'n', type = 'n')
    for (k in 1:length(platenames)){
      name = platenames[k]
      x = na.omit(log10(conc.s[which(conc.s[, "plate.id"] == name), assayid]))
      if (length(x) == 0) {
        k = k + 1
      } else {
        vioplot(x, names = name, add = T, at = k, col = col.pal[k])
        label = paste0(name, " (n=", length(x), ")")
        axis(1, at = k, labels = label, las = 2, cex.axis = 0.6)
      }
    }
    
    #Output progress
    print(paste0(signif(i/n.assay*100, digits = 4),"% Complete!"))
  }
  graphics.off()
  
  write.table(x = table.report, file = "QC Step by Step Report.csv", sep = ",", na = "", row.names = FALSE)
  
  pdf(file = "QC Undetected Rate.pdf", onefile = T, paper = "letter", height = 10)
  #display the bar plot with missing rates
  par(mfrow = c(1,1))
  x = as.numeric(table.report[9,-1])*100
  col = rep(x = "green", length(x))
  col[which(x > 70)] <- "red"
  legend <- c(paste0("N of passed = ", sum(x<70)),
              paste0("N of failed = ", sum(x>70)))
  barplot(height = x, names.arg = colnames(table.report)[-1], cex.names = 0.7, las = 2, main = "Percentage below LOD",
          ylab = "% Undetected", col = col)
  legend("topright", legend = legend, fill = c("green", "red"))
  abline(h = 60, col = "orange")
  abline(h = 70, col = "red", lwd = 2)
  
  #automatically generate QC report for each analyte
  for (p in colnames(table.report)[-1]){
    #Step 1
    step.1 <- table.report[c(2:3), c("QC Steps", p)]
    step.1.table <- data.frame(matrix(0, nrow = 2, ncol = (length(platenames) + 1)))
    colnames(step.1.table) <- c("Std QC Step", platenames)
    step.1.table[1, ] <- unlist(c(step.1[1, 1], strsplit(step.1[1, 2], split = ","))) 
    step.1.table[2, ] <- unlist(c(step.1[2, 1], strsplit(step.1[2, 2], split = ",")))  
    theme <- ttheme_default(base_size = 10)
    step.1.grid <- tableGrob(step.1.table, rows = NULL, cols = gsub("\\ ", "\\\n", colnames(step.1.table)), theme = theme)
    
    #Step 2
    step.2 <- table.report[c(5:6), c("QC Steps", p)]
    step.2.table <- data.frame(matrix(0, nrow = 1, ncol = (length(platenames) + 1)))
    colnames(step.2.table) <- c("Accuracy QC Step", platenames)
    step.2.table[1, ] <- unlist(c(step.2[1, 1], strsplit(step.2[1, 2], split = ","))) 
    step.2.grid <- tableGrob(step.2.table, rows = NULL, cols = gsub("\\ ", "\\\n", colnames(step.2.table)), theme = theme)
    
    #Step 3
    step.3 <- table.report[c(8:9), c("QC Steps", p)]
    step.3[, 2] <- signif(as.numeric(step.3[, 2]), 3) 
    step.3.grid <- tableGrob(step.3, rows = NULL, theme = theme)
    
    grid.arrange(step.1.grid, step.2.grid, step.3.grid, newpage = T, nrow = 3, top = p)
  }
  graphics.off()
  
  write.table(x = conc.s, file = "Calculated Concentration.csv", sep = ",", na = "", row.names = FALSE)
  
}

######################################################################################
#filepattern.samples = "xxx FI.csv"
#filepattern.knownValues = "xxx Expected Conc.csv"
#dir <- "/Users/Rufei/Desktop/OMRF Analysis/OLC 200"
#sampleid = "Barcode"
#pattern.s = "Std"

#pre.analysis.qc(filepattern.samples = "xxx FI.csv", filepattern.knownValues = "xxx Expected Conc.csv", dir = "/Users/rufeilu/Desktop/OMRF Analysis/OLC 200", sampleid = "Barcode", pattern.s = "Std", pattern.blank = "Blank.*$", lob = "yes")

