########################################
getPlates <- function (pattern, x){
  ret <- sub(pattern, "\\1", x, ignore.case = TRUE)
  return(ret)
}
########################################

########################################
#x <- expectedConc
SmartOmit <- function(x, threshold = 2){
  rowDel <- c()
  for (i in 1:nrow(x)){
    countNA <- sum(is.na(x[i, ]))
    if (countNA > threshold){
      rowDel <- c(rowDel, i) 
    }
  }
  
  x <- x[-rowDel,]
  return(x)
}
########################################

########################################
WellToColumns <- function(x, testName){
  no.well <- (ncol(x) - 1) * (nrow(x) - 1)
  tempColumns <- data.frame(c(1:no.well), c(1:no.well)) #Two columns with 96 wells
  colnames(tempColumns) <- c("Well", testName)
  for (k in 1:nrow(x)){
    for (l in 2:ncol(x)){
      row <- (k - 1) * 12 + (l - 1)
      tempColumns[row, 1] <- paste0(x[k, 1], l - 1)
      tempColumns[row, 2] <- x[k, l]
    }
  }
  
  return(tempColumns)
}
########################################

########################################
#Reformatter function to translate 96 well plate into pre-stat QC usable files

Reformatter <- function (dir, wellPattern, analytePattern){
  require(XLConnectJars)
  require(XLConnect)
  
  setwd(dir) #Set working directory 
  fileNames <- dir(dir) #Get file names in the working directory
  
  #Automatically generate file finder based on the the 96-well plate format file pattern
  wellFiles_p <- paste0(sub("(.*?)xxx(.*?)$", "\\1", wellPattern),"(.*?)",
         sub("(.*?)xxx(.*?)$", "\\2", wellPattern),"$")
  wellFiles <- grep(pattern = wellFiles_p, x = fileNames, ignore.case = TRUE)
  #Get plate names automatically based on file pattern
  wellPlates <- lapply(X = fileNames[wellFiles], FUN = function(x) getPlates(pattern = wellFiles_p, x))
  wellPlates <- unlist(wellPlates)
  
  #Automatically generate file finder based on the multi-analytes format file pattern
  analyteFiles_p <- paste0(sub("(.*?)xxx(.*?)$", "\\1", analytePattern),"(.*?)",
                      sub("(.*?)xxx(.*?)$", "\\2", analytePattern),"$")
  analyteFiles <- grep(pattern = analyteFiles_p, x = fileNames, ignore.case = TRUE)
  #Get plate names automatically based on file pattern
  analytePlates <- lapply(X = fileNames[analyteFiles], FUN = function(x) getPlates(pattern = analyteFiles_p, x))
  analytePlates <- unlist(analytePlates)
  
  #Automatically find matching plates
  matchedPlates = c() #Initialize the variable 
  for (i in 1:length(wellPlates)){
    count = length(grep(pattern = wellPlates[i], x = analytePlates, ignore.case = TRUE))
    if (count > 0) {
      matchedPlates = c(matchedPlates, wellPlates[i])
    }
  }
  
  for (i in 1:length(matchedPlates)){
    #Read in xlsx file that contains 96-Well Plate format multi-sheets
    tempWell <- loadWorkbook(filename = paste0(matchedPlates[i], sub("(.*?)xxx(.*?)$", "\\2", wellPattern)),
                              create = FALSE)
    analyteNames <- getSheets(tempWell) #Get name of all analytes
    exclude.id1 <- grep(pattern = "Assaychex", analyteNames)
    exclude.id2 <- grep(pattern = "Standard Curve", analyteNames)
    exclude <- c(exclude.id1, exclude.id2)
    analyteNames <- analyteNames[-c(exclude)]
    
    #Parse out the Assaychex and Standard curve tabs
    
    
    #Read in xlsx file that contains multi-analytes files
    tempAnalyte <- loadWorkbook(filename = paste0(matchedPlates[i], sub("(.*?)xxx(.*?)$", "\\2", analytePattern)),
                                create = FALSE)
    
    #Grab the sheet called FI, which contains the template
    template <- na.omit(readWorksheet(object = tempAnalyte, sheet = "FI", startCol = 1, endCol = 2))
    colnames(template) <- template[1, ]
    template <- template[-1, ]
    
    #Grab the expected concentration
    expectedConc <- readWorksheet(object = tempAnalyte, sheet = "Exp Conc")
    expectedConc <- SmartOmit(x = expectedConc)
    
    #Be carefully translating here!!!!! Make sure format is correct, the code is semi-SMART
    #This part of the code assume a standard format from the machine
    colNames <- unlist(c(expectedConc[2, c(1, 2)], expectedConc[1, c(3:ncol(expectedConc))]))
    colnames(expectedConc) <- colNames
    colnames(expectedConc)[1] <- "Barcode" #Rename first column to Barcode to standardize
    expectedConc <- expectedConc[-c(1:2), -2]
    if ("Assaychex" %in% colnames(expectedConc)){
      expectedConc <- expectedConc[, -c(grep(pattern = "Assaychex", colnames(expectedConc)))]
    }
   
    
    #Write csv file into the expected concentration
    ExpConcName <- paste0(matchedPlates[i], " Expected Conc.csv")
    write.csv(x = expectedConc, file = ExpConcName, row.names = FALSE)
    
    #Sort the template out based on the first two columns in the FI tab
    for (j in 1:nrow(template)){
      temp <- unlist(strsplit(x = template[j, 2], split = ","))
      if (length(temp) > 1){
        template[j, 2] <- temp[1]
        for (k in 2:length(temp)){
          tempRow <- c(template[j, 1], temp[k])
          template <- rbind(template, tempRow)  
        }
      }
    }
    template <- template[order(template[, 2]), ] #Sort the template based on teh 96-well plate
    
    tempData <- template
    for (j in 1:length(analyteNames)){
      analyteTemp <- readWorksheet(object = tempWell, sheet = analyteNames[j])
      analyteTemp <- SmartOmit(analyteTemp)
      
      #Sort 96-Well plate into column forms
      tempColumns <- WellToColumns(x = analyteTemp, testName = analyteNames[j])
      
      tempData <- cbind(tempData, c(1:96)) 
      colnames(tempData)[ncol(tempData)] <- analyteNames[j]
      for (o in 1:nrow(tempData)){
        temp <- which(tempColumns[, 1] == tempData[o, 2])
        temp <- tempColumns[temp, 2]
        tempData[o, ncol(tempData)] <- temp
      }
    }
    tempData <- tempData[, -2]
    colnames(tempData)[1] <- "Barcode"
    
    FIname <- paste0(matchedPlates[i], " FI.csv")
    write.csv(x = tempData, file = FIname, row.names = FALSE)
    
  }

}
#########################################
#Example of running the code
#dir <- "/Users/Rufei/Desktop/OMRF Analysis/BOLD cytokine data/M016 BOLD Challenge Flu 50-plex poly plasma/Data/Raw Data Processing"
#wellPattern <- "xxx FI Output 96-well format.xlsx"
#analytePattern <- "xxx opt stds multianalyte labeled.xlsx"
#Reformatter(dir = dir, wellPattern = wellPattern, analyatePattern = analyatePattern)
