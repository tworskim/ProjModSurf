####
# 1- Creating a noised file from an .off file
####
args = commandArgs(trailingOnly=TRUE)
print(args)
#file = args[1]
#noisestd = as.integer(args[2])
file = "../Data//eight.off"
noisestd = 0.1

#Reading OFF files
data <-readLines(file)
numberOfElements <- strsplit(data[2], "[ |\t]+")
Nvertices <- as.integer(numberOfElements[[1]][1])
Nfaces <- as.integer(numberOfElements[[1]][2])

#Noising data
for (i in 1:Nvertices){
  splitedString <- strsplit(data[2 + i], "[ |\t]+")
  splitedString <- as.double(splitedString[[1]])
  splitedString <- splitedString + rnorm(3, 0, noisestd)
  data[2 + i] <- paste(splitedString[1], splitedString[2], splitedString[3])
}

for (i in 1:Nfaces){
  splitedString <- strsplit(data[2 + i + Nvertices], "[ |\t]+")
  splitedString <- as.integer(splitedString[[1]])
  data[2 + i + Nvertices] <- paste(splitedString[1], splitedString[2], splitedString[3], splitedString[4])
}

#Saving noised file
writeLines(data, con = "input.off")