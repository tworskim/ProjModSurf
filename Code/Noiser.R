####
# 1- Creating a noised file from an .off file
####
args = commandArgs(trailingOnly=TRUE)
print(args)
dir = "../Meshs/"
listmesh = setdiff(list.files(path = dir), list.dirs(path = dir, recursive = FALSE, full.names = FALSE))
meshs = strsplit(listmesh, ".off")

#file = args[1]
#noisestd = as.integer(args[2])
for (file in meshs){
for (noisestd in c(0.001, 0.002,0.005,0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)){
filepath = paste(dir, file,".off", sep = "")
#dir.create(paste(dir, file, sep = ""))
#dir.create(paste(dir, file, "/", "input",  sep = ""))


#Reading OFF files
data <-readLines(filepath)
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
  data[2 + i + Nvertices] <- paste(splitedString[1:length(splitedString)], collapse = " ")
}

#Saving noised file
dirwl = paste(dir,file, "/input/", sep = "")
filwl = paste(noisestd, ".off", sep = "")
writeLines(data, con = paste(dirwl, filwl, sep = ""))
print(file)
print(noisestd)
}}
list.dirs("../Inputs/")