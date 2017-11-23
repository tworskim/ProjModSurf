####
# 1- Creating a noised file from an .off file
####
file = "../Data/chat.off"
noisestd = 0.5

#Reading OFF files
data <-readLines(file)
numberOfElements <- strsplit(data[2], " ")
Nvertices <- as.integer(numberOfElements[[1]][1])
Nfaces <- as.integer(numberOfElements[[1]][2])

#Noising data
for (i in 1:Nvertices){
  splitedString <- strsplit(data[2 + i], " ")
  splitedString <- as.double(splitedString[[1]])
  splitedString <- splitedString + rnorm(3, 0, noisestd)
  data[2 + i] <- paste(splitedString[1], splitedString[2], splitedString[3])
}

#Saving noised file
writeLines(data, con = "testnex.off")

