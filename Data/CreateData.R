#This file export "drive.csv", "fitness.csv", "initial.csv" and "mutation.csv" for gene drive simultaions
#Model parameters:
nloci <- 2
locus1 <- c("Aw","Ad","Ar")
locus2 <- c("Bw","Bd","Br")

#Gametes
length(locus1)
for (i in 1:length(locus1)){
    for (j in 1:length(locus2)){
        print(j)
    }
}

0.5*0.94+0.7*0.06