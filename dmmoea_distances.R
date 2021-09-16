load.gene.distance <- function(path){
    
    # Load data matrix
    data.matrix<-as.matrix(read.table(paste(path,"/BD.csv",sep=""), header=T, sep=","))
    n.genes <-nrow(data.matrix)
    
    # Load biology-based distance matrix
    bio.dist <-read.table(paste(path,"/matriz_distancia_biologia_wang.csv",sep=""), header=T, sep=",")
    
    # Load expression-based distance matrix
    exp.dist <-read.table(paste(path,"/matriz_distancia_expresion_pearson.csv",sep=""), header=T, sep=",")
    
    return(list(data.matrix=data.matrix, exp.dist=exp.dist, bio.dist=bio.dist, n.genes=n.genes))  
}