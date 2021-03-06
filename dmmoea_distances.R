load.gene.distance <- function(path, alpha=0.5){
    
    # Load data matrix
    data.matrix<-as.matrix(read.table(file.path(path,"BD.csv"), header=T, sep=",", row.names=1, as.is = TRUE))
    n.genes <-nrow(data.matrix)
    
    # Load biology-based distance matrix
    bio.dist <-read.table(paste(path,"/matriz_distancia_biologia_wang.csv",sep=""), header=T, sep=",")
    
    # Load expression-based distance matrix
    exp.dist <-read.table(paste(path,"/matriz_distancia_expresion_pearson.csv",sep=""), header=T, sep=",")
    
    comp.dist <- alpha*exp.dist + (1-alpha)*exp.dist
    #print(comp.dist[1:20, 1:20])
    
    return(list(data.matrix=data.matrix, exp.dist=exp.dist, bio.dist=bio.dist, comp.dist=comp.dist, n.genes=n.genes))  
}