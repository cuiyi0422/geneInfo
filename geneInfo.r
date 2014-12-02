geneIdentify = function(DBfile, id){
  
  ## Building the gene information database
  gdbBuilder = function(DBfile){
    
    # reading unformed database file
    geneDB <- readLines(DBfile)
    geneDB <- as.matrix(geneDB)
    
    geneRowName <- geneDB[1]
    geneRowName = unlist(strsplit(geneRowName,'\t'))
    database = matrix(ncol = length(geneRowName),nrow = length(geneDB)-1,data = "NA")
    colnames(database) = geneRowName
    rm(geneRowName)
    
    for (i in 2:length(geneDB)){
      for (j in 1:15){
        ids=geneDB[i]
        ids=unlist(strsplit(ids,'\t'))
        database[i-1,j]=ids[j]
      }
    }
    rm(i,j,ids)
    return(database)
    # write.table(database, 'formedGeneDB.txt', quote = FALSE, row.names = FALSE, sep='\t')
  }
  
  geneInfoDB <- gdbBuilder(DBfile)
  
  geneInfoDB=as.data.frame(geneInfoDB)
  geneInfoDB$tax_id=NULL
  geneInfoDB$chromosome=NULL
  geneInfoDB=as.matrix(geneInfoDB)

  id = 'TPP1'

  index <- which(geneInfoDB == id)
  colInd <- NULL
  rowInd <- NULL

  for (k in 1:length(index)){
    i <- 0
    while (i <= dim(geneInfoDB)[2]){
      i <- i+1
      if (i*dim(geneInfoDB)[1] >= index[k] && (i-1)*dim(geneInfoDB)[1] < index[k]){
        j = index[k] - (i-1)*dim(geneInfoDB)[1]
        colInd = c(colInd, j)
        rowInd = c(rowInd, i)
        }
      }
    }
  rm(k, i, j)

  geneInfoDB[colInd[1], rowInd[1]]
}
