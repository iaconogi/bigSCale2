
#' homogenize.networks
#' @param input.networks A list with in which each element is the output of a previous compute.network()
#' @param tick Pearson correlation steps sued to homogenize the numner of edges. The smallest the larger the computational time but the better the homogenization.
#' @examples
#' homogenize.networks(list(output1,output2))
#' @export
homogenize.networks<-function(input.networks){
    
    input.networks=homogenize.networks.internal(input.networks,tick = 0.025)
    input.networks=homogenize.networks.internal(input.networks,tick = 0.01)
    input.networks=homogenize.networks.internal(input.networks,tick = 0.0025)
    input.networks=homogenize.networks.internal(input.networks,tick = 0.001)
    return(input.networks)
  }




homogenize.networks.internal<-function(input.networks,tick){
  
  edges=c()
  cutoff=c()
  for (k in 1:length(input.networks))
  {
    edges[k]=igraph::gsize(input.networks[[k]]$graph)
    cutoff[k]=input.networks[[k]]$cutoff.p
  }
  
  print('Number of edges for each network')
  print(edges)
  median.edges=median(edges)+1
  print(sprintf('Aiming at %g edges for each network',median.edges))
  
  if (sum(edges==median.edges)>0)
    error('One network has exactly the amount of edges we are aiming to. This will cause an infinite loop. In the very rare case plese contact developer at gio.iacono.work@gmail.com')
  
  for (k in 1:length(input.networks))
  {
    dynamic.tick=tick
    
    gap=edges[k]-median.edges
    tick.sign=sign(gap)
    dynamic.tick=tick.sign*tick
    
    while(1)
    {
      dummy=compute.network(previous.output = input.networks[[k]],quantile.p = (cutoff[k]+dynamic.tick))
      new.gap=igraph::gsize(dummy$graph)-median.edges
      if (abs(new.gap)>abs(gap))
        break
      else
      {
        input.networks[[k]]=dummy
        gap=new.gap
        dynamic.tick=dynamic.tick+tick.sign*tick
      }
    }
  }
  
  return(input.networks)
  
}



batch.normalize<-function(mat,sizeF){
  
  blockS=min(5000,ncol(mat))
  avg.library.size=mean(sizeF)
  A=1
  B=A+blockS-1
  last.round=0
  
  pb <- progress::progress_bar$new(format = "Normalizing cells [:bar] :current/:total (:percent) eta: :eta", total = ncol(mat))
  
  for (k in 1:ceiling(ncol(mat)/blockS))
  {
    dummy=as.matrix(mat[,A:B])*avg.library.size
    for (h in 1:ncol(dummy))
      dummy[,h]=dummy[,h]/sizeF[h]
    
    #start_time <- Sys.time()
    if (k==1)
      output.mat=Matrix::Matrix(dummy)
    else
      output.mat=cbind(output.mat,Matrix::Matrix(dummy))
    #end_time <- Sys.time()
    #print(end_time - start_time)
    A=A+blockS
    B=B+blockS
    pb$tick(B-A+1)
    if (A > ncol(mat)) break
    if (B > ncol(mat)) B=ncol(mat)
  }
  return(output.mat)
}



#' bigscale.toCytoscape
#' @param G The graph in igraph format previoulsy calculated with compute.network
#' @param file.name Directory and/or name of the Cytoscape file 
#' @examples
#' bigscale.toCytoscape(G,'graph.json')
#' @export

toCytoscape <- function (G, file.name) {
  # Extract graph attributes
  library(igraph)
  graph_attr = graph.attributes(G)
  
  # Extract nodes
  node_count = length(V(G))
  if('name' %in% list.vertex.attributes(G)) {
    V(G)$id <- V(G)$name
  } else {
    V(G)$id <- as.character(c(1:node_count))
  }
  
  nodes <- V(G)
  v_attr = vertex.attributes(G)
  v_names = list.vertex.attributes(G)
  
  nds <- array(0, dim=c(node_count))
  for(i in 1:node_count) {
    if(i %% 1000 == 0) {
      print(i)
    }
    nds[[i]] = list(data = mapAttributes(v_names, v_attr, i))
  }
  
  edges <- get.edgelist(G)
  edge_count = ecount(G)
  e_attr <- edge.attributes(G)
  e_names = list.edge.attributes(G)
  
  attr_exists = FALSE
  e_names_len = 0
  if(identical(e_names, character(0)) == FALSE) {
    attr_exists = TRUE
    e_names_len = length(e_names)
  }
  e_names_len <- length(e_names)
  
  eds <- array(0, dim=c(edge_count))
  for(i in 1:edge_count) {
    st = list(source=toString(edges[i,1]), target=toString(edges[i,2]))
    
    # Extract attributes
    if(attr_exists) {
      eds[[i]] = list(data=c(st, mapAttributes(e_names, e_attr, i)))
    } else {
      eds[[i]] = list(data=st)
    }
    
    if(i %% 1000 == 0) {
      print(i)
    }
  }
  
  el = list(nodes=nds, edges=eds)
  
  x <- list(data = graph_attr, elements = el)
  print("Done.  Writing Json...")
  #return (toJSON(x))
  
  #merda<<-jsonlite::toJSON(x)
  fileConn<-file(file.name)
  writeLines(RJSONIO::toJSON(x), fileConn)
  close(fileConn)
  
  
}


mapAttributes <- function(attr.names, all.attr, i) {
  attr = list()
  cur.attr.names = attr.names
  attr.names.length = length(attr.names)
  
  for(j in 1:attr.names.length) {
    if(is.na(all.attr[[j]][i]) == FALSE) {
      #       attr[j] = all.attr[[j]][i]
      attr <- c(attr, all.attr[[j]][i])
    } else {
      cur.attr.names <- cur.attr.names[cur.attr.names != attr.names[j]]
    }
  }
  names(attr) = cur.attr.names
  return (attr)
}





#' bigscale.generate.report
#' @param sce The single cell object for which you want to generate the report
#' @param file.name Directory and/or name of the excel file containing the formatted report
#' @examples
#' bigscale.generate.report('report.xlsx')
#' @export


bigscale.generate.report = function(sce,file.name)
{
  export.data=as.data.frame(as.numeric(table(getClusters(sce))))
  
  colnames(export.data)='Cell number'
  
  avg.lib.size=c()
  dummy=sizeFactors(sce)
  for (k in 1:max(getClusters(sce)))
    avg.lib.size[k]=mean(dummy[which(getClusters(sce)==k)])
  
  avg.det.genes=c()
  for (k in 1:max(getClusters(sce)))
    avg.det.genes[k]=mean(Rfast::colsums(as.matrix(normcounts(sce)[,which(getClusters(sce)==k)]>0)))
  
  export.data$Average_Library_Size=avg.lib.size
  export.data$Average_Detected_Genes=avg.det.genes
  export.data$Complexity=avg.lib.size/avg.det.genes
  
  export.data$Markers_Lv1=sce@int_metadata$Mlist.counts$Markers_LV1
  export.data$Markers_Lv2=sce@int_metadata$Mlist.counts$Markers_LV2
  
  #R.utils::copyFile(paste(system.file(package="bigSCale"),'/data/template.xlsx',sep = ''),file.name)
  
  
  #Level 1 markers
  Mlist=getMarkers(sce)
  to.write=Mlist[[1,1]][1,]
  for (k in 2:max(getClusters(sce)))
    to.write=rbind(to.write,Mlist[[k,1]][1,])
 export.data=cbind(export.data,to.write[,2:4])
  
  
  #Level 1 markers
  to.write=Mlist[[1,2]][1,]
  for (k in 2:max(getClusters(sce)))
    to.write=rbind(to.write,Mlist[[k,2]][1,])
  export.data=cbind(export.data,to.write[,2])
  
  names=colnames(export.data)
  names[7]='Best Lv1'
  names[10]='Best Lv2'
  
  colnames(export.data)=names
  xlsx::write.xlsx(export.data,file.name)
}



pharse.10x.peaks = function(file,keep='promoter',reject='distal')
{
  peaks = read.delim(file, header=FALSE, stringsAsFactors=FALSE)
  peaks=peaks[-1,]
  good.ix=rep(0,nrow(peaks))
  feature.names=c()
  for (k in 1:nrow(peaks))
  {
    if (peaks[k,4]!=reject)
    {
      types=unlist(strsplit(peaks[k,4], ";"))
      genes=unlist(strsplit(peaks[k,2], ";"))
      for ( h in 1:length(types))
      {
        if (types[h]==keep)
        {
          good.ix[k]=1
          if (h==1) feature.names[k]=genes[h]
          else { feature.names[k]=paste(feature.names[k],';', genes[h]) }
        }
      }
    }
  }
  return(list(feature.names=feature.names[which(good.ix==1)],good.ix=which(good.ix==1)))
}



pick.sequences = function(file, numbers=NA)
{
  library("BSgenome.Hsapiens.UCSC.hg19")
  peaks = read.delim(file, header=FALSE, stringsAsFactors=FALSE) 
  if (is.na(numbers[1]))
    numbers=c(1:nrow(peaks))
  sequences=getSeq(Hsapiens, as.character(peaks[numbers,1]), start=peaks[numbers,2],end=peaks[numbers,3])
}

sub.communities <- function(G,dim,module)
{
  mycl.new=rep(0,length(igraph::V(G)))
  unclusterable=mycl.new  
  
  #node.names=V(G)$names
  #node.ix=c(1:length(igraph::V(G)))
  
  mycl=igraph::cluster_louvain(G,weights = NULL)$membership
  cat(sprintf('\n Starting from %g clusters',max(mycl)))
  #current.nodes=node.ix
  round=1
  while (1)
  {
    action.taken=0
    for (k in 1:max(mycl))
    {
      ix=which(mycl==k)  
      if (length(ix)>dim & sum(unclusterable[ix])==0 )
      {
        G.sub=igraph::induced.subgraph(graph = G,vids = which(mycl == k))
        out=igraph::cluster_louvain(G.sub,weights = NULL)$membership
        mycl.new[ix]=out+max(mycl.new)
        action.taken=1
      }
      else
      {
        mycl.new[ix]=1+max(mycl.new)
        unclusterable[ix]=1 
      }
    }
    round=round+1
    cat(sprintf('\n Round %g) Obtained %g clusters',round,max(mycl.new)))  
    if (action.taken==0) 
      break
    else
    {
      mycl=mycl.new
      mycl.new=rep(0,length(igraph::V(G)))
    }
  }
  
  
  return(mycl)
  
}



Pdistance <- function(counts) {

  counts=counts>0
  pop.size=nrow(counts)
  pvalues=matrix(NA,ncol(counts),ncol(counts))
  samp.size=Rfast::colsums(counts)
  
  for (h in 1:ncol(counts))
  {
    samp.hits=Rfast::colsums(  counts[which(counts[,h]>0),]  )
    pop.hits=rep(samp.size[h],ncol(counts))
    pvalues[h,]=phyper(samp.hits, pop.hits, pop.size-pop.hits, samp.size,lower.tail = FALSE)+dhyper(samp.hits, pop.hits, pop.size-pop.hits, samp.size)
    #pvalues[h,]=(samp.hits/samp.size)/(pop.hits/pop.size)
  }
  pvalues=1/(-log10(pvalues))
  diag(pvalues)=0
  pvalues.out<<-pvalues
  samp.size.out<<-samp.size
  return(pvalues)
}


jaccard_dist_text2vec_04 <- function(x, y = NULL, format = 'dgCMatrix') {
  if (!inherits(x, 'sparseMatrix'))
    stop("at the moment jaccard distance defined only for sparse matrices")
  # union x
  rs_x = Matrix::rowSums(x)
  if (is.null(y)) {
    # intersect x
    RESULT = Matrix::tcrossprod(x)
    rs_y = rs_x
  } else {
    if (!inherits(y, 'sparseMatrix'))
      stop("at the moment jaccard distance defined only for sparse matrices")
    # intersect x y
    RESULT = Matrix::tcrossprod(x, y)
    # union y
    rs_y = Matrix::rowSums(y)
  }
  RESULT = as(RESULT, 'dgTMatrix')
  # add 1 to indices because of zero-based indices in sparse matrices
  # 1 - (...) because we calculate distance, not similarity
  RESULT@x <- 1 - RESULT@x / (rs_x[RESULT@i + 1L] + rs_y[RESULT@j + 1L] - RESULT@x)
  if (!inherits(RESULT, format))
    RESULT = as(RESULT, format)
  RESULT
}

bigscale.classifier = function(expr.counts,gene.names,selected.genes,stop.at){
  
  
if (length(selected.genes)==1 | length(selected.genes)>3)
  stop('Accepts only two or three markers. Contact the developer if you more.')
  
pvalues=all.coexp(expr.counts = expr.counts,gene.names = gene.names,selected.genes = selected.genes)
    
if (is.na(stop.at))
  tot.markers=20
else
  tot.markers=stop.at


if (length(selected.genes)==2){

ix=order(pvalues[1,])
markers1=gene.names[setdiff(ix,which(pvalues[2,]<0.999999))]
ix=order(pvalues[2,])
markers2=gene.names[setdiff(ix,which(pvalues[1,]<0.999999))]

final.counts=matrix(0,4,tot.markers)
for (k in 1:tot.markers)
{
  
  testing1=which(is.element(gene.names,markers1[1:k]))
  testing2=which(is.element(gene.names,markers2[1:k]))
  
  if (k==1)
  {
    counts1=as.matrix(expr.counts[testing1,]>0)
    counts2=as.matrix(expr.counts[testing2,]>0)
  }    
  else
  {
    counts1=Rfast::colsums(as.matrix(expr.counts[testing1,]>0))
    counts2=Rfast::colsums(as.matrix(expr.counts[testing2,]>0)) 
  }
  
  final.counts[,k]=c(sum(counts1>0 & counts2==0),sum(counts2>0 & counts1==0),sum(counts1==0 & counts2==0),sum(counts1>0 & counts2>0))
  
}

plotting.text=c('none')
for (k in 1:tot.markers)
  plotting.text=c(plotting.text,sprintf('%s / %s',markers1[k],markers2[k]))
  
clusters=rep(0,length(counts1))
clusters[which(counts1>0 & counts2==0)]=1
clusters[which(counts2>0 & counts1==0)]=2
clusters[which(counts1==0 & counts2==0)]=3
clusters[which(counts1>0 & counts2>0)]=4


final.counts=cbind(c(0,0,ncol(expr.counts),0),final.counts)/ncol(expr.counts)
rownames(final.counts)=c(markers1[1],markers2[1],'others','doublets')
final.counts=t(final.counts)
final.counts=as.data.frame(final.counts)


p <- plotly::plot_ly(data = final.counts, x = c(0:tot.markers), y = final.counts[,1], name = selected.genes[1], type='scatter', mode = 'lines+markers',text=plotting.text)
p=plotly::add_trace(p,y = final.counts[,2], name = selected.genes[2], mode = 'lines+markers',text=plotting.text)
p=plotly::add_trace(p,y = final.counts[,3], name = 'others', mode = 'lines+markers',text=plotting.text)
p=plotly::add_trace(p,y = final.counts[,4], name = 'doublets', mode = 'lines+markers',text=plotting.text)
p=plotly::layout(p,xaxis = list(title='Number of genes used'), yaxis = list(title='Percentage of cells'),title=sprintf('%g cells %s+, %g cells %s+, %g others and %g Doublets',sum(clusters==1),selected.genes[1],sum(clusters==2),selected.genes[2],sum(clusters==3),sum(clusters==4)))
print(p)


return(clusters)
}


if (length(selected.genes)==3){
  
  ix=order(pvalues[1,])
  markers1=gene.names[setdiff(ix,which(pvalues[2,]<0.999999 | pvalues[3,]<0.999999))]
  
  ix=order(pvalues[2,])
  markers2=gene.names[setdiff(ix,which(pvalues[1,]<0.999999 | pvalues[3,]<0.999999))]
  
  ix=order(pvalues[3,])
  markers3=gene.names[setdiff(ix,which(pvalues[1,]<0.999999 | pvalues[2,]<0.999999))]
  
  final.counts=matrix(0,5,tot.markers)
  for (k in 1:tot.markers)
  {
    
    testing1=which(is.element(gene.names,markers1[1:k]))
    testing2=which(is.element(gene.names,markers2[1:k]))
    testing3=which(is.element(gene.names,markers3[1:k]))
    
    if (k==1)
    {
      counts1=as.matrix(expr.counts[testing1,]>0)
      counts2=as.matrix(expr.counts[testing2,]>0)
      counts3=as.matrix(expr.counts[testing3,]>0)
    }    
    else
    {
      counts1=Rfast::colsums(as.matrix(expr.counts[testing1,]>0))
      counts2=Rfast::colsums(as.matrix(expr.counts[testing2,]>0)) 
      counts3=Rfast::colsums(as.matrix(expr.counts[testing3,]>0)) 
    }
    
    final.counts[1:4,k]=c(sum(counts1>0 & counts2==0 & counts3==0),sum(counts2>0 & counts1==0 & counts3==0),sum(counts3>0 & counts1==0 & counts2==0),sum(counts1==0 & counts2==0 & counts3==0))
    final.counts[5,k]=ncol(expr.counts)-sum(final.counts[1:4,k])
  }
  
  plotting.text=c('none')
  for (k in 1:tot.markers)
    plotting.text=c(plotting.text,sprintf('%s / %s / %s',markers1[k],markers2[k],markers3[k]))
  
  clusters=rep(0,length(counts1))
  clusters[which(counts1>0 & counts2==0 & counts3==0)]=1
  clusters[which(counts2>0 & counts1==0 & counts3==0)]=2
  clusters[which(counts3>0 & counts1==0 & counts2==0)]=3
  clusters[which(counts1==0 & counts2==0 & counts3==0)]=4
  clusters[which(clusters==0)]=5
  
  
  final.counts=cbind(c(0,0,0,ncol(expr.counts),0),final.counts)/ncol(expr.counts)
  rownames(final.counts)=c(markers1[1],markers2[1],markers3[1],'others','doublets')
  final.counts=t(final.counts)
  final.counts=as.data.frame(final.counts)
  
  
  p <- plotly::plot_ly(data = final.counts, x = c(0:tot.markers), y = final.counts[,1], name = selected.genes[1], type='scatter', mode = 'lines+markers',text=plotting.text)
  p=plotly::add_trace(p,y = final.counts[,2], name = selected.genes[2], mode = 'lines+markers',text=plotting.text)
  p=plotly::add_trace(p,y = final.counts[,3], name = selected.genes[3], mode = 'lines+markers',text=plotting.text)
  p=plotly::add_trace(p,y = final.counts[,4], name = 'others', mode = 'lines+markers',text=plotting.text)
  p=plotly::add_trace(p,y = final.counts[,5], name = 'doublets', mode = 'lines+markers',text=plotting.text)
  p=plotly::layout(p,xaxis = list(title='Number of genes used'), yaxis = list(title='Percentage of cells'),title=sprintf('%g cells %s+, %g cells %s+, %g cells %s+, %g others and %g Doublets',sum(clusters==1),selected.genes[1],sum(clusters==2),selected.genes[2],sum(clusters==3),selected.genes[3],sum(clusters==4),sum(clusters==5)))
  print(p)
  
  
  return(clusters)
}

}




all.coexp = function (expr.counts,gene.names,selected.genes){
  
  expr.counts=as.matrix(expr.counts)
  pop.size=ncol(expr.counts)
  pvalues=matrix(NA,length(selected.genes),length(gene.names))
  #fe=pvalues
  #counts=pvalues


  for (h in 1:length(selected.genes))
  {
    # v1=expr.counts[which(gene.names==selected.genes[h]),]
    # pop.hits=sum(v1>0)
    # print(sprintf('Computing p-values for gene %s, it should take few minutes ...',selected.genes[h]))
    # for (k in 1:length(gene.names))  
    # {
    #   samp.size=sum(expr.counts[k,]>0)
    #   samp.hits=sum(v1>0 & expr.counts[k,]>0)
    #   pvalues[h,k]=phyper(samp.hits, pop.hits, pop.size-pop.hits, samp.size,lower.tail = FALSE)+dhyper(samp.hits, pop.hits, pop.size-pop.hits, samp.size)
    #   #fe[h,k]=(samp.hits/samp.size)/(pop.hits/pop.size)
    #   #counts[h,k]=samp.hits
    # }
    print(sprintf('Computing p-values for gene %s',selected.genes[h]))
    
    v1=expr.counts[which(gene.names==selected.genes[h]),]
    
    samp.size=Rfast::rowsums(expr.counts>0)
    samp.hits=Rfast::rowsums(expr.counts[,which(v1>0)]>0)
    pop.hits=rep(sum(v1>0),nrow(expr.counts))
    pvalues[h,]=phyper(samp.hits, pop.hits, pop.size-pop.hits, samp.size,lower.tail = FALSE)+dhyper(samp.hits, pop.hits, pop.size-pop.hits, samp.size)
  }
  return(pvalues)
}



#' deconvolute
#'
#' @export
#' 
#' 
deconvolute = function(icells.end,icells.start){
  
  
  icells=matrix(NA,nrow(icells.end),ncol(icells.end)*ncol(icells.start))
  
  for (k in 1:nrow(icells.end))
  {
    dummy=icells.end[k,]
    dummy=dummy[!is.na(dummy)]
    icells.temp=c()
    for (j in 1:length(dummy))
      icells.temp=c(icells.temp,icells.start[dummy[j],])
    icells.temp=icells.temp[!is.na(icells.temp)]
    icells[k,1:(length(icells.temp))]=icells.temp
  }

  return(icells)
  
}

#' extract.cells
#'
#' Works with the .mex files. It subsets to the specified set of cells and write them into the output .mex file.
#' 
#' @param input.mex File name of the input dataset
#' @param input.mex File name of the output dataset
#' @param cells Indices of the cells you want to extract
#' @param chunks A techincal parameter, you should not really touch this. It represent how large is the chunk of cells read when sequentially pharsing the input.mex
#' 
#' @return  Writes to disk the subsetted dataset. Please note that if \code{cells} will be automatially sorted. For example, is \code{cells=c(10,1)}, meaning that you want to extract cells 10 and 1, the oupput matrix will contain, in order, cells 1 and 10.
#'
#'
#' @examples
#' check the online tutorial at Github
#'
#' @export

extract.cells = function(input.mex,output.mex,cells,chuncks=50000){

  
  # Reading header for sample info
  out=bigscale.readMM(file=input.mex,size.only=TRUE)
  tot.cells=out[1]
  tot.genes=out[2]
  tot.numbers=out[3]
  avg.genes.cell=tot.numbers/tot.cells
  chuncks.nz=chuncks*avg.genes.cell
  
  # Counting nz
  pointer=0
  tot.cells.read=0
  nz=0
  while (TRUE)
    {
    
    out=bigscale.readMM(file=file.dir,max.vals = chuncks.nz,skip.vals = pointer)
    pointer=out$pointer
    
    if (length(out$dgTMatrix)==0) break

    cell.references=c(1:ncol(out$dgTMatrix))+tot.cells.read
    tot.cells.read=tot.cells.read+ncol(out$dgTMatrix)
    selected.cells=which(is.element(cell.references,cells))
    print(sprintf('Read %g cells, selected %g cells',tot.cells.read,length(selected.cells)))
    
    if (length(selected.cells)>0)
         nz=nz+sum(out$dgTMatrix[,selected.cells]>0) 

    }
      
      
  # Writing the beginning of the output file
  temp.file=gzfile(output.mex)
  open(temp.file,"w")
  writeLines("%%MatrixMarket matrix coordinate integer general", temp.file)
  write.table(c(tot.genes,length(cells),nz),temp.file,row.names = FALSE,col.names = FALSE)    
  
  
  
  # Reading and writing
  pointer=0
  tot.cells.read=0
  tot.cells.selected=0
  while (TRUE)
    {
    out=bigscale.readMM(file=file.dir,max.vals = chuncks.nz,skip.vals = pointer)
    pointer=out$pointer
    if (length(out$dgTMatrix)==0) break    
 
    cell.references=c(1:ncol(out$dgTMatrix))+tot.cells.read
    tot.cells.read=tot.cells.read+ncol(out$dgTMatrix)
    selected.cells=which(is.element(cell.references,cells))
    
    print(sprintf('Read %g cells, of which selected %g',tot.cells.read,length(selected.cells)))
    
    if (length(selected.cells)>0)
      {
      temp.mat=Matrix::summary(out$dgTMatrix[,selected.cells])
      temp.mat[,2]=temp.mat[,2]+as.integer(tot.cells.selected)
      tot.cells.selected=tot.cells.selected+length(selected.cells)
      write.table(temp.mat,temp.file,row.names = FALSE,col.names = FALSE)
      }
    
  }
  
  close(temp.file)
}


merge.chunks = function(verbose){
  
  file.names=list.files(pattern = "Temp_")
  if (verbose==TRUE) print(sprintf('Detected %g chunks to merge',length(file.names)))
  
  nc=0
  nz=0
  for (k in 1:length(file.names))
  {
    numbers=bigscale.readMM(file.names[k],size.only = TRUE)
    nc=nc+numbers[1]
    nz=nz+numbers[3]
    nr=numbers[2] # this actually enough to be done only once
  }
  
  
  temp.file=gzfile("iCells.mtx.gz")
  open(temp.file,"w")
  
  writeLines("%%MatrixMarket matrix coordinate integer general", temp.file)
  write.table(c(nr,nc,nz),temp.file,row.names = FALSE,col.names = FALSE)
  
  current.cell.num=0
  for (k in 1:length(file.names))
  {
    if (verbose==TRUE) print(sprintf('Appending chunk %g/%g',k,length(file.names)))
    temp.mat=Matrix::readMM(file.names[k])
    temp.mat=Matrix::summary(temp.mat)
    
    if (verbose==TRUE) print(sprintf("Input cells from %g to %g",min(temp.mat[,2]),max(temp.mat[,2])))
    current.cells=max(temp.mat[,2])
    temp.mat[,2]=temp.mat[,2]+as.integer(current.cell.num)
    current.cell.num=current.cell.num+current.cells
    if (verbose==TRUE) print(sprintf("Writing cells from %g to %g",min(temp.mat[,2]),max(temp.mat[,2])))
    
    write.table(temp.mat,temp.file,row.names = FALSE,col.names = FALSE)
  }
  
  
  
  file.remove(file.names)
  close(temp.file)
  
}


#' bigscale.writeMM
#'
#' @export
#' 
bigscale.writeMM <- function(data,file.name)
{
  nc=ncol(data)
  nr=nrow(data)
  data=Matrix::summary(data)
  nz=nrow(data)
 
  temp.file=gzfile(file.name)
  open(temp.file,"w")
  
  writeLines("%%MatrixMarket matrix coordinate integer general", temp.file)
  write.table(c(nr,nc,nz),temp.file,row.names = FALSE,col.names = FALSE)
  # data[,1]=as.integer(data[,1])
  # data[,2]=as.integer(data[,2])
  #data[,3]=as.integer(data[,3])
  write.table(data,temp.file,row.names = FALSE,col.names = FALSE)
  close(temp.file)
}


#' bigscale.convert.h5
#'
#' Converts and .HDF5 format to a .MEX format
#' @export
#' 
bigscale.convert.h5 <- function(input.file,output.file,counts.field,filter.cells=0)
{
library(Matrix)
library(rhdf5)
  
dims = h5read(input.file, paste(counts.field, "/shape",sep=""))
chunks=50000
current.cell=1

indices = h5read(input.file,paste(counts.field, "/indices",sep=""))
if (length(indices)>2000000000)
  {
  print('We have a very large dataset, storing inptr as double instead of integer')
  indptr = h5read(input.file,paste(counts.field, "/indptr",sep=""),bit64conversion='double')
  }
else
  indptr = h5read(input.file,paste(counts.field, "/indptr",sep=""))

data = h5read(input.file,paste(counts.field, "/data",sep=""))
det.genes=diff(indptr)



# Writing the beginning of the output file
temp.file=gzfile(output.file)
open(temp.file,"w")
writeLines("%%MatrixMarket matrix coordinate integer general", temp.file)
write.table(c(dims[1],sum(det.genes>=filter.cells),length(indices)),temp.file,row.names = FALSE,col.names = FALSE)  
print(sprintf("Writing nc=%g",sum(det.genes>=filter.cells)))
pb <- progress::progress_bar$new(format = "Converting cells [:bar] :current/:total (:percent) eta: :eta", dims[2])
pb$tick(0)

current.cell=1
current.cell.written=1

while (TRUE)
  {
  
  if ((current.cell+chunks)>length(indptr))
    chunks=length(indptr)-current.cell  
  ix=c( (indptr[current.cell]+1) : indptr[current.cell+chunks] )
  
  i_local=as.integer(indices[ix])
  x_local=as.numeric(data[ix])
  
  j_local=(rep(0,length(i_local)))
  indptr_local=indptr[current.cell:(current.cell+chunks)]
  indptr_local=as.integer(indptr_local-indptr_local[1])
  for (k in 1:chunks)
    j_local[(indptr_local[k]+1):indptr_local[k+1]]=k
  j_local=as.integer(j_local)
  
  det.genes.local=det.genes[current.cell:(current.cell+chunks-1)]
  cells.okay=which(det.genes.local>=filter.cells)
  
  temp.matrix=new("dgTMatrix", Dim = as.integer(c(dims[1], chunks)), i = i_local,j = j_local - 1L, x = x_local)
  temp.matrix=temp.matrix[,cells.okay]
  temp.matrix=Matrix::summary(temp.matrix)
  
  
  print(sprintf('Kept %g/%g cells with at least %g espressed genes,current cell written=%g',length(cells.okay),chunks,filter.cells,current.cell.written))
  temp.matrix[,2]=temp.matrix[,2]+as.integer(current.cell.written-1)

  current.cell=current.cell+chunks
  current.cell.written=current.cell.written+length(cells.okay)
  write.table(temp.matrix,temp.file,row.names = FALSE,col.names = FALSE)
  pb$tick(chunks)
  
  if (current.cell==length(indptr)) break
  }
close(temp.file)
h5closeAll()
return(list(filtered.cells=which(det.genes>=filter.cells),det.genes=det.genes))
}


  
  # for ( k in cells.to.read)
  #   {
  #   print(k)
  #   ix=c((indptr[k]+1):(indptr[k+1]))
  #   temp.matrix[indices[ix],k]=data[ix]
  #   }


#' bigscale.readMM
#'
#' @export
#' 
#' 
bigscale.readMM <- function(file,max.vals,skip.vals=0,size.only=FALSE,current.condition=NA)
{
  library(Matrix)
  #print('Starting bigscale.readMM')
  #gc.out<-gc()
  #print(gc.out)
  if (is.character(file))
    file <- if(file == "") stdin() else file(file)
  if (!inherits(file, "connection"))
    stop("'file' must be a character string or connection")
  if (!isOpen(file)) {
    open(file)
    on.exit(close(file))
  }
  scan1 <- function(what, ...)
    scan(file, nmax = 1, what = what, quiet = TRUE, ...)
  
  if (scan1(character()) != "%%MatrixMarket")# hdr
    stop("file is not a MatrixMarket file")
  
  if (!(typ <- tolower(scan1(character()))) %in% "matrix")
    stop(gettextf("type '%s' not recognized", typ), domain = NA)
  
  if (!(repr <- tolower(scan1(character()))) %in% c("coordinate", "array"))
    stop(gettextf("representation '%s' not recognized", repr), domain = NA)
  
  elt <- tolower(scan1(character()))
  if (!elt %in% c("real", "complex", "integer", "pattern"))
    stop(gettextf("element type '%s' not recognized", elt), domain = NA)
  
  sym <- tolower(scan1(character()))
  if (!sym %in% c("general", "symmetric", "skew-symmetric", "hermitian"))
    stop(gettextf("symmetry form '%s' not recognized", sym), domain = NA)
  
  nr <- scan1(integer(), comment.char = "%")
  nc <- scan1(integer())
  nz <- scan1(double())
  
  if (size.only) return(c(nc,nr,nz))
  
  checkIJ <- function(els) 
  {
    if(any(els$i < 1 | els$i > nr))
      stop("readMM(): row	 values 'i' are not in 1:nr", call.=FALSE)
    if (any(els$j < 1 | els$j > nc))
      stop(sprintf("readMM(): column values 'j' are in %g:%g and not in 1:nc (1:%g)",min(els$j),max(els$j),nc), call.=FALSE)
  }
  
  
  
  if (repr == "coordinate") {
    switch(elt,
           "real" = ,
           "integer" = {
             ## TODO: the "integer" element type should be returned as
             ##       an object of an "iMatrix" subclass--once there are
             #print(sprintf('max.vals = %g, nz=%g',max.vals,nz))
              if (max.vals>(nz-skip.vals))
              {
              #print('Reading all file at once') 
              els <- scan(file, skip = skip.vals, quiet = TRUE,what= list(i= integer(), j= integer(), x= numeric()))
              }
              else
              els <- scan(file, nmax = min(max.vals,nz-skip.vals), skip = skip.vals, quiet = TRUE,what= list(i= integer(), j= integer(), x= numeric()))
              
              #checkIJ(els) 
              
             if (length(els$i)==0) return(c())
             
            
            if (length(intersect(els$j,current.condition))>0) # we have read cell of the current condition
              pointer=max(which(els$j==current.condition))
              else
                 if (length(unique(els$j))==1 | max.vals>(nz-skip.vals)) 
                   pointer=length(els$i) # keeping all read values
                 else   
                   pointer=max(which(els$j==(max(els$j)-1))) # keeping all values up to j-1 cell
             
             els$i=els$i[1:pointer]
             els$j=els$j[1:pointer]
             els$j=as.integer(els$j-min(els$j)+1)
             els$x=els$x[1:pointer]
             nc=max(els$j)
             #print(nc)
             #pointer=read.to/3
             
             checkIJ(els)
             
             switch(sym,
                    "general" = {
                      return(list(dgTMatrix=new("dgTMatrix", Dim = c(nr, nc), i = els$i - 1L,j = els$j - 1L, x = els$x),pointer=(pointer+skip.vals)))
                    },
                    "symmetric" = {
                      stop("general symmetry form 'symmetric' not yet implemented for reading")
                      #new("dsTMatrix", uplo = "L", Dim = c(nr, nc),
                      #    i = els$i - 1L, j = els$j - 1L, x = els$x)
                    },
                    "skew-symmetric" = {
                      stop("general symmetry form 'skew-symmetric' not yet implemented for reading")
                      ## FIXME: use dgT... but must expand the (i,j,x) slots!
                      new("dgTMatrix", uplo = "L", Dim = c(nr, nc),
                          i = els$i - 1L, j = els$j - 1L, x = els$x)
                      
                    },
                    "hermitian" = {
                      stop("general symmetry form 'hermitian' not yet implemented for reading")
                    },
                    ## otherwise (not possible; just defensive programming):
                    stop(gettextf("symmetry form '%s' is not yet implemented",
                                  sym), domain = NA)
             )
           },
           "pattern" = {
             els <- scan(file, nmax = nz, quiet = TRUE,
                         what = list(i = integer(), j = integer()))
             checkIJ(els)
             switch(sym,
                    "general" = {
                      stop("pattern form 'symmetric' not yet implemented for reading")
                      #new("ngTMatrix", Dim = c(nr, nc),
                      #    i = els$i - 1L, j = els$j - 1L)
                    },
                    "symmetric" = {
                      stop("pattern form 'symmetric' not yet implemented for reading")
                      #new("nsTMatrix", uplo = "L", Dim = c(nr, nc),
                      #    i = els$i - 1L, j = els$j - 1L)
                    },
                    "skew-symmetric" = {
                      stop("pattern form 'skew-symmetric' not yet implemented for reading")
                      ## FIXME: use dgT... but must expand the (i,j,x) slots!
                      new("ngTMatrix", uplo = "L", Dim = c(nr, nc),
                          i = els$i - 1L, j = els$j - 1L)
                      
                    },
                    "hermitian" = {
                      stop("pattern form 'hermitian' not yet implemented for reading")
                    },
                    ## otherwise (not possible; just defensive programming):
                    stop(gettextf("symmetry form '%s' is not yet implemented",
                                  sym), domain = NA)
             )
           },
           "complex" = {
             stop("element type 'complex' not yet implemented")
           },
           ## otherwise (not possible currently):
           stop(gettextf("'%s()' is not yet implemented for element type '%s'",
                         "readMM", elt), domain = NA))
  }
  else
    stop(gettextf("'%s()' is not yet implemented for  representation '%s'",
                  "readMM", repr), domain = NA)
}





pool.icell.fast = function (all.distances,all.neighbours,pooling.factor,cutoff,verbose){

  icells=c()

  #all.distances.out<<-all.distances
  #all.neighbours.out<<-all.neighbours
  
  
  #initializing min neighbours
  if (pooling.factor>=4) min.neighbours=2
  else min.neighbours=1
  
  if (verbose==TRUE)
    {
    print(sprintf('USING CUTOFF %.2f',cutoff))
    print(sprintf('Actual 0.1 quantile of local distances %.2f',quantile(all.distances,0.1)))
    print(sprintf('Sorting %g cells with %g neighbours...',nrow(all.distances),ncol(all.distances)))
    }
  
  for (k in 1:nrow(all.distances))# sorting cells for distances
      {
      ix=order(all.distances[k,])
      all.distances[k,]=all.distances[k,ix]
      all.neighbours[k,]=all.neighbours[k,ix]
      }
  
    
  #cleaning step to reduce the columns and fasten the code!!!!
  good.neighbours=Rfast::rowsums(all.distances<=cutoff)
  best.cell=max(good.neighbours)
    if (best.cell>1)
        {
        if (verbose==TRUE) print(sprintf('Best cell has %g good neighbours, restricting the data',best.cell))
        all.distances=all.distances[,1:best.cell]
        all.neighbours=all.neighbours[,1:best.cell]
        }
    
    
  if (verbose==TRUE) print('Starting to pool...')

  
  
  #ordering by average distance
  avg.dist=Rfast::rowmeans(all.distances)
  ix=order(avg.dist,decreasing = FALSE)
  all.distances=all.distances[ix,]
  all.neighbours=all.neighbours[ix,]
  good.neighbours=good.neighbours[ix]
  column.reference=ix
  rm(ix)
  gc()
  
  #marking useless distances
  ix=which(all.distances>cutoff)
  all.distances[ix]=NA
  all.neighbours[ix]=NA
  
  cell.available=rep(1,nrow(all.distances))
  used.cells=rep(NA,nrow(all.distances))
 
  
  counts.used=0
  
  for (k in 1:nrow(all.distances)) #main cycle
    {
    available=setdiff(all.neighbours[k,1:good.neighbours[k]],used.cells) #neighbours already sorted from clostest to fartest
    if (cell.available[column.reference[k]]>0)
      if (length(setdiff(available,NA))>=min.neighbours)
        {
        selected=available[1:min(pooling.factor,length(available))]
        dummy.pool=c(column.reference[k],selected,rep(NA,pooling.factor-length(selected)))
        used.cells[(counts.used+1):(counts.used+length(dummy.pool))]=dummy.pool
        counts.used=counts.used+length(dummy.pool)
        icells=rbind(icells,dummy.pool)
        cell.available[dummy.pool]=0
        }
    
    }
  return(icells)
}







pool.icell = function (all.distances,all.neighbours,pooling.factor,cutoff,verbose){
  
pooling.factor.start=pooling.factor  
icells=c()

column.reference=c(1:nrow(all.distances)) 
starting.cells=nrow(all.distances)
round.count=1

if (verbose==TRUE)  print(sprintf('USING CUTOFF %.2f',cutoff))

if (verbose==TRUE) print(sprintf('Actual 0.1 quantile of local distances %.2f',quantile(all.distances,0.1)))

while (TRUE)# Main loop
  
  {
  
  if (verbose==TRUE) print(sprintf('Sorting %g cells with %g neighbours...',nrow(all.distances),ncol(all.distances)))
  

  for (k in 1:nrow(all.distances))# sorting cells for distances
    {
      ix=order(all.distances[k,])
      all.distances[k,]=all.distances[k,ix]
      all.neighbours[k,]=all.neighbours[k,ix]
    }
  
  if (round.count==1)
    {
    #cleaning step to reduce the columns and fasten the code!!!!
    good.neighbours=Rfast::rowsums(all.distances<=cutoff)
    best.cell=max(good.neighbours)
    if (best.cell>1)
      {
      if (verbose==TRUE) print(sprintf('Best cell has %g good neighbours, restricting the data',best.cell))
      all.distances=all.distances[,1:best.cell]
      all.neighbours=all.neighbours[,1:best.cell]
      }
    }
  
  
  
  if (verbose==TRUE) print('Starting to pool...')
  sorted=FALSE
  while (ncol(all.distances)>=pooling.factor) # straight assignment
  {
   

    if (pooling.factor>1)
      available=which( Rfast::rowsums(is.na(all.distances[,1:pooling.factor])) == 0)
    else
      available=which( is.numeric(all.distances[,1:pooling.factor]) )  
    
    closest.pool=available[which.min(all.distances[available,pooling.factor])]  
    
    # print(closest.pool)
    # print(all.neighbours[closest.pool,1:pooling.factor])
    # print(all.distances[closest.pool,1:pooling.factor] ) 
    # cat(sprintf('\n'))
    
    if (length(closest.pool)>0)
      
      if (all.distances[closest.pool,pooling.factor]<cutoff )
        {
          dummy.pool=c(column.reference[closest.pool],all.neighbours[closest.pool,1:pooling.factor]) # dummy.pool in cell numbers
          dummy.pool=c(dummy.pool,rep(NA,pooling.factor.start-pooling.factor))

          icells=rbind(icells,dummy.pool)
          
          all.distances[which(is.element(column.reference,dummy.pool)),]=NA
          all.distances[which(is.element(all.neighbours,dummy.pool))]=NA
          sorted=TRUE
          # all.distances=all.distances[-pooled.cells,] # faster if I do it here?????
          # all.neighbours=all.neighbours[-pooled.cells,]
          # column.reference=column.reference[-pooled.cells]
        }
      else
         break # if closest distances have some NAs or they exceed limit
    else
      break # if there is no closets pool (all NAs)

  }
  
  unique.icells=setdiff(icells,NA)
  
  if (verbose==TRUE) print(sprintf('After pooling reached %g iCells containing %g original cells',nrow(icells),length(unique.icells)))
  
  if (length(unique.icells)>(starting.cells-pooling.factor))
    return(icells)
  
  if (sorted=='FALSE')# Could not pool any cell
  {
    pooling.factor=pooling.factor-1
    if(pooling.factor==0)
      break
    else
      if (verbose==TRUE) print(sprintf('Decreasing pooling factor to %g',pooling.factor))
  }
    
  # removing pooled cells
  
  pooled.cells=which(is.element(column.reference,unique.icells))
  if (length(pooled.cells)>0)
    {
    if (verbose==TRUE) print(sprintf('Removing %g used cells...',length(pooled.cells)))
    all.distances=all.distances[-pooled.cells,]
    all.neighbours=all.neighbours[-pooled.cells,]
    column.reference=column.reference[-pooled.cells]
    gc()
    }

  
  round.count=round.count+1
}

return(icells)
}








check.conditions = function (sample.conditions,verbose){
 
if (length(sample.conditions)>1)   
  {
  all.sample.conditions=unique(sample.conditions) # keeps in order of fisrt appearance
    
  avg.pos=c()
  for (k in 1:length(all.sample.conditions))
    avg.pos[k]=mean(which(sample.conditions==all.sample.conditions[k]))
  
  incremental=unique(diff(avg.pos)>0)
  
  if (length(incremental)>1) error('Conditions must be in incremental order, contact bigSCale developer gio.iacono.work@gmail.com')
  if (incremental==FALSE) error('Conditions must be in incremental order, contact bigSCale developer gio.iacono.work@gmail.com')
  
  
  sample.conditions=cumsum(as.numeric(table(sample.conditions)))
  
  if (verbose==TRUE) print('Detected conditions for splitting the data:')
  if (verbose==TRUE) print(sample.conditions)
}
  
return(sample.conditions)
  
}








#' iCells
#' 
#' Creates an iCell dataset starting from the typical outputs of CellRanger (.mex or .h5)
#' 
#' @param file.dir input dataset, file name
#' @param target.cells Number of iCells you would like to obtain
#' @param sample.conditions optional,  a factor indicating your sample conditions (for example, stages, treatments). This will prevent cells from different conditions to be pooled into the same iCell.
#' @param neighbours How many neighbours are used to search for a mate. Increasing it will increase fidelity of iCells, but only slighlty. In fact, multiple, iterative searches for neighbours are perfomed anyway.
#' @param verbose Whether to print on screen all the processing information

#' @param pooling Advanced use only. A technical parameter, you should not touch this.
#' @param q.cutoffs Advanced use only. The cutoff for pooling cells, default \bold{0.05}, meaning only cells closer than \bold{0.05} percentile are pooled. Deacresing the values yields better quality iCells but longer times, and viceversa.
#' @param preproc.cells Advanced use only. how many cells to use for the initial creation of the model. Reduce it if you have memory issues. Useless to increase it.
#' @param icells.chuncks Advanced use only. Size of the chunks of the original cells. Reduce it if you have memory issues. Probably useless to increase it.
#' @param preproc.chuncks Advanced use only. Reduce it if you have memory issues. Increase to fasten the initial step where the model is calculated.
#' @param min_ODscore Advanced use only. Increasing the value will result in using less highly variable genes for the creation of iCells.
#' @return  A list with two elements: icells.data (for you) and debugging (for me, if there are problems).In addition, it writes to current working directory the icells matrix automatically under the name icells.mtx.gz.
#'
#' \itemize{i
#' \item {\bold{icell.mat}} {icell Expression counts in the Matricx format}
#' \item {\bold{iCells}} {indices of the original cells}
#' \item {\bold{output.conditions}} {If you forced a pooling by condition (parameter \code{sample.conditions}) this vector contains the condition of each icell.}
#' \item {\bold{model}} {The numerical model used to compute cell to cell distances. Its meaning is explained in the first part (bigSCale core, advanced use) of the GitHub tutorial}
#' \item {\bold{driving.genes}} {The highly variable genes used to caculate cell to cell distances}
#' }
#' 
#' @examples
#' check the online tutorial at Github
#'
#' @export

iCells= function (file.dir,target.cells,sample.conditions=NA,neighbours=500,verbose=FALSE,pooling='relative',q.cutoffs=0.05,preproc.cells=10000,icells.chuncks=100000,preproc.chuncks=200000,min_ODscore=3){
  
  

  
# Reading header for sample info
out=bigscale.readMM(file=file.dir,size.only=TRUE)
tot.cells=out[1]

estimated.pooling=tot.cells/target.cells
print(sprintf('Attempting to reduce the size approximately %g times',estimated.pooling))

if (estimated.pooling<=5)
  {
  result.icells=iCells.simple(file.dir,pooling.factor=round(estimated.pooling-1),sample.conditions,pooling,q.cutoffs,verbose,preproc.cells,icells.chuncks,neighbours,preproc.chuncks,min_ODscore)
  return(result.icells)
  }

if (estimated.pooling>5 & estimated.pooling<=25)
  {
  number=sqrt(estimated.pooling)
  
  pooling1=ceiling(number)-1
  result1=iCells.simple(file.dir = file.dir,pooling.factor=pooling1,sample.conditions = sample.conditions,pooling = pooling,q.cutoffs = q.cutoffs,verbose=verbose,preproc.cells=preproc.cells,icells.chuncks=icells.chuncks,neighbours=neighbours,preproc.chuncks=preproc.chuncks,min_ODscore = min_ODscore,intermediate = TRUE)           
  if (pooling1>2)
   {
   preproc.cells=round(preproc.cells/(pooling1-1))
   icells.chuncks=round(icells.chuncks/(pooling1-1))
   print(sprintf('Reducing preproc.cells to %g and icells.chuncks to %g',preproc.cells,icells.chuncks))
   }
  
  pooling2=floor(number)-1
  result2=iCells.simple(file.dir = 'iCells.mtx.gz',pooling.factor=pooling2,sample.conditions = result1$output.conditions,pooling = pooling,q.cutoffs = q.cutoffs,verbose=verbose,preproc.cells=preproc.cells,icells.chuncks=icells.chuncks,neighbours=neighbours,preproc.chuncks=preproc.chuncks,min_ODscore = min_ODscore)
  
  
  iCells.final=deconvolute(result2$iCells,result1$iCells)
  result.final=result2
  result.final$iCells=iCells.final
  return(list(icell.data=result.final,debugging=list(intermediate.raw1=result1,intermediate.raw2=result2)))
  }
  
  
  
  
if (estimated.pooling>25)
{
  number=estimated.pooling^(1/3)
  
  pooling1=min(ceiling(number)-1,4)
  pooling2=min(ceiling(number)-1,4)
  pooling3=min(floor(number)-1,4)
  
  result1=iCells.simple(file.dir = file.dir,pooling.factor=pooling1,sample.conditions = sample.conditions,pooling = pooling,q.cutoffs = q.cutoffs,verbose=verbose,preproc.cells=preproc.cells,icells.chuncks=icells.chuncks,neighbours=neighbours,preproc.chuncks=preproc.chuncks,min_ODscore = min_ODscore,intermediate = TRUE)                  
  if (pooling1>2)
  {
    preproc.cells=round(preproc.cells/(pooling1-1))
    icells.chuncks=round(icells.chuncks/(pooling1-1))
    print(sprintf('Reducing preproc.cells to %g and icells.chuncks to %g',preproc.cells,icells.chuncks))
  }
  
  result2=iCells.simple(file.dir = 'iCells.mtx.gz',pooling.factor=pooling2,sample.conditions = result1$output.conditions,pooling = pooling,q.cutoffs = q.cutoffs,verbose=verbose,preproc.cells=preproc.cells,icells.chuncks=icells.chuncks,neighbours=neighbours,preproc.chuncks=preproc.chuncks,min_ODscore=min_ODscore,intermediate = TRUE) 
  
  result3=iCells.simple(file.dir = 'iCells.mtx.gz',pooling.factor=pooling3,sample.conditions = result2$output.conditions,pooling = pooling,q.cutoffs = q.cutoffs,verbose=verbose,preproc.cells=preproc.cells,icells.chuncks=icells.chuncks,neighbours=neighbours,preproc.chuncks=preproc.chuncks,min_ODscore=min_ODscore)
  
  
  iCells.final=deconvolute(result2$iCells,result1$iCells)
  iCells.final=deconvolute(result3$iCells,iCells.final)
  result.final=result3
  result.final$iCells=iCells.final
  return(list(icell.data=result.final,debugging=list(intermediate.raw1=result1,intermediate.raw2=result2,intermediate.raw3=result3)))
}

  

}


iCells.simple = function (file.dir,pooling.factor,sample.conditions,pooling,q.cutoffs,verbose,preproc.cells,icells.chuncks,neighbours,preproc.chuncks,min_ODscore,intermediate=FALSE){
 
library(SingleCellExperiment)
  


# Reading header for sample info
out=bigscale.readMM(file=file.dir,size.only=TRUE)

tot.cells=out[1]
tot.numbers=out[3]
avg.genes.cell=tot.numbers/tot.cells

print(sprintf('Total of %g cells, to be reduced with pooling factor %g (=%g+1)',tot.cells,pooling.factor+1,pooling.factor))

icells.chuncks=round(icells.chuncks*pooling.factor/4)
print(sprintf('Adjusting  icells.chuncks to %g cells',icells.chuncks))

# Initializing pre-processing parameters
proproc.chuncks.nz=preproc.chuncks*avg.genes.cell
downsamp=preproc.cells/tot.cells
if (downsamp>1) downsamp=1

print(sprintf('For the proprocessing, downsampling of %f',downsamp))

icells.chuncks.nz=icells.chuncks*avg.genes.cell

reps=round(tot.cells/100000)
if(reps==0) reps=1
if(reps>5) reps=5


if (verbose==TRUE) print(sprintf('reps=%g',reps))

driving.genes=c()

for (propro.rep in 1:reps)
      {
      # Reading dataset chunk by chunk to create a subset for pre-processing
      pointer=0
      round=1
      
      while (TRUE)
        {
          out=bigscale.readMM(file=file.dir,max.vals = proproc.chuncks.nz,skip.vals = pointer)
          if (length(out$dgTMatrix)==0) break
          
          ix.rand=sample( ncol(out$dgTMatrix),round(ncol(out$dgTMatrix)*downsamp))
          out$dgTMatrix=out$dgTMatrix[,ix.rand]
          
          if (round==1)
            tot.mat=out$dgTMatrix
          else
            tot.mat=cbind(tot.mat,out$dgTMatrix)
          
          print(sprintf('Incorporated %g cells for pre-processing',ncol(tot.mat)))
          pointer=out$pointer
          round=round+1
          rm(out)
          gc()
        }
      

        gene.names=c()
        for (k in 1:nrow(tot.mat))
        gene.names[k]=sprintf('Gene_%g',k)
        sce = SingleCellExperiment(assays = list(counts = tot.mat))
        rm(tot.mat)
        gc()
        rownames(sce)=gene.names
        
        
        sce=preProcess(sce)
        if (propro.rep>1) sce@int_metadata$edges=edges # overwriting edges with the initial ones
        sce = storeNormalized(sce,memory.save=FALSE)
        sce=setModel(sce)
        sce=setODgenes(sce,min_ODscore = min_ODscore)
        if (propro.rep==reps & pooling=='absolute') sce=setDistances(sce) # only in the last round and if absolute pooling
        
        if (propro.rep==1) 
          N=sce@int_metadata$N
        else
          N=N+sce@int_metadata$N
        edges=sce@int_metadata$edges
        dummy=which(sce@int_elementMetadata$ODgenes==1)
        driving.genes=union(driving.genes,sce@int_metadata$express.filtered[dummy])
        if (verbose==TRUE) print(sprintf('Pre-processing round %g, cumulated %g OD genes, sum(N)=%g',propro.rep,length(driving.genes),sum(N)))
      }

# Calculating percentiles
N_pct=matrix(NA,length(edges)-1,length(edges)-1)
for (k in 1:nrow(N))
  for (j in 1:nrow(N))
    N_pct[k,j] = sum(N[k,j:ncol(N)])/sum(N[k,])

N_pct[is.na(N_pct)]=1 # fix 0/0

if (verbose==TRUE) print(sprintf('Total sum of cases in N %g',sum(N)))

if (pooling=='absolute')
  {
  q.cutoff.narrow=quantile(sce@int_metadata$D,q.cutoffs[1])
  q.cutoff.large=quantile(sce@int_metadata$D,q.cutoffs[2])
  if (verbose==TRUE) print(sprintf('Estimated cutoffs for cell pooling: narrow %.2f (quantile %.2f), large %.2f (quantile %.2f)',q.cutoff.narrow,q.cutoffs[1],q.cutoff.large,q.cutoffs[2]))
  }
else
  {
  q.cutoff.narrow=q.cutoffs[1]
  q.cutoff.large=NA
  }
  



# all pre-proc data is now available. proceeding to convolute one piece at a time
rm(list=setdiff(setdiff(ls(), lsf.str()), c('file.dir','N_pct','edges','driving.genes','library.size','icells.chuncks.nz','avg.genes.cell','tot.cells','icells.chuncks','q.cutoff.narrow','q.cutoff.large','sample.conditions','pooling.factor','verbose','neighbours','intermediate')))

pointer=0
round=1
current.condition=1 # can be NA

condition.names=as.character(unique(sample.conditions))
print(condition.names)
sample.conditions=check.conditions(sample.conditions,verbose)

condition.pos=1
if (length(sample.conditions)>1)
  current.condition=sample.conditions[condition.pos]
else
  current.condition=NA

tot.cells.read=0
output.conditions=c()

pb <- progress::progress_bar$new(format = "Analyzing cells [:bar] :current/:total (:percent) eta: :eta", tot.cells)
pb$tick(0)

while (TRUE)
  {
  cells.to.read=round(tot.cells-pointer/avg.genes.cell)

  if (cells.to.read<(1.5*icells.chuncks)) icells.chuncks.nz=Inf

  out=compute.distances.icell(file.dir=file.dir , max.vals=icells.chuncks.nz , pointer = pointer , driving.genes = driving.genes , N_pct = N_pct, edges = edges,q.cutoff.narrow = q.cutoff.narrow,q.cutoff.large=q.cutoff.large,chunk.number = round,current.condition=current.condition,pooling.factor=pooling.factor,verbose=verbose,neighbours=neighbours)

  if(length(out)==0) break
  
  pointer=out$pointer
  
  if (round==1)
    iCells=out$iCells
  else
    iCells=rbind(iCells,(out$iCells+tot.cells.read))

  if (verbose==TRUE) print(range(out$iCells,na.rm=TRUE))
  if (verbose==TRUE) print(range(out$iCells+tot.cells.read,na.rm=TRUE))
  
  tot.cells.read=tot.cells.read+out$read.cells
  
  if (length(sample.conditions)>1)
    if (tot.cells.read>=current.condition) 
      {
      if (verbose==TRUE) print('Found a change of condition')
      output.conditions=c(output.conditions,rep(condition.names[condition.pos],nrow(out$iCells)))
      condition.pos=condition.pos+1
      current.condition=sample.conditions[condition.pos]
      }
    
  if (verbose==TRUE) print(sprintf('Total cells read (cumulative): %g, current group of iCells ranges from %g to %g, total iCells range from %g to %g',tot.cells.read,min(out$iCells,na.rm=TRUE),max(out$iCells,na.rm=TRUE),min(iCells,na.rm=TRUE),max(iCells,na.rm=TRUE)))
  
  pb$tick(out$read.cells)
  
  round=round+1
  gc()
  
  }

merge.chunks(verbose)

if (intermediate==FALSE)
  {
  icell.mat=Matrix::readMM("iCells.mtx.gz")
  return(list(icell.mat=icell.mat,iCells=iCells,output.conditions=output.conditions,model=N_pct,driving.genes=driving.genes))
  }
else
  {
  return(list(iCells=iCells,output.conditions=output.conditions,model=N_pct,driving.genes=driving.genes))
  }

}





#' Compute iCell distances
#'
#' @export

compute.distances.icell = function (file.dir,max.vals,pointer,driving.genes,N_pct,edges,q.cutoff.narrow,q.cutoff.large,chunk.number,current.condition,pooling.factor,verbose,neighbours){

 
out=bigscale.readMM(file=file.dir,max.vals = max.vals,skip.vals = pointer,current.condition=current.condition)
if (length(out$dgTMatrix)==0) return(c())
pointer.new=out$pointer
library.size=Matrix::colSums(out$dgTMatrix)
read.cells=ncol(out$dgTMatrix)

# normalize expression data for library size without scaling to the overall average depth
expr.driving.norm=as.matrix(out$dgTMatrix[driving.genes,])
rm(out)
gc()
for (k in 1:ncol(expr.driving.norm))
  expr.driving.norm[,k]=expr.driving.norm[,k]/library.size[k]
gc()


log.scores=get.log.scores(N_pct)

# .....................................................................................
# Allocating vector  .................................................................
indA.size=1000000
  
critical.value=max(library.size)*max(expr.driving.norm)
if (critical.value>indA.size/10)
  {
  indA.size=indA.size*50
    warning(sprintf('Critical value very high (%g): Increased memory usage!!',critical.value))
  if (critical.value>indA.size/10)
    stop(sprintf('Critical value way too high (%g): Stopping the analysis!!',critical.value))
  }
  
#print('Proceding to allocate large vector')
vector=c(0:indA.size)/10 # increase if code blocks, It can assign a gene exprssion level up to 10000000
ind.A=as.integer(cut(vector,edges,include.lowest = TRUE))
rm(vector)
ind.A=ind.A-1
# .....................................................................................
# ......................................................................................


#  ......................................................................................
# CORE PART
#  ......................................................................................

iCells.tot=c()
cell.reference=c(1:ncol(expr.driving.norm))
round=1
q.cutoff=q.cutoff.narrow # remains like that only if "absolute" pooling is used


while (TRUE) 
  
  {
  if (verbose==TRUE) print(sprintf('Starting from %g unpooled cells....',ncol(expr.driving.norm)))
  if (verbose==TRUE) print('')
  
  if (ncol(expr.driving.norm)<10) 
    {
    print('Less than 10 cells remaining, quitting')
    break
    }
  
  tot.cells=ncol(expr.driving.norm)
  sample.size=min(neighbours,(tot.cells-1))
  all.distances=matrix(0,tot.cells,sample.size)
  all.neighbours=matrix(0,tot.cells,sample.size)
  
  random.neighbours=sample(tot.cells,sample.size)
  
  if (verbose==TRUE) print('Computing distances ...')
  for ( k in 1:tot.cells)
    {
    random.neighbours=setdiff(sample(tot.cells,(sample.size+1)),k)
    random.neighbours=random.neighbours[1:sample.size]
    all.distances[k,]=distances_icells(expr.driving.norm[,c(k,random.neighbours)],log.scores,ind.A,library.size[c(k,random.neighbours)])
    all.neighbours[k,]=random.neighbours
    }
  
  if (is.na(q.cutoff.large)) # then we are using "relative" pooling
    {
    q.cutoff=quantile(all.distances,q.cutoff.narrow)
    if (verbose==TRUE) print(sprintf('Relative pooling: Estimated relative distance cutoff %.2f (quantile %.2f)',q.cutoff,q.cutoff.narrow))
    }
    
    
  if (verbose==TRUE) print('Launching iCells pooling ...')
  # Computing iCells ...............................................
  icells=pool.icell.fast(all.distances = all.distances,all.neighbours = all.neighbours,pooling.factor = pooling.factor , cutoff = q.cutoff,verbose)
  icells.local=icells
  if (length(icells)>0)
      {
      for (k in 1:nrow(icells)) # remapping icells to reference safe numbers
        icells[k,]=cell.reference[icells[k,]]
      iCells.tot=rbind(iCells.tot,icells)
      pooling.done=TRUE
      }
  else
      pooling.done=FALSE
  
  # Unpooled cells .................................................
  unpooled.cells=setdiff(1:tot.cells,icells.local)
  if (verbose==TRUE) print(sprintf('%g cells are still unpooled .... (pooling.factor=%g)',length(unpooled.cells),pooling.factor))
  
  
  # if (length(unpooled.cells)<=(pooling.factor-1)) # -1 because if we have 4 cells left, we only have 3 possible neighbours,for example. we are done, pooed enough
  #   break
  if (pooling.done==FALSE) # could not pool more even after looking for more neighbours
    if (q.cutoff==q.cutoff.large | is.na(q.cutoff.large)) # again we are done
      break
    else
      q.cutoff=q.cutoff.large
  
  expr.driving.norm=expr.driving.norm[,unpooled.cells]
  cell.reference=cell.reference[unpooled.cells]
  library.size=library.size[unpooled.cells]
  gc()
  round=round+1
  }

rm(list=setdiff(setdiff(ls(), lsf.str()), c('file.dir','max.vals','pointer','pointer.new','iCells.tot','chunk.number','read.cells','verbose')))
gc()


# reading again same part as beginning
if (verbose==TRUE) print('Reading again from source')
out=bigscale.readMM(file=file.dir,max.vals = max.vals,skip.vals = pointer)

#adding a fake all zero column
out$dgTMatrix=cbind(out$dgTMatrix,rep(0,nrow(out$dgTMatrix)))
fake.column=ncol(out$dgTMatrix)
# initializing other stuff
size.icells=ncol(iCells.tot)
iCells.mat=matrix(0,nrow(out$dgTMatrix),nrow(iCells.tot))

icells.at.once=500
starting.icell=1

while (TRUE)
  {
  gc()
  end.icell=min((starting.icell+icells.at.once-1),nrow(iCells.tot))
  
  iCells.tot.in.use=iCells.tot[starting.icell:end.icell, ]
  
  if (verbose==TRUE) print(sprintf('Processing iCells from %g to %g',starting.icell,end.icell))
  
  # initializing icells.vector
  icells.vector=as.vector(t(iCells.tot.in.use))
  icells.vector[is.na(icells.vector)]=fake.column
  jumps=c(0,c(1:nrow(iCells.tot.in.use))*size.icells)
  
  temp.full.reshuff=as.matrix(out$dgTMatrix[,icells.vector])

  
  #print(sprintf('k runs from %g to %g',1,(length(jumps)-1)))
  for (k in 1:(length(jumps)-1))
    {
    #print(sprintf('k=%g (out of %g),using %g -%g out of %g(%g)',k,length(jumps),jumps[k]+1,jumps[k+1],length(icells.vector),ncol(temp.full.reshuff)))
    iCells.mat[,k+starting.icell-1]=Rfast::rowsums(temp.full.reshuff[,(jumps[k]+1):jumps[k+1]]) 
    }
  starting.icell=end.icell
  
  if (starting.icell==nrow(iCells.tot))
    break
  
  }

# print('Creating iCells matrix ...')
# out$dgTMatrix=as.matrix(out$dgTMatrix)
# gc()
# for (k in 1:nrow(iCells.tot))
# {
#   icell.to.use=iCells.tot[k,!is.na(iCells.tot[k,])]
#   iCells.mat[,k]=Rfast::rowsums(out$dgTMatrix[,icell.to.use]) # fix NA problem
# }
  
iCells.mat=Matrix::Matrix(iCells.mat)
gc()

bigscale.writeMM(data = iCells.mat,file.name = sprintf('Temp_Part_%g.mtx.gz',chunk.number))


if (verbose==TRUE) print(dim(iCells.tot))

return(list(iCells=iCells.tot,pointer=pointer.new,read.cells=read.cells))
}








compute.max.inter = function (D,clusters) {

all.inter=c()
for (k in 1:max(clusters))
  for (h in 1:max(clusters))
    all.inter=c(all.inter,mean(D[which(clusters==h),which(clusters==k)]))

print(sprintf('Computed maximum inter-cluster distance %g',max(all.inter)))

max.all.inter=max(all.inter)
return(max.all.inter)
    
}





bigscale.recursive.clustering = function (expr.data.norm,model,edges,lib.size,fragment=FALSE,create.distances=FALSE,modality) {
  
 
gc()
  

num.samples=ncol(expr.data.norm)
  
if (fragment==FALSE)
  {
  # Adjusting max_group_size according to cell number
  if (num.samples<1000) dim.cutoff=10
  if (num.samples>=1000 & num.samples<5000) dim.cutoff=50
  if (num.samples>=5000 & num.samples<10000) dim.cutoff=100
  if (num.samples>=10000) dim.cutoff=150
  }
else
{
  if (fragment==TRUE)
      {
      if (num.samples<1000) dim.cutoff=10
      else dim.cutoff=50
      }
  else
    {
    dim.cutoff=fragment*ncol(expr.data.norm)
    cat(sprintf('\nClustering to groups of at most %g cells (%g %%)',dim.cutoff,fragment*100))
    }
    
}
    

print(sprintf('Clustering cells down to groups of approximately %g-%g cells',dim.cutoff,dim.cutoff*5))   
#dim.cutoff=ncol(expr.data.norm)*min.group.size  

mycl=rep(1,ncol(expr.data.norm))
tot.recursive=1
#current.cutting=40
unclusterable=rep(0,length(mycl))

if (create.distances==TRUE) D.final=matrix(0,ncol(expr.data.norm),ncol(expr.data.norm))

while(1){

  cat(sprintf('\nRecursive clustering, beginning round %g ....',tot.recursive))
  action.taken=0
  mycl.new=rep(0,length(mycl))
  
  
  for (k in 1:max(mycl))
  {
    #print(sprintf('Checking cluster %g/%g, %g cells',k,max(mycl),length(which(mycl==k))))
    if (length(which(mycl==k))>dim.cutoff & sum(unclusterable[which(mycl==k)])==0 ) # then it must be re-clustered
    {
      #print('Computing Overdispersed genes ...')
      ODgenes=calculate.ODgenes(expr.data.norm[,which(mycl==k)],verbose = FALSE,min_ODscore = 2)
      
      if (is.null(ODgenes))
        {
        D=NA
        unclusterable[which(mycl==k)]=1   
        }
      else
        {
        dummy=as.matrix(ODgenes[[1]])
        ODgenes=which(dummy[,1]==1)
        #print('Computing distances ...')
        D=compute.distances(expr.norm = expr.data.norm[,which(mycl==k)],N_pct = model,edges = edges,driving.genes = ODgenes,lib.size = lib.size[which(mycl==k)],modality=modality)
        temp.clusters=bigscale.cluster(D,plot.clusters = FALSE,clustering.method = 'low.granularity',granularity.size=dim.cutoff,verbose=FALSE)$clusters #cut.depth=current.cutting,method.treshold = 0.2
        if (max(temp.clusters)>1) 
          action.taken=1
        else
          unclusterable[which(mycl==k)]=1  
        }
      
      
      if (create.distances==TRUE)
        {
        D.final[which(mycl==k),which(mycl==k)]=D # assigning distances
        max.inter=compute.max.inter(D,temp.clusters)
        D.final[which(mycl==k),which(mycl!=k)]=D.final[which(mycl==k),which(mycl!=k)]+max.inter
        D.final[which(mycl!=k),which(mycl==k)]=D.final[which(mycl!=k),which(mycl==k)]+max.inter
        }
      #print(sprintf('Partitioned cluster %g/%g in %g sub clusters ...',k,max(mycl),max(temp.clusters)))
      mycl.new[which(mycl==k)]=temp.clusters+max(mycl.new)
      
    }
    else
    {
      #print(sprintf('Cluster %g/%g is of size %g and cannot be further partitioned',k,max(mycl),length(which(mycl==k))))
      mycl.new[which(mycl==k)]=1+max(mycl.new)
    }
  }
  
  tot.recursive=tot.recursive+1  
  
  cat(sprintf('\nRecursive clustering, after round %g obtained %g clusters',tot.recursive,max(mycl.new)))
  # if (modality=='jaccard' & tot.recursive==5)
  # {
  #   cat(sprintf('\n\nYou are analyzing ATAC-seq data, I stop the clustering to avoid creating too many clusters. If you want to customize the analysis with more/less rounds of recursive contact the developer at gio.iacono.work@gmail.com'))
  #   break
  # }
  #current.cutting=current.cutting+10
  if (action.taken==0) 
    break
  else
    mycl=mycl.new 
  
  
  
  #if(tot.recursive==5) break
}

# Preventing small clusters to be used, with a ugly trick

bad.cells=c()
bad.clusters=c()
indexes=list()
for (k in 1:max(mycl))
    {
    indexes[[k]]=which(mycl==k)
    if (length(which(mycl==k))<5)
      {
      bad.clusters=c(bad.clusters,k)
      bad.cells=c(bad.cells,which(mycl==k))
      }
    }
if (length(bad.clusters)>0)
  {
  indexes=indexes[-bad.clusters]
  mycl=rep(0,length(mycl))
  for (k in 1:length(indexes))
    mycl[indexes[[k]]]=k
  warning(sprintf('Hidden %g clusters of %g cells total to be used: they are too small',length(bad.clusters),length(bad.cells)))
}

if (create.distances==FALSE)
  D.final=NA
return(list(mycl=mycl,D=D.final))
}


#' Compare gene centralities 
#'
#' Works with any given number of networks (N).  The centralities previously calculated with \code{compute.network()} for N networks are given as input.
#' The script sorts the genes accoring to their change in centrality between the first network (first element of the input list) and the other N-1 networks (the rest of the list)
#'
#' @param centralities List of at least two elements. Each element must be the (\code{data.frames}) of the centralities previously calculated by \code{compute.network()}. 
#' @param names.conditions character of the names of the input networks, same length of the \bold{compute.network()}
#' 
#' @return  A list with a (\code{data.frame}) for each centrality. In each (\code{data.frame}) the genes are ranked for thier change in centrality.
#'
#'
#' @examples
#' check the online tutorial at Github
#'
#' @export
 
compare.centrality <- function(centralities,names.conditions)
{
  
  
  centrality.names=colnames(centralities[[1]])
  
  total.genes=c()
  gene.set=list()
  pos=list()
  
  for (k in 1:length(centralities))
  {
    gene.set[[k]]=rownames(centralities[[k]])
    total.genes=union(total.genes,gene.set[[k]])
  }
  
  for (k in 1:length(centralities))    
    pos[[k]]=id.map(gene.set[[k]],total.genes)
  
  output=list()  
  
  for (k in 1:length(centrality.names))
  {
    result=matrix(0,length(total.genes),length(centralities)+2)
    
    for (j in 1:length(centralities))
    {
      dummy.vector=centralities[[j]][,k]
      result[pos[[j]],j]=dummy.vector
    }
    
    
    if (length(centralities)>2)
    {
      result[,length(centralities)+1]=result[,1]-(apply(X = result[,2:length(centralities)],MARGIN = 1,FUN = max))
      ranking=order(result[,length(centralities)+1],decreasing = TRUE)
      #dummy=which(ranking>(length(ranking)/2))
      #ranking[dummy]=length(ranking)-ranking[dummy]
      result[,length(centralities)+2]=ranking
    }
    else
    {
      result[,3]=(result[,1]-result[,2])
      ranking=rank(result[,3])
      # dummy=which(ranking>(length(ranking)/2))
      # ranking[dummy]=length(ranking)-ranking[dummy]
      result[,4]=ranking
    }
    
    

    
    table.title=c()
    for (j in 1:length(centralities))
      table.title[j]=sprintf('%s.%s',centrality.names[k],names.conditions[j])
    table.title[length(table.title)+1]='DELTA'
    table.title[length(table.title)+1]='Ranking'
    
    result=data.frame(result)
    colnames(result)=table.title
    rownames(result)=total.genes
    

    result=result[order(result[,length(centralities)+1]),]
    
    output[[k]]=result
  }
  
  names(output)=centrality.names
  
  return(output)
  
}   



polish.graph = function (G)
{
  # if(organism=='human')
  #   {
  #   org.ann.1=org.Hs.egALIAS2EG
  #   org.ann.2=org.Hs.egENSEMBL
  #   }
  #   else 
  #     if (organism=='mouse')
  #       {
  #       org.ann.1=org.Mm.egALIAS2EG 
  #       org.ann.2=org.Mm.egENSEMBL
  #       }
  #     else
  #       stop('Ask the programmer to add other organisms')
  
  gene.names=igraph::vertex_attr(graph = G,name = 'name')
  
  mapped=list()
  org.ann=list()
  org.ann[[1]]=as.list(org.Hs.eg.db::org.Hs.egALIAS2EG)
  org.ann[[2]]=as.list(org.Hs.eg.db::org.Hs.egENSEMBL2EG)
  org.ann[[3]]=as.list(org.Mm.eg.db::org.Mm.egALIAS2EG)
  org.ann[[4]]=as.list(org.Mm.eg.db::org.Mm.egENSEMBL2EG)
  class.names=c('Human, Gene Symbol','Human, ENSEMBL','Mouse, Gene Symbol','Mouse, ENSEMBL')
  mapped$map1=which(is.element(gene.names,names(org.ann[[1]])))
  mapped$map2=which(is.element(gene.names,names(org.ann[[2]])))
  mapped$map3=which(is.element(gene.names,names(org.ann[[3]])))
  mapped$map4=which(is.element(gene.names,names(org.ann[[4]])))
  
  hits=unlist(lapply(mapped, length))
  best.hit=which(hits==max(hits))
  rm(mapped,org.ann)
  gc()
  print(sprintf('Recognized %g/%g (%.2f%%) as %s',max(hits),length(gene.names),max(hits)/length(gene.names)*100,class.names[best.hit]))
  
  if (best.hit==1 | best.hit==2) organism.detected='human'
  if (best.hit==3 | best.hit==4) organism.detected='mouse'
  if (best.hit==1 | best.hit==3) code.detected='gene.name'
  if (best.hit==2 | best.hit==4) code.detected='ensembl'
  
  if (organism.detected=='human')  GO.ann = as.list(org.Hs.eg.db::org.Hs.egGO2ALLEGS)
  if (organism.detected=='mouse')  GO.ann = as.list(org.Mm.eg.db::org.Mm.egGO2ALLEGS)

  regulators.entrez <- GO.ann$`GO:0010468`
  
  if (organism.detected=='human' & code.detected=='gene.name')  org.ann=as.list(org.Hs.eg.db::org.Hs.egSYMBOL)
  if (organism.detected=='human' & code.detected=='ensembl')  org.ann=as.list(org.Hs.eg.db::org.Hs.egENSEMBL)
  if (organism.detected=='mouse' & code.detected=='gene.name') org.ann=as.list(org.Mm.eg.db::org.Mm.egSYMBOL)
  if (organism.detected=='mouse' & code.detected=='ensembl') org.ann=as.list(org.Mm.eg.db::org.Mm.egENSEMBL)
  
  
  regulators=unique(unlist(org.ann[regulators.entrez]))
  
  end.nodes=igraph::ends(G,igraph::E(G))
  testing=Rfast::rowsums(cbind(is.element(end.nodes[,1],regulators),is.element(end.nodes[,2],regulators)))
  G=igraph::delete_edges(G, which(testing==0))
  G=igraph::delete_vertices(G, which(igraph::degree(G)==0))
  
  return(G)
  
}
  

#' ATAC-seq Gene regulatory network
#'
#' Infers the gene regulatory network from single cell ATAC-seq data
#'
#' @param expr.data matrix of expression counts. Works also with sparse matrices of the \pkg{Matrix} package.
#' @param gene.names character of gene names, now it supports Gene Symbols or Ensembl, Mouse and Human.
#' @param clustering type of clustering and correlations computed to infer the network.
#' \itemize{
#'  \item {\bold{recursive}} Best quality at the expenses of computational time. If the dataset is larger than 10-15K cells and is highly heterogeneous this can lead to very long computational times (24-48 hours depending of the hardware).
#'  \item {\bold{direct}} Best trade-off between quality and computational time. If you want to get a quick output not much dissimilar from the top quality of \bold{recursive} one use this option. Can handle quickly also large datasets (>15-20K cells in 30m-2hours depending on hardware)
#'  \item {\bold{normal}} To be used if the correlations (the output value \bold{cutoff.p}) detected with either \bold{direct} or \bold{recursive} are too low. At the moment, bigSCale displays a warning if the correlation curoff is lower than 0.8 and suggests to eithe use \bold{normal} clustering or increase the input parameter \bold{quantile.p}
#' }
#' @param quantile.p only the first \eqn{1 - quantile.p} correlations are used to create the edges of the network. If the networ is too sparse(dense) decrease(increase) \eqn{quantile.p}
#' @param speed.preset Used only if  \code{clustering='recursive'} . It regulates the speed vs. accuracy of the Zscores calculations. To have a better network quality it is reccomended to use the default \bold{slow}.
##' \itemize{
#'   \item {\bold{slow}} {Highly reccomended, the best network quality but the slowest computational time.} 
#'   \item {\bold{normal}} {A balance between network quality and computational time. }
#'   \item {\bold{fast}} {Fastest computational time, worste network quality.}
#' }
#' @param previous.output previous output of \code{compute.network()} can be passed as input to evaluate networks with a different quantile.p without re-running the code. Check the online tutorial at https://github.com/iaconogi/bigSCale2.
#' 
#' @return  A list with the following items:
#' \itemize{
#' \item {\bold{centrality}} {Main output: a Data-frame with the network centrality (Degree,Betweenness,Closeness,PAGErank) for each gene(node) of the network}
#' \item {\bold{graph}} {The regulatory network in iGraph object}
#' \item {\bold{correlations}} {All pairwise correlations between genes. The correlation is an average between \emph{Pearson} and \emph{Spearman}. Note that it is stored in single precision format (to save memory space) using the package \pkg{float32}.To make any operation or plot on the correlations first transform it to the standard double precisione by running \code{correlations=dbl(correlations)} }
#' \item {\bold{cutoff.p}} {The adptive cutoff used to select significant correlations}
#' \item {\bold{tot.scores}} {The Z-scores over which the correlations are computed. The visually check the correlation between to genes \emph{i} and \emph{j} run  \code{plot(tot.scores[,i],tot.scores[,j])} }
#' \item {\bold{clusters}} {The clusters in which the cells have been partitioned}
#' \item {\bold{model}} {Bigscale numerical model of the noise}
#' }
#'
#'
#' @examples
#' out=compute.network(expr.data,gene.names)
#'
#' @export


compute.atac.network = function (expr.data,feature.file,quantile.p=0.998){
  
  clustering='direct'
  print('Pharsing the list of 10x annotated peaks')
  out=pharse.10x.peaks(feature.file)
  print(sprintf('Found %g peaks in promoters',length(out$good.ix)))
  expr.data=expr.data[out$good.ix,]
  feature.names=out$feature.names
  peaks.in.use=out$good.ix
  rm(out)
  gc()
  
  if (ncol(expr.data)>20000)
      warning('It seems you are running compute.network on a kind of large dataset, it it failed for memory issues you should try the compute.network for large datasets')

  print('1) Pre-processing) Removing null features ')
  exp.genes=which(Matrix::rowSums(expr.data)>0)
  if ((nrow(expr.data)-length(exp.genes))>0)
      print(sprintf("Discarding %g peaks with all zero values",nrow(expr.data)-length(exp.genes)))
  expr.data=expr.data[exp.genes,]
  gc()
  feature.names=feature.names[exp.genes]
  expr.data=Matrix::Matrix(expr.data)
  tot.cells=Matrix::rowSums(expr.data)
  
  print('2) Clustering ...')
  
  if (clustering=="direct")
      mycl=bigscale.recursive.clustering(expr.data.norm = expr.data,fragment=TRUE,modality = 'jaccard')$mycl

  tot.clusters=max(mycl)
    
  
  pass.cutoff=which(tot.cells>(max(15,ncol(expr.data)*0.005)))
  feature.names=feature.names[pass.cutoff]
    
  print(sprintf('Assembling cluster average expression for %g features detected in at least %g cells',length(pass.cutoff),max(15,ncol(expr.data)*0.005)))
  tot.scores=matrix(0,length(pass.cutoff),tot.clusters)
  for (k in 1:(tot.clusters))
        tot.scores[,k]=Rfast::rowmeans(as.matrix(expr.data[pass.cutoff,which(mycl==k)]>0))
  #tot.scores=log2(tot.scores+1)
  tot.scores=t(tot.scores)

  o=gc()
  print(o)
  print('Calculating Pearson ...')
  rm(list=setdiff(setdiff(ls(), lsf.str()), c('tot.scores','quantile.p','feature.names','tot.scores','mycl','peaks.in.use')))
  o=gc()
  print(o)
  
  Dp=Rfast::cora(tot.scores)
  o=gc()
  print(o)
  print('Calculating quantile ...')
  cutoff.p=quantile(Dp,quantile.p,na.rm=TRUE)
  print(sprintf('Using %f as cutoff for pearson correlation',cutoff.p))
  

  Dp=float::fl(Dp)
  gc()
  
  print('Calculating Spearman ...')
  Ds=cor(x = tot.scores,method = "spearman")
  Ds=float::fl(Ds)
  gc()
  
  print('Calculating the significant links ...')
  rm(list=setdiff(setdiff(ls(), lsf.str()), c("Dp","Ds","cutoff.p",'feature.names','tot.scores','mycl','peaks.in.use')))  
  gc()
  network=((Dp>cutoff.p & Ds>0.7) | (Dp<(-cutoff.p) & Ds<(-0.7)))
  diag(network)=FALSE
  network[is.na(network)]=FALSE
  degree=Rfast::rowsums(network)
  to.include=which(degree>0)
  
  G=igraph::graph_from_adjacency_matrix(adjmatrix = network[to.include,to.include],mode = 'undirected')
  G=igraph::set_vertex_attr(graph = G,name = "name", value = feature.names[to.include])
  rm(network)
  gc()
  
  
  print('Calculating the final score ...')
  Df=(Ds+Dp)/float::fl(2)
  rm(Dp)
  rm(Ds)
  gc()
  rownames(Df)=feature.names
  colnames(Df)=feature.names
  
  
  print(sprintf('Inferred the raw regulatory network: %g nodes and %g edges (ratio E/N)=%f',length(igraph::V(G)),length(igraph::E(G)),length(igraph::E(G))/length(igraph::V(G))))
  
  print('Computing the centralities')
  Betweenness=igraph::betweenness(graph = G,directed=FALSE,normalized = TRUE)
  Degree=igraph::degree(graph = G)
  PAGErank=igraph::page_rank(graph = G,directed = FALSE)$vector
  Closeness=igraph::closeness(graph = G,normalized = TRUE)

  if (cutoff.p<0.7)
    warning('bigSCale: the cutoff for the correlations seems very low. You should either increase the parameter quantile.p or select clustering=normal (you need to run the whole code again in both options,sorry!). For more information check the quick guide online')
  
  return(list(graph=G,correlations=Df,tot.scores=tot.scores,clusters=mycl,centrality=as.data.frame(cbind(Degree,Betweenness,Closeness,PAGErank)),cutoff.p=cutoff.p,peaks.in.use=peaks.in.use[to.include]))
}




#' Compute network model
#'
#' Compute network model
#'
#' @param expr.data matrix of expression counts. Works also with sparse matrices of the \pkg{Matrix} package.
#' @return  Model
#'
#'
#' @examples
#' out=compute.network.model(expr.data,gene.names)
#'
#' @export


compute.network.model = function (expr.data)
{
  sce = SingleCellExperiment(assays = list(counts = expr.data))  
  sce = preProcess(sce)
  sce = storeNormalized(sce)
  sce=setModel(sce)
  return(list(model=sce@int_metadata$model,edges=sce@int_metadata$edges))
}
  
  
#' Gene regulatory network
#'
#' Infers the gene regulatory network from single cell data
#'
#' @param expr.data matrix of expression counts. Works also with sparse matrices of the \pkg{Matrix} package.
#' @param gene.names character of gene names, now it supports Gene Symbols or Ensembl, Mouse and Human.
#' @param clustering type of clustering and correlations computed to infer the network.
#' \itemize{
#'  \item {\bold{recursive}} Best quality at the expenses of computational time. If the dataset is larger than 10-15K cells and is highly heterogeneous this can lead to very long computational times (24-48 hours depending of the hardware).
#'  \item {\bold{direct}} Best trade-off between quality and computational time. If you want to get a quick output not much dissimilar from the top quality of \bold{recursive} one use this option. Can handle quickly also large datasets (>15-20K cells in 30m-2hours depending on hardware)
#'  \item {\bold{normal}} To be used if the correlations (the output value \bold{cutoff.p}) detected with either \bold{direct} or \bold{recursive} are too low. At the moment, bigSCale displays a warning if the correlation curoff is lower than 0.8 and suggests to eithe use \bold{normal} clustering or increase the input parameter \bold{quantile.p}
#' }
#' @param quantile.p To select how many correlatoins are retained. By default \eqn{quantile.p=0.9} meaning that only Pearson Correlatins>0.9 are retained to create the network. If the number is higher than 1 then it is interpreted as a quantile, meaning tha only the first \eqn{1 - quantile.p} correlations are used to create the edges of the network. 
#' If \eqn{0<quantile.p<1} ten it is interprested directly as the minimum Pearson correlation. If the networ is too sparse(dense) decrease(increase) \eqn{quantile.p}
#' @param speed.preset Used only if  \code{clustering='recursive'} . It regulates the speed vs. accuracy of the Zscores calculations. To have a better network quality it is reccomended to use the default \bold{slow}.
##' \itemize{
#'   \item {\bold{slow}} {Highly reccomended, the best network quality but the slowest computational time.} 
#'   \item {\bold{normal}} {A balance between network quality and computational time. }
#'   \item {\bold{fast}} {Fastest computational time, worste network quality.}
#' }
#' @param previous.output previous output of \code{compute.network()} can be passed as input to evaluate networks with a different quantile.p without re-running the code. Check the online tutorial at https://github.com/iaconogi/bigSCale2.
#' 
#' @return  A list with the following items:
#' \itemize{
#' \item {\bold{centrality}} {Main output: a Data-frame with the network centrality (Degree,Betweenness,Closeness,PAGErank) for each gene(node) of the network}
#' \item {\bold{graph}} {The regulatory network in iGraph object}
#' \item {\bold{correlations}} {All pairwise correlations between genes. The correlation is an average between \emph{Pearson} and \emph{Spearman}. Note that it is stored in single precision format (to save memory space) using the package \pkg{float32}.To make any operation or plot on the correlations first transform it to the standard double precisione by running \code{correlations=dbl(correlations)} }
#' \item {\bold{cutoff.p}} {The adptive cutoff used to select significant correlations}
#' \item {\bold{tot.scores}} {The Z-scores over which the correlations are computed. The visually check the correlation between to genes \emph{i} and \emph{j} run  \code{plot(tot.scores[,i],tot.scores[,j])} }
#' \item {\bold{clusters}} {The clusters in which the cells have been partitioned}
#' \item {\bold{model}} {Bigscale numerical model of the noise}
#' }
#'
#'
#' @examples
#' out=compute.network(expr.data,gene.names)
#'
#' @export

  
compute.network = function (expr.data,gene.names,modality='pca',model=NA,clustering='recursive',quantile.p=0.9,speed.preset='slow',previous.output=NA){

if (is.na(previous.output))  
  {
  
  if (ncol(expr.data)>20000)
    warning('It seems you are running compute.network on a kind of large dataset, it it failed for memory issues you should try the compute.network for large datasets')
  expr.data=as.matrix(expr.data)
  
  
  # Part 1) Initial processing of the dataset  ************************************************************
  print('Pre-processing) Removing null rows ')
  exp.genes=which(Rfast::rowsums(expr.data)>0)
  if ((nrow(expr.data)-length(exp.genes))>0)
    print(sprintf("Discarding %g genes with all zero values",nrow(expr.data)-length(exp.genes)))
  expr.data=expr.data[exp.genes,]
  gc()
  gene.names=gene.names[exp.genes]
  tot.cells=Rfast::rowsums(expr.data>0)
  
  print('PASSAGE 1) Setting the size factors ....')
  lib.size = Rfast::colsums(expr.data)
    
    
  
  if (is.list(model))
  {
  print('Using edges given as input')
  edges=model$edges
  }
  else
  {
    print('PASSAGE 2) Setting the bins for the expression data ....')
    edges=generate.edges(expr.data)
  }
    
  print('PASSAGE 3) Storing in the single cell object the Normalized data ....')
  avg.library.size=mean(lib.size)
  for (k in 1:ncol(expr.data)) expr.data[,k]=expr.data[,k]/lib.size[k]*avg.library.size
  expr.data.norm=expr.data
  rm(expr.data)
  expr.data.norm=Matrix::Matrix(expr.data.norm)
  gc()
  
  if (is.list(model))
    {
    print('Using Model given in the input')
    model=model$model
    }
  else
    if (!(modality=='pca' & clustering=='direct'))
    {  
    print('PASSAGE 4) Computing the numerical model (can take from a few minutes to 30 mins) ....')
    model=fit.model(expr.data.norm,edges,lib.size)$N_pct
    gc()
    }
    

    
  print('PASSAGE 5) Clustering ...')
  if (clustering=="direct" | clustering=="recursive")
    {
    if (clustering=="direct")
      mycl=bigscale.recursive.clustering(expr.data.norm = expr.data.norm,model = model,edges = edges,lib.size = lib.size,fragment=TRUE,modality=modality)$mycl
    if (clustering=="recursive")
      mycl=bigscale.recursive.clustering(expr.data.norm = expr.data.norm,model = model,edges = edges,lib.size = lib.size,modality=modality)$mycl
    }
  else
    {
    ODgenes=calculate.ODgenes(expr.data.norm)
    dummy=as.matrix(ODgenes[[1]])
    ODgenes=which(dummy[,1]==1)
    print('Computing distances ...')
    D=compute.distances(expr.norm = expr.data.norm,N_pct = model,edges = edges,driving.genes = ODgenes,lib.size = lib.size)
    mycl=bigscale.cluster(D,plot.clusters = TRUE,method.treshold = 0.5)$clusters
    }
  tot.clusters=max(mycl)
  
  
  
  
  #filtering the number of genes
  pass.cutoff=which(tot.cells>(max(15,ncol(expr.data.norm)*0.005)))
  #pass.cutoff=which(tot.cells>0)
  gene.names=gene.names[pass.cutoff]
  
  
  if (clustering=="direct")
      {
      print(sprintf('Assembling cluster average expression for %g genes expressed in at least %g cells',length(pass.cutoff),max(15,ncol(expr.data.norm)*0.005)))
      tot.scores=matrix(0,length(pass.cutoff),tot.clusters)
      for (k in 1:(tot.clusters))
          tot.scores[,k]=Rfast::rowmeans(as.matrix(expr.data.norm[pass.cutoff,which(mycl==k)]))
      tot.scores=log2(tot.scores+1)
      }
    else
      {
      print(sprintf('Calculating Zscores for %g genes expressed in at least %g cells',length(pass.cutoff),max(15,ncol(expr.data.norm)*0.005)))  
      Mscores=calculate.marker.scores(expr.norm = expr.data.norm[pass.cutoff,], clusters=mycl, N_pct=model, edges = edges, lib.size = lib.size, speed.preset = speed.preset)$Mscores
  
      rm(expr.data.norm)
      gc()
      
      tot.scores=matrix(0,length(Mscores[[1,2]]),(tot.clusters*(tot.clusters-1)/2))
      
      # assembling the matrix of Mscores
      counts=1
      for (k in 1:(tot.clusters-1))
        for (h in (k+1):tot.clusters)
        {
          tot.scores[,counts]=Mscores[[k,h]]
          counts=counts+1 
        }
      }
  tot.scores=t(tot.scores)
  }
  
  else
  {
    print('It appears you want to tweak previously created networks with a different quantile.p, proceeding ....')
    tot.scores=previous.output$tot.scores
    gene.names=colnames(previous.output$correlations)
    mycl=previous.output$clusters
    model=previous.output$model
    rm(previous.output)
    gc()
  }
  
  




o=gc()
print(o)
print('Calculating Pearson ...')
rm(list=setdiff(setdiff(ls(), lsf.str()), c('tot.scores','quantile.p','model','gene.names','tot.scores','mycl')))
o=gc()
print(o)

Dp=Rfast::cora(tot.scores)
o=gc()
print(o)

if (quantile.p>1)
  {
  print('Calculating correlation quantile associated with your quantile.p ...')
  cutoff.p=quantile(Dp,quantile.p/100,na.rm=TRUE)
  }
else
  {
  print('Detacted quantile.p<1, using it directly as correlation treshold ...')
  cutoff.p=quantile.p
  }

print(sprintf('Using %f as cutoff for pearson correlation',cutoff.p))


Dp=float::fl(Dp)
gc()

print('Calculating Spearman ...')
Ds=cor(x = tot.scores,method = "spearman")
Ds=float::fl(Ds)
gc()

print('Calculating the significant links ...')
rm(list=setdiff(setdiff(ls(), lsf.str()), c("Dp","Ds","cutoff.p",'gene.names','tot.scores','mycl','model')))  
gc()
network=((Dp>cutoff.p & Ds>0.7) | (Dp<(-cutoff.p) & Ds<(-0.7)))
diag(network)=FALSE
network[is.na(network)]=FALSE
degree=Rfast::rowsums(network)
to.include=which(degree>0)

G=igraph::graph_from_adjacency_matrix(adjmatrix = network[to.include,to.include],mode = 'undirected')
G=igraph::set_vertex_attr(graph = G,name = "name", value = gene.names[to.include])
rm(network)
gc()


print('Calculating the final score ...')
Df=(Ds+Dp)/float::fl(2)
rm(Dp)
rm(Ds)
gc()
rownames(Df)=gene.names
colnames(Df)=gene.names


print(sprintf('Inferred the raw regulatory network: %g nodes and %g edges (ratio E/N)=%f',length(igraph::V(G)),length(igraph::E(G)),length(igraph::E(G))/length(igraph::V(G))))

G=polish.graph(G)
comp=igraph::components(G)
small.comp=which(comp$csize<0.01*sum(comp$csize))
to.remove=which(is.element(comp$membership,small.comp))
G=igraph::delete_vertices(G, to.remove)
comp=igraph::components(G)
cat('\n')
print(sprintf('Final network after GO filtering: %g nodes and %g edges (ratio E/N)=%f and %g components',length(igraph::V(G)),length(igraph::E(G)),length(igraph::E(G))/length(igraph::V(G)), comp$no))


print('Computing the centralities')
Betweenness=igraph::betweenness(graph = G,directed=FALSE,normalized = TRUE)
Degree=igraph::degree(graph = G)
PAGErank=igraph::page_rank(graph = G,directed = FALSE)$vector
Closeness=igraph::closeness(graph = G,normalized = TRUE)

G=igraph::set.vertex.attribute(graph = G,name = 'bigSCale.Degree',value = Degree)
G=igraph::set.vertex.attribute(graph = G,name = 'bigSCale.PAGErank',value = PAGErank)
G=igraph::set.vertex.attribute(graph = G,name = 'bigSCale.Closeness',value = Closeness)
G=igraph::set.vertex.attribute(graph = G,name = 'bigSCale.Betweenness',value = Betweenness)

if (cutoff.p<0.7)
  warning('bigSCale: the cutoff for the correlations seems very low. You should either increase the parameter quantile.p or select clustering=normal (you need to run the whole code again in both options,sorry!). For more information check the quick guide online')

return(list(graph=G,correlations=Df,tot.scores=tot.scores,clusters=mycl,centrality=as.data.frame(cbind(Degree,Betweenness,Closeness,PAGErank)),cutoff.p=cutoff.p,model=model))
}










compute.pseudotime = function (minST){
  
  
  minST=igraph::set_vertex_attr(graph = minST,name = "name", value = c(1:max(igraph::V(minST))))
  
  minST.core=minST
  
  #estimation of largest shortest path
  #random.vertices=sample(c(1:length(V(minST))),100)
  #dummy=distances(graph = minST,v = random.vertices,weights = NA)
  #network.diameter=max(dummy)
    
  # At every cycle I remove the tail nodes with degree = 1
  cycles=round(length(igraph::V(minST.core))/100)
  #cycles=round(network.diameter/5)

  print('Searching for tail clusters ...')
  for (k in 1:cycles)
  {
    tails=which(igraph::degree(minST.core)==1)
    minST.core=igraph::delete_vertices(minST.core, tails)
    core.nodes=igraph::vertex.attributes(minST.core)$name
    minST.tails=igraph::delete_vertices(minST, core.nodes)
    if ((igraph::components(graph = minST.tails)$no)<20)
    {
      if ((igraph::components(graph = minST.tails)$no)==1)
        {stop('BigSCale author needs to fix a problem here, in the computation of the pseudotime')}
      break
    }
  }
  

  
  #plot(minST.tails,vertex.size=2,vertex.label=NA,layout=layout_with_kk(graph = minST.tails,weights = 1/E(minST.tails)$weight))
  #plot(minST.tails,vertex.size=2,vertex.label=NA,layout=layout_with_kk(graph = minST.tails,weights = NA))
       
  graph.comp=igraph::components(minST.tails)
  
  
  end.tails=list()
  count.comp=1
  for (k in 1:graph.comp$no)
    if (graph.comp$csize[k]>10)
    {
      end.tails[[count.comp]]=igraph::vertex.attributes(minST.tails)$name[which(graph.comp$membership==k)]
      count.comp=count.comp+1
    }
  
  print(sprintf('Computing distances of %g tail clusters ...',length(end.tails)))
  out=c()
  for (k in 1:length(end.tails))
  {
    # minST.out<<-minST
    # end.tails.out<<-end.tails
    out[k]=mean(igraph::distances(graph = minST,v = end.tails[[k]],to = unlist(end.tails[setdiff.Vector(1:length(end.tails),k)])))
  }
  
  
  #plot(minST,vertex.size=2,vertex.label=NA,layout=layout_with_kk(graph = minST,weights = NA),mark.groups=end.tails[[1]])
  
  
  # FIX THIS USING MORE THAN ONE NODE AND ONLY COUNTING THE DISTANCES AGINST OUTER CLUSTERS
  startORend=end.tails[[which.max(out)]]
  dist.to.all=Rfast::rowmeans(igraph::distances(graph = minST,v = startORend))
  end.node=startORend[which.max(dist.to.all)]
  pseudotime=igraph::distances(graph = minST,v = end.node)
  
  #removing outliers
  outliers.dn=which(pseudotime<(median(pseudotime)-3*sd(pseudotime)))
  outliers.up=which(pseudotime>(median(pseudotime)+3*sd(pseudotime)))
  lower=min(pseudotime[-c(outliers.up,outliers.dn)])
  upper=max(pseudotime[-c(outliers.up,outliers.dn)])
  pseudotime[outliers.up]=upper
  pseudotime[outliers.dn]=lower
  
  pseudotime=shift.values(x = pseudotime,A = lower,B = upper)
  
  return(pseudotime)
  
}

bigSCale.violin = function (gene.expr,groups,gene.name){
  
  tot.clusters=length(unique(groups))
  df=adjust.expression(gene.expr = gene.expr, groups = groups)

  p=ggplot2::ggplot(df, ggplot2::aes(x=groups, y=log2(gene.expr+1))) + ggplot2::geom_violin(fill = RgoogleMaps::AddAlpha("grey90",0.5), colour = RgoogleMaps::AddAlpha("grey90",0.5), trim=TRUE,scale='width') + ggbeeswarm::geom_quasirandom(mapping = ggplot2::aes(color=groups),varwidth = TRUE) + ggplot2::scale_colour_manual(name="color", values=set.quantitative.palette(tot.clusters)) + ggplot2::theme_bw() + ggplot2::ggtitle(gene.name)+ ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  print(p)
  return(p)
  
}


adjust.expression = function (gene.expr,groups){
  
  unique.groups=unique(groups,fromLast = FALSE)
  result=c()
  for (k in 1:length(unique.groups))
    result[k]=sum(gene.expr[which(groups==unique.groups[k])]>0)/sum(groups==unique.groups[k])
  
  cut.to=max(result)
  print(result)
  print(cut.to)
  
  TOTretained=c()
  
  for (k in 1:length(unique.groups))
    {
    group.pos=which(groups==unique.groups[k])
    order=sort(gene.expr[group.pos],index.return=T,decreasing = TRUE)
    retained=order$ix[1:round(length(order$ix)*cut.to)]
    retained=group.pos[retained]
    TOTretained=c(TOTretained,retained)
    }
  
  gene.expr=gene.expr[TOTretained]
  groups=groups[TOTretained]
  #groups = factor(groups, levels = unique.groups)
  groups = factor(groups)
  return(data.frame(groups=groups,gene.expr =gene.expr))
  
}


bigSCale.barplot = function (ht,clusters,gene.expr,gene.name){
  
 # STA MERDA DI R NON MI LASCIA PLOTTARE DENDROFRAMMA + GENE, RIMANDO A PIU' TARDI, PROBLEMA E´CHE DENDRO E´r BASE OBJECT E DEVE ESSERE CONVERITTO IN GRID FORMAT O CHE SO IO
   
 # setting up and coloring the dendrogram
 #d=as.dendrogram(ht)
 tot.clusters=max(clusters)
 palette=set.quantitative.palette(tot.clusters)
 #d = dendextend::color_branches(d, k = tot.clusters,col = palette)
 #test <- as.grob(plot(d,leaflab = 'none',axes = FALSE))

  
 df=data.frame(X=c(1:length(gene.expr)),Y=gene.expr[ht$order],clusters=as.factor(clusters[ht$order]))
 
 palette=set.quantitative.palette(tot.clusters)
 palette=RgoogleMaps::AddAlpha(palette,0.4)
 temp=unique(clusters[ht$order])
 temp=sort(temp,index.return=T)
 palette.sorted=palette[temp$ix]
 

 
 p=ggpubr::ggbarplot(df, x = "X", y = "Y", fill = "clusters",color = "clusters", palette = palette.sorted,sort.by.groups = FALSE,x.text.angle = 90,width=1)
 
 cluster.order=unique(df$clusters)
 print(cluster.order)
 cluster.ycoord=-max(gene.expr)/40
 
 for ( k in 1:tot.clusters)
    {
    cluster.xcoord=mean(which(df$clusters==cluster.order[k]))
    p = p + ggplot2::annotate("text", x = cluster.xcoord, y = cluster.ycoord, label = sprintf("C%g",k))
    }
 p=ggpubr::ggpar(p,legend='none')
 p = p +  ggplot2::theme(axis.title.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(), axis.line=ggplot2::element_blank(),axis.text.x=ggplot2::element_blank())+ggplot2::ylab('EXPRESSION COUNTS') + ggplot2::ggtitle(gene.name)+ ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
 
 return(p)
}
  
  
  
bigSCale.signature.plot = function (ht,clusters,colData,data.matrix,signatures,size.factors,type){
  
  gc()
  
  if (is.na(ht[1]))
  {
    print('You are using ViewSignatures with ATAC-seq data or with recursive clustering with RNA-seq')
    print('In these cases it works diffrently than normal')
    print('calling ViewSignatures will return a data frame in which the i-th column is the average expression of the i-th signature')
    print('You can pass these colmuns one by one to ViewGeneViolin() or ViewReduced() to plot the signature expression')
    
    data.matrix=Matrix::as.matrix(data.matrix)
    plotting.data=c()
    row.labels=c()
    for (k in 1:length(signatures))
      {
      plotting.data=cbind(plotting.data,shift.values(Rfast::colmeans(data.matrix[signatures[[k]]$GENE_NUM,]),0,1))
      row.labels[k]=sprintf('Signature %g',k)
      }
    
    plotting.data=as.data.frame(plotting.data)
    colnames(plotting.data)=row.labels
    return(plotting.data)
  }
  
  
  data.matrix=Matrix::as.matrix(data.matrix)
  
  # setting up and coloring the dendrogram
  d=as.dendrogram(ht)
  tot.clusters=max(clusters)
  palette=set.quantitative.palette(tot.clusters)
  d = dendextend::color_branches(d, k = tot.clusters,col = palette)  
  
  # sorting the palette to match the oder of the dendrorgam
  temp=unique(clusters[ht$order])
  temp=sort(temp,index.return=T)
  palette.sorted=palette[temp$ix]

  # 1) Plotting data of the signatures
  plotting.data=c()
  row.labels=c()
   for (k in 1:length(signatures))
   {
    plotting.data=cbind(plotting.data,shift.values(Rfast::colmeans(data.matrix[signatures[[k]]$GENE_NUM,]),0,1))
    if (type=='signatures')
      row.labels[k]=sprintf('Signature %g (%g genes)',k, length(signatures[[k]]$GENE_NUM))
    else
      row.labels[k]=sprintf('Markers Level %g (%g genes)',k, length(signatures[[k]]$GENE_NUM))
   }
    
  plotting.data=t(plotting.data)
  
  
  
  # 2) Plotting data of the size factors
  palette.sizefactors=set.graded.palette(size.factors)
  ColSideColors = cbind(palette.sizefactors,palette.sorted[clusters])
  colnames(ColSideColors)=c('TRANSCRIPTOME SIZE','CLUSTERS')
  
  
  # 3) Plotting user custom colData of single cell object
  if (ncol(colData)>0 & length(colData)>0)
  {
    names.user.metadata=colnames(colData)
    for (k in 1:ncol(colData))
      {
        vector=colData[,k]
        if (is.factor(vector))
          {
          palette=set.quantitative.palette(length(unique(vector)))
          ColSideColors = cbind(palette,ColSideColors)      
          }
        else
          {
          palette=set.graded.palette(vector)
          ColSideColors = cbind(palette,ColSideColors)      
          }         
      }
    colnames(ColSideColors)=c(names.user.metadata[length(names.user.metadata):1],'TRANSCRIPTOME SIZE','CLUSTERS')
  }
  

  if (type=='signatures') 
    heatmap3::heatmap3(plotting.data,Colv = d, showRowDendro = T , ColSideColors = ColSideColors,labCol = c(''),labRow  = row.labels,cexRow =1,col = colorRampPalette(c('black','yellow'))(1024),scale='none',useRaster = TRUE)#,ColSideAnn=colData,ColSideFun=function(x) plot(as.matrix(x)),ColSideWidth=)
  else
    heatmap3::heatmap3(plotting.data[c(nrow(plotting.data):1),],Colv = d, Rowv = NA, ColSideColors = ColSideColors,labCol = c(''),labRow  = row.labels[c(nrow(plotting.data):1)],cexRow =1,col = colorRampPalette(c('black','yellow'))(1024),scale='none',useRaster = TRUE)#,ColSideAnn=colData,ColSideFun=function(x) plot(as.matrix(x)),ColSideWidth=)
  
  gc()
}

set.graded.palette = function (x)
  {
  # mask color of outliers
  value.up=mean(x)+3*sd(x)
  value.dn=mean(x)-3*sd(x)
  outliers.up=which(x>value.up)
  outliers.dn=which(x<value.dn)
  x[outliers.up]=value.up
  x[outliers.dn]=value.dn
  # assign palette
  palette=colorRampPalette(c('ivory','black'))(1024)
  x=shift.values(x,0,1)
  x=as.integer(cut(x,c(0:1024)/1024,include.lowest = TRUE))
  return(palette[x])
  }


assign.color = function (x, scale.type=NA){
  

  if (is.na(scale.type))
    scale.type='expression'
   
  if (scale.type=='expression')
      {
      P1=colorRampPalette(c(rgb(255,255,255, maxColorValue=255),rgb(255,255,0, maxColorValue=255)),space='rgb')(10)
      P2=colorRampPalette(c(rgb(255,255,0, maxColorValue=255),rgb(255,153,0, maxColorValue=255)),space='rgb')(507)
      P3=colorRampPalette(c(rgb(255,153,0, maxColorValue=255),rgb(204,0,0, maxColorValue=255)),space='rgb')(507)
      P=c(P1,P2,P3)  
      }
  if (scale.type=='pseudo')
      {
      P=colorRampPalette(c(rgb(255,255,255, maxColorValue=255),rgb(0,0,0, maxColorValue=255)),space='rgb')(1024)
      }
  
  
  limits=NULL
    
  if(is.null(limits)) 
    limits=range(x)
  
  assigned.color=P[findInterval(x,seq(limits[1],limits[2],length.out=length(P)+1), all.inside=TRUE)]
    
  
  
  
  return(assigned.color)
}


# bigSCale.metadata.plot = function (ht,clusters,colData){
#   
#   gc()
#   
#   # setting up and coloring the dendrogram
#   class(colData)
#   d=as.dendrogram(ht)
#   tot.clusters=max(clusters)
#   palette=set.quantitative.palette(tot.clusters)
#   d = dendextend::color_branches(d, k = tot.clusters,col = palette)  
#   
#   # sorting the palette to match the oder of the dendrorgam
#   temp=unique(clusters[ht$order])
#   temp=sort(temp,index.return=T)
#   palette.sorted=palette[temp$ix]
# 
#   
#   if (missing(colData)) # then plotting only the clusters
#     myplot = plot(d,leaflab = 'none',axes = FALSE) %>% colored_bars(colors = palette.sorted[clusters], dend = d, order_clusters_as_data = TRUE, rowLabels = c('CLUSTERS'),y_shift=0)
#   else
#     #if (is.object(colData))  # then plotting also the user cell metadata
#     #  {
#       myplot = plot(d,leaflab = 'none',axes = FALSE) %>% colored_bars(colors = colData, dend = d, order_clusters_as_data = TRUE, rowLabels = colnames(colData)) %>% colored_bars(colors = palette.sorted[clusters], dend = d, order_clusters_as_data = TRUE, rowLabels = c('CLUSTERS'),y_shift=0)
#     #  }
# 
#   gc()
# }


bigSCale.tsne.plot = function (tsne.data,color.by,fig.title,colorbar.title,cluster_label){
  
  gc()
  #setting cell names for hovering
  cell.names=c()
  for (k in 1:nrow(tsne.data))
    cell.names[k]=sprintf('Cell_%g  Cluster_%g',k,cluster_label[k]) 
  
  if (is.factor(color.by))
    {
    p=interactive.plot(x=tsne.data[,1],y=tsne.data[,2],color=color.by,colors = set.quantitative.palette(length(unique(color.by))),text=cell.names)#) %>% layout(title = sprintf("%s",fig.title))
    p=plotly::layout(p,title = sprintf("%s",fig.title))
    print(p)    
    }
  else
  {
   custom.colorscale=list(c(0, 'rgb(225,225,225)'), c(0.01, 'rgb(255,255,0)'),c(0.5, 'rgb(255,153,0)'), c(1, 'rgb(204,0,0)'))
   Expression=color.by
   p=interactive.plot(x=tsne.data[,1],y=tsne.data[,2],marker=list(color=~Expression,colorbar=list(title=colorbar.title),colorscale=custom.colorscale,reversescale =F),text=cell.names)
   p=plotly::layout(p,title = sprintf("%s",fig.title),xaxis = list(zeroline = FALSE),yaxis = list(zeroline = FALSE))#,cauto=F,cmin=0,cmax=5))   #cauto=F,cmin=0,cmax=1,
    print(p)
   }

  gc()
}

set.quantitative.palette = function(tot.el){
  
  gc()
  if (tot.el<=11)
      {
      palette=c('#ffe119', '#4363d8', '#f58231', '#e6beff', '#800000', '#000075', '#a9a9a9', '#000000','#FF0000','#fffac8','#f032e6')
      palette=palette[1:tot.el]
      #https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
      }
    else
    {
      #palette=c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabebe', '#469990', '#e6beff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff', '#000000')
      #if (tot.el>length(palette))
      #  stop('You have more than 22 clusters and I cannot find enough colors for them. Contact the developer at gio.iacono.work@gmail.com to fix this issue')
      #palette=palette[1:tot.el]
      #https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
      color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
      palette=sample(color, tot.el, replace = FALSE)
    }

    
  
    gc()                    
return(palette)
  
}

organize.markers = function (Mscores,Fscores,cutoff=NA,gene.names){

# PART1 : calculating the marker list for all the clusters and all the levels using the given cutoff
gc()
  


if (is.na(cutoff)) cutoff=3
tot.clusters=nrow(Mscores)

Mlist = matrix(vector('list',tot.clusters*(tot.clusters-1)),tot.clusters,(tot.clusters-1))

  

for ( j in 1:tot.clusters)
      {
      # concatenate the DEscores of the kth line
      block.M=c()
      block.F=c()
      for ( k in 1:tot.clusters)
        if (k!=j)
          {
          block.M=cbind(block.M,-Mscores[[j,k]])
          block.F=cbind(block.F,-Fscores[[j,k]])
          }
        
      if (tot.clusters>2)
        {
        result=Rfast::rowsums(block.M>cutoff)
        order.by.M=Rfast::rowOrder(block.M)
        block.M.sorted=block.M # only to initialize the momory space
        block.F.sorted=block.F # only to initialize the momory space
        for (row.ix in 1:nrow(order.by.M))
          {
          block.M.sorted[row.ix,] = block.M[row.ix,order.by.M[row.ix,]]
          block.F.sorted[row.ix,] = block.F[row.ix,order.by.M[row.ix,]]
          }
        }
      else
        {
        result=as.numeric(block.M>cutoff)
        block.M.sorted = block.M
        block.F.sorted = block.F
        }
      
      for (level in 1:(tot.clusters-1))
        {
        tresh.min=tot.clusters-level # So for example, let's have 5 clusters, then level 1 markers must be up-regulated 4 times
        detected.markers=which(result>=tresh.min)
        if (tot.clusters>2)
          temp=data.frame(GENE_NUM=detected.markers,GENE_NAME=gene.names[detected.markers],Z_SCORE=block.M.sorted[detected.markers,level],LOG_2=block.F.sorted[detected.markers,level] )
        else
          temp=data.frame(GENE_NUM=detected.markers,GENE_NAME=gene.names[detected.markers],Z_SCORE=block.M.sorted[detected.markers],LOG_2=block.F.sorted[detected.markers] )
        ix=sort(temp$Z_SCORE, index.return=TRUE,decreasing = TRUE)$ix
        temp=temp[ix,]
        Mlist[[j,level]]=temp
        rm(temp)
        }
    }
gc()
return(Mlist)
  
}


calculate.signatures = function (Mscores,gene.names,cutoff=3){
  
gc()

tot.clusters=nrow(Mscores) 
tot.scores=matrix(0,length(Mscores[[1,2]]),(tot.clusters*(tot.clusters-1)/2))

# assembling the matrix of Mscores
counts=1
for (k in 1:(tot.clusters-1))
  for (h in (k+1):tot.clusters)
    {
      tot.scores[,counts]=Mscores[[k,h]]
      counts=counts+1 
    }


if (is.vector(tot.scores) | ncol(tot.scores)==1)
  {
  signatures=NA
  return(signatures)
  }
  

# removing genes without significant DE
expressed=which(Rfast::rowsums(abs(tot.scores)>cutoff)>0)
if (length(expressed)>15000)
{
  print('Found too many DE genes/features, limiting to top 15000')
  ranking=Rfast::rowsums(abs(tot.scores))
  expressed=order(ranking,decreasing = T)[1:15000]
}
  

tot.scores=tot.scores[expressed,]
gene.names=gene.names[expressed]

scores=Rfast::rowsums(abs(tot.scores))/ncol(tot.scores)


print(sprintf("Clustering %g genes differentially expressed with a cutoff of %.2f...",nrow(tot.scores),cutoff))

# clustering and creating the list of markers
clusters=bigscale.cluster(1-Rfast::cora(t(tot.scores)))$clusters
signatures=list()
for (k in 1:max(clusters))
  {
  dummy=data.frame(GENE_NUM=expressed[which(clusters==k)],GENE_NAME=gene.names[which(clusters==k)],SCORE=scores[which(clusters==k)])
  signatures[[k]]=dummy[order(dummy$SCORE,decreasing = TRUE),]  
  }
  

gc()

return(signatures)
  
}





calculate.marker.scores = function (expr.norm, clusters, N_pct, edges, lib.size, speed.preset,cap.ones=F){
  
  gc()  
  tot.clusters=max(clusters)
  
  # Initializing the 2D list (a matrix list)
  Mscores = matrix(vector('list',tot.clusters*tot.clusters),tot.clusters,tot.clusters)
  Fscores = matrix(vector('list',tot.clusters*tot.clusters),tot.clusters,tot.clusters)
    
 
  to.calculate=tot.clusters*(tot.clusters-1)/2
  computed=0
  # Rolling trought the pairwise cluster comparisons
  
  pb <- progress::progress_bar$new(format = "Running Diff.Expr. [:bar] :current/:total (:percent) eta: :eta", total = tot.clusters*(tot.clusters-1)/2)
  pb$tick(0)

  if (cap.ones==T)
    expr.norm[expr.norm>0]=1
  
  for (j in 1:(tot.clusters-1))
    for (k in (j+1):tot.clusters)
    {
      out=bigscale.DE(expr.norm = expr.norm, N_pct = N_pct, edges = edges, lib.size = lib.size, group1 = which(clusters==j),group2 = which(clusters==k),speed.preset = speed.preset)
      Mscores[[j,k]] = out[,1] # DEscores.final
      Fscores[[j,k]] = out[,2] # Log2(Fold change)
      Mscores[[k,j]] = -out[,1] # DEscores.final
      Fscores[[k,j]] = -out[,2] # Log2(Fold change)
      computed=computed+1
      pb$tick(1)
    }

out=list(Mscores=Mscores,Fscores=Fscores)
gc()
return(out)
  
}


get.log.scores = function (N_pct,p.cutoff=0.01){
  

  tot.el=nrow(N_pct)
  
  # Making N_pct symmetric
  for (k in 1:tot.el)
    for (h in 1:tot.el)
      if (k>h) N_pct[k,h]=N_pct[h,k]
  
  # 2) Finding deregulated
  dereg=N_pct<p.cutoff
  
  # 3) Calculating log_scores and removing infinite values
  log.scores=-log10(N_pct)
  for (k in 1:tot.el)
  {
    max.val=max(log.scores[k,is.finite(log.scores[k,])])
    infinite.el=which(is.infinite(log.scores[k,]))
    log.scores[k,infinite.el]=max.val
    log.scores[infinite.el,k]=max.val
  }
  
  # zeroing the diagonal and setting pval<dereg zone
  diag(log.scores)=0
  log.scores=abs(log.scores*dereg)
  
}



bigscale.DE = function (expr.norm, N_pct, edges, lib.size, group1,group2,speed.preset='slow',plot.graphic=FALSE){
  
#print('Starting DE')  
gc.out=gc()

  


num.genes.initial=nrow(expr.norm)

if (class(expr.norm)=='big.matrix')
  expr.norm=bigmemory::as.matrix(expr.norm[,c(group1,group2)])
else
  expr.norm=as.matrix(expr.norm[,c(group1,group2)])


tot.lib.size=lib.size
lib.size=lib.size[c(group1,group2)]
group1=c(1:length(group1))
group2=c((length(group1)+1):(length(group1)+length(group2)))

#remove the genes with zero reads
#genes.expr =which(Rfast::rowsums(expr.norm[,c(group1,group2)])>0)
genes.expr =which(Rfast::rowsums(expr.norm)>0)
#print(sprintf('I remove %g genes not expressed enough', (nrow(expr.norm)-length(genes.expr))))
expr.norm=expr.norm[genes.expr,]
gc()
num.genes=nrow(expr.norm)
num.samples=ncol(expr.norm)
tot.el=nrow(N_pct)  


# normalize expression data for library size without scaling to the overall average depth
if (speed.preset!='fast')
  expr.norm=expr.norm/mean(tot.lib.size)



# saving the FOLD changes
F.change=log2(Rfast::rowmeans(expr.norm[,group2])/Rfast::rowmeans(expr.norm[,group1]))
dummy=matrix(0,num.genes.initial,1)
dummy[genes.expr]=F.change
F.change=dummy
rm(dummy)

# if ( (length(group1)+length(group2))==1774 )
# {
#   # expr.norm.out<<-expr.norm
#   # group2.out<<-group2
#   # group1.out<<-group1
#   # num.genes.initial.out<<-num.genes.initial
#   F.change.out<<-F.change
#   #print('Fold change of positin 40764')
#   #print(F.change[40764])
# }
#   


# Calculating Wilcoxon
dummy=rep(0,num.genes)
DE.scores.wc=rep(0,num.genes.initial)
info.groups=c(rep(TRUE,length(group1)),rep(FALSE,length(group2)))
input.matrix=expr.norm[,c(group1,group2)]


#print(sprintf('Launching the Wilcoxon, %g VS %g cells',length(group1),length(group2)))
dummy=BioQC::wmwTest(t(input.matrix), info.groups , valType = 'p.two.sided' )
# fixing the zeroes in wilcoxon
zeroes=which(dummy==0)
if (length(zeroes)>0)
  dummy[zeroes]=min(dummy[-zeroes])
#print('Computed Wilcoxon')
dummy[dummy>0.5]=0.5
DE.scores.wc[genes.expr]=-qnorm(dummy)
DE.scores.wc=abs(DE.scores.wc)*sign(F.change)
rm(dummy,input.matrix)

if (speed.preset=='fast')
  return(cbind(DE.scores.wc,F.change))

problem.size=mean(length(group1),length(group2))

if (speed.preset=='normal')
  {
  max.size=2000
  if (problem.size<=200) 
    max.rep =  5
  else if (problem.size>200 & problem.size<=500)
    max.rep =  3
  else if (problem.size>500 & problem.size<=1000)
    max.rep =  2
  else if (problem.size>1000)
   max.rep =  1
  }
  
  
if (speed.preset=='slow')
  {
  max.size=5000
  if (problem.size<=1000) 
    max.rep =  4
  else if (problem.size>1000 & problem.size<=2500)
    max.rep =  3
  else if (problem.size>2500 & problem.size<=4000)
    max.rep =  2
  else if (problem.size>4000)
    max.rep =  1
  } 


# Downsampling according to MAX-SIZE
if (length(group1)>max.size)
{
  group1=sample(group1)
  group1=group1[1:max.size]
}
if (length(group2)>max.size)
{
  group2=sample(group2)
  group2=group2[1:max.size]
} 


# Making N_pct symmetric
for (k in 1:tot.el)
  for (h in 1:tot.el)
    if (k>h) N_pct[k,h]=N_pct[h,k]


# 3) Calculating log_scores and removing infinite values
log.scores=-log10(N_pct)
for (k in 1:tot.el)
{
  max.val=max(log.scores[k,is.finite(log.scores[k,])])
  infinite.el=which(is.infinite(log.scores[k,]))
  log.scores[k,infinite.el]=max.val
  log.scores[infinite.el,k]=max.val
}

# zeroing the diagonal and setting pval<dereg zone
diag(log.scores)=0
for (k in 1:tot.el)
  for (h in 1:tot.el)
    if (k>h) log.scores[k,h]=-log.scores[h,k]

# Vector is a trick to increase speed in the next C++ part
indA.size=1000000
critical.value=max(lib.size)*max(expr.norm)

if (critical.value>indA.size/10)
{
  indA.size=indA.size*50
  warning(sprintf('Critical value very high (%g): Increased memory usage!!',critical.value))
  if (critical.value>indA.size/10)
    stop(sprintf('Critical value way too high (%g): Stopping the analysis!!',critical.value))
}

#print('Proceding to allocate large vector')
vector=c(0:indA.size)/10 # increase if code blocks, It can assign a gene exprssion level up to 10000000
gc()
ind.A=as.integer(cut(vector,edges,include.lowest = TRUE))
rm(vector)
gc()
ind.A=ind.A-1
gc()



# Calculating scores of real DE
#print('Rcpp computing real DE')
results.DE.real=DE_Rcpp(expr.norm[,group1],expr.norm[,group2], log.scores,  ind.A, lib.size[group1], lib.size[group2])
gc()
DE.scores.real=results.DE.real[1:num.genes]
DE.counts.real=results.DE.real[(num.genes+1):length(results.DE.real)]



# Calculating scores of random permutations
#print('Rcpp computing randomly reshuffled DEs')
idx=c(group1,group2)
Total.scores.fake=c()
Total.counts.fake=c()
for (repetitions in 1:max.rep)
  {
  #print(sprintf('Random repetition %d',repetitions))
  idx=sample(idx)
  fake.A=idx[1:length(group1)]
  fake.B=idx[(length(group1)+1):length(idx)]
  results.random.DE=DE_Rcpp(expr.norm[,fake.A],expr.norm[,fake.B], log.scores,  ind.A, lib.size[fake.A], lib.size[fake.B])
  gc()
  Total.scores.fake=c(Total.scores.fake,results.random.DE[1:num.genes])
  Total.counts.fake=c(Total.counts.fake, results.random.DE[(num.genes+1):length(results.random.DE)])
  gc()
  }

  
# Fitting the random DEs: Moving standard deviation DE_scores against DE_counts
#print('Fitting the random DEs')
sa=sort(Total.counts.fake, index.return = TRUE)
sd_width=200
movSD=zoo::rollapply(Total.scores.fake[sa$ix], width = sd_width,FUN=sd, fill = NA)
last.value.y=movSD[(length(movSD)-sd_width/2)] # From now a quick fix to the two ends of movSD
before.last.value.y=movSD[(length(movSD)-sd_width)]
last.value.x=sa$x[(length(movSD)-sd_width/2)]
before.last.value.x=sa$x[(length(movSD)-sd_width)]
estimated.slope=(last.value.y-before.last.value.y)/(last.value.x-before.last.value.x)
estimated.slope=max(0,estimated.slope)
estimated.increase=estimated.slope*(sa$x[length(sa$x)]-last.value.x)
movSD[length(movSD)]=last.value.y+estimated.increase
movSD[1:sd_width/2]=min(movSD[(sd_width/2+1):sd_width])
f=smooth.spline(x=sa$x[!is.na(movSD)],y=movSD[!is.na(movSD)],df = 32)


yy=approx(f$x,f$y,DE.counts.real)$y

  # f.out<<-f
  #yy.out<<-yy
  #DE.counts.real.out<<-DE.counts.real
  # x.out<<-sa$x[!is.na(movSD)]
  # y.out<<-movSD[!is.na(movSD)]


if (any(is.na(yy))) # fixing the NA in the yy, when they coincide with the highest DE.counts.real
  {
    null.fits=which(is.na(yy))
    for (scroll.null.fits in 1:length(null.fits))
    {
      result=unique(DE.counts.real[-null.fits]<DE.counts.real[null.fits[scroll.null.fits]])
      if (length(result)>1) stop('Fix this bug, programmer!')
      if (result=='FALSE') 
        yy[null.fits[scroll.null.fits]]=min(yy[-null.fits])
      else
        yy[null.fits[scroll.null.fits]]=max(yy[-null.fits])
    }
  }
 

# I want at least 5 non zeros cell in one group if the other is full empty
treshold=min( 5*max(length(group1),length(group2)),(length(group1)*length(group2))/2) # values fitted below the 400 comparisons (eg.20cell*20cells) are considered noisy and replace with the lower value of next interval
if (sum(DE.counts.real>treshold)==0)
  error('Clusters too small for DE')
yy[DE.counts.real<=treshold]=min(yy[DE.counts.real>treshold])


# Setting to zero the scores of all gene with less than 5 cells.
#DE.scores.real.out<<-DE.scores.real
#DE.scores.wc.out<<-DE.scores.wc
#yy.out<<-yy

DE.scores=DE.scores.real/yy
## debugging


# saving the DE scores
DE.scores=rep(0,num.genes.initial)
DE.scores[genes.expr]=DE.scores.real/yy
DE.scores=abs(DE.scores)*sign(F.change)

#DE.scores.out<<-DE.scores


if (plot.graphic)
  {
  yy2=approx(f$x,f$y,Total.counts.fake)$y
  df=as.data.frame(t(rbind(DE.counts.real,DE.scores.real,yy,Total.scores.fake,Total.counts.fake,yy2)))
  g1=ggplot2::ggplot(df) + ggplot2::ggtitle('Fitting over random resampling')  + ggplot2::geom_point(ggplot2::aes(x=Total.counts.fake, y=Total.scores.fake,color='gray'),alpha=0.5)  + ggplot2::geom_line(ggplot2::aes(x=Total.counts.fake, y=yy2),colour="black") +
    ggplot2::xlab('DE_counts') + ggplot2::ylab('DE_score')
  print(g1)
  g2=ggplot2::ggplot(df) + ggplot2::ggtitle('DE results')  + ggplot2::geom_point(ggplot2::aes(x=DE.counts.real, y=DE.scores.real,color='gray'),alpha=0.5)  + ggplot2::geom_line(ggplot2::aes(x=DE.counts.real, y=yy),colour="black") + ggplot2::xlab('DE_counts') + ggplot2::ylab('DE_score')
  print(g2)
  }


     

# Merging Wilcoxon and bigscale pvalues to increase DE quality

factor1=max(DE.scores.wc[!is.infinite(DE.scores.wc)])/max(DE.scores[!is.na(DE.scores)])
factor2=min(DE.scores.wc[!is.infinite(DE.scores.wc)])/min(DE.scores[!is.na(DE.scores)])
factor=mean(c(factor1,factor2))
# print(max(DE.scores.wc[!is.infinite(DE.scores.wc)]))
# print(max(DE.scores[!is.na(DE.scores)]))
# print(factor1)
# print(min(DE.scores.wc[!is.infinite(DE.scores.wc)]))
# print(min(DE.scores[!is.na(DE.scores)]))
# print(factor2)
if (is.infinite(factor) | factor<0)
{
  stop('Problems with the factor')
}
#print(sprintf('Factor1 %.2f, factor2 %.2f, average factor %.2f',factor1,factor2,factor))

DE.scores.wc=DE.scores.wc/factor
DE.scores.final=sqrt(DE.scores.wc^2+DE.scores^2)*sign(F.change)

#print(Sys.time()-start)

gc()

return(cbind(DE.scores.final,F.change))

  
}


bigscale.cluster = function (D,plot.clusters=FALSE,cut.depth=NA,method.treshold=0.5,clustering.method='high.granularity',granularity.size,verbose=TRUE){

gc()
  
  if (class(D)=='big.matrix')
    D=bigmemory::as.matrix(D)
  if (class(D)=='dgCMatrix')
    error('Error of the stupid biSCale2 programmer!')
  
ht=hclust(as.dist(D),method='ward.D')
ht$height=round(ht$height,6)


if (is.na(cut.depth) & clustering.method=='high.granularity')
  {
  # Trying a variation of the elbow method to automatically set the cluster number
  result=rep(0,100)
  for (k in 1:100)
    {
    mycl <- cutree(ht, h=max(ht$height)*k/100)
    result[k]=max(mycl)
    }
  movAVG=zoo::rollapply(data = -diff(result),10,mean, fill = NA)

  cut.depth=max(which(movAVG>method.treshold)) #1

  if (is.infinite(cut.depth)) cut.depth=50
  
  }

if (is.na(cut.depth) & clustering.method=='low.granularity')
{

  
  result=c()
  progressive.depth=c(90, 80, 70, 60, 50, 40, 30, 20, 15, 10)
  
  for (k in 1:length(progressive.depth))
  {
    mycl <- cutree(ht, h=max(ht$height)*progressive.depth[k]/100)
    result[k]=max(mycl)
  }
  #print(result)
  result= ( diff(result)>0 )
  
  #print(result)
  
  detected.positions=c()
  for (k in 1:(length(result)-2))
    if (result[k]==TRUE & result[k+1]==TRUE)
      detected.positions=c(detected.positions,k)
      
  
  if (sum(result)==0)
    cut.depth=30
  else  
  {
    if (length(detected.positions)==0) 
      cut.depth=progressive.depth[min(which(result==TRUE))]
    else
      cut.depth=progressive.depth[min(detected.positions)]
  }

  
  
  #print(sprintf('Cut depth before checking granularity: %g percent',cut.depth))
  if (cut.depth>=30 & length(ht$order)<granularity.size) cut.depth=100 # do not cluster at all if it starts fragemnting so high alreay
  if (cut.depth>=40 & length(ht$order)>=granularity.size & length(ht$order)<(granularity.size*2)) cut.depth=100 # do not cluster at all if it starts fragemnting so high alreay
  if (cut.depth>=50 & length(ht$order)>=(granularity.size*2) & length(ht$order)<(granularity.size*3)) cut.depth=100 # do not cluster at all if it starts fragemnting so high alreay
  if (cut.depth>=60 & length(ht$order)>=(granularity.size*3) & length(ht$order)<(granularity.size*4)) cut.depth=100 # do not cluster at all if it starts fragemnting so high alreay
  if (cut.depth>=70 & length(ht$order)>=(granularity.size*4) & length(ht$order)<(granularity.size*5)) cut.depth=100 # do not cluster at all if it starts fragemnting so high alreay
  if (cut.depth>=80 & length(ht$order)>=(granularity.size*5) & length(ht$order)<(granularity.size*6)) cut.depth=100 # do not cluster at all if it starts fragemnting so high alreay
}


if (cut.depth<100)
  mycl = cutree(ht, h=max(ht$height)*cut.depth/100)
else
  mycl = rep(1,length(ht$order)) 

if (verbose)
print(sprintf('Automatically cutting the tree at %g percent, with %g resulting clusters',cut.depth,max(mycl)))

#if (length(unique(mycl))>1 & cut.depth==100) stop('Error in the clustering, wake up!')

# Fixing the ordering of the clusters
vicious.order=unique(mycl[ht$order])
clusters=rep(0,length(mycl))
for (k in 1:max(mycl))
  clusters[which(mycl==k)]=which(vicious.order==k)

if (plot.clusters) # plotting the dendrogram
  {
  #print('We are here')
  d=as.dendrogram(ht)
  tot.clusters=max(clusters)
  palette=set.quantitative.palette(tot.clusters)
  d = dendextend::color_branches(d, k = tot.clusters,col = palette)  
  plot(d)
  }

gc()
return(list(clusters=clusters,ht=ht))

}


compute.distances = function (expr.norm, N_pct , edges, driving.genes , genes.discarded,lib.size,modality='pca',pca.components=25){
  
  
  print(sprintf('Proceeding to calculated cell-cell distances with %s modality',modality))
  
  if (modality=='jaccard')
      {
      print("Calculating Jaccard distances ...")
      D=as.matrix(jaccard_dist_text2vec_04(x = Matrix::Matrix( t(as.matrix(expr.norm[driving.genes,]>0)), sparse = T  )))
      return(D)
      }
 
  if (modality=='pca')
  {
    
    
    print(sprintf('Using %g PCA components for %g genes and %g cells',pca.components,length(driving.genes),ncol(expr.norm)))
    if(max(expr.norm)==1)
      {
      print('Detecting ATAC-seq data...')
      dummy=svd(expr.norm[driving.genes,],0,pca.components)  
      }
      else
      dummy=svd(log10(expr.norm[driving.genes,]+1),0,pca.components)
    
    for (k in 1:ncol(dummy$v))
      dummy$v[,k]=dummy$v[,k]*dummy$d[k]
    print('Computing distance from PCA data...')
    D=dist(dummy$v,method = 'euclidean')
    return(D) # is a distance object
  }  
  
  # normalize expression data for library size without scaling to the overall average depth
  if (class(expr.norm)=='big.matrix')
    expr.driving.norm=bigmemory::as.matrix(expr.norm[driving.genes,])/mean(lib.size)
  else
    expr.driving.norm=as.matrix(expr.norm[driving.genes,])/mean(lib.size)
  gc()
  
  if (modality=='correlation')
    {
    print("Calculating normalized-transformed matrix ...")
    expr.norm.transformed = transform.matrix(expr.driving.norm , 2 )
    gc()
    print("Calculating Pearson correlations ...")
    D=1-Rfast::cora(expr.norm.transformed)
    rm(expr.norm.transformed)
    gc()
    return(D)
  }
  
  
  
  
  if (!hasArg(genes.discarded)) genes.discarded =c()

  
  log.scores=get.log.scores(N_pct)
  
  # Consumes several Gb of memory this step!
  # Vector is a trick to increase speed in the next C++ part

  indA.size=1000000

  critical.value=max(lib.size)*max(expr.driving.norm)
  if (critical.value>indA.size/10)
    {
    indA.size=indA.size*50
      warning(sprintf('Critical value very high (%g): Increased memory usage!!',critical.value))
    if (critical.value>indA.size/10)
      stop(sprintf('Critical value way too high (%g): Stopping the analysis!!',critical.value))
    }
  
  #print('Proceding to allocate large vector')
  vector=c(0:indA.size)/10 # increase if code blocks, It can assign a gene exprssion level up to 10000000
  gc()
  ind.A=as.integer(cut(vector,edges,include.lowest = TRUE))
  rm(vector)
  gc()
  ind.A=ind.A-1
  gc()
  
  

  
  if (length(genes.discarded)>0)
  
    print('Detect that you want to remove genes')
  # vector.weights=ones(num_genes,1)
  # vector_weights(remove_genes)=0;
  # tic; D_ls=SC_dist_MEX_double_rv(total_data_norm,log_scores,ind_A,sum_ex,0,vector_weights); t=toc
  else
    {
    start.time=Sys.time() 
    D=C_compute_distances(expr.driving.norm,log.scores,ind.A,lib.size)
    #print("Time elapsed to calculate distances")
    #print(Sys.time()-start.time)
    }

  #D=Dist(Rfast::squareform(D), method = "euclidean", square = FALSE,vector=FALSE)
  #D=(1-fastCor(Rfast::squareform(D)))
  #diag(D)=0
  D=Rfast::squareform(D)
  
  gc()
  
  return(D)
  

}










# fit.model.v2 = function(expr.data,edges) {
#   
#   gc() 
#   
#   # % Performing down-sampling if ds>0
#   # if (ds>0)
#   #   num_samples=length(total_data(1,:))
#   # disp('Performing downsampling, as you requested');
#   # selected=randi(num_samples,ds,1);
#   # if issparse(total_data)
#   # total_data=full(total_data(:,selected));        
#   # else
#   #   total_data=total_data(:,selected);
#   # end
#   # 
#   # end
#   
#   
#   
#   #Remove genes with 1 or 0 umis/reads in the whole dataset. 
#   genes.low.exp =which(Rfast::rowsums(expr.data)<=1)
#   if (length(genes.low.exp))
#   {
#     print(sprintf('I remove %g genes not expressed enough', length(genes.low.exp)))
#     expr.data=expr.data[-genes.low.exp,]
#   }
#   num.genes=nrow(expr.data)
#   num.samples=ncol(expr.data)
#   
#   
#   # normalize expression data for library size
#   lib.size=colsums(expr.data)
#   expr.norm=expr.data/rep_row(lib.size, nrow(expr.data))*mean(lib.size)
#   expr.norm = transform.matrix(expr.norm , 2 )
#   D=fastCor(expr.norm) #improvable with BLAS, or cor(expr.norm,method = 'pearson')
#   D.sorted=t(apply(D,1,order,decreasing = TRUE)) # decreasing because this is a correlatio matri, high correlation more similar
#   couples=cbind(1:nrow(D),D.sorted[,2]) # taking the second column becase first column is the same cell twice (best correlation = 1 from same cell)
#   couples=unique(t(apply(couples,1,sort)))
#   N=enumerate.v2(expr.data,edges,couples)
# 
# 
#   
#   # Calculating percentiles
#   N_pct=matrix(NA,length(edges)-1,length(edges)-1)
#   for (k in 1:nrow(N))
#     for (j in 1:nrow(N))
#       N_pct[k,j] = sum(N[k,1:j])/sum(N[k,])
#   
#   gc()
#   
#   return(N)
#   
#   
# }

fit.model = function(expr.norm,edges,lib.size,plot.pre.clusters=TRUE){
  
  # Performing down-sampling, model does not require more than 5000 cells
  if (ncol(expr.norm)>5000)
    {
    selected=sample(1:ncol(expr.norm),5000)
    if (class(expr.norm)=='big.matrix')
      expr.norm=bigmemory::as.matrix(expr.norm[,selected])
    else
      expr.norm=as.matrix(expr.norm[,selected])
    }
  else
    if (class(expr.norm)=='big.matrix')
      expr.norm=bigmemory::as.matrix(expr.norm)
    else
      expr.norm=as.matrix(expr.norm)
  
  print(dim(expr.norm))
  gc()

  # if (pipeline=='atac')
  #   cutoff=0
  # else
    cutoff=1
  #Remove genes with 1 or 0 UMIs/reads in the whole dataset.
  genes.low.exp =which(Rfast::rowsums(expr.norm>0)<=cutoff)
  if (length(genes.low.exp))
    {
    print(sprintf('I remove %g genes not expressed enough', length(genes.low.exp)))
    expr.norm=expr.norm[-genes.low.exp,]
    }
  
  num.genes=nrow(expr.norm)
  num.samples=ncol(expr.norm)
  
  
  # if (pipeline=='atac')
  # {
  #   print("Calculating distances for ATAC-seq pre-clustering...")
  #   D=as.matrix(jaccard_dist_text2vec_04(x = Matrix::Matrix(t(expr.norm>0))))
  #   gc()
  #   }  
  # else
  #   {
    print("Calculating normalized-transformed matrix ...")
    expr.norm.transformed = transform.matrix(expr.norm , 2 )
    gc()
    print("Calculating Pearson correlations ...")
    #D=fastCor(expr.norm) #improvable with BLAS, or cor(expr.norm,method = 'pearson')
    D=Rfast::cora(expr.norm.transformed)
    rm(expr.norm.transformed)
    gc()
    
    
  
  D=as.dist(1-D)
  print("Clustering  ...")
  ht=hclust(D,method='ward.D')
  
  
  if (plot.pre.clusters) # plotting the dendrogram
    plot(as.dendrogram(ht))
  
  
  # Adjusting max_group_size according to cell number
  if (num.samples<1250) max.group.size=0.10
  if (num.samples>=1250 & num.samples<=2000) max.group.size=0.08
  if (num.samples>=2000 & num.samples<=3000) max.group.size=0.06
  if (num.samples>=3000 & num.samples<=4000) max.group.size=0.05
  if (num.samples>=4000) max.group.size=0.04
  
  # Calculating optimal cutting depth for the pre-clustering
  print('Calculating optimal cut of dendrogram for pre-clustering')
  MAX_CUT_LEVEL=0.9#0.6
  cut.level=MAX_CUT_LEVEL
  while (TRUE)
    {
    while (TRUE)
      {
      if (cut.level<0) break
      T = cutree(ht, h=max(ht$height)*cut.level)
      if ( max(table(as.factor(T)))<(max.group.size*num.samples) & length(unique(T))<(0.05*num.samples) ) break
      cut.level=cut.level-0.01
      }
    if (cut.level>0) break
    else
      {
      max.group.size=max.group.size+0.01;
      cut.level=MAX_CUT_LEVEL;
      }
  }
  
   #g.dendro = fviz_dend(ht, h=max(ht$height)*cut.level)
   mycl <- cutree(ht, h=max(ht$height)*cut.level)
   
   if (plot.pre.clusters) # plotting the dendrogram
   {
     print('We are here')
     d=as.dendrogram(ht)
     tot.clusters=max(mycl)
     palette=set.quantitative.palette(tot.clusters)
     d = dendextend::color_branches(d, k = tot.clusters,col = palette)
     plot(d)
   }
   
   
   print(sprintf('Pre-clustering: cutting the tree at %.2f %%: %g pre-clusters of median(mean) size %g (%g)',cut.level*100,max(mycl),median(as.numeric(table(mycl))),mean(as.numeric(table(mycl))) ))
   # For loop over the clusters

   pb <- progress::progress_bar$new(format = "Analyzing cells [:bar] :current/:total (:percent) eta: :eta", total = length(mycl))
   
   N=matrix(0,length(edges)-1,length(edges)-1)
   for (i in 1:max(mycl))
      {
       pos=which(mycl==i)
      if (length(pos)>3) # just to aviod bugs
        {
        #print(sprintf("Pre-cluster %g/%g",i,max(mycl)))
        N=N+enumerate(expr.norm[,pos],edges,lib.size)
        }
       pb$tick(length(pos))
      }
     
   
   # Calculating percentiles
   N_pct=matrix(NA,length(edges)-1,length(edges)-1)
   for (k in 1:nrow(N))
     for (j in 1:nrow(N))
       N_pct[k,j] = sum(N[k,j:ncol(N)])/sum(N[k,])
  
  N_pct[is.na(N_pct)]=1 # fix 0/0
  print(sprintf("Computed Numerical Model. Enumerated a total of %g cases",sum(N)))
  gc()
  

  
  return(list(N_pct=N_pct,N=N))
   
   
}



enumerate = function(expr.norm,edges,lib.size) {
  gc()
  num.genes=nrow(expr.norm)
  num.samples=ncol(expr.norm) 
  min.cells = min ( round(num.samples/40) , 10 )

  
  # normalize expression data for library size without scaling to the overall average depth
  expr.norm=expr.norm/mean(lib.size)
  gc()
  
  # Removing lowly expressed genes
  genes.low.exp =which(Rfast::rowsums(expr.norm>0)<=min.cells)
  #print(sprintf('Enumeration of combinations, removing  %g genes not expressed enough', length(genes.low.exp)))
  expr.norm=expr.norm[-genes.low.exp,]
  
  # Initiating variable N
  #print(sprintf('Calculating %g couples over %g cells',(num.samples*num.samples - num.samples)/2,num.samples))
  N=matrix(0,length(edges)-1,length(edges)-1)
  
  # Enumerating  combinations
  for (i in 1:num.samples)
    for (j in i:num.samples)
      if (i!=j)
        {
        factor = mean ( lib.size[c(i,j)] )
        xgroup=cut(expr.norm[,i]*factor,edges,include.lowest = TRUE)
        ygroup=cut(expr.norm[,j]*factor,edges,include.lowest = TRUE)
        xy = table(xgroup, ygroup)
        N=N+xy
      }
  gc()
  return(N)
}

# enumerate.v2 = function(expr.data,edges,positions) {
#   gc()
#   num.genes=nrow(expr.data)
#   num.samples=ncol(expr.data) 
#   min.cells = min ( round(num.samples/40) , 10 )
#   
#   # normalize expression data for library size
#   lib.size=colsums(expr.data)
#   expr.norm=expr.data/rep_row(lib.size, nrow(expr.data))
#   rm(expr.data)
#   
#   # Removing lowly expressed genes
#   genes.low.exp =which(Rfast::rowsums(expr.norm>0)<=min.cells)
#   print(sprintf('Enumeration of combinatons, removing  %g genes not expressed enough', length(genes.low.exp)))
#   expr.norm=expr.norm[-genes.low.exp,] 
#   
#   # Initiating variable N
#   print(sprintf('Calculating %g couples',nrow(positions)))
#   N=matrix(0,length(edges)-1,length(edges)-1)
# 
#   
#   # Enumerating  combinations
#   for (i in 1:nrow(positions))
#       {
#         #print(i)
#         factor = mean ( lib.size[c(positions[i,1],positions[i,2])] )
#         xgroup=cut(expr.norm[,positions[i,1]]*factor,edges,include.lowest = TRUE)
#         ygroup=cut(expr.norm[,positions[i,2]]*factor,edges,include.lowest = TRUE)
#         xy = table(xgroup, ygroup)
#         #print(sum(xy))
#         N=N+xy
#   }
#   gc()
#   return(N)
# }


calculate.ODgenes = function(expr.norm,min_ODscore=2.33,verbose=TRUE,use.exp=c(0,1)) {

  if (ncol(expr.norm)<15000)
    downsample.od.genes = c(1:ncol(expr.norm))
  else
    downsample.od.genes = sample(ncol(expr.norm),15000)
  
  if (class(expr.norm)=='big.matrix')
    expr.norm=bigmemory::as.matrix(expr.norm[,downsample.od.genes])
  else
    expr.norm=as.matrix(expr.norm[,downsample.od.genes])
  
  #start.time <- Sys.time() 
  
    
  num.samples=ncol(expr.norm) 
  num.genes=nrow(expr.norm) 
  min.cells=max( 15,  round(0.002*length(expr.norm[1,]))) 
  skwed.cells=max( 5,  round(0.002*length(expr.norm[1,]))) 
  
  #if (verbose)
    print(sprintf('Analyzing %g cells for ODgenes, min_ODscore=%.2f',ncol(expr.norm),min_ODscore))
  
  
  # Discarding skewed genes
  if (verbose)
    print('Discarding skewed genes')
  
  if (max(expr.norm>1))
      {
      expr.row.sorted=Rfast::rowSort(expr.norm) #MEMORY ALERT with Rfast::sort_mat
      a=Rfast::rowmeans(expr.row.sorted[,(num.samples-skwed.cells):num.samples])
      La=log2(a)
      B=Rfast::rowVars(expr.row.sorted[,(num.samples-skwed.cells):num.samples], suma = NULL, std = TRUE)/a
      rm(expr.row.sorted)
      gc()
      f=smooth.spline(x=La[La>0],y=B[La>0],df = 12) # smoothing sline approximatinf the La/B relationship
      yy=approx(f$x,f$y,La)$y # smoothing spline calculate on each exact value of La
    
      skewed=rep(0,length(yy))
      skewed_genes=which(B/yy>4)
      skewed[skewed_genes]=1
       df=as.data.frame(t(rbind(La,B,yy,skewed)))
      df$skewed=as.factor(df$skewed)
      g1=ggplot2::ggplot(df) + ggplot2::ggtitle('Discarding skewed genes')  + ggplot2::geom_point(ggplot2::aes(x=La, y=B,color=skewed),alpha=0.5)  + ggplot2::scale_color_manual(breaks = c("0", "1"),values=c("gray", "red")) + ggplot2::geom_line(ggplot2::aes(x=La, y=yy),colour="black")
      }
    else
      {
      skewed_genes=c()
      g1=plot(1,1)
      }
 
  
  
  
  
  # Fitting OverDispersed genes
  okay=which( Rfast::rowsums(expr.norm>0)>min.cells )
  if (verbose)
    print(sprintf('Using %g genes detected in at least >%g cells',length(okay),min.cells))

  okay=setdiff(okay,skewed_genes)
  if (verbose)
    print(sprintf('Further reducing to %g geni after discarding skewed genes', length(okay)))

  if (length(okay)<200)
    {
    print('Returning no highly variable genes, too few genes, too much noise!')
    return(NULL)
    }
  
  # STEP1: local fit of expression and standard deviation
  expr.norm=expr.norm[okay,]
  #expr.norm.odgenes<<-expr.norm
  gc()
  a=Rfast::rowmeans(expr.norm)
  La=log2(a)
  B=Rfast::rowVars(expr.norm,std = TRUE)/a
  f=smooth.spline(x=La,y=B,df = 8) # smoothing sline approximatinf the La/B relationship
  yy=approx(f$x,f$y,La)$y # smoothing spline calculate on each exact value of La
  
  df=as.data.frame(t(rbind(La,B,yy)))
  g2=ggplot2::ggplot(df) + ggplot2::ggtitle('STEP1: local fit of expression and standard deviation')  + ggplot2::geom_point(ggplot2::aes(x=La, y=B),colour="brown",alpha=0.5)  + ggplot2::geom_line(ggplot2::aes(x=La, y=yy),colour="black") +
    ggplot2::xlab('Log2(expression)') + ggplot2::ylab('Standard Deviation')
  
  
  
  # STEP2: moving standard deviation of fit1
  B_corr1=B-yy
  sLa=sort(La, index.return = TRUE)
  sd_width=100
  movSD=zoo::rollapply(B_corr1[sLa$ix], width = sd_width,FUN=trim_sd, fill = NA)
  f=smooth.spline(x=sLa$x[!is.na(movSD)],y=movSD[!is.na(movSD)],df = 16)
  yy=approx(f$x,f$y,La)$y
  
  pos=max(which(!is.na(yy[sLa$ix])))
  yy[sLa$ix[pos:length(yy)]]=yy[sLa$ix[pos]]

  
  # movSD2=runSD(B_corr1[sLa$ix],sd_width)
  # f2=smooth.spline(x=sLa$x[!is.na(movSD)],y=movSD2[!is.na(movSD)],df = 16)
  # yy2=approx(f$x,f$y,La)$y
  
  
  ODscore=B_corr1/yy
  od_genes=which(ODscore>min_ODscore)

  
  od_genes=intersect(od_genes,which( La>=quantile(La,use.exp[1]) &  La<=quantile(La,use.exp[2]) ))
      
  ODgenes=rep(0,length(La))
  ODgenes[od_genes]=1
  
  df=as.data.frame(t(rbind(La,B_corr1,yy,ODgenes,ODscore)))
  df$ODgenes=as.factor(df$ODgenes)
  
  g3 = ggplot2::ggplot(df) + ggplot2::ggtitle('Determining overdispersed genes (OGgenes)')  + 
    ggplot2::geom_point(ggplot2::aes(x=La, y=B_corr1,color=ODgenes),alpha=0.5)  + ggplot2::scale_color_manual(breaks = c("0", "1"),values=c("gray", "red")) + ggplot2::geom_line(ggplot2::aes(x=La, y=yy),colour="black") +
  ggplot2::xlab('Log2(expression)') + ggplot2::ylab('Adjusted standard deviation ')
  
  # B_corr1.out=B_corr1
  # yy.out<<-yy
  # La.out<<-La
  
  g4=ggplot2::ggplot(df) + ggplot2::ggtitle('Determining overdispersed genes (OGgenes)')  + ggplot2::geom_point(ggplot2::aes(x=La, y=ODscore,color=ODgenes),alpha=0.5)  + ggplot2::scale_color_manual(breaks = c("0", "1"),values=c("gray", "red"))
  
  # generatnig a data.frame for the export
  dummy=matrix(0,num.genes,2)
  #print(dim(dummy))
  #print(length(okay))
  #print(dim(t(rbind(ODgenes,ODscore))))
  dummy[okay,]=t(rbind(ODgenes,ODscore))
  DFout=as.data.frame(dummy)
  colnames(DFout)=c('ODgenes','ODscore')
  
  
  if (verbose)
    print(sprintf('Determined  %g overdispersed genes',length(od_genes)))
  
  #end.time <- Sys.time()
  #print(end.time - start.time)
  
  gc()
  return(list(DFout,g1,g2,g3,g4))
  
  
}









remove.batches = function(expr.data,samples.f,batches.f) {
# Input expression counts,conditions and batches. Output, expression counts corrected for batch effect.  
  
  gc()
  start.time <- Sys.time()
  
  if (missing(sample.f)) samples.f=rep(0,length(batches.f))
  
  # some pre-processing of the variables (change factors to list of positions)
  num_genes=length(expr.data[,1])
  samples=as.numeric(samples.f)
  batches=as.numeric(batches.f)
  percentiles=0:100/100
  
  samples_idx= vector('list',nlevels(samples.f)) 
  batches_idx = vector('list',nlevels(batches.f)) 
  for (k in 1:nlevels(samples.f)){
    samples_idx[[k]]=which(samples == k) 
  }
  for (k in 1:nlevels(batches.f)){
    batches_idx[[k]]=which(batches == k) 
  }   
    
    
  #set random seed
  

  
  #correction for unbalanced batch es
  unbalanced_correction=0
  
  #normalize expression data
  tot_counts=Rfast::colsums(expr.data)
  expr.norm=expr.data/rep_row(tot_counts, nrow(expr.data))*mean(tot_counts)
  
  
  # Check that each pools contains all the samples (using table command on Factors)
  result=table(data.frame(samples.f,batches.f))
  positions=which(Rfast::colsums(result>0)==max(Rfast::colsums(result>0)))
  if (!length(positions)==length(result[1,])) {
    stop('Unbalanced batch design!') # cambialo normalizzando nella maniera corretta
  }
  
  
  # MAIN CYCLE *****************************************
  
  
  for (condition in 1:nlevels(samples.f))
  {
    
    print(condition)
    indexes= vector('list',nlevels(batches.f))
    all_indexes=c()
    
    for (k in 1:nlevels(batches.f))
    {
    indexes[[k]]=intersect(samples_idx[[condition]],batches_idx[[k]])
    all_indexes=c(all_indexes,indexes[[k]])
    } 
    
    
    
    # removing empty pools
    bad_batches=which(table(batches.f)<8);
    if (length(bad_batches)>0)
      {
      indexes(bad_batches)=c();
      printf('Condition %g) Not considering %g bad batches',condition,length(bad_batches))
      }

    
    # If we have more than one pool with at least 8 cells
    if ( length(indexes) > 1 )
    {
      
      

      
      for (h in 1:num_genes)
        {
        if ( nnz(expr.norm[h,all_indexes]) > 1 )
          {
          vq=matrix(0,length(indexes),length(percentiles))
          tot=matrix(0,length(indexes),1)
          order= vector('list',length(indexes))

          
          for (k in 1:length(indexes)) #test batch by batch
            {
              # Calculating percentiles for each combination condition/batch
              data_in_use=expr.norm[h,indexes[[k]]]
              vq[k,] = quantile(data_in_use,percentiles);
              tot[k] = length(indexes[[k]]);
              dummy=sort(data_in_use, index.return = TRUE)
              
              # Ensuring that the sorting is random for the positions with the same expression (not by first-last)
              nums=unique(dummy$x)
              for (n in 1:length(nums))
                {
                pos=which(dummy$x==nums[n])
                if (length(pos)>1)
                  {
                  values=dummy$ix[pos]
                  dummy$ix[pos]=sample(values)
                  }
                }
              order[[k]]=dummy$ix
          }
          


          
          #if (unbalanced_correction==1 & any(intersect(conditions,bad_conditions)) )
          #  cleaned=sum(vq(good_pools,:).*repmat(tot(good_pools),1,length(vq(1,:))))/sum(tot(good_pools));
          #else
          cleaned = colsums( vq*mio_repvector_col(tot,length(percentiles)) ) / colsums(tot) # did not put Rfast::colsums to gain time
          #end
          
          for (k in 1:length(indexes))
            {
            dummy=indexes[[k]]
            expr.norm[h,dummy[order[[k]]]] = approx(percentiles, cleaned,  seq(0,1,1/(tot[k]-1)) )$y
            }
          } # if
        

          
        } # num_genes

    }
  }# MAIN CYCLE *****************************************
  
  end.time <- Sys.time()
  print(end.time - start.time)
      
  plot(Rfast::colsums(expr.norm),type='h',ylim=c(0,max(Rfast::colsums(expr.norm))))
  expr_batch= ( expr.norm/mean(Rfast::colsums(expr.data)) ) * mio_repvector_row(Rfast::colsums(expr.data),nrow(expr.norm))
  expr_batch=round(expr_batch)  

  
  gc()
  return(expr_batch)


  
}


nnz<-function(x){
  sum(x != 0)
}

mio_repvector_row<-function(x,n){
  matrix(data = rep(x,n), nrow = n, ncol = length(x), byrow = TRUE)
}

mio_repvector_col<-function(x,n){
  matrix(data = rep(x,n), nrow = length(x), ncol = n, byrow = FALSE)
}


trim_sd<-function(x){
  x=x[!is.na(x)]
  
  tmp=quantile(abs(x),c(0.05,0.95))
  out=sd(x[x>tmp[1] & x<tmp[2]])
}




generate.edges<-function(expr.data){
  
  edges=c(0,0.0000001)
  

  # if (pipeline=='atac')
  #   {
  #   step=0.75
  #   for (k in 1:40)
  #     edges[length(edges)+1]=edges[length(edges)]+step
  #   edges[length(edges)+1]=Inf
  #   return(edges)
  #   }
  if (ncol(expr.data)>6000)
  {
  print('Subsetting dataset...')
  ix=sample(c(1:ncol(expr.data)),5000)
  expr.data=expr.data[,ix]    
  }
  print('Creating edges...')
  percentage.exp=(sum(expr.data<10 & expr.data>0)/sum(expr.data>0))
  # Testing how many nonzeros values are below 70, if number is large than we have UMIs like distribution
  if (percentage.exp>0.9)
      { 
      modality='UMIs'
      num.edges=80
      period=10
      print(sprintf('%.1f %% of elements < 10 counts, therefore Using a UMIs compatible binning',percentage.exp*100))
      }#'UMIs'
    else
    { 
      modality='reads'
      num.edges=80
      period=5
      print(sprintf('%.1f %% of elements < 10 counts, therefore Using a reads compatible binning',percentage.exp*100))
    }
  
  if (modality=='UMIs')
      {
      step=0.75
      period.reset=period
      for (k in 3:num.edges)
        {
        edges[k]=edges[k-1]+step
        period=period-1
        if (period==0) 
          {
          period=period.reset
          step=step*2
          }
        }
      }
  else
      {
        count.periods=0
        step=5
        period.reset=period
        for (k in 3:num.edges)
        {
          edges[k]=edges[k-1]+step
          period=period-1
          if (period==0) 
          {
            period=period.reset
            step=step+10
          }
        }
      }
  edges[length(edges)+1]=Inf
  return(edges)
}

transform.matrix<-function(expr.norm,case){

  print('Computing transformed matrix ...')
  print('Converting sparse to full Matrix ...')
  
  expr.norm=as.matrix(expr.norm)
  
  
  if (case==2)# model=2. Log(x+1), 
    expr.norm=log2(expr.norm+1)
  
  if (case==4)# capped to 95% expression
    {
    print('Capping expression gene by gene ...') 
    for (k in 1:nrow(expr.norm))
      expr.norm[k,]=cap.expression(expr.norm[k,])
    }
  gc()
  print('Normalizing expression gene by gene ...')  
  # each row (gene) normalized between [0:1]  
  
  for (k in 1:nrow(expr.norm))
    expr.norm[k,]=shift.values(expr.norm[k,],0,1)
  #t(apply(expr.norm, 1,shift.values,A=0,B=1))
  
  
  gc()
  return(expr.norm)
  
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

cap.expression<-function(expr.vector){
  
  cutoff=quantile(expr.vector[expr.vector>0],0.95)
  expr.vector[expr.vector>cutoff]=cutoff
  return(expr.vector)
}

interactive.plot<-function(x,y,...)
  {
  if (hasArg(x)) 
    {
        if (is.matrix(x))
        {
         X = 1:nrow(x)
         data=as.data.frame(cbind(X,x))
         colnames(data)=c("X","V1","V2")
         p=plotly::plot_ly(data, x = ~X, y = ~V1, type = 'scatter', mode = 'lines+markers') %>% add_trace(y = ~V2, type = 'scatter', mode = 'lines+markers')
        }
        else
        {
        data = data.frame(x, y)
        p=plotly::plot_ly(data, x = ~x, y = ~y, type = 'scatter', mode = 'markers',...) 
        #p <- plot_ly(type = 'scatter',mode='markers',x = ~x, y=~y,...)
        }
    }
  else
    {
    x = 1:length(y)
    data = data.frame(x, y)
    p=plotly::plot_ly(data, x = ~x, y = ~y, type = 'scatter', mode = 'lines+markers') 
    }
  return(p)
  }




bigscale.showLegend<-function(legend=c("Group A","Group B"),lwd=3,cex=1.1,col=c("red","blue"),...) {
  plot(0,xaxt="n",bty="n",yaxt="n",type="n",xlab="",ylab="")
  legend("topleft",legend=legend,lwd=lwd,col=col,bty="n",cex=cex,...)
}


shift.values<-function(x,A,B){
  a=min(x)
  b=max(x)
  x_out=(x-a)/(b-a)*(B-A)+A
}


id.map<-function(gene.list,all.genes){
  pos=rep(0,length(gene.list))
  for (k in 1:length(gene.list))
  {
    dummy=which(all.genes==gene.list[k])
    if (length(dummy)>0)
      pos[k]=dummy[1]
  }
  return(pos)
}



#' bigSCale2.0
#'
#' Compute cell clusters, markers and pseudotime
#'
#' @param sce object of the SingleCellExperiment class. The required elements are \code{counts(sce)} and \code{rownames(sce)}. Optionally, you can also fill \code{colData(sce)} with any annotation such as batch, condition: they will be displayed in the plots. 
#' @param speed.preset regulates the speed vs. accuracy in the computation of the marker and differentially expressed genes.
##' \itemize{
#'   \item {\bold{slow}} { Reccomended for most datasets, provides best marker accuracy but slowest computational time.} 
#'   \item {\bold{normal}} {A balance between marker accuracy and computational time. }
#'   \item {\bold{fast}} {Fastest computational time, if you are in a hurry and you have lots of cell (>15K) you can use this}
#' }
#' @param memory.save enables a series of tricks to reduce RAM memory usage. Aware of one case (in linux) in which this option causes irreversible error.
#' @param clustering setting \code{clustering='recursive'} forces the immediate and accurate detection of cell subtypes. 
#' 
#' @return  An sce object storing the markers, pseudotime, cluster and other results. To access the results you can use several S4 methods liste below. Also check the online quick start tutorial over
#'
#'
#' @examples
#' sce=bigscale(sce)
#'
#' @export
#'   
#'   
#' @seealso    
#' [ViewSignatures()]  

bigscale = function (sce,modality='pca',speed.preset='slow',compute.pseudo=TRUE, clustering='normal'){
  
if ('counts' %in% assayNames(sce))
  {
   print('PASSAGE 1) Setting the bins for the expression data ....')
   sce=preProcess(sce)
  
   print('PASSAGE 2) Storing the Normalized data ....')
   sce = storeNormalized(sce)

   print('PASSAGE 3) Computing the numerical model (can take from a few minutes to 30 mins) ....')
   if (!(modality=='pca' & speed.preset=='fast'))
   sce=setModel(sce)
  }
  
  print('Computing Overdispersed genes ...')
  sce=setODgenes(sce)#favour='high'

 if (modality=='pca')
   sce=setDCT(sce) 
  
 else # bigscale pipeline
  {
  print('Storing the Normalized-Transformed data (needed for some plots) ....')
  sce = storeTransformed(sce)
    
  if (clustering=='normal')
    {
    print('Computing cell to cell distances ...')
    sce=setDistances(sce,modality = modality)
    }
  else
    {
    print('Going for recursive clustering')
    sce=RecursiveClustering(sce)
    }
  sce=storeTsne(sce)
  }
  
  sce=setClusters(sce)
  print('Computing the markers (slowest part) ...')
  sce=computeMarkers(sce,speed.preset=speed.preset)
  print('Organizing the markers ...')
  sce=setMarkers(sce)
 

 return(sce)
 
}










#' bigSCale ATAC 
#' 
#' 
#' Compute cell clusters, markers and pseudotime
#'
#' @param sce object of the SingleCellExperiment class. The required elements are \code{counts(sce)} and \code{rownames(sce)}.  
#' @param memory.save If the code fails due to memory problems switch this option on to try to overcome them.
#' @param fragment ATAC-seq clusters are created with a recurive apporach. The clustersing stops whan the size of the smallest clusters is around \code{fragment=0.1} (ie 10 percent) of the total cell number. If you want to partition more/less then decrease/increase the value.
#' @return  An sce object storing the markers, pseudotime, cluster and other results. To access the results you can use several S4 methods. Check the online quick start tutorial for more info.
#'
#'
#' @examples
#' sce=bigscale.atac(sce)
#'
#' @export


bigscale.atac = function (sce,pca.components=25){
  
if ('counts' %in% assayNames(sce))
  {
  print('PASSAGE 1) Per-processing the dataset')
  sce=preProcess(sce,pipeline='atac')
  sce = storeNormalized(sce,memory.save = FALSE)
  }
  sce=setODgenes(sce,min_ODscore = 2)
  
  normcounts(sce)[normcounts(sce)>0]=1
  sce=setDistances(sce,modality='pca',pca.components=pca.components)
  sce=setClusters(sce)

  print('PASSAGE 7) Computing TSNE ...')
  sce=storeTsne(sce)
  
  print('PASSAGE 10) Computing the markers (slowest part) ...')
  sce=computeMarkers(sce,speed.preset='fast',cap.ones=T)
  
  
  print('PASSAGE 11) Organizing the markers ...')
  sce=setMarkers(sce)
  
  return(sce)
  
}

#' bigSCale ATAC Pseudo  (IN DEVELOPMENT, DO NOT USE)
#'
#' Pseudotime for ATAC-seq data
#'
#' @param sce object of the SingleCellExperiment class. The required elements are \code{counts(sce)} and \code{rownames(sce)}.
#' @return  An sce object storing the markers cluster and other results. Check the online quick start tutorial \url{https://github.com/iaconogi/bigSCale2#atac-seq-data} to learn more. 
#'
#'
#' @examples
#' sce=bigscale.atac(sce)
#'
#' @export


bigscale.atac.pseudo = function (sce, memory.save=FALSE){
  
  if ('counts' %in% assayNames(sce))
  {
    print('PASSAGE 1) Setting the bins for the expression data ....')
    sce=preProcess(sce,pipeline='atac') # so does not create edges
    
    print('PASSAGE 2) Storing the Normalized data ....')
    sce = storeNormalized(sce,memory.save)
  }
  
  else
    sce = storeNormalized(sce,memory.save)
  
  sce=setODgenes(sce,min_ODscore = 4)
  sce=setDistances(sce,modality='jaccard')
  
  sce=storePseudo(sce)
  
  if (memory.save==TRUE)
  { 
    print('PASSAGE 12) Restoring full matrices of normalized counts and transformed counts...')
    sce=restoreData(sce)
  }
  
  return(sce)
  
}

