bigscale.recursive.clustering = function (expr.data.norm,model,edges,lib.size,fragment=FALSE) {
  
num.samples=ncol(expr.data.norm)
  
if (fragment==FALSE)
  {
  # Adjusting max_group_size according to cell number
  if (num.samples<5000) dim.cutoff=50
  if (num.samples>=5000 & num.samples<10000) dim.cutoff=100
  if (num.samples>=10000) dim.cutoff=150
  }
else
  dim.cutoff=50  

print(sprintf('Clustering cells down to groups of approximately %g-%g cells',dim.cutoff,dim.cutoff*5))   
#dim.cutoff=ncol(expr.data.norm)*min.group.size  

mycl=rep(1,ncol(expr.data.norm))
tot.recursive=1
#current.cutting=40
unclusterable=rep(0,length(mycl))

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
      ODgenes=calculate.ODgenes(expr.data.norm[,which(mycl==k)],verbose = FALSE)
      dummy=as.matrix(ODgenes[[1]])
      ODgenes=which(dummy[,1]==1)
      #print('Computing distances ...')
      D=compute.distances(expr.norm = expr.data.norm[,which(mycl==k)],N_pct = model,edges = edges,driving.genes = ODgenes,lib.size = lib.size[which(mycl==k)])
      temp.clusters=bigscale.cluster(D,plot.clusters = FALSE,clustering.method = 'low.granularity',granularity.size=dim.cutoff,verbose=FALSE)$clusters #cut.depth=current.cutting,method.treshold = 0.2
      if (max(temp.clusters)>1) 
        action.taken=1
      else
        unclusterable[which(mycl==k)]=1
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


return(mycl)
}


#' Compare gene centralities 
#'
#' Works with any given number of networks (N).  The centralities previously calculated with \code{compute.network()} for N networks are given as input.
#' The script sorts the genes accoring to their change in centrality between the first network (first element of the input list) and the other N-1 networks (the rest of the list)
#'
#' @param centralities List of at least two elements. Each elemnt must be the (\code{data.frames}) of the centralities previously calculated by \code{compute.network()}. 
#' @param names.conditions character of the names of the input networks, same length of the \bold{compute.network()}
#' 
#' @return  A list with a (\code{data.frame}) for each centrality. In each (\code{data.frame}) the genes are ranked for thier change in centrality.
#'
#'
#' @examples
#' out=compute.network(expr.data,gene.names)
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
      ranking=rank(result[,length(centralities)+1])
      dummy=which(ranking>(length(ranking)/2))
      ranking[dummy]=length(ranking)-ranking[dummy]
      result[,length(centralities)+2]=ranking
    }
    else
    {
      result[,3]=(result[,1]-result[,2])
      ranking=rank(result[,3])
      dummy=which(ranking>(length(ranking)/2))
      ranking[dummy]=length(ranking)-ranking[dummy]
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
#' @param quantile.p only the first \eqn{1 - quantile.p} correlations are used to create the edges of the network. If the networ is too sparse(dense) decrease(increase) \eqn{quantile.p}
#' @param speed.preset Used only if  \code{clustering='recursive'} . It regulates the speed vs. accuracy of the Zscores calculations. To have a better network quality it is reccomended to use the default \bold{slow}.
##' \itemize{
#'   \item {\bold{slow}} {Highly reccomended, the best network quality but the slowest computational time.} 
#'   \item {\bold{normal}} {A balance between network quality and computational time. }
#'   \item {\bold{fast}} {Fastest computational time, worste network quality.}
#' }
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

  
compute.network = function (expr.data,gene.names,clustering='recursive',quantile.p=0.998,speed.preset='slow'){

  
  
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
  
  
print('PASSAGE 2) Setting the bins for the expression data ....')
edges=generate.edges(expr.data)
  
print('PASSAGE 3) Storing in the single cell object the Normalized data ....')
avg.library.size=mean(lib.size)
for (k in 1:ncol(expr.data)) expr.data[,k]=expr.data[,k]/lib.size[k]*avg.library.size
expr.data.norm=expr.data
rm(expr.data)
expr.data.norm=Matrix::Matrix(expr.data.norm)
gc()
  
print('PASSAGE 4) Computing the numerical model (can take from a few minutes to 30 mins) ....')
model=fit.model(expr.data.norm,edges,lib.size)
gc()
  
print('PASSAGE 5) Clustering ...')
if (clustering=="direct" | clustering=="recursive")
  #mycl=sample(real.clusters)
  mycl=bigscale.recursive.clustering(expr.data.norm = expr.data.norm,model = model,edges = edges,lib.size = lib.size,fragment=TRUE)
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


o=gc()
print(o)
print('Calculating Pearson ...')
rm(list=setdiff(setdiff(ls(), lsf.str()), c('tot.scores','quantile.p','model','gene.names','tot.scores','mycl')))
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
  
  tot.clusters=max(groups)
  df=adjust.expression(gene.expr = gene.expr, groups = groups)

  p=ggplot2::ggplot(df, ggplot2::aes(x=groups, y=log2(gene.expr+1))) + ggplot2::geom_violin(fill = RgoogleMaps::AddAlpha("grey90",0.5), colour = RgoogleMaps::AddAlpha("grey90",0.5), trim=TRUE,scale='width') + ggbeeswarm::geom_quasirandom(mapping = ggplot2::aes(color=groups),varwidth = TRUE) + ggplot2::scale_colour_manual(name="color", values=set.quantitative.palette(tot.clusters)) + ggplot2::theme_bw() + ggplot2::ggtitle(gene.name)+ ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  print(p)
  return(p)
  
}


adjust.expression = function (gene.expr,groups){
  
  result=c()
  for (k in 1:max(groups))
    result[k]=sum(gene.expr[which(groups==k)]>0)/sum(groups==k)
  
  cut.to=max(result)
  print(result)
  print(cut.to)
  
  TOTretained=c()
  
  for (k in 1:max(groups))
    {
    group.pos=which(groups==k)
    order=sort(gene.expr[group.pos],index.return=T,decreasing = TRUE)
    retained=order$ix[1:round(length(order$ix)*cut.to)]
    retained=group.pos[retained]
    TOTretained=c(TOTretained,retained)
    }
  
  gene.expr=gene.expr[TOTretained]
  groups=groups[TOTretained]
  
  return(data.frame(groups=as.factor(groups),gene.expr =gene.expr))
  
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
  
  
  
bigSCale.signature.plot = function (ht,clusters,colData=NA,data.matrix,signatures,size.factors,type){
  
  gc()
  

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
  if (!is.na(colData) & ncol(colData)>0 & length(colData)>0)
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


bigSCale.tsne.plot = function (tsne.data,color.by,fig.title,colorbar.title){
  
  gc()
  #setting cell names for hovering
  cell.names=c()
  for (k in 1:nrow(tsne.data))
    cell.names[k]=sprintf('Cell_%g',k) 
  
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
  
  if (tot.el<3)
    palette=randomcoloR::distinctColorPalette(tot.el)
  else
    if (tot.el<=12)
      {
      palette=c('#ffe119', '#4363d8', '#f58231', '#e6beff', '#800000', '#000075', '#a9a9a9', '#000000','#FF0000','#fffac8','#f032e6')
      palette=palette[1:tot.el]
      #https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
      }
    else
      palette=randomcoloR::distinctColorPalette(tot.el)

    
  
    gc()                    
return(palette)
  
}

organize.markers = function (Mscores,Fscores,cutoff=NA,gene.names){

# PART1 : calculating the marker list for all the clusters and all the levels using the given cutoff
gc()
  


if (is.na(cutoff)) cutoff=3
tot.clusters=nrow(Mscores)

Mlist = matrix(vector('list',tot.clusters*(tot.clusters-1)),tot.clusters,(tot.clusters-1))

  if (tot.clusters<=2)
    return(0)
  
  else
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
      result=Rfast::rowsums(block.M>cutoff)
      block.M.sorted = Rfast::sort_mat(block.M,by.row = TRUE)
      block.F.sorted = Rfast::sort_mat(block.F,by.row = TRUE)
      
      for (level in 1:(tot.clusters-1))
        {
        tresh.min=tot.clusters-level # So for example, let's have 5 clusters, then level 1 markers must be up-regulated 4 times
        detected.markers=which(result>=tresh.min)
        temp=data.frame(GENE_NUM=detected.markers,GENE_NAME=gene.names[detected.markers],Z_SCORE=block.M.sorted[detected.markers,level],LOG_2=block.F.sorted[detected.markers,level] )
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

# removing genes without significant DE
expressed=which(Rfast::rowsums(abs(tot.scores)>cutoff)>0)



tot.scores=tot.scores[expressed,]
gene.names=gene.names[expressed]
scores=Rfast::rowsums(abs(tot.scores))/ncol(tot.scores)
print(sprintf("Clustering %g genes differentially expressed...",nrow(tot.scores)))

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





calculate.marker.scores = function (expr.norm, clusters, N_pct, edges, lib.size, speed.preset){
  
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



bigscale.DE = function (expr.norm, N_pct, edges, lib.size, group1,group2,speed.preset='slow',plot.graphic=FALSE){
  
#print('Starting DE')  
gc.out=gc()
#print(gc.out)
  

#POSSIBLY FASTER PASSING DIRECTLY THE TWO MATRICES WITHOUT GROUP1/GROUP2
start=Sys.time()


num.genes.initial=nrow(expr.norm)

expr.norm=bigmemory::as.matrix(expr.norm[,c(group1,group2)])
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
expr.norm=expr.norm/mean(lib.size)

# saving the FOLD changes
F.change=log2(Rfast::rowmeans(expr.norm[,group2])/Rfast::rowmeans(expr.norm[,group1]))
dummy=matrix(0,num.genes.initial,1)
dummy[genes.expr]=F.change
F.change=dummy
rm(dummy)

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
  stop(sprintf('Critical value too high (%g): Note from bigSCale author, you have to increase indA.size',critical.value))
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
 

treshold=min(400,(length(group1)*length(group2))/2) # values fitted below the 400 comparisons (eg.20cell*20cells) are considered noisy and replace with the lower value of next interval
if (sum(DE.counts.real>treshold)==0)
  error('Clusters too small for DE')
yy[DE.counts.real<=treshold]=min(yy[DE.counts.real>treshold])

DE.scores=DE.scores.real/yy


## debugging
# sa_inside <<- sa$x
# movSD_inside <<- movSD 
# yy_inside <<- yy 
# DE.scores_inside <<- DE.scores
# f_inside<<-f
# DE.counts.real_inside<<-DE.counts.real

# saving the DE scores
DE.scores=rep(0,num.genes.initial)
DE.scores[genes.expr]=DE.scores.real/yy
DE.scores=abs(DE.scores)*sign(F.change)


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


if (is.infinite(factor) | factor<0)
  stop('Problems with the factor')

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
  
ht=hclust(as.dist(D),method='ward.D')



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


compute.distances = function (expr.norm, N_pct , edges, driving.genes , genes.discarded,lib.size){
  
  #expr.norm=as.matrix(expr.norm)
  gc()
  

  
  if (!hasArg(genes.discarded)) genes.discarded =c()
  
  num.genes=nrow(expr.norm)
  num.samples=ncol(expr.norm)
  tot.el=nrow(N_pct)
  #if issparse(total_data)

  
  # Configure deregulated treshold 
  p.cutoff=0.01
  
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
                 
  # normalize expression data for library size without scaling to the overall average depth
  expr.driving.norm=bigmemory::as.matrix(expr.norm[driving.genes,])/mean(lib.size)
  rm(expr.norm)
  gc()
  # Consumes several Gb of memory this step!
  # Vector is a trick to increase speed in the next C++ part

  indA.size=1000000
  expr.driving.norm.out<<-expr.driving.norm
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
    #expr.norm=as.matrix(expr.norm[,selected])
    expr.norm=bigmemory::as.matrix(expr.norm[,selected])
    }
  else
    expr.norm=bigmemory::as.matrix(expr.norm)
  
  print(dim(expr.norm))
  gc()


  
  #Remove genes with 1 or 0 UMIs/reads in the whole dataset.
  genes.low.exp =which(Rfast::rowsums(expr.norm>0)<=1)
  if (length(genes.low.exp))
    {
    print(sprintf('I remove %g genes not expressed enough', length(genes.low.exp)))
    expr.norm=expr.norm[-genes.low.exp,]
    }
  
  num.genes=nrow(expr.norm)
  num.samples=ncol(expr.norm)
  
  

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
  
  
  # if (plot.pre.clusters) # plotting the dendrogram
  # {
  #   d=as.dendrogram(ht)
  #   plot(d)
  # }
  
  
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
  
  print(sprintf("Computed Numerical Model. Enumerated a total of %g cases",sum(N)))
  gc()
  

  
  return(N_pct)
   
   
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


calculate.ODgenes = function(expr.norm,min_ODscore=2.33,verbose=TRUE,favour='none') {

  #print(dim(expr.norm))
  gc()
  expr.norm=bigmemory::as.matrix(expr.norm)
  gc()

  #start.time <- Sys.time() 
  
  num.samples=ncol(expr.norm) 
  num.genes=nrow(expr.norm) 
  min.cells=max( 15,  round(0.002*length(expr.norm[1,]))) 
  skwed.cells=max( 5,  round(0.002*length(expr.norm[1,]))) 

  
  # Discarding skewed genes
  if (verbose)
    print('Discarding skewed genes')
  expr.row.sorted=Rfast::sort_mat(expr.norm, by.row = TRUE) #MEMORY ALERT with Rfast::sort_mat
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
  
  
  
  
  
  # Fitting OverDispersed genes
  okay=which( Rfast::rowsums(expr.norm>0)>min.cells )
  if (verbose)
    print(sprintf('Using %g genes detected in at least >%g cells',length(okay),min.cells))

  okay=setdiff(okay,skewed_genes)
  if (verbose)
    print(sprintf('Further reducing to %g geni after discarding skewed genes', length(okay)))

  # STEP1: local fit of expression and standard deviation
  expr.norm=expr.norm[okay,]
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
  
  if (favour=='high') od_genes=intersect(od_genes,which(La>quantile(La,0.75)))

  if (favour=='low')  od_genes=intersect(od_genes,which(La<quantile(La,0.25)))
      
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
  
  
  if (case==4)
    expr.norm=bigmemory::as.matrix(expr.norm)
  
  
  print('Computing transformed matrix ...')
  
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
  
  if (case==4)
    {
    print('saving to swap transformed matrix ...')  
    expr.norm=bigmemory::as.big.matrix(expr.norm)
    }
    
  
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

bigscale = function (sce,speed.preset='slow',compute.pseudo=TRUE, memory.save=TRUE){
  
 
 # Generate the edges for the binning
 print('PASSAGE 1) Setting the bins for the expression data ....')
 sce=preProcess(sce)

  # Calculate and store the normalized expression data ()
 # !!!!!!! ADD AUTOMATIC USE OF BATCH.CORRECTED IF PRESENT
 print('PASSAGE 3) Storing the Normalized data ....')
 sce = storeNormalized(sce)
 
 # Compute the empirical model of the noise
 print('PASSAGE 4) Computing the numerical model (can take from a few minutes to 30 mins) ....')
 sce=setModel(sce)
 
 # Remove the bacth effect !! To be automated
 #sce=remove.batch.effect(sce, batches=sce$batches, conditions=sce$conditions)
 
 # # Calculate and store the matrix with transformation = 4 (capped expression) for the signature plots
  print('PASSAGE 5) Storing the Normalized-Transformed data (needed for some plots) ....')
 sce = storeTransformed(sce)
 
 # Compute the overdispersed genes
 print('PASSAGE 5) Computing Overdispersed genes ...')
 sce=setODgenes(sce)#favour='high'
 
 # Compute cell to cell distances
 print('PASSAGE 6) Computing cell to cell distances ...')
 sce=setDistances(sce)
 
 print('PASSAGE 7) Computing TSNE ...')
 sce=storeTsne(sce)
 
 # Cluster the cells
 print('PASSAGE 8) Computing the clusters ...')
 sce=setClusters(sce)
 
  # Store the Pseudotime information and removes distances (not used anymore)
 if (compute.pseudo)
    {
    print('PASSAGE 9) Storing the pseudotime order ...')
    sce=storePseudo(sce)
    }
   
 
# Use the gene Zscore information to organize them in groups of markers
 print('PASSAGE 10) Computing the markers (slowest part) ...')
 sce=computeMarkers(sce,speed.preset=speed.preset)
 
 
 print('PASSAGE 11) Organizing the markers ...')
 sce=setMarkers(sce)
 
 print('PASSAGE 12) Restoring full matrices of normalized counts and transformed counts...')
 sce=restoreData(sce)
 
 
 # #viewStuff
 # ViewPseudo(sce,color.by)
 # viewGeneViolin(object = ,gene.name = ,groups = )
 # viewGeneBarPlot(object = ,gene.list = )
 # viewSignatures
 # viewModel
 # View(sce@int_metadata$Mlist.counts)
 # datatable(sce@int_metadata$Mlist[[5,1]])
}



