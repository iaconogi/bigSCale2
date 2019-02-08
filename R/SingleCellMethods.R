#' getEdges
#'
#' Retrives the values used to bin expressoin data-
#'
#' @param sce object of the SingleCellExperiment class.
#' 
#' @return  A numeric vector containing the bin limits
#'
#'
#' @examples
#' edges=getEdges(sce)
#'
#' @export


setGeneric(name="getEdges",
           def=function(object)
           {
             standardGeneric("getEdges")
           }
)

setMethod(f="getEdges",
          signature="SingleCellExperiment",
          definition=function(object)
          {
            return(object@int_metadata$edges)
          }
)

setGeneric(name="setEdges",
           def=function(object)
           {
             standardGeneric("setEdges")
           }
)

setMethod(f="setEdges",
          signature="SingleCellExperiment",
          definition=function(object)
          {
            object@int_metadata$edges=generate.edges(counts(object))
            #validObject(object)
            return(object)
          }
)


# create a method to set the Model
setGeneric(name="setModel",
           def=function(object)
           {
             standardGeneric("setModel")
           }
)

setMethod(f="setModel",
          signature="SingleCellExperiment",
          definition=function(object)
          {
            object@int_metadata$model=fit.model(expr.norm = object@int_metadata$normcounts,edges = getEdges(object),lib.size=sizeFactors(object))
            gc()
            #validObject(object)
            return(object)
          }
)

# create a method to view the Model
setGeneric(name="viewModel",
           def=function(object)
           {
             standardGeneric("viewModel")
           }
)

setMethod(f="viewModel",
          signature="SingleCellExperiment",
          definition=function(object)
          {
            plot_ly(z=object@int_metadata$model) %>% add_surface()
          }
)

# create a method to set overdispersed genes
setGeneric(name="setODgenes",
           def=function(object,...)
           {
             standardGeneric("setODgenes")
           }
)

setMethod(f="setODgenes",
          signature="SingleCellExperiment",
          definition=function(object,...)
          {
            #if ('counts.batch.removed' %in% assayNames(object))
            #    {
                print('Calculating ODgenes using normalized expression counts')
                out=calculate.ODgenes(expr.norm = object@int_metadata$normcounts,...)
            #    }
            # else
            #     {
            #     print('Calculating ODgenes using expression counts NOT batch corrected')
            #     out=calculate.ODgenes(expr.data = counts(object),...)
            #     }  
            object@int_metadata$ODgenes.plots=out[2:5]
            dummy=as.matrix(out[[1]])
            object@int_elementMetadata$ODgenes=dummy[,1]
            object@int_elementMetadata$ODscore=dummy[,2]
            rm(dummy)
            gc()
            #validObject(object)
            return(object)
          }
)



# create a method to view processing overdispersed genes
setGeneric(name="viewODgenes",
           def=function(object)
           {
             standardGeneric("viewODgenes")
           }
)

setMethod(f="viewODgenes",
          signature="SingleCellExperiment",
          definition=function(object)
          {
            print(object@int_metadata$ODgenes.plots)
            #validObject(object)
            return(object)
          }
)


# create a method to remove bacth effect
setGeneric(name="remove.batch.effect",
           def=function(object,batches,conditions)
           {
             standardGeneric("remove.batch.effect")
           }
)

setMethod(f="remove.batch.effect",
          signature="SingleCellExperiment",
          definition=function(object,batches,conditions)
          {
            if (missing(batches)) stop('You must provide the batches if you want to remove them ....')
            assay(object, "counts.batch.removed")==remove.batches(expr.data = counts(object), batches.f = batches, samples.f = conditions)
            #validObject(object)
            return(object)
          }
)



# create a method to set cell to cell distances
setGeneric(name="setDistances",
           def=function(object)
           {
             standardGeneric("setDistances")
           }
)

setMethod(f="setDistances",
          signature="SingleCellExperiment",
          definition=function(object)
          {
            object@int_metadata$D=compute.distances(expr.norm = object@int_metadata$normcounts,N_pct = object@int_metadata$model, edges = getEdges(object), driving.genes = which(object@int_elementMetadata$ODgenes==1),lib.size = sizeFactors(object))
            gc()
            #validObject(object)
            return(object)
          }
)


# create a method to get the distances
setGeneric(name="getDistances",
           def=function(object)
           {
             standardGeneric("getDistances")
           }
)

setMethod(f="getDistances",
          signature="SingleCellExperiment",
          definition=function(object)
          {
            return(object@int_metadata$D)
          }
)

# create a method to set clusters
setGeneric(name="setClusters",
           def=function(object,customClust,...)
           {
             standardGeneric("setClusters")
           }
)

setMethod(f="setClusters",
          signature="SingleCellExperiment",
          definition=function(object,customClust,...)
          {
            if (missing(customClust))
              {
              print("Calculating the clusters")
              out=bigscale.cluster(object@int_metadata$D,...)
              object@int_colData$clusters=out$clusters
              object@int_metadata$htree=out$ht
              }
            else
              {
              print("Setting custom defined clusters") 
              object@int_colData$clusters=customClust
              }
            gc()
            #validObject(object)
            return(object)
          }
)

# create a method to get clusters
setGeneric(name="getClusters",
           def=function(object,customClust)
           {
             standardGeneric("getClusters")
           }
)

setMethod(f="getClusters",
          signature="SingleCellExperiment",
          definition=function(object)
          {
              return(object@int_colData$clusters)
          }
)

# create a method to store reduced dimension Tsne object with Rtsne
setGeneric(name="storeTsne",
           def=function(object,...)
           {
             standardGeneric("storeTsne")
           }
)

setMethod(f="storeTsne",
          signature="SingleCellExperiment",
          definition=function(object,...)
          {
            tsne.data=Rtsne::Rtsne(X=as.dist(object@int_metadata$D),is_distance = TRUE, ...)
            reducedDims(object)$TSNE=tsne.data$Y
            rm(tsne.data)
            gc()
            return(object)
          }
)

# create a method to store reduced dimension Tsne object with Rtsne
setGeneric(name="storeUMAP",
           def=function(object,...)
           {
             standardGeneric("storeUMAP")
           }
)

setMethod(f="storeUMAP",
          signature="SingleCellExperiment",
          definition=function(object,...)
          {
            umap.data=umap::umap(object@int_metadata$D)
            reducedDims(object)$UMAP=umap.data$layout
            rm(umap.data)
            gc()
            return(object)
          }
)


# create a method to view TSNE (+set same colors )
setGeneric(name="viewReduced",
           def=function(object,method ='UMAP',color.by,transform.in='trim')
           {
             standardGeneric("viewReduced")
           }
)

setMethod(f="viewReduced",
          signature="SingleCellExperiment",
          definition=function(object,method ='UMAP',color.by,transform.in='trim')
          {
            #if (missing(transform.in)) transform.in='trim'
            if (color.by[1]=='clusters') # we color by cluster
              bigSCale.tsne.plot(tsne.data = reducedDim(object, method),color.by = as.factor(getClusters(object)),fig.title='CLUSTERS')
            if (is.numeric(color.by) | is.factor(color.by)) # we color by custom used defined condition
              bigSCale.tsne.plot(tsne.data = reducedDim(object, method),color.by = as.factor(color.by),fig.title='Custom Condition')
            if (is.character(color.by[1]) & color.by[1]!='clusters')
              {
              gene.name=color.by
              gene.pos=grep(pattern = sprintf('\\b%s\\b',gene.name),x=rownames(object))
              gene.pos=gene.pos[1]
              print(gene.pos)
              if (transform.in=='log') 
                {
                color.by=log2((object@int_metadata$normcounts[gene.pos,]+1))
                colorbar.title='LOG2(COUNTS)'
                }
              if (transform.in=='trim') 
                {
                color.by=cap.expression((object@int_metadata$normcounts[gene.pos,]))
                colorbar.title='COUNTS'
                }
              bigSCale.tsne.plot(tsne.data = reducedDim(object, method),color.by,fig.title=sprintf('%s expression',gene.name),colorbar.title)
              }
              

          }
)


#' viewSignatures
#'
#' Pots and heatmap with dendrogram, clusters, pseudotime, transcriptome complexity, custom user colData and expression values for markers and signtures.
#'
#' @param sce object of the SingleCellExperiment class.
#' 
#' @return  A numeric vector containing the bin limits
#'
#'
#' @examples
#' viewSignatures(sce) #plots all the signatures of correlated markers
#' viewSignatures(sce,selected.cluster=1) # plots the markers (of all levels) of cluster 1
#' 
#' @export


setGeneric(name="viewSignatures",
           def=function(object,selected.cluster)
           {
             standardGeneric("viewSignatures")
           }
)

setMethod(f="viewSignatures",
          signature="SingleCellExperiment",
          definition=function(object,selected.cluster)
          {
            print('Computing normalized - transformed matrix ... ') 
            data.dummy=transform.matrix(object@int_metadata$normcounts,case = 4)
            if (missing(selected.cluster))
              bigSCale.signature.plot(ht = object@int_metadata$htree, clusters=getClusters(object),colData = as.data.frame(colData(object)),data.matrix = data.dummy, signatures = object@int_metadata$Signatures,size.factors = sizeFactors(object),type='signatures') 
            else
            {
              Mlist=object@int_metadata$Mlist
              signatures=list()
              for (k in 1:ncol(Mlist))
              {
                dummy=Mlist[[selected.cluster,k]]
                signatures[[k]]=dummy$GENE_NUM
              }
              bigSCale.signature.plot(ht = object@int_metadata$htree, clusters=getClusters(object),colData = as.data.frame(colData(object)),data.matrix = data.dummy, signatures = signatures,size.factors = sizeFactors(object),type='levels')  
            }
          }
)


# create a method to store reduced dimension Tsne objct with Rtsne
setGeneric(name="storeNormalized",
           def=function(object,...)
           {
             standardGeneric("storeNormalized")
           }
)

setMethod(f="storeNormalized",
          signature="SingleCellExperiment",
          definition=function(object)
          {
            dummy=counts(object)
            counts(object)=c()
            gc()
            library.size=sizeFactors(object)
            avg.library.size=mean(library.size)
            for (k in 1:ncol(dummy)) dummy[,k]=dummy[,k]/library.size[k]*avg.library.size
            #dummy=dummy/rep_row(library.size, nrow(dummy))*mean(library.size)
            object@int_metadata$normcounts=Matrix::Matrix(dummy)
            rm(dummy)
            gc()
            return(object)
          }
)

# create a method to store matrix normalized and transformed for the signature plots
setGeneric(name="storeTransformed",
           def=function(object,...)
           {
             standardGeneric("storeTransformed")
           }
)

setMethod(f="storeTransformed",
          signature="SingleCellExperiment",
          definition=function(object)
          {
            assay(object,'transformed4')=transform.matrix(normcounts(object),case = 4)
            return(object)
          }
)


# create a method to store reduced dimension Tsne objct with Rtsne
setGeneric(name="computeMarkers",
           def=function(object,...)
           {
             standardGeneric("computeMarkers")
           }
)

setMethod(f="computeMarkers",
          signature="SingleCellExperiment",
          definition=function(object,...)
          {
            out=calculate.marker.scores(expr.norm = object@int_metadata$normcounts, clusters=getClusters(object), N_pct=object@int_metadata$model, edges = getEdges(object), lib.size = sizeFactors(object), ...)
            object@int_metadata$Mscores=out$Mscores
            object@int_metadata$Fscores=out$Fscores
            return(object)
          }
)

# create a method to store reduced dimension Tsne objct with Rtsne
setGeneric(name="setMarkers",
           def=function(object, ...) 
           {
             standardGeneric("setMarkers")
           }
)

setMethod(f="setMarkers",
          signature="SingleCellExperiment",
          definition=function(object, ...)
          {
            Mlist=organize.markers( Mscores = object@int_metadata$Mscores , Fscores = object@int_metadata$Fscores , gene.names = rownames(object), ... )
            object@int_metadata$Signatures=calculate.signatures( Mscores = object@int_metadata$Mscores , ... )
            # retriving numer of markers for each cluster
            n.markers=matrix(0,nrow(Mlist),ncol(Mlist))
            row.names=c()
            col.names=c()
            for (i in 1:nrow(Mlist))
              for (j in 1:ncol(Mlist))
                n.markers[i,j]=nrow(Mlist[[i,j]])
            for (i in 1:nrow(Mlist))
              row.names[i]=sprintf('C%g',i)
            for (j in 1:ncol(Mlist))
              col.names[j]=sprintf('Markers_LV%g',j)
            
            colnames(n.markers)=col.names
            rownames(n.markers)=row.names
            
            object@int_metadata$Mlist.counts=as.data.frame(n.markers)
            
            
            object@int_metadata$Mlist=Mlist
            return(object)
          }
)

# create a method to view TSNE (+set same colors )
setGeneric(name="viewGeneBarPlot",
           def=function(object,gene.list)
           {
             standardGeneric("viewGeneBarPlot")
           }
)

setMethod(f="viewGeneBarPlot",
          signature="SingleCellExperiment",
          definition=function(object,gene.list)
          {
          p=list()
          for ( k in 1:length(gene.list))
            {
            pos=grep(pattern = sprintf('\\b%s\\b',gene.list[k]),x=rownames(object))
            print(pos)
            p[[k]]=bigSCale.barplot(ht = object@int_metadata$htree,clusters = getClusters(object),gene.expr = object@int_metadata$normcounts[pos,],gene.name=gene.list[k])
          }
          if (length(gene.list)>1)
            grid.arrange(grobs=p)
          else
            { 
            plot(p[[1]])
            return(p[[1]])
            }

          }
)

# create a method to view TSNE (+set same colors )
setGeneric(name="viewGeneViolin",
           def=function(object,gene.name,groups)
           {
             standardGeneric("viewGeneViolin")
           }
)

setMethod(f="viewGeneViolin",
          signature="SingleCellExperiment",
          definition=function(object,gene.name,groups)
          {
            p=list()
            for ( k in 1:length(gene.name))
              {
              pos=grep(pattern = sprintf('\\b%s\\b',gene.name[k]),x=rownames(object))
              print(pos)
              }
            if (length(gene.name)>1)
              print('Violin plots expects only one gene at the time. Multiple gene names supperted by viewGeneBarPlot')
            else
             if (missing(groups))
                  bigSCale.violin(gene.expr = normcounts(object)[pos,],groups = getClusters(object),gene.name=gene.name)
               else
                  bigSCale.violin(gene.expr = normcounts(object)[pos,],groups = groups,gene.name=gene.name)
              
            
          }
)


# create a method to store matrix normalized and transformed for the signature plots
setGeneric(name="storePseudo",
           def=function(object,...)
           {
             standardGeneric("storePseudo")
           }
)

setMethod(f="storePseudo",
          signature="SingleCellExperiment",
          definition=function(object)
          {
            gc()
            print('Generating the graph ...') # CRUSHING POINT!!!!
            G=graph_from_adjacency_matrix(adjmatrix = object@int_metadata$D,mode = 'undirected',weighted = TRUE)
            print('Deleting cell-cell distances (not needed anymore, to save memory) ...')
            object@int_metadata$D=c()
            gc()
            print('Computing MST ...')
            object@int_metadata$minST=mst(G, weights = E(G)$weight,algorithm = 'prim')
            print('Computing the layout ...')
            #correcting for outliers in the weights
            #weights.corrected=E(object@int_metadata$minST)$weight
            #outliers=which(weights.corrected>quantile(weights.corrected,0.99))
            #weights.corrected[outliers]=quantile(weights.corrected,0.99)
            #object@int_metadata$minST.layout=layout_with_kk(graph = object@int_metadata$minST,weights = weights.corrected)
            object@int_metadata$minST.layout=layout_with_kk(graph = object@int_metadata$minST,weights = NA)
            print('Computing the pseudotime ...')
            object$pseudotime=compute.pseudotime(object@int_metadata$minST)
            return(object)
          }
)


# create a method to store matrix normalized and transformed for the signature plots
setGeneric(name="ViewPseudo",
           def=function(object,color.by,transform.in)
           {
             standardGeneric("ViewPseudo")
           }
)

setMethod(f="ViewPseudo",
          signature="SingleCellExperiment",
          definition=function(object,color.by,transform.in)
          {
            if (is.character(color.by)) 
                if (color.by=='Clusters' | color.by=='clusters')
                {
                mycl=getClusters(object)  
                my.palette=set.quantitative.palette(max(mycl))
                plot(object@int_metadata$minST,vertex.size=2,vertex.label=NA,vertex.color=my.palette[mycl],layout=object@int_metadata$minST.layout)
                }
                else if (color.by=='Pseudo' | color.by=='pseudo')
                  {
                  to.plot=assign.color(object$pseudotime, scale.type = 'pseudo')     
                  plot(object@int_metadata$minST,vertex.size=2,vertex.label=NA,vertex.color=to.plot, layout=object@int_metadata$minST.layout)
                  }
                  else # then gene name
                  {
                  gene.name=color.by
                  gene.pos=grep(pattern = sprintf('\\b%s\\b',gene.name),x=rownames(object))
                  gene.pos=gene.pos[1]
                  print(gene.pos)
                  if (missing(transform.in)) transform.in='trim'
                  if (transform.in=='log') 
                    to.plot=assign.color(x = log2(normcounts(object)[gene.pos,]+1), scale.type = 'expression')
                  if (transform.in=='trim') 
                    to.plot=assign.color(cap.expression(normcounts(object)[gene.pos,]), scale.type = 'expression')         
                  plot(object@int_metadata$minST,vertex.size=2,vertex.label=NA,vertex.color=to.plot,layout=object@int_metadata$minST.layout)
                  }
            else # then groups
                plot(object@int_metadata$minST,vertex.size=2,vertex.label=NA,vertex.color=color.by,layout=object@int_metadata$minST.layout)
          }
)
