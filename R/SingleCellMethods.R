#' restoreData
#'
#' To be run as the last step of the analysis. Retrives data termporarily stored in  virtual memory and stores if safely in the single cell object.
#' The data was put in virtual meomory in the first place to reduce RAM memory usage.
#'
#' @param sce object of the SingleCellExperiment class.
#' 
#' @return  object of the SingleCellExperiment class, with the assays \code{normcounts()} and \code{transcounts()} restored.
#'
#'
#' @examples
#' sce=restoreData(sce)
#'
#' @export


setGeneric(name="restoreData",
           def=function(object)
           {
             standardGeneric("restoreData")
           }
)

setMethod(f="restoreData",
          signature="SingleCellExperiment",
          definition=function(object)
          {
            gc()
            normcounts(object)=Matrix::Matrix(bigmemory::as.matrix(object@int_metadata$expr.norm.big))
            X=object@int_metadata$expr.norm.big
            rm(X)
            gc()
            
            gc()
            assay(object, "transcounts")=Matrix::Matrix(bigmemory::as.matrix(object@int_metadata$transformed.big))
            X=object@int_metadata$transformed.big
            rm(X)
            gc()

            return(object)
          }
)

#' getMarkers
#'
#' Retrives a 2D list containing the markers of different specificity for each cluster.
#' For more information check the online tutorial at www.github.com
#' 
#' @param sce object of the SingleCellExperiment class.
#' 
#' @return  2D list, each row is a cluster, each column is a level.
#'
#' @examples
#' Mlist=getMarkers(sce)
#'
#' @export


setGeneric(name="getMarkers",
           def=function(object)
           {
             standardGeneric("getMarkers")
           }
)

setMethod(f="getMarkers",
          signature="SingleCellExperiment",
          definition=function(object)
          {
            return(object@int_metadata$Mlist)
          }
)

#' getSignatures
#'
#' Retrives the signatures of co-expressed genes.
#' For more information check the online tutorial at www.github.com
#' 
#' @param sce object of the SingleCellExperiment class.
#' 
#' @return  A list of co-expressed genes.
#' 
#' @examples
#' Signatures=getSignatures(sce)
#' 
#' @export



setGeneric(name="getSignatures",
           def=function(object)
           {
             standardGeneric("getSignatures")
           }
)

setMethod(f="getSignatures",
          signature="SingleCellExperiment",
          definition=function(object)
          {
            return(object@int_metadata$Signatures)
          }
)


#' preProcess
#'
#' Pre-porocessing of the dataset. Removes genes with all zero values, sets the size factors for normalization and other interal things you do not need to know about.
#' In the future I will add the possibility to add custom size factors here.
#' 
#' @param sce object of the SingleCellExperiment class.
#' 
#' @return  object of the SingleCellExperiment class, pre-processed.
#' 
#' @examples
#' sce=preProcess(sce)
#' 
#' @export

 
setGeneric(name="preProcess",
           def=function(object,pipeline='rna')
           {
             standardGeneric("preProcess")
           }
)

setMethod(f="preProcess",
          signature="SingleCellExperiment",
          definition=function(object,pipeline='rna')
          {
            
            # Removing zeros rows
            print('Pre-processing) Removing null rows ')
            gene.names=rownames(object)
            exp.genes=which(Matrix::rowSums(counts(object))>0)
            
            if ((nrow(counts(object))-length(exp.genes))>0)
            {
              print(sprintf("Discarding %g genes with all zero values",nrow(counts(object))-length(exp.genes)))
              object=object[exp.genes,]
              #object <- SingleCellExperiment(assays = list(counts = counts(object)[exp.genes,]))
              #rownames(object)=as.matrix(gene.names[exp.genes]) 
              gc()
            }

            object@int_metadata$express.filtered=exp.genes
            
            # Assign the size factors
            print('Setting the size factors ....')
            sizeFactors(object) = Matrix::colSums(counts(object))
            
            if (pipeline=='atac')
              counts(object)[counts(object)>0]=1
            
            if (pipeline=='rna')
              {
              print('Generating the edges ....')
              object@int_metadata$edges=generate.edges(counts(object))
              }

            return(object)
          }
)


#' setModel
#'
#' Computes the numerical model at the core of all bigSCale2 analysis.
#' 
#' @param sce object of the SingleCellExperiment class.
#' 
#' @return  object of the SingleCellExperiment class, with the model stored inside.
#' 
#' @examples
#' sce=setModel(sce)
#' 
#' @export

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
            model.output=fit.model(expr.norm = normcounts(object),edges = object@int_metadata$edges,lib.size=sizeFactors(object))
            object@int_metadata$model=model.output$N_pct
            object@int_metadata$N=model.output$N
            return(object)
          }
)

#' viewModel
#'
#' Plots the numerical model calculated by bigSCale2. The numerical model is the core for all the analysis performed by bigSCale2.
#'
#' @param sce object of the SingleCellExperiment class.
#' 
#' @examples
#' viewModel(sce)
#'
#' @export


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
            p=plotly::plot_ly(z=object@int_metadata$model)
            p=plotly::add_surface(p)
            p
          }
)

#' setODgenes
#'
#' Computes the numerical model at the core of all bigSCale2 analysis.
#' 
#' @param sce object of the SingleCellExperiment class.
#' @param min_ODscore the treshold (Z-score) used to select highly variable genes. Increasing(decreasing) it results in less(more) highly variable genes.
#' @param use.exp A vactor of two elements \code{c(a,b)}: lower and upper quantile of expression for restricting the selection of highly variable genes. 
#' Only the genes whose expression quantile is higher than \code{c(a)} and lower than \code{c(b)} will be used for highyl variable genes.
#' For example, if your clustering is driven by lots of highly expressed genes (Ribosomal, Mitochondrial, ..) you can set something like \code{use.exp=c(0,0.9)} to discard them.
#' @param custom.genes A character vector with a custo set of highly variable genes. Overrides any of the previous arguments
#' 
#' @return  object of the SingleCellExperiment class, with the highly variable genes stored inside.
#' 
#' @examples
#' sce=setODgenes(sce)
#' sce=setODgenes(sce,min_ODscore=2) # I want to use more highly variable geness so I lower the threshold.
#' 
#' @export

setGeneric(name="setODgenes",
           def=function(object,min_ODscore=1,use.exp=c(0,1),custom.genes=NA)
           {
             standardGeneric("setODgenes")
           }
)

setMethod(f="setODgenes",
          signature="SingleCellExperiment",
          definition=function(object,min_ODscore=1,use.exp=c(0,1),custom.genes=NA)
          {
            
            #if ('counts.batch.removed' %in% assayNames(object))
            #    {
            
            if (!is.na(custom.genes[1])){
              pos=which(is.element(rownames(object),custom.genes))
              print(sprintf("Using %g/%g of your custom list of genes",sum(pos>0),length(custom.genes)))
              ODgenes=rep(0,length(rownames(object)))
              ODgenes[pos]=1
              object@int_elementMetadata$ODgenes=ODgenes
              return(object)
            }
            
            rounds=1
            if (ncol(normcounts(object))>30000)
              rounds=round(ncol(normcounts(object))/15000)
            if (rounds>5) rounds=5
            print(sprintf('Proceeding with %g rounds of search for driving genes',rounds))

       
            for (k in 1:rounds) 
            {
              out=calculate.ODgenes(expr.norm = normcounts(object),min_ODscore=min_ODscore,use.exp=use.exp)
              dummy=as.matrix(out[[1]])
              if (k==1)
                {
                ODgenes=dummy[,1]
                ODscores=dummy[,2]
                }
                 else
                {
                ODgenes=cbind(ODgenes,dummy[,1])
                ODscores=cbind(ODscores,dummy[,2])
                }                 
            }
            
            #ODgenes.out<<-ODgenes
            #ODscores.out<<-ODscores
            
            if(rounds>1)
            {
              ODgenes=Rfast::rowsums(ODgenes)
              ODgenes[which(ODgenes>0)]=1
              ODscores=Rfast::rowmeans(ODscores) 
              print(sprintf('Union results in %g driving genes',sum(ODgenes)))
              object@int_elementMetadata$ODgenes=ODgenes
              object@int_elementMetadata$ODscore=ODscores
              return(object)
            }

            

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



#' viewODgenes
#'
#' Visualizes a seris of plots showing how the highly variable genes were selected.
#' 
#' @param sce object of the SingleCellExperiment class.
#' 
#' 
#' @examples
#' viewODgenes(sce)
#' @export


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
            assay(object, "counts.batch.removed")=remove.batches(expr.data = counts(object), batches.f = batches, samples.f = conditions)
            #validObject(object)
            return(object)
          }
)



#' setDistances
#'
#' Computes cell to cell distances using the highly variable genes.
#' 
#' @param sce object of the SingleCellExperiment class.
#' @param modality If \code{classic} uses the default bigSCale algorythm. If \code{alternative}, uses an alternative approach which gives a different clustering/TNSE in case you are not satisfied with the default one.
#' @return  object of the SingleCellExperiment class, with cell to cell distances stored inside (in the virtual memory)
#' 
#' @examples
#' sce=setDistances(sce)
#' @export

setGeneric(name="setDistances",
           def=function(object,modality='pca',pca.components=25)
           {
             standardGeneric("setDistances")
           }
)

setMethod(f="setDistances",
          signature="SingleCellExperiment",
          definition=function(object,modality='pca',pca.components=25)
          {
            object@int_metadata$D=compute.distances(expr.norm = normcounts(object),N_pct = object@int_metadata$model, edges = object@int_metadata$edges, driving.genes = which(object@int_elementMetadata$ODgenes==1),lib.size = sizeFactors(object),modality=modality,pca.components=pca.components)
            gc()
            #validObject(object)
            return(object)
          }
)



#' setDCT
#'
#' computes t-SNE, UMAP and clusters
#' 
#' @param sce object of the SingleCellExperiment class.
#' @param modality If \code{classic} uses the default bigSCale algorythm. If \code{alternative}, uses an alternative approach which gives a different clustering/TNSE in case you are not satisfied with the default one.
#' @return  object of the SingleCellExperiment class, with cell to cell distances stored inside (in the virtual memory)
#' 
#' @examples
#' sce=setDCT(sce)
#' @export

setGeneric(name="setDCT",
           def=function(object,pca.components=25)
           {
             standardGeneric("setDCT")
           }
)

setMethod(f="setDCT",
          signature="SingleCellExperiment",
          definition=function(object,pca.components=25)
          {  
            
            driving.genes=which(object@int_elementMetadata$ODgenes==1)
            print(sprintf('Using %g PCA components for %g genes and %g cells',pca.components,length(driving.genes),ncol(normcounts(object))))
            if(max(normcounts(object))==1)
            {
              print('Detecting ATAC-seq data...')
              dummy=svd(normcounts(object)[driving.genes,],0,pca.components)  
            }
            else
              dummy=svd(log10(normcounts(object)[driving.genes,]+1),0,pca.components)
            
            for (k in 1:ncol(dummy$v))
              dummy$v[,k]=dummy$v[,k]*dummy$d[k]
            
            print('Computing t-SNE and UMAP...')
            tsne=Rtsne::Rtsne(X = dummy$v,normalize = FALSE,pca = FALSE)# transforming to full matrix
            umap=umap::umap(d = dummy$v)
            reducedDims(object)$TSNE=tsne$Y
            reducedDims(object)$UMAP=umap$layout
            gc()
            return(object)
          }  


          
)




#' getDistances
#'
#' Retrives the cell to cell distances.
#' 
#' @param sce object of the SingleCellExperiment class.
#' @return  matrix with all cell to cell distances, recovered from virtual memory.
#' 
#' @examples
#' sce=getDistances(sce)
#' @export


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
            return(as.matrix(object@int_metadata$D))
          }
)


#' setClusters
#'
#' Calculates the clusters of cell types.
#' 
#' @param object object of the SingleCellExperiment class.
#' @param customClust A numeric vector containg your custom cluster assignment, overrides all previous settings. 
#' @param plot.clusters By default \code{plot.clusters=FALSE}. If \code{plot.clusters=TRUE} plots a dendrogram of the clusters while making the analysis.
#' @param k.f Inversely proportinal to the number of clusters
#' @param cut.depth By default not used. It overrides the internal decisions of bigSCale2 and forces it to cut the dendrogram at cut.depth (0-100 percent). Works ONLY if using bigSCale modality.
#' @param classifier New option which allows to cluster the cells according to two or three genes given as input.
#' @param num.classifiers How many markers should be used for each group. Check online tutorial for further help https://github.com/iaconogi/bigSCale2#classifier
#' 
#' @return  SingleCellExperiment object with the clusters stored inside.
#' 
#' @examples
#' sce=setClusters(sce)
#' sce=setClusters(sce,classifier='Cd4','Cd8')
#' @export

setGeneric(name="setClusters",
           def=function(object,customClust=NA,classifier=NA,num.classifiers=NA,k.f=100,...)
           {
             standardGeneric("setClusters")
           }
)

setMethod(f="setClusters",
          signature="SingleCellExperiment",
          definition=function(object,customClust=NA,classifier=NA,num.classifiers=NA,k.f=100,...)
          {
            if (is.na(customClust[1]) & is.na(classifier[1]))
              {
              print("Calculating the clusters")
              
              
              if (is.null(object@int_metadata$D))
                {
                #round(nrow(reducedDim(object, 'UMAP'))/k.num)
                k.num=round(nrow(reducedDim(object, 'UMAP'))/k.f)
                knn.norm = FNN::get.knn(reducedDim(object, 'UMAP'), k = k.num )
                knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index),k.num), to = as.vector(knn.norm$nn.index), weight = 1/(1 + as.vector(knn.norm$nn.dist)))
                knn.norm = igraph::graph_from_data_frame(knn.norm, directed = FALSE)
                knn.norm = igraph::simplify(knn.norm)
                object@int_colData$clusters=igraph::cluster_louvain(knn.norm)$membership
                }
              else
                {
                out=bigscale.cluster(object@int_metadata$D,...)
                object@int_colData$clusters=out$clusters
                object@int_metadata$htree=out$ht
                }
              }
            if (length(customClust)>1)
              {
              print("Setting custom defined clusters") 
              object@int_colData$clusters=customClust
              }
            if (is.na(classifier[1])==0 )
            {
              print(sprintf("Using the classifier ...s")) 
              object@int_colData$clusters=bigscale.classifier(expr.counts=normcounts(object),gene.names=rownames(object),selected.genes=classifier,stop.at=num.classifiers)
            }      
            gc()
            #validObject(object)
            print(sprintf('Determined %g clusters',max(object@int_colData$clusters)))
            return(object)
          }
)

#' getClusters
#'
#' Retrives the clusters of cell types.
#' 
#' @param sce object of the SingleCellExperiment class.
#' 
#' @return  numeric vector with cluster assignemnts for each cell
#' 
#' @examples
#' sce=getClusters(sce)
#' @export

setGeneric(name="getClusters",
           def=function(object)
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

#' storeTsne
#'
#' Calculates and stores the TSNE of your dataset. It uses the package \pkg{Rtsne}
#' 
#' @param sce object of the SingleCellExperiment class.
#' @param ... passed to \code{Rtsne()}
#' 
#' @return   object of the SingleCellExperiment class with TSNE stored inside
#' 
#' @examples
#' sce=storeTsne(sce)
#' @export

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
            gc.out=gc()
            print(gc.out)
            if ('normcounts' %in% assayNames(object))
              tsne.data=Rtsne::Rtsne(X=object@int_metadata$D,is_distance = TRUE,normalize=FALSE, ...)
            else
              tsne.data=Rtsne::Rtsne(X=bigmemory::as.matrix(object@int_metadata$D),is_distance = TRUE,normalize=FALSE, ...)
            
            reducedDims(object)$TSNE=tsne.data$Y
            rm(tsne.data)
            gc()
            return(object)
          }
)


#' storeUMAP
#'
#' Calculates and stores the UMAP of your dataset. It uses the package \pkg{UMAP}
#' 
#' @param sce object of the SingleCellExperiment class.
#' @param ... passed to \code{UMAP()}
#' 
#' @return   object of the SingleCellExperiment class with UMAP stored inside
#' 
#' @examples
#' sce=storeUMAP(sce)
#' @export

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
            umap.config=umap::umap.defaults
            umap.config$input="dist"
            umap.config$verbose=TRUE
            
            if (class(object@int_metadata$D)=='big.matrix')
              umap.data=umap::umap(bigmemory::as.matrix(object@int_metadata$D))  
            else
              umap.data=umap::umap(as.matrix(object@int_metadata$D),config = umap.config)
            reducedDims(object)$UMAP=umap.data$layout
            rm(umap.data)
            gc()
            return(object)
          }
)


#' viewReduced
#'
#' View 2D reduced dimentions. For the moment only supports TSNE.
#'
#' @param sce object of the SingleCellExperiment class.
#' @param color.by \bold{clusters} to color by clusters, or a gene name to color by gene names ( must be the names used in \code{rownames()}, or a factor to color by custom grouping.
#' @param transform.in adjustmet of gene expression for better visualuzation. Default is \bold{trim}, you can use \bold{log} to visualize \emph{log2(counts)+1}
#' @param method can only be set to \bold{TSNE} for the moment
#' 
#' @examples
#' viewReduced(sce,"clusters") # colors with clusters
#' viewReduced(sce,"Arx") # colors with Arx gene expression. The gene name must be same as in rownames()
#' 
#' @export
#' 
setGeneric(name="viewReduced",
           def=function(object,color.by='clusters',transform.in='trim',method ='TSNE')
           {
             standardGeneric("viewReduced")
           }
)

setMethod(f="viewReduced",
          signature="SingleCellExperiment",
          definition=function(object,color.by='clusters',transform.in='trim',method ='TSNE')
          {
            #if (missing(transform.in)) transform.in='trim'
            if (color.by[1]=='clusters') # we color by cluster
              bigSCale.tsne.plot(tsne.data = reducedDim(object, method),color.by = as.factor(getClusters(object)),fig.title='CLUSTERS',cluster_label=getClusters(object))
            if (is.numeric(color.by) | is.factor(color.by)) # we color by custom used defined condition
              bigSCale.tsne.plot(tsne.data = reducedDim(object, method),color.by = color.by,fig.title='Custom Values',colorbar.title='Custom values',cluster_label=getClusters(object))
            if (is.character(color.by[1]) & color.by[1]!='clusters')
              {
              gene.name=color.by
              gene.pos=grep(pattern = sprintf('\\b%s\\b',gene.name),x=rownames(object))
              gene.pos=gene.pos[1]
              print(gene.pos)
              if (transform.in=='log') 
                {
                color.by=log2((normcounts(object)[gene.pos,]+1))
                colorbar.title='LOG2(COUNTS)'
                }
              if (transform.in=='trim') 
                {
                color.by=cap.expression((normcounts(object)[gene.pos,]))
                colorbar.title='COUNTS'
                }
              bigSCale.tsne.plot(tsne.data = reducedDim(object, method),color.by,fig.title=sprintf('%s expression',gene.name),colorbar.title,cluster_label=getClusters(object))
              }
              

          }
)


#' atures
#'
#' Plots and heatmap with dendrogram, clusters, pseudotime, transcriptome complexity, custom user colData and expression values for markers and signtures.
#'
#' @param sce object of the SingleCellExperiment class.
#' 
#' @examples
#' viewSignatures(sce) # plots all the signatures of correlated markers
#' viewSignatures(sce,selected.cluster=1) # plots the markers (of all specificity levels) of cluster 1
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
            if (is.na(object@int_metadata$htree[1]))
              return(bigSCale.signature.plot(ht = object@int_metadata$htree, clusters=getClusters(object),colData = as.data.frame(colData(object)),data.matrix = assay(object,'transcounts'), signatures = object@int_metadata$Signatures,size.factors = sizeFactors(object),type='signatures') )
              
            if (missing(selected.cluster))
              bigSCale.signature.plot(ht = object@int_metadata$htree, clusters=getClusters(object),colData = as.data.frame(colData(object)),data.matrix = assay(object,'transcounts'), signatures = object@int_metadata$Signatures,size.factors = sizeFactors(object),type='signatures') 
            else
            {
              Mlist=object@int_metadata$Mlist
              signatures=list()
              for (k in 1:ncol(Mlist))
              {
                #dummy=Mlist[[selected.cluster,k]]
                signatures[[k]]=Mlist[[selected.cluster,k]]#dummy$GENE_NUM
              }
              bigSCale.signature.plot(ht = object@int_metadata$htree, clusters=getClusters(object),colData = as.data.frame(colData(object)),data.matrix = assay(object,'transcounts'), signatures = signatures,size.factors = sizeFactors(object),type='levels')  
            }
          }
)


#' storeNormalized
#'
#' Calculates and stores the normalized counts
#' 
#' @param sce object of the SingleCellExperiment class.
#' @param memory.save enables a series of tricks to reduce RAM memory usage. Aware of one case (in linux) in which this option causes irreversible error.
#' @return   object of the SingleCellExperiment class with normalized counts saved inside (in the virtual memory)
#' 
#' @examples
#' sce=storeNormalized(sce)
#' @export

setGeneric(name="storeNormalized",
           def=function(object, memory.save=TRUE)
           {
             standardGeneric("storeNormalized")
           }
)

setMethod(f="storeNormalized",
          signature="SingleCellExperiment",
          definition=function(object) #dummy=dummy/Rfast::rep_row(library.size, nrow(dummy))*mean(library.size)
          {
            if ('normcounts' %in% assayNames(object))
              return(object)
              
            dummy=counts(object)
            counts(object)=c()
            gc()
            normcounts(object)=batch.normalize(dummy,sizeFactors(object))
            gc()
            return(object)
          }
)

#' storeTransformed
#'
#' Calculates and stores a matrix of transformed counts used by bigSCale to make some plots.
#' 
#' @param sce object of the SingleCellExperiment class.
#' 
#' @return   object of the SingleCellExperiment class with transformed counts saved inside (in the virtual memory)
#' 
#' @examples
#' sce=storeTransformed(sce)
#' @export

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
          assay(object,'transcounts')=transform.matrix(normcounts(object),case = 4)
          return(object)
          }
)


#' computeMarkers
#'
#' Calculates markers and differentially expressed genes.
#' 
#' @param sce object of the SingleCellExperiment class.
#' @param speed.preset by default \code{speed.preset='slow'}. It regulates the speed vs. accuracy in the computation of the marker and differentially expressed genes.
#' \itemize{
#'   \item {\bold{slow}} { Reccomended for most datasets, provides best marker accuracy but slowest computational time.} 
#'   \item {\bold{normal}} {A balance between marker accuracy and computational time. }
#'   \item {\bold{fast}} {Fastest computational time, if you are in a hurry and you have lots of cell (>15K) you can use this}
#' }
#' 
#' 
#' @return   object of the SingleCellExperiment class with raw data of differential expression safely stored inside in an hidden place.
#' 
#' @examples
#' sce=computeMarkers(sce)
#' @export

setGeneric(name="computeMarkers",
           def=function(object,cap.ones=F,...)
           {
             standardGeneric("computeMarkers")
           }
)

setMethod(f="computeMarkers",
          signature="SingleCellExperiment",
          definition=function(object,cap.ones=F,...)
          {
            if ('normcounts' %in% assayNames(object))
              out=calculate.marker.scores(expr.norm = normcounts(object), clusters=getClusters(object), N_pct=object@int_metadata$model, edges = object@int_metadata$edges, lib.size = sizeFactors(object),cap.ones=cap.ones,...)
            else
              out=calculate.marker.scores(expr.norm = object@int_metadata$expr.norm.big, clusters=getClusters(object), N_pct=object@int_metadata$model, edges = object@int_metadata$edges, lib.size = sizeFactors(object),cap.ones=cap.ones, ...)
            
            object@int_metadata$Mscores=out$Mscores
            object@int_metadata$Fscores=out$Fscores
            return(object)
          }
)

#' setMarkers
#'
#' Organizes the differentially expressed genes in markers of different levels and signatures
#'
#' @param sce object of the SingleCellExperiment class.
#' @param cutoff Z-score cutoff to retain only genes with significant changes of expression. By default \code{cutoff=3}. If you feel you have too many(not enough) markers, decrease its value.
#' 
#' @return   object of the SingleCellExperiment class with markers and signatures safely stored inside.
#' 
#' @examples
#' sce=setMarkers(sce)
#' 
#' @export

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
            object@int_metadata$Signatures=calculate.signatures( Mscores = object@int_metadata$Mscores, gene.names = rownames(object), ... )
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
            
            # 
            Mmatrix=matrix(0,length(rownames(object)),ncol(object@int_metadata$Mscores))
            for (k in 1:ncol(object@int_metadata$Mscores))
            {
              dummy.vec=c()
              for (h in 1:ncol(object@int_metadata$Mscores))
                dummy.vec=cbind(dummy.vec,object@int_metadata$Mscores[[k,h]])
              Mmatrix[,k]=-Rfast::rowmeans(dummy.vec)
            }
              
            
            object@int_metadata$Mmatrix=Mmatrix
            object@int_metadata$Mlist=Mlist
            return(object)
          }
)

#' viewGeneBarPlot
#'
#' View 2D reduced dimentions. For the moment only supports TSNE.
#'
#' @param sce object of the SingleCellExperiment class.
#' @param gene.list a character vector with a list of gene names
#' @examples
#' viewGeneBarPlot(sce,c("Arx","NeuroD1","Aqp4")) # colors with clusters
#' @export
#' 
#' 
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
            pos=pos[1]
            gene.name=rownames(object)[pos]
            p[[k]]=bigSCale.barplot(ht = object@int_metadata$htree,clusters = getClusters(object),gene.expr = normcounts(object)[pos,],gene.name=gene.name)
          }
          if (length(gene.list)>1)
            gridExtra::grid.arrange(grobs=p)
          else
            { 
            plot(p[[1]])
            return(p[[1]])
            }

          }
)

#' viewGeneViolin
#'
#' View violin plots
#'
#' @param sce object of the SingleCellExperiment class.
#' @param gene.name character name of the gene to plot. Must be same names use in \code{rownames()}
#' @param groups optional grouping  \code{factor} to be used in place of clusters
#' @examples
#' viewGeneViolin(sce,"Penk")
#' viewGeneViolin(sce,"Penk",groups=custom.groups) # custom.groups must be a factor variable
#' @export
#' 
#' 
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
          if (length(gene.name)>1)
            stop('Violin plots expects only one gene at the time. Multiple gene names supperted by viewGeneBarPlot')

          pos=grep(pattern = sprintf('\\b%s\\b',gene.name),x=rownames(object))
          print(pos)
          pos=pos[1]
          gene.name=rownames(object)[pos]
          
          if (missing(groups))
              {
              p.out=bigSCale.violin(gene.expr = normcounts(object)[pos,],groups = getClusters(object),gene.name=gene.name)
              groups = getClusters(object)
              gene.expr = normcounts(object)[pos,]
              }
          else
              {
              p.out=bigSCale.violin(gene.expr = normcounts(object)[pos,],groups = groups,gene.name=gene.name)
              gene.expr = normcounts(object)[pos,]
              }
            
          all.groups=sort(unique(groups))       
          avg.exp=c()
            for (k in 1:length(all.groups))
              avg.exp[k]=mean(gene.expr[which(groups==all.groups[k])])
          return(avg.exp)  
          }
)

#' storePseudo
#'
#' Computes and stores the pseudotime
#'
#' @param sce object of the SingleCellExperiment class.
#' 
#' @return   object of the SingleCellExperiment class with pasudotime data stored inside.
#' 
#' @examples
#' sce=storePseudo(sce)
#' @export


setGeneric(name="storePseudo",
           def=function(object)
           {
             standardGeneric("storePseudo")
           }
)

setMethod(f="storePseudo",
          signature="SingleCellExperiment",
          definition=function(object)
          {
            
            if (class(object@int_metadata$D)=='dgCMatrix')
              error('Error of the stupid biSCale2 programmer!')
            print('Generating the graph ...') # CRUSHING POINT!!!!
            G=igraph::graph_from_adjacency_matrix(adjmatrix = as.matrix(object@int_metadata$D),mode = 'undirected',weighted = TRUE)
            gc()

            
            # gc()
            # 
            # ix=sample(ncol(D)^2,max(2000000,ncol(D)^2))
            # cutoff=quantile(D[ix],0.25)
            # print(sprintf('PSEUDOTIME: Using %.2f as cutoff for cleaning the distances',cutoff))
            # 
            # for (k in 1:nrow(D))
            #   {
            #   ix=which(D[k,]>cutoff)
            #   if ( length(ix)<(nrow(D)-10) ) # at leats 10 close neighbours, otherwise leave all numbers
            #       D[k,ix]=0
            #   }
            # gc()
            # D=Matrix::Matrix(D)
            # gc()
            
            #G.out<<-G
            print('Computing MST ...')
            object@int_metadata$minST=igraph::mst(G, weights = igraph::E(G)$weight,algorithm = 'prim')
            print('Computing the layout ...')

            #correcting for outliers in the weights
            weights.corrected=igraph::E(object@int_metadata$minST)$weight
            outliers=which(weights.corrected>quantile(weights.corrected,0.99))
            weights.corrected[outliers]=quantile(weights.corrected,0.99)
            
            
            layout=igraph::layout_with_fr(graph = object@int_metadata$minST,weights = 1/weights.corrected,niter = 50000)
            layout=igraph::layout_with_kk(graph = object@int_metadata$minST,weights = NA,coords = layout)
            object@int_metadata$minST.layout=layout

            #object@int_metadata$minST.layout=igraph::layout_with_kk(graph = object@int_metadata$minST,weights = NA)
            print('Computing the pseudotime ...')
            object$pseudotime=compute.pseudotime(object@int_metadata$minST)
            
            gc.out=gc()
            print(gc.out)
            return(object)
          }
)


#' ViewPseudo
#'
#' View the Pseudotime and colors the cell with different options
#'
#' @param sce object of the SingleCellExperiment class.
#' @param color.by \bold{clusters} to color cell by the cluster of belonging. \bold{pseudo} to color cell by the pseudotime. A \bold{gene name} to color cell by its expression.
#' @param transform.in adjustmet of gene expression for better visualuzation. Default is \bold{trim}, you can use \bold{log} to visualize \emph{log2(counts)+1}
#' 
#' 
#' @examples
#' ViewPseudo(sce,'clusters')
#' ViewPseudo(sce,'pseudo')
#' ViewPseudo(sce,'Olig1')
#' @export


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




#' RecursiveClustering
#'
#' Clusters in an accurate recursive way to force immediate identification of cell-subtypes
#' 
#' @param sce object of the SingleCellExperiment class.
#' @param fragment describe fragment
#' 
#' @return  object of the SingleCellExperiment class, with cell-cell distances stored inside
#' 
#' @examples
#' sce=RecursiveClustering(sce)
#' 
#' @export

setGeneric(name="RecursiveClustering",
           def=function(object,modality='bigscale',fragment=FALSE)
           {
             standardGeneric("RecursiveClustering")
           }
)

setMethod(f="RecursiveClustering",
          signature="SingleCellExperiment",
          definition=function(object,modality='bigscale',fragment=FALSE)
          {
            if ('normcounts' %in% assayNames(object))
              {
              out=bigscale.recursive.clustering(expr.data.norm = normcounts(object),model= object@int_metadata$model,edges = object@int_metadata$edges,lib.size = sizeFactors(object),fragment=fragment,create.distances=TRUE,modality=modality)
              object@int_metadata$D=out$D
              }
            else
              {
              out=bigscale.recursive.clustering(expr.data.norm = object@int_metadata$expr.norm.big,model= object@int_metadata$model,edges = object@int_metadata$edges,lib.size = sizeFactors(object),fragment=fragment,create.distances=TRUE,modality=modality)
              object@int_metadata$D=bigmemory::as.big.matrix(out$D)
              }
            gc()
            object@int_colData$clusters=out$mycl
            object@int_metadata$htree=NA
            return(object)
          }
)



#' Differential Expression
#'
#' Performs DE over two groups of cells
#'
#' @param object object of the SingleCellExperiment class.
#' @param group1 a numeric vector with the indices of cells of group1
#' @param group2 a numeric vector with the indices of cells of group2
#' @param speed.preset by default \code{speed.preset='slow'}. It regulates the speed vs. accuracy in the computation of the marker and differentially expressed genes.
#' \itemize{
#'   \item {\bold{slow}} { Reccomended for most datasets, provides best marker accuracy but slowest computational time.} 
#'   \item {\bold{normal}} {A balance between marker accuracy and computational time. }
#'   \item {\bold{fast}} {Fastest computational time, if you are in a hurry and you have lots of cell (>15K) you can use this}
#' }
#' 
#' @examples
#' DE=bigscale.DE(sce,group1=c(1:100),group1=c(101:200))
#' @export
#' 
#' 


setGeneric(name="bigscaleDE",
           def=function(object,group1,group2,speed.preset='slow',cap.ones=FALSE)
           {
             standardGeneric("bigscaleDE")
           }
)

setMethod(f="bigscaleDE",
          signature="SingleCellExperiment",
          definition=function(object,group1,group2,speed.preset='slow',cap.ones=FALSE)
          {
            expr.norm = normcounts(object)
            
            if (cap.ones==T)
              expr.norm[expr.norm>0]=1
            
            out=bigscale.DE(expr.norm, N_pct = object@int_metadata$model, edges = object@int_metadata$edges, lib.size = sizeFactors(object), group1 = group1,group2 = group2,speed.preset = speed.preset)
            out=as.data.frame(out)
            gene.names=rownames(object)
            if(length(unique(gene.names)) < length(gene.names))
            {
              print('You have some duplicated gene names. I will append a number to every gene name')
              for (k in 1:length(gene.names))
                gene.names[k]=paste(gene.names[k],k)
            }
            rownames(out)=gene.names
            colnames(out)=c('Z-score','Fold-change')
            return(out)
          }
)






setGeneric(name="ATACimpute",
           def=function(object,tot.cells=5)
           {
             standardGeneric("ATACimpute")
           }
)

setMethod(f="ATACimpute",
          signature="SingleCellExperiment",
          definition=function(object,tot.cells=5)
          {
            D=as.matrix(jaccard_dist_text2vec_04(x = Matrix::t(counts(object)>0)))
            ix=matrix(0,nrow(D),tot.cells)
            for (k in 1:nrow(D))
              ix[k,]=order(D[k,])[1:tot.cells]
            rm(D)
            gc()
            
            full.counts=as.matrix(counts(object)>0)
            imputed.counts=matrix(0,nrow(full.counts),ncol(full.counts))
            gc()
            for (k in 1:nrow(ix))
              imputed.counts[,k]=Rfast::rowsums(full.counts[,ix[k,]])

            dummy=Matrix::colSums(counts(sce)>0)
            avg.before=mean(dummy)
            sd.before=sd(dummy)
            dummy=Rfast::colsums(imputed.counts>0)
            avg.after=mean(dummy)
            sd.after=sd(dummy)
            
            
            print(sprintf('Before Imputation: %g average open sites per cell, %.2f standard deviation over mean',avg.before,sd.before/avg.before*100))
            print(sprintf('After Imputation: %g average open sites per cell, %.2f standard deviation over mean',avg.after,sd.after/avg.after*100))
            counts(object)=Matrix::Matrix(imputed.counts>0)
            return(object)
          }
)
