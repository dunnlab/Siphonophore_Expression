  QCplot <- function(object) {

	#plot of counts. Note that these are raw, not normalized counts. 
  Dmds<-tibble(
			count=colSums(object@edgeR$counts), 
			library_id=object@library_id, 
			individual=object@individual, 
			treatment=object@treatment
		) %>% dplyr::arrange( treatment ) %>% tibble::add_column(species=rep(object@species, length(object@treatment)))
  
	total <- ggplot( Dmds, aes( factor( library_id,levels=unique( library_id ) ) ,count ) ) + 
		geom_bar( stat = "identity", aes( fill=treatment ) ) + 
		xlab( "Samples" ) + 
		ylab( "Total read counts" ) + 
		theme( axis.text.x=element_text( angle=40, hjust=1 )) + labs(title=Dmds$species[[1]])

	# plot of 0 reads
	ZeroReads <- tibble(
			count=colSums(object@edgeR$counts==0), 
			library_id=object@library_id, 
			individual=object@individual, 
			treatment=object@treatment
		) %>%
		dplyr::arrange( treatment ) %>%
	  tibble::add_column(species=rep(object@species, length(object@treatment)))

	zero <- ggplot(ZeroReads, aes(factor(library_id, levels=unique(library_id)),count)) + 
		geom_bar(stat = "identity", aes(fill=treatment)) + 
		xlab("Samples") + 
		ylab("Total 0 counts") + 
		theme(axis.text.x=element_text(angle=40, hjust=1)) + labs(title=ZeroReads$species[[1]])

	#plot of rRNA
	rRNA <- tibble(
			count=object@rRNA, 
			library_id=object@library_id, 
			individual=object@individual, 
			treatment=object@treatment
		) %>%
	  dplyr::arrange( treatment ) %>%
	  tibble::add_column(species=rep(object@species, length(object@treatment)))

	rRNA <- ggplot( rRNA, aes( factor( library_id, levels=unique( library_id ) ), count ) ) + 
		geom_bar( stat = "identity", aes( fill=treatment ) ) + 
		xlab( "Samples" ) + 
		ylab( "rRNA counts" ) + 
		theme( axis.text.x=element_text( angle=40, hjust=1 ) ) + labs(title=rRNA$species[[1]])
	
	df <- 
		as.data.frame(log2(object@edgeR$counts + 1), stringsAsFactors=FALSE) %>% 
		tidyr::gather( everything(), key="Samples", value="value" ) %>%
		dplyr::left_join( 
			Dmds %>% 
				dplyr::select( library_id, individual, treatment, species ) %>% 
				dplyr::mutate( library_id=as.character( library_id ) ), 
			by=c( "Samples" = "library_id" ) 
		) %>%
		dplyr::arrange( treatment ) %>%
		dplyr::mutate( Samples = as.factor( Samples ) )

	density <- ggplot(df, aes(x = value, colour = Samples, fill=Samples)) + 
		ylim(c(0, 0.25)) + 
		geom_density(alpha = 0.05, size = 0.5) + 
		facet_wrap(~treatment) + 
		theme(legend.position = "top") + 
		xlab(expression(log[2](count + 1))) 
	
	density_all <- ggplot(df, aes(x = value, colour = Samples, fill=Samples)) + 
	  ylim(c(0, 0.25)) + geom_density(alpha = 0.05, size = 0.5) + theme(legend.position = "top") + xlab(expression(log[2](count + 1))) 
	
	setClass("gg")
	  setClass(
	  Class="QCPlot",
	  representation = representation(
	    total="gg",
	    zero="gg",
	    rRNA="gg",
	    density="gg",
	    density_all="gg",
	    Dmds="data.frame",
	    df="data.frame"
	  )
	  )

	
	 QCPlot <- methods::new( "QCPlot" )

	 QCPlot@total<-total
	 QCPlot@zero <- zero
	 QCPlot@rRNA<- rRNA
	 QCPlot@density<- density
	 QCPlot@Dmds<-as.data.frame(Dmds)
	 QCPlot@df<-df
	 QCPlot@density_all<- density_all
	 
	 return(QCPlot)
}


####dge results function#####

dgeresults <- function( object ){
	if( length( levels( object@individual ) ) >=2 ) {
		
		ddsMF <- agalmar::create_DESeq2( object, design = ~individual + treatment )
		
		###collapse technical replicates
		
		run=paste(ddsMF$individual,ddsMF$treatment,sep="_")
		
		if(length(unique(run))!=length(run)){
		ddsMF$run = ddsMF$run=as.numeric(as.factor(run))
		
		ddsMF <- DESeq2::collapseReplicates( ddsMF,ddsMF$run, renameCols=FALSE )
		}
		
		ddsMF <- DESeq2::DESeq( ddsMF )
		
		#plots
		
		plotMA1 <- function(ddsMF){
		  
			pairs <- combn(as.character(unique(ddsMF$treatment)),2)
			plotMA1<-list()
			for(i in seq(from=1, to=ncol(pairs))){
				resulttable <- DESeq2::results( ddsMF,contrast=c( "treatment", as.character( pairs[,i] ) ), pAdjustMethod="bonferroni", alpha=0.05 )
				resulttable <- as.data.frame(resulttable)
				resulttable$padj[is.na(resulttable$padj)]<-1
				plotMA1[[i]]<-ggplot(resulttable)+geom_point(alpha=0.8,aes(x=log(baseMean),y=log2FoldChange,colour= padj < 0.05))+scale_colour_manual(values =c("TRUE"="red","FALSE"="grey40"))+ggtitle(paste( as.character( pairs[,i]),collapse =" vs " ))
				plotMA1<<-plotMA1
			}
		}

	  plotMA1(ddsMF)
		
		#Transforming data
		rld <- DESeq2::rlog( ddsMF ) #rlog transform
		vsd <- DESeq2::varianceStabilizingTransformation( ddsMF ) #variance Stabilizing Transformation
		
		notAllZero <- ( rowSums( DESeq2::counts( ddsMF ) )>0 )

		#adding heatmap
		select <- order( rowMeans( counts( ddsMF,normalized=TRUE ) ), decreasing=TRUE )[1:20]
		nt <- DESeq2::normTransform( ddsMF )
		log2.norm.counts <- assay( nt )[ select, ]
		df <- as.data.frame( ddsMF@colData[ ,c( "treatment","individual" ) ] )
		
		count_heatmap <- pheatmap::pheatmap( log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df )
		
		rld_heatmap<- pheatmap::pheatmap( assay( rld )[select,], cluster_rows=FALSE, show_rownames=FALSE,
													cluster_cols=FALSE, annotation_col=df )
		
		vsd_heatmap<- pheatmap::pheatmap( assay( vsd )[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df )
		
		#sample distances
		sampleDists <- stats::dist( t( assay( rld ) ) )
		
		sampleDistMatrix <- as.matrix( sampleDists )
		rownames( sampleDistMatrix ) <- paste( rld$treatment, rld$individual, sep="-" )
		
		colnames( sampleDistMatrix ) <- NULL
		colors <- grDevices::colorRampPalette( rev( brewer.pal( 9, "Greens" ) ) )( 255 )
		distance_heatmap<- pheatmap::pheatmap( sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors )
		
		
		#PCA plot
		data <- DESeq2::plotPCA( rld, intgroup=c("treatment", "individual"), returnData=TRUE )
		percentVar <- round( 100 * attr(data, "percentVar") )
		
		pca <- ggplot( data, aes( PC1, PC2, color=individual, shape=treatment ) ) + scale_shape_manual(values=0:length(data$treatment)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +  ylab(paste0("PC2:",percentVar[2],"% variance"))
		
		setClass("gg")
		setClass(
		  Class="DGEplotresults",
		  representation = representation(
		    ddsMF="DESeqDataSet",
		    plotMA="list",
		    count_heatmap="pheatmap",
		    rld_heatmap="pheatmap",
		    vsd_heatmap="pheatmap",
		    distance_heatmap="pheatmap",
		    pca="gg",
		    rld="DESeqTransform",
		    vsd="DESeqTransform"
		  )
		)
		
		DGEPlots <- methods::new( "DGEplotresults" )
		  DGEPlots@ddsMF<-ddsMF
		  DGEPlots@plotMA<-plotMA1
		  DGEPlots@count_heatmap<-pheatmap( log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df )
		  DGEPlots@rld_heatmap<-pheatmap( assay( rld )[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df )
		  DGEPlots@vsd_heatmap<-pheatmap( assay( vsd )[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df )
		  DGEPlots@distance_heatmap<-pheatmap( sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors )
		  DGEPlots@pca<-pca
		  DGEPlots@rld<-rld
		  DGEPlots@vsd<-vsd
		  
	 return(DGEPlots)
		
	}
}

#######tree manipulation functions ###########################

# Clone gene tree edge lengths from species tree
clone_edge_lengths <- function( nhx, node_ages, ... ) {
  tags <- nhx@data
  
  # Only consider internal nodes
  tags <- tags[ is.na( tags$species ), ]
  
  # only consider speciation events
  tags <- tags[ ( tags$Ev == "S" ), ]
  
  # Only consider nodes for which age is available
  tags <- tags[ tags$S %in% names( node_ages ) , ]
  
  ages <- node_ages[ match( tags$S, names( node_ages ) ) ]
  
  calibration<-data.frame(
    node=tags$node,
    age.min=ages,
    age.max=ages,
    soft.bounds=rep( TRUE, nrow( tags ) )
  )
  
  # Verify the coherence of node ages, ie make sure that all calibration points are older than 
  # calibrations on younger nodes.
  oldest_descendant <- sapply( 1:max( nhx@phylo$edge ), function( x ) {
    progeny <- hutan::descendants( nhx@phylo, x )
    oldest <- 0
    if ( length( progeny ) > 0 ) {
      indices <- match( progeny, tags$node )
      if ( ! all( is.na( indices ) ) ) {
        oldest <- max( ages[ indices ], na.rm = TRUE )
      }
    }
    return( oldest )
  })
  
  ages_all <- rep( NA,max( nhx@phylo$edge ) ) 
  ages_all[ calibration$node ] <- calibration$age.min
  
  stopifnot( all( ages_all > oldest_descendant, na.rm=TRUE ) )
  
  # ggtree(nhx, branch.length="none") + geom_text(aes(label=round(oldest_descendant,3)), vjust=-.5, hjust=1, size=2.5, col="darkslategray4")+ geom_text(aes(label=S), vjust=1.5, hjust=2, size=2.5, col="blue") + geom_point(aes(color=Ev)) + geom_text(aes(label=round(ages_all,3)), vjust=-.5, hjust=2, size=2.5, col="red") + geom_text(aes(label=1:max(nhx@phylo$edge)), vjust=1.5, hjust=3, size=2.5, col="green")
  
  calibrated <- try( ape::chronos( nhx@phylo, calibration=calibration) )
  
  # Return the try error if failed, otherwise return updated nhx object
  if( "try-error" %in% class( calibrated ) ) {
    return( calibrated )
  }
  else{
    class( calibrated ) = "phylo"
    nhx@phylo=calibrated
    return( nhx )
  }
  
}


# Returns a logical vector indicating which rows in nhx_tags are for tips
is.tip.nhx <- function( nhx ) {
  is.tip <- rep( FALSE, nrow( nhx@data ) )
  is.tip[ 1:length( nhx@phylo$tip.label ) ] = TRUE
  is.tip
}

# Drops tips without expression data.

drop_empty_tips <- function( nhx, col_values ){
  if ( class( nhx ) != "treedata" ) {
    return( NA )
  }
  
  if (length(col_values)==1){
    to_drop <-which(is.na(nhx@data[ is.tip.nhx( nhx ), ][,names(nhx@data[is.tip.nhx( nhx ), ]) == col_values]))
    return( drop.tip.mod( nhx, to_drop ) )
  }
  
  if(length(col_values)>1){
    if(length(which(names(nhx@data[is.tip.nhx( nhx ), ]) %in% col_values))==1){
      to_drop <-which(is.na(nhx@data[ is.tip.nhx( nhx ), ][,names(nhx@data[is.tip.nhx( nhx ), ]) %in% col_values]))
      return( drop.tip.mod( nhx, to_drop ) )
    }
    else{ to_drop <- nhx@data[ is.tip.nhx( nhx ), ][rowSums(nhx@data[ is.tip.nhx( nhx ), ][,names(nhx@data[is.tip.nhx( nhx ), ]) %in% col_values],na.rm=TRUE)==0,]$node
  return( drop.tip.mod( nhx, to_drop ))
    }
    }
  }


# Run pic and add results to nhx object
pic.nhx <- function( nhx, col_value ) {
  
  if ( class( nhx ) != "treedata" ) {
    return( NA )
  }
  
  if(col_value %in% names(nhx@data)){
  
  p <- ape::pic( nhx@data[,names(nhx@data)==col_value][ is.tip.nhx( nhx ) ], nhx@phylo, var.contrasts=TRUE )
  
  # Verify that the pic names match the internal node names
  stopifnot( all( ( rownames( p ) == nhx@data$phy_node_names[ ! is.tip.nhx( nhx ) ] ) ) )
  
  m<-data.frame(pic=c( rep( NA, length( nhx@phylo$tip.label ) ), p[ ,1 ] ), var_exp=c( rep( NA, length( nhx@phylo$tip.label ) ), p[ ,2 ] ))
  
  names(m)<-c(paste("pic",col_value,sep="_"),paste("var_exp",col_value,sep="_"))
  
  nhx@data<-dplyr::bind_cols(nhx@data,m)

  return( nhx )
  }
}

# Transfers pic and var_exp results from nhx1 to nhx2, where nhx1 is a subtree of nhx2
merge_pic <- function( nhx1, nhx2 ) {
  if ( class( nhx2 )!="treedata" ) {
    return( nhx2 )
  }
  
  sub1 <- nhx1@data[ ,names( nhx1@data ) %in% c( "ND", "pic", "var_exp" ) ]
  
  merged <- merge( nhx2@data, sub1, all=TRUE, by="ND" )
  merged <- dplyr::arrange( merged, node )
  nhx2@data <- merged
  
  return( nhx2 )
}

# put node ages in @data
store_node_age <- function( nhx ) {
  
  if ( class( nhx ) != "treedata" ) {
    return( nhx )
  }
  
  node_age <- hutan::distance_from_tip( nhx@phylo )
  
  # make sure the dataframe is ordered by consecutive nodes
  stopifnot( all( nhx@data$node == 1:length( nhx@data$node ) ) )
  
  nhx@data$node_age <- node_age 
  
  return( nhx )
}

add_de_to_nhx <- function( nhx, combined_de ){
  
  # Get the de rows that apply to these tips
  de <- combined_de[ match( nhx@data$sequence_ids, combined_de$sequence_ids ) , ]
  de$species <- NULL
  de$gene_tree <- NULL
  
  nhx@data <- dplyr::bind_cols( nhx@data, de[,2:ncol(de)] )
  
  return(nhx)
}

# Return the first value in a vector that isn't NA
first_not_na <- function ( x ) {
  x <- x[ ! is.na( x ) ]
  if ( length( x ) == 0 ) {
    return( NA )
  }
  else{
    return ( x[1] )
  }
}

# Calculate the phylogenetic signal on an nhx tree for the variable in a specified nhx tag
phylosignal.nhx <- function ( nhx, tag, ... ) {
  if ( class( nhx ) != "treedata" ) {
    return( NA )
  }
  
  #exclude tags where number of tips <3
  if ( length(nhx@data[[ tag ]][!is.na(nhx@data[[ tag ]])]) < 3 ) {
    return( NA )
  }
  
  x <- nhx@data[[ tag ]]
  x <- x[ 1:length( nhx@phylo$tip.label ) ]
  names( x ) <- nhx@phylo$tip.label
  return ( phytools::phylosig( nhx@phylo, x, ... ) )
}

##from https://bitbucket.org/caseywdunn/sicb2013/
regularization_by_thresholding <- function ( cor_matrix, n ) {
  # Regularization by thresholding
  # Bickel, P. J. & Levina, E. Covariance regularization by thresholding. 
  # Ann. Statist. 36, 2577–2604 (2008). http://dx.doi.org/10.1214/08-AOS600
  
  p <- ncol( cor_matrix )
  
  regularized <- cor_matrix * ( abs(cor_matrix) > sqrt(log(p)/n) )
  
  return( regularized )
}



short_pairs<-function(pair){ 
  short<- pair%>%
    stringr::str_split(.," ")%>%
    unlist() %>%
    stringr::str_sub(.,1,3) %>%
    stringr::str_c(.,collapse="")
  return(short)
}

short_pairs_tpm<-function(pair){ 
  short<- pair%>%
    stringr::str_c(.,collapse="")
  return(short)
}

#Function to get DGE results for every pair
resultsdge <- function (object,pair) {
  treat1<-pair[1]
  treat2<-pair[2]
  
  if("CMN0040" %in% colnames(object@ddsMF) & treat1=="Gastrozooid mature")
  {
    treat1<-"Gastrozooid white mature"
  }
  
  else if("CMN0040" %in% colnames(object@ddsMF) & treat2=="Gastrozooid mature"){
    treat2<-"Gastrozooid white mature"
  }
  
  #check that the two treatments are present, and that there are more than two replicates
  if( (length(unique(object@ddsMF@colData$treatment)[unique(object@ddsMF@colData$treatment) %in% c(treat1,treat2) ==TRUE])==2)&((length(object@ddsMF@colData$treatment[object@ddsMF@colData$treatment==treat1])>=2)&(length(object@ddsMF@colData$treatment[object@ddsMF@colData$treatment==treat2])>=2)))  {
    
    
    short<-short_pairs(c(treat1,treat2))
    
    result <- data.frame() 
    res<-DESeq2::results(object@ddsMF, contrast=c("treatment", treat1, treat2)) %>% #write results to res
      as.data.frame
    names(res)<- paste(names(res), short ,sep="_")
    
    res<-tibble::rownames_to_column(res, var="sequence_ids")
    
    res$sequence_ids<-as.numeric(res$sequence_ids)
    
    result<-res
    
    return(result)
  }
  if((length(unique(object@ddsMF@colData$treatment)[unique(object@ddsMF@colData$treatment) %in% c(treat1,treat2) ==TRUE])<=1)){
    return(NA)
  }
  
  if((length(object@ddsMF@colData$treatment[object@ddsMF@colData$treatment==treat1])<=1)|(length(object@ddsMF@colData$treatment[object@ddsMF@colData$treatment==treat2])<=1)){
    return(NA)
  }
  
}

#pull out DGE results for every pair

pairsdge<- function(object){
  y<-list()
  for(i in 1:ncol(pairs)){
    res<-resultsdge(object,pairs[,i])
    y[[i]] <- res
  }
  return(y)
}

sort_by_pair<-function(object){
  m<-list()
    for(i in 1:length(pairs[1,])){
        m[[i]]<-lapply(object, extract, i) %>%  bind_rows(.)  
    }
      return(m)
  }
    
extract<-function(object,i){
  if(!is.data.frame(object[[i]])){
    return(NULL)
  }
  else{
    return(object[[i]])
  }
}

    
#pull out log tpm values

tpm_output<-function(object, treatment){
  
  if(object@species=="Bargmannia elongata" & treatment=="Gastrozooid mature"){
    treatment="Gastrozooid white mature"
  }
  
  if(length(object@treatment[object@treatment==treatment])>=2)  {
    
    Short_treat<- treatment%>%
      stringr::str_split(.," ")%>%
      unlist() %>%
      stringr::str_sub(.,1,3) %>%
      stringr::str_c(.,collapse="")
    
    count<-data.frame(rowMeans(object@tpm[,object@treatment ==treatment]))
    count<-log(count+1)
    names(count)<-Short_treat
    count<-tibble::rownames_to_column(count, var="sequence_ids")
    count$sequence_ids<-as.numeric(count$sequence_ids)
    return(count)
  }
  else if(length(object@treatment[object@treatment==treatment])<=1){
    return(NA)
  }
}

#pull out tpm values for each treatment

tpmtreat<-function(object){
  y<-list()
  for(i in 1:length(focal_treatments)){
    res<-tpm_output(object,focal_treatments[i])
    y[[i]] <- res
  }
  return(y)
}


sort_tpm<-function(object){
  m<-list()
  for(i in 1:length(focal_treatments)){
    m[[i]]<-lapply(object, extract, i) %>%  bind_rows(.)  
  }
  return(m)
}

#This code is modified from the function drop.tip from, https://github.com/GuangchuangYu/treeio. This is used to keep the original node labels associated with the data even after subsetting

add_node_name<-function(object){
    if ( class( object ) != "treedata" ) {
    return( NA )
  }
  
	if (is.null(object@phylo$node.label)) {
	  object@phylo$node.label <- ape::Ntip(object@phylo) + (1:treeio::Nnode(object@phylo))
	}
		 object@data$node.label_original <- c(object@phylo$tip.label, as.character(object@phylo$node.label))
	return(object)
}


###This is modified code from https://github.com/GuangchuangYu/treeio . This code is temporary, as there are issues calling drop.tip (see https://github.com/GuangchuangYu/treeio/issues/3). This will be fixed in the next release. 
drop.tip.mod <- function (object, tip) {
  
  node_label_name = "cd8128f329f72c167a8028cf8"
  
  if (!is.null(object@phylo$node.label)) {
    # Tree has node labels. Put these in data 
    # for safe keeping and remove them from tree
    # for now
    
    labels = c( 
      rep( NA, length( object@phylo$tip.label ) ), 
      object@phylo$node.label
    )
    object@data[[ node_label_name ]] <- labels
    
    object@phylo$node.label <- NULL
  }
  
  ## label the internal tree nodes by their number
  object@phylo$node.label <- ape::Ntip(object@phylo) + (1:treeio::Nnode(object@phylo))
  
  ## Prepare the nhx object for subsampling
  object@data$node <- as.numeric(object@data$node)
  object@data <- object@data[order(object@data$node),]
  
  ## add a colmn that has labels for both tips and internal nodes
  object@data$node.label <- c(object@phylo$tip.label, as.character(object@phylo$node.label))
  
  ## Will need to take different approaches for subsampling tips
  ## and internal nodes, add a column to make it easy to tell them apart
  object@data$is_tip <- object@data$node <= ape::Ntip(object@phylo)
  
  ## Remove tips
  object@phylo = ape::drop.tip( object@phylo, tip )
  
  ## Subsample the tags
  object@data = object@data[object@data$node.label %in% (c(object@phylo$tip.label, as.character(object@phylo$node.label))),]
  
  ## Update tip node numbers
  tip_nodes <- object@data$node.label[ object@data$is_tip ]
  object@data$node[ object@data$is_tip ] = match(object@phylo$tip.label, tip_nodes)
  
  internal_nodes <- object@data$node.label[ !object@data$is_tip ]
  object@data$node[ !object@data$is_tip ] = match(object@phylo$node.label, internal_nodes) + length(object@phylo$tip.label)
  
  ## Clean up
  object@data$node.label = NULL
  object@data$is_tip = NULL
  
  ## Add node labels back to tree, if there were any
  if (node_label_name %in% names( object@data ) ) {
    labels = object@data[[ node_label_name ]]
    ntips = length ( object@phylo$tip.label )
    labels = labels[ (ntips+1):nrow( object@data ) ]
    object@phylo$node.label = labels
    object@data[[ node_label_name ]] <- NULL
  }
  
  return(object)
}


anc.recon.new<-function(x,phy){
  #make new node labels that can be put back in later
  phy$node.label <- paste("node",ape::Ntip(phy)+(1:treeio::Nnode(phy)))   
  x=as.data.frame(x)
  #drop tips with missing data
  phy=ape::drop.tip(phy,which(is.na(x)))
  #Run ancestral reconstruction
  y= Rphylopars::anc.recon(x[!is.na(x)],phy,CI=TRUE)
  #add the old rownames back in
  rownames(y$Yhat)=phy$node.label
  rownames(y$lowerCI)=phy$node.label
  rownames(y$upperCI)=phy$node.label
  return(y)
}

add_ace_to_tree<- function(nhx,col_value) {
  
  if ( class( nhx ) != "treedata" ) {
    return( NA )
  }
  
  #Node storage methods used in this function adapted from treeio::drop.tip
  
  node_label_name = "cd8128f329f72c167a8028cf8"
  
  if (!is.null(nhx@phylo$node.label)) {
    # Tree has node labels. Put these in data 
    # for safe keeping and remove them from tree
    # for now
    
    labels = c( 
      rep( NA, length( nhx@phylo$tip.label ) ), 
      nhx@phylo$node.label
    )
    nhx@data[[ node_label_name ]] <- labels
    
    nhx@phylo$node.label <- NULL
  }
  
 # y=anc.recon.new(nhx@data[,names(nhx@data)==col_value][is.tip.nhx(nhx),],nhx@phylo)
  y=anc.recon.new(dplyr::pull(nhx@data, col_value)[is.tip.nhx(nhx)],nhx@phylo)
  
  
  nhx@phylo$node.label <- paste("node",Ntip(nhx@phylo)+(1:treeio::Nnode(nhx@phylo)))
  
  nhx@data$node.label <- c(nhx@phylo$tip.label, as.character(nhx@phylo$node.label))
  
  #add the ace values into a new column
  nhx@data[,paste("ace",col_value,sep="_")]<-rep(NA,nrow(nhx@data)) %>% as.numeric()
  nhx@data[,paste("ace_CI_1",col_value,sep="_")]<-rep(NA,nrow(nhx@data)) %>% as.numeric()
  nhx@data[,paste("ace_CI_2",col_value,sep="_")]<-rep(NA,nrow(nhx@data))%>% as.numeric()
  nhx@data[,paste("ace",col_value,sep="_")][which(nhx@data$node.label %in% rownames(y$Yhat)),]<-y$Yhat %>% as.numeric()
  nhx@data[,paste("ace_CI_1",col_value,sep="_")][which(nhx@data$node.label %in% rownames(y$Yhat)),]<-y$lowerCI %>% as.numeric()
  nhx@data[,paste("ace_CI_2",col_value,sep="_")][which(nhx@data$node.label %in% rownames(y$Yhat)),]<-y$upperCI %>% as.numeric()
  
  ## Add node labels back to tree, if there were any
  if (node_label_name %in% names(nhx@data ) ) {
    labels = nhx@data[[ node_label_name ]]
    ntips = length ( nhx@phylo$tip.label )
    labels = labels[ (ntips+1):nrow( nhx@data ) ]
    nhx@phylo$node.label = labels
    nhx@data[[ node_label_name ]] <- NULL
  }
  
  return(nhx)
}

ace_trees<- function(nhx,col_value){
  n_focal_tips <-function( nhx, col_value ){
    if ( class( nhx ) != "treedata" ) {
      return( NA )
    }
    x<-sum(!is.na(nhx@data[,names(nhx@data)==col_value]))
    return(x)
  }
  y<-lapply(nhx,n_focal_tips,col_value) %>%
    unlist()
  
  subtree <- nhx[ y >= min_tips ]
  
  subtree_ace<-lapply(subtree,add_ace_to_tree, col_value)
  
  return(subtree_ace)	
}

merge_ace <- function( nhx1, nhx2) {
  #Where sub is a subtree that is reduced down to only the added ace and node.label columns.
  if ( class( nhx2 )!="treedata" ) {
    return( nhx2 )
  }
  
  if (any(grep("ace",names(nhx1@data)))==FALSE) {
    return( nhx2 )
  }
  
  else{
  sub1 <- nhx1@data[ ,c(paste(names(nhx1@data)[grep("ace",names(nhx1@data))]),"ND")]
  
  merged <- merge( nhx2@data, sub1, all=TRUE, by="ND" )
  merged <- dplyr::arrange( merged, node )
  nhx2@data <- as_tibble(merged)
  
  return( nhx2 )
  }
}


###modified from agalmar::summarize_edges
summarize_ace_edges = function ( nhx, default_length_val=NA ) {
  if ( class( nhx ) == "treedata" ) {
  # Create a data frame of internal node annotations
  tags = nhx@data
  tags$node =  tags$node 
  tags$S = tags$S 
  tags$ND = tags$ND 
  tags = tags[order( tags$node ),]
  
  parents = nhx@phylo$edge[,1]
  children = nhx@phylo$edge[,2]
  
  terminal = rep( FALSE, nrow( nhx@phylo$edge ) )
  terminal[ children <= length( nhx@phylo$tip.label ) ] = TRUE
  
  #add ace values to one column
  
  edge_length=nhx@phylo$edge.length
  
  df = data.frame( 
    gene_tree = rep( digest::digest( nhx ), nrow( nhx@phylo$edge ) ),
    length = edge_length, 
    Ev_parent = tags$Ev[parents],
    S_parent  = tags$S[parents] ,
    ND_parent = tags$ND[parents] ,
    node_parent = tags$node[parents],
    node_depth_parent=tags$node_depth[parents] ,
    node_age_parent=tags$node_age[parents],
    Gasdev_parent= tags$Gasdev[parents] ,
    Palmat_parent= tags$Palmat[parents] ,
    Gasmat_parent= tags$Gasmat[parents] ,
    Necdev_parent= tags$Necdev[parents] ,
    Pne_parent= tags$Pne[parents] ,
    Bradev_parent= tags$Bradev[parents] ,
    Gonmal_parent= tags$Gonmal[parents] ,
    Gonfem_parent= tags$Gonfem[parents] ,
    PalmatGasmat_parent=(tags$Palmat[parents]+1)/(tags$Gasmat[parents]+1) ,
    GasdevGasmat_parent=(tags$Gasdev[parents]+1)/(tags$Gasmat[parents]+1) ,
    NecdevGasmat_parent=(tags$Necdev[parents]+1)/(tags$Gasmat[parents]+1) ,
    PneGasmat_parent=(tags$Pne[parents]+1)/(tags$Gasmat[parents]+1) ,
    tau_parent=tags$tau[parents],
    Ev_child = tags$Ev[children],
    S_child  = tags$S[children] ,
    ND_child = tags$ND[children] ,
    node_child = tags$node[children],
    node_depth_child=tags$node_depth[children] ,
    node_age_child=tags$node_age[children],
    Gasdev_child= tags$Gasdev[children] ,
    Palmat_child= tags$Palmat[children] ,
    Gasmat_child= tags$Gasmat[children] ,
    Necdev_child= tags$Necdev[children] ,
    Pne_child= tags$Pne[children] ,
    Bradev_child= tags$Bradev[children] ,
    Gonmal_child= tags$Gonmal[children] ,
    Gonfem_child= tags$Gonfem[children] ,
    PalmatGasmat_child=(tags$Palmat[children]+1)/(tags$Gasmat[children]+1) ,
    GasdevGasmat_child=(tags$Gasdev[children]+1)/(tags$Gasmat[children]+1) ,
    NecdevGasmat_child=(tags$Necdev[children]+1)/(tags$Gasmat[children]+1) ,
    PneGasmat_child=(tags$Pne[children]+1)/(tags$Gasmat[children]+1) ,
    tau_child=tags$tau[children],
    Gasdev_scaled=(tags$Gasdev[children] -tags$Gasdev[parents] )/edge_length,
    Palmat_scaled=(tags$Palmat[children] -tags$Palmat[parents] )/edge_length,
    Gasmat_scaled=(tags$Gasmat[children] -tags$Gasmat[parents] )/edge_length,
    Necdev_scaled=(tags$Necdev[children] -tags$Necdev[parents] )/edge_length,
    Pne_scaled=(tags$Pne[children] -tags$Pne[parents] )/edge_length,
    Bradev_scaled=(tags$Bradev[children] -tags$Bradev[parents] )/edge_length,
    Gonmal_scaled=(tags$Gonmal[children] -tags$Gonmal[parents] )/edge_length,
    Gonfem_scaled=(tags$Gonfem[children] -tags$Gonfem[parents] )/edge_length,
    PalmatGasmat_scaled = ( ((tags$Palmat[children]+1)/(tags$Gasmat[children]+1)) - ((tags$Palmat[parents]+1)/(tags$Gasmat[parents]+1)) )/edge_length ,
    GasdevGasmat_scaled = ( ((tags$Gasdev[children]+1)/(tags$Gasmat[children]+1)) - ((tags$Gasdev[parents]+1)/(tags$Gasmat[parents]+1)) )/edge_length ,
    NecdevGasmat_scaled = ( ((tags$Necdev[children]+1)/(tags$Gasmat[children]+1)) - ((tags$Necdev[parents]+1)/(tags$Gasmat[parents]+1)) )/edge_length  ,
    PneGasmat_scaled = ( ((tags$Pne[children]+1)/(tags$Gasmat[children]+1)) - ((tags$Pne[parents]+1)/(tags$Gasmat[parents]+1)) )/edge_length ,
    tau_scaled=(tags$tau[children]-tags$tau[parents])/edge_length,
    terminal = terminal,
    default_length = FALSE,
    stringsAsFactors = FALSE
  )
  
  if ( ! is.na( default_length_val ) ){
    df$default_length = dplyr::near( nhx@phylo$edge.length, default_length_val )
  }
  
  return( df )
  }
}

summarize_trees = function( gene_trees, col_value ) {
  
  
  tree_summary = 
    lapply( 
      gene_trees, 
      function( nhx ){
        if ( class( nhx ) == "treedata" ) {
        tags = nhx@data
        phy = nhx@phylo
        x = nhx@data[,names(nhx@data)==col_value][1:length(phy$tip.label),] %>% unlist()
        names( x ) = phy$tip.label
        if(length(x[!is.na(x)])>3){
        species<-tags$species[!is.na(x)][1:length(phy$tip.label)]
        phy=ape::drop.tip(phy,which(is.na(x)))
        event<-nhx@data$Ev[nhx@data$phy_node_names %in% phy$node.label]
        x<-x[!is.na(x)]
        
        tibble(
          gene = digest( nhx ),
          zooid = col_value,
          n_tips = length( x ),
          mean = mean( x ),
          var = var( x ),
          K = phytools::phylosig( phy, x, method="K" ),
          n_species = length(unique(species)),
          n_dup = length(event[event=="D"]),
          n_spec = length(event[event=="S"])
        )
        
        }
        }
      }
    ) %>%
    bind_rows() 
  
  return( tree_summary ) 
}

summarize_nodes.mod = function ( nhx, default_length_val=NA ) {
  if ( class( nhx ) == "treedata" ) {
  # Create a data frame of internal node annotations
  tags = cbind( 
    gene_tree= digest::digest( nhx ), 
    nhx@data
  )
  
  # Add a boolean column that indicates if nodes are parents to
  # edges with default length 
  tags$default_length = FALSE
  if ( ! is.na( default_length_val ) ){
    default_edges = dplyr::near( nhx@phylo$edge.length, default_length_val )
    parent_nodes = nhx@phylo$edge[ , 1 ]
    default_nodes = parent_nodes[ default_edges ]
    tags$default_length[ default_nodes ] = TRUE
  }
  
  tags %<>% 
    dplyr::mutate_if( is.factor, as.character )
  
  return( tags )
  }
}

tau=function(nhx){
  if ( class( nhx ) != "treedata" ) {
    return( NA )
  }
  
  x=nhx@data
  x=x[,names(nhx@data) %in% c("Gasdev","Palmat","Gasmat","Necdev","Pne","Gonmal","Gonfem")]
  tau=apply(x, 1, function(y) sum(y/max(y))/length(y))
  nhx@data$tau<-tau
  return(nhx)
}

sim_values = function(nhx, col_value, dup_adjust=1, a=NA ) {
  
  if ( class( nhx ) != "treedata" ) {
    return( NA )
  }
  
  phy = nhx@phylo
  
  stopifnot(nhx@data$phy_node_names[is.tip.nhx(nhx)]==phy$tip.label)
  
  col_value_original = dplyr::pull(nhx@data,col_value)[is.tip.nhx(nhx)]
  names( col_value_original ) = phy$tip.label
  
  phy <- ape::drop.tip(phy,names(which(is.na(col_value_original))))
  
  col_value_mod = col_value_original[!is.na(col_value_original)]
  
  brownian_model = tryCatch(fitCont( 
    phy, 
    col_value_mod,
    model="BM"), error=function(e) NULL)
  
  if(!is.null(brownian_model)){
  # Simulate trait given the tree and parameter estimates - bounding between 0 to 9 as this matches empirical observations
    
  x = phytools::fastBM( 
    phy, 
    a = brownian_model$a, 
    sig2 = brownian_model$sig2,
    bounds=c(0,9)
  ) %>% abs()
  
  #remove ancestral trait reconstructions
  nhx@data[,names(nhx@data)==col_value][! is.tip.nhx(nhx), ]<-rep(NA,length(nhx@data[,names(nhx@data)==col_value][! is.tip.nhx(nhx),]))
  
  names( x ) = NULL
  nhx@data[,names(nhx@data)==col_value][which(nhx@data$phy_node_names %in% phy$tip.label),] <- x
  
  #subset data to only include only values of interest
  nhx@data<-nhx@data[,names(nhx@data) %in% c("ND",col_value)]
  
  rm(phy)
  
  return( nhx )
  }
  else{
    return(NA)
  }
}


#must subset down to ensure that n=3 in the trees for BM & OU reconstruction
add_model_parameters_trees<- function(nhx,col_value){
  
  n_focal_tips <-function( nhx, col_value ){
    if ( class( nhx ) != "treedata" ) {
      return( 0 )
    }
    
    x<-sum(!is.na(nhx@data[,names(nhx@data)==col_value]))
    return(x)
  }
  
  y<-lapply(nhx,n_focal_tips,col_value) %>%
    unlist()
  
  subtree <- nhx[ y >= 3 ]
  
  #simulate BM models on the subtree
  subtree_sim<-lapply(subtree,sim_values, col_value)
  
  #generate ace values for the simulated data
  subtree_ace<-parallel::mclapply(subtree_sim,add_ace_to_tree, col_value, mc.cores=cores)
  
  return(subtree_ace)	
}

merge_sim <- function( nhx1, nhx2) {
  if ( class( nhx1 ) != "treedata" ) {
    return( nhx2 )
  }
  
  #Where sub is a subtree that is reduced down to only the added ace and node.label columns.
  if ( class( nhx2 )!="treedata" ) {
    return( nhx2 )
  }
    
  if (any(grep("ace",names(nhx1@data)))==FALSE) {
    return( nhx2 )
  }
    
    sub1<- nhx1@data[,!(names(nhx1@data) %in% c("Ev","S","node","phy_node_names","species","sequence_ids","node_depth","blast_hit","node_age","node.label_original","node.label"))]
    
    merged <- merge( nhx2@data, sub1, all=TRUE, by="ND" )
    merged <- dplyr::arrange( merged, node )
    nhx2@data <- as_tibble(merged)
    
    return( nhx2 )
}


####code below from the very helpful Liam Revell http://blog.phytools.org/2014/10/alternative-implementations-of.html######

## helper function to get the root node number
getRoot<-function(tree) ape::Ntip(tree)+1

## here is our fitContinuous lite function
fitCont<-function(tree,x,model="BM",interval=NULL){
  
  if(model!="BM"){
    lk<-function(par,tree,x,model) phytools::brownie.lite(phytools::paintSubTree(geiger::rescale(tree,model,par),getRoot(tree),"1"),x)$logL1
    oFit<-stats::optimize(lk,interval,tree=tree,x=x,model=model,maximum=TRUE)
    cFit<-phytools::brownie.lite(phytools::paintSubTree(geiger::rescale(tree,model,oFit$maximum),getRoot(tree),"1"),x)
    obj<-list(par=oFit$maximum,model=model,sig2=cFit$sig2.single,a=cFit$a.single,logLik=oFit$objective)
  } else {
    fit<-phytools::brownie.lite(phytools::paintSubTree(tree,getRoot(tree),"1"),x)
    obj<-list(model=model,sig2=fit$sig2.single,a=fit$a.single,logLik=fit$logL1)
  }
  obj
}

####code above from the very helpful Liam Revell http://blog.phytools.org/2014/10/alternative-implementations-of.html######

branch_changes<-function(branch){
  
setClass(Class="branchchanges",
      representation = representation(
      filtered_edge="data.frame",
      positive="data.frame",
      negative="data.frame",
      neutral="data.frame" 
    )
)
  
filtered_edges<- edges_tidy %>% dplyr::filter(change>change_deciles[2] & change<change_deciles[10] & S_child==branch) %>% group_by(gene_tree,treatment) %>% summarize(positive=length(change[change> 0.5]),negative=length(change[change< -0.5]),neutral=length(change[change> -0.5 & change< 0.5]),blast_hit=blast_hit[[1]]) %>% dplyr::arrange(desc(positive))
  
positive <- filtered_edges %>% group_by(gene_tree) %>% summarize(n_positive=length(positive[positive>0]),positive_gasmat=length(positive[treatment=="Gasmat" & positive>0]),positive_palmat=length(positive[treatment=="Palmat" & positive>0]),positive_necdev=length(positive[treatment=="Necdev" & positive>0]), positive_pne=length(positive[treatment=="Pne" & positive>0]),positive_gasdev=length(positive[treatment=="Gasdev" & positive>0]),positive_gonfem=length(positive[treatment=="Gonfem" & positive>0]),positive_gonmat=length(positive[treatment=="Gonmat" & positive>0]),blast_hit=blast_hit[[1]]) %>% dplyr::arrange(desc(n_positive))
  
negative <- filtered_edges %>% group_by(gene_tree) %>% summarize(n_neutral=length(neutral[neutral>0]),neutral_gasmat=length(neutral[treatment=="Gasmat" & neutral>0]),neutral_palmat=length(neutral[treatment=="Palmat" & neutral>0]),neutral_necdev=length(neutral[treatment=="Necdev" & neutral>0]), neutral_pne=length(neutral[treatment=="Pne" & neutral>0]),neutral_gasdev=length(neutral[treatment=="Gasdev" & neutral>0]),neutral_gonfem=length(neutral[treatment=="Gonfem" & neutral>0]),neutral_gonmal=length(neutral[treatment=="Gonmal" & neutral>0]),blast_hit=blast_hit[[1]]) %>% dplyr::arrange(desc(n_neutral))
  
neutral <- filtered_edges %>% group_by(gene_tree) %>% summarize(n_neutral=length(neutral[neutral>0]),neutral_gasmat=length(neutral[treatment=="Gasmat" & neutral>0]),neutral_palmat=length(neutral[treatment=="Palmat" & neutral>0]),neutral_necdev=length(neutral[treatment=="Necdev" & neutral>0]), neutral_pne=length(neutral[treatment=="Pne" & neutral>0]),neutral_gasdev=length(neutral[treatment=="Gasdev" & neutral>0]),neutral_gonfem=length(neutral[treatment=="Gonfem" & neutral>0]),neutral_gonmal=length(neutral[treatment=="Gonmal" & neutral>0]),blast_hit=blast_hit[[1]]) %>% dplyr::arrange(desc(n_neutral))
  
branchchanges <- methods::new( "branchchanges" )
branchchanges@filtered_edge <-filtered_edges
branchchanges@positive <- positive
branchchanges@negative <- negative
branchchanges@neutral <- neutral
  
  return(branchchanges)
}

#### Read in gene trees from orthofinder -- this is project specfic code

read_gene_trees<-function(file){
  
  treetext <- readLines(file, warn=FALSE)
  treetext <- treetext[treetext != ""]
  treetext <- treetext[treetext != " "]
  
  if (length(treetext) > 1) {
    treetext <- paste0(treetext, collapse = '')
  }
  treetext %<>% gsub(" ", "",. )
  
  phylo <- read.tree(text=treetext)
  nnode <- treeio::Nnode(phylo, internal.only=FALSE)
  
  Ev<-rep(NA,nnode) %>% as.character()
  
  
  S<-rep(NA,nnode) %>% as.character()
  
  ND <- rep(NA,nnode) %>% as.character()
  node <- 1:nnode %>% as.numeric()
  
  phy_node_names <- c(phylo$tip.label, phylo$node.label)
  
  # Parse sequence id, the integer after @, from the tip names
  sequence_ids <- phy_node_names
  sequence_ids[ !grepl( '@', sequence_ids ) ] <- NA
  sequence_ids <- as.numeric( sub( '^.+@', '', sequence_ids, perl=TRUE ) )
  
  # Parse species, the character string before @, from the tip names
  species_names <- phy_node_names
  species_names[ !grepl( '@', species_names ) ] <- NA
  species_names <- sub( '@.+$', '', species_names, perl=TRUE )
  species_names <- sub( '_', ' ', species_names, perl=TRUE )
  
  node_depth = ape::node.depth( phylo )
  
  tags<-tibble(Ev=Ev,
               ND=ND,
               node=node,
               phy_node_names=phy_node_names, 
               species=species_names,
               sequence_ids=sequence_ids,
               node_depth=node_depth,
               stringsAsFactors=FALSE
  )
  
  
  Orthogroup_specific <- sub(".*(OG\\w+)_tree.txt","\\1",file, perl=TRUE)
  Orthodf<- Duplication_nodes %>% dplyr::filter(Orthogroup == Orthogroup_specific) %>% dplyr::rename("S" = `Species Tree Node`,"ND" = `Gene Tree Node`) %>% dplyr::select(-c(Type,`Genes 1`,`Genes 2`))
  
  #Annotate with D or S
  
  tags[tags$phy_node_names %in% Orthodf$ND,] %<>% dplyr::mutate(Ev="D")
  tags[!tags$phy_node_names %in% Orthodf$ND,] %<>% dplyr::mutate(Ev="S")
  
  #Add Species and gene tree node
  tags$ND <- c(seq(1:length(phylo$tip.label)), phylo$node.label)
  
  stopifnot(
    all(tags[is.na(tags$species),]$ND == tags[is.na(tags$species),]$phy_node_names)
    
  )
  
  tags<- dplyr::left_join(tags,Orthodf,by=c("phy_node_names"="ND")) %>% dplyr::select(Ev, ND, S, node, phy_node_names,species,sequence_ids,node_depth,Support)
  
  #there's probably a cleaner way to do this... adding tip numbers to S in tags. Reads species name and checks correspondence
  tags[which(tags$node_depth=="1"& tags$Ev=="S"),]$S <- node_comparison$numbered[match(tags[which(tags$node_depth=="1"& tags$Ev=="S"),]$species,node_comparison$species)] 
  
  getspeciesnode = function(node_n){
    descendant=tags$species[phytools::getDescendants(phylo, node_n)] %>% unique %>% .[!is.na(.)]
    MRCA=phytools::findMRCA(phyldog_species_numbered_tree,descendant)
    x=as.character(node_comparison$species[MRCA])
    return(x)
  }
  
  tags[which(tags$node_depth!="1"& tags$Ev=="S"),] %<>% dplyr::mutate(S=sapply(node,getspeciesnode))
  
  #now go in and correct the instances where the node value is a species name
  
  if(any(tags[which(tags$node_depth!="1"& tags$Ev=="D"),]$S %in% phyldog_species_numbered_tree$tip.label ==TRUE)){
    tags[which(tags$node_depth!="1"& tags$Ev=="D"),] %<>% dplyr::mutate(S=node_comparison$numbered[match(tags[which(tags$node_depth!="1"& tags$Ev=="D"),]$S,node_comparison$species)])
  }
  
  #now go in and correct errors at S nodes arising when topologies are out of order and descendents are the same. 
  
  descendant_check<-function(node_n){
    node_S<-tags$S[tags$node==node_n]
    descendant_stats=data.frame(
    descendant_S=tags$S[phytools::getDescendants(phylo, node_n)],
    descendant_n=tags$node[phytools::getDescendants(phylo, node_n)], 
    descendant_Ev=tags$Ev[phytools::getDescendants(phylo, node_n)]
    )
    if(any(descendant_stats$descendant_S %in% node_S)){
      substitute_node<-descendant_stats %>% dplyr::filter(descendant_S %in% node_S & descendant_Ev=="S") %>% .$descendant_n
      tags$S[tags$node%in%substitute_node]<-"XX"
    }
    
    return(tags)
  }
  
  key_nodes <-tags[which(tags$node_depth!="1"& tags$Ev=="S"),]$node
  if(length(key_nodes>=1)){
  for(key_node in 1:length(key_nodes)){
    tags$S[descendant_check(key_nodes[[key_node]])$S%in%"XX"]<-"XX"
  }
  }
  
tree<- new("treedata",
             phylo = phylo,
             data = tags
  )
  
  return(tree)
  
}


summarize_edges_new = function ( nhx, default_length_val=NA ) {
  
  # Create a data frame of internal node annotations
  tags = nhx@data
  tags$node = tags$node
  tags$S = tags$S 
  tags$ND = tags$ND
  tags = tags[order( tags$node ),]
  
  parents = nhx@phylo$edge[,1]
  children = nhx@phylo$edge[,2]
  
  terminal = rep( FALSE, nrow( nhx@phylo$edge ) )
  terminal[ children <= length( nhx@phylo$tip.label ) ] = TRUE
  
  df = data.frame(
    gene_tree = rep( digest::digest( nhx ), nrow( nhx@phylo$edge ) ),
    length = nhx@phylo$edge.length,
    Ev_parent = tags$Ev[parents],
    S_parent  = tags$S[parents],
    ND_parent = tags$ND[parents],
    Ev_child = tags$Ev[children],
    S_child  = tags$S[children],
    ND_child = tags$ND[children],
    terminal = terminal,
    default_length = FALSE,
    stringsAsFactors = FALSE
  )
  
  if ( ! is.na( default_length_val ) ){
    df$default_length = dplyr::near( nhx@phylo$edge.length, default_length_val )
  }
  
  return( df )
}

clean_single_copy<-function(expression){	
  keep<-rep(NA, length(expression$treatment))
  for(treatment in 1:length(expression$treatment)){
    keep[[treatment]]= length(expression$treatment[expression$treatment %in% expression$treatment[[treatment]]])>=2
  }
  expression$library_id <- expression$library_id[keep]
  expression$treatment <- expression$treatment[keep]
  expression$individual <- expression$individual[keep]
  expression$sample_prep <- expression$sample_prep[keep]
  expression$read_count <- expression$read_count[keep]
  expression$count <- expression$count[,keep]
  expression$fpkm <- expression$fpkm[,keep]
  expression$tpm <- expression$tpm[,keep]
  
  return(expression)
  
}

GOseqresults<-function(geneList,lengthlist, annot){
  pwf<-goseq::nullp(geneList,bias.data = lengthlist)
  GO.wall<-goseq::goseq(pwf,gene2cat = annot)
  k <- as.data.frame(GO.wall)
  k$bh_adjust <-  p.adjust(k$over_represented_pvalue,method="BH")
  enr <- subset(k, k$bh_adjust <.05)
  return(enr)
}

Go_lookup<-function(species, zooid_list){
  geneList <- factor(as.integer(names(get(paste0(species,"_GO",sep=""))) %in% zooid_list))
  names(geneList) <- names(get(paste0(species,"_GO",sep="")))
  Go_terms<-GOseqresults(geneList,get(paste0(species,"_length",sep="")),get(paste0(species,"_annot",sep="")))
  return(Go_terms)
}
