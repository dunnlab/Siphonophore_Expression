 # Load libraries
	#biocLite libraries
	#source("https://bioconductor.org/biocLite.R")
		# biocLite("BiocUpgrade")
  #biocLite()
	library(DESeq2)
	library(vsn)
	library(geneLenDataBase)
	library(ggtree)
  library(edgeR)

	#devtools
  # library( devtools )
	library(hutan)
		# install_github( 'caseywdunn/hutan' )
	library("agalmar")
    #install_github( 'caseywdunn/agalmar')
 
	#CRAN libraries
#	library(gridExtra)
  library(ggplot2)
	library(knitr)
	library(fields)
  library(ape)
	library(picante)
	library(jsonlite)
	library(RColorBrewer)
	library(pheatmap)
	library(readr)
	library(SparseM)
	library(tidyverse)
  library(ggrepel)
  library(GGally)
  library(digest)
  library(parallel)
  library(magrittr)
	library(phytools)

  #this is a biocLite library
  library(treeio)

	source( "functions.R" )

  	## The minimum number of counts to pass gene sampling criteria
	min_count <- 1 
	
	## p value cutoff for evaluating differential expression significance
	p_value_threshold <- 0.05 
	
	focal_species <- c( "Agalma elegans", "Nanomia bijuga", "Bargmannia elongata", "Frillagalma vityazi", "Diphyes dispar", "Physalia physalis", "Apolemia lanosa" )
	
	### The default branch length pace holder used in the gene tree inference software
	default_length_val <- 1e-06 
	
	### Exclude trees with branches that exceed length threshold
	edge_length_max_threshold <- 2

	### Exclude trees that exceed specified fraction of branches with default values
	fraction_default_max_threshold <- 0.25
	
	### The minimum number of tips with expression data for a tree to be considered 
	min_tips <- 3
	
	### The maximum root depth of gene trees, where the root of the species tree is 1
	### Gene trees that exceed this threshold are removed
	max_root_depth <- 5
	
	### Number of speciation nodes that don't match topology of species tree that can be tolerated 
	### A pecularity of Orthofinder's species overlap method (a clade of genes which have a single copy of each gene will be considered a speciation event regardless of topology)
	### I used a method in the function read_gene_trees to identify instances where this occurs and assign X to descendant nodes (otherwise descendant nodes will have the same node ID 'S' as the parent node and this affects calibration). 
	### Too many Xs in a tree may indicate a poorly constructed tree, so I need to filter this out. 
	max_x_per_node<-0.3
	
	# Set system computational parameters
	cores <- detectCores() - 1
	if ( cores < 1 ) {
		cores <- 1
	}
	set.seed( 23456 )
	
	cores <- detectCores() - 1
	
	unlink( ".RData" )

	data_path <- "export.json"

	agalma.data <- jsonlite::fromJSON( data_path )
	
	tpm<-function(count, length){
	  constant<- sum(count/length )
	  tpm_new<- (1000000 * count) / (length *constant)
	  return(tpm_new)
	}
	
	#adjust tpm values by a scalar to adjust for library size
	adjusttpm<-function(tpm_new){
	  tpm_new<- (tpm_new * nrow(tpm_new))/10000
	  return(tpm_new)
	}
	
	applytpm<-function(object){
	
	length<-object@lengths %>% as.numeric() 
	#get the effective length, in this case 50bp-1
	length<-length-49
	tpm_new<-apply(object@x, 2, tpm,length)
	tpm_new<-adjusttpm(tpm_new)
	
	return(tpm_new)
	
	}
	
#remove expression values where only one replicate is present
	
agalma.data$expression <-	lapply(agalma.data$expression,clean_single_copy)

#clean up Physalia names

agalma.data$expression$`HWI-ST625-73-C0JUVACXX-7-TTAGGC`$treatment[agalma.data$expression$`HWI-ST625-73-C0JUVACXX-7-TTAGGC`$treatment=="Ampule mature"] <-"Tentacular palpon mature"
agalma.data$expression$`HWI-ST625-73-C0JUVACXX-7-TTAGGC`$treatment[agalma.data$expression$`HWI-ST625-73-C0JUVACXX-7-TTAGGC`$treatment=="Ampule developing"] <-"Tentacular palpon developing"
	
e <- parallel::mclapply( agalma.data$expression, Expression, mc.cores=cores )
	
e<-lapply(e, function(object){
  object@tpm<-applytpm(object)
  return(object)
	})
	
	phyldog_species_numbered_tree <- ape::read.tree( text= agalma.data$speciestree_numbered)
	
	species_tree <- ape::read.tree( text=agalma.data$speciestree )
  
	gene_trees_path <- "Results_Feb16/Resolved_Gene_Trees"
	gene_trees_paths <- list.files(gene_trees_path, "*_tree.txt", full.names =TRUE)
	gene_trees_names <- dir(gene_trees_path, pattern= "*_tree.txt") %>% str_remove(., "_tree.txt")
	
	Duplication_nodes<-readr::read_tsv("Results_Feb16/Gene_Duplication_Events/Duplications.tsv")
	
	#Cleanup names for later
	Duplication_nodes$`Species Tree Node`<- Duplication_nodes$`Species Tree Node` %>% sub( "_\\d+", "", . )%>%
	  sub( '_', ' ', . ) %>%  sub( '_', '', . )
	
	rm( agalma.data )

#generate DGE results
	
dge <- lapply( e, dgeresults )


#library_summary<- mclapply(e,summarize_libraries)

##QC plots


QCPlot<-lapply( e, QCplot )

all<-	rbind(QCPlot$`K00162-189-HJTYGBBXX-7-NCAGTG`@df,QCPlot$`HWI-ST625-159-C4MVCACXX-5-CCGTCC`@df,QCPlot$`HWI-ST625-73-C0JUVACXX-7-AGALMA2`@df,QCPlot$`HWI-ST625-51-C02UNACXX-8-BARGMANNIA`@df,QCPlot$`HWI-ST625-51-C02UNACXX-6-FRILLAGALMA`@df,QCPlot$`HWI-ST625-51-C02UNACXX-7-NANOMIA`@df,QCPlot$`HWI-ST625-73-C0JUVACXX-7-TTAGGC`@df)

all_dist<-ggplot(all, aes(x = value, colour = treatment, fill=treatment)) + 
		ylim(c(0, 0.25)) + 
		geom_density(alpha = 0.05, size = 0.5) + 
		facet_wrap(~species) + 
		theme(legend.position = "top") + 
		xlab(expression(log[2](count + 1))) 

combined_dist<-ggplot(all, aes(x = value, colour = treatment, fill=treatment)) + 	ylim(c(0, 0.25)) +		geom_density(alpha = 0.05, size = 0.5) + 		theme(legend.position = "top") + 	xlab(expression(log[2](count + 1))) 

  # Parse the blast results into a data frame

	blast_lookup <- 
	  parallel::mclapply( 
			e, 
			function( x ) data.frame( 
				sequence_ids=rownames( x@x ), 
				blast_hit=x@blast_hit, 
				stringsAsFactors=FALSE 
			), 
			mc.cores=cores 
		) %>%
		dplyr::bind_rows()

	blast_lookup %<>% dplyr::filter( !is.na( blast_hit ) )
#  blast_lookup$sequence_ids %<>% as.integer( )
	
	# Shorten blast hit descriptions be removing redundant information
	blast_lookup$blast_hit %<>% stringr::str_replace_all( "swissprot\\|", "" )
	blast_lookup$blast_hit %<>% stringr::str_replace_all( " OS=.+", "" )
	
	#e.g usage blast_lookup$blast_hit[blast_lookup$sequence_ids==4130517]


	# Repartition length of edges descended from root so they are non-zero
	species_tree <- hutan::slide_root_edges( species_tree )

	# species_tree = ape::root( species_tree, c( "Nematostella_vectensis", "Aiptasia_pallida" ), resolve.root=TRUE )
	# species_tree = ape::unroot(species_tree)
	
	stopifnot( species_tree$Nnode == length( species_tree$node.label ) )
	
	species_tree <- ape::ladderize( species_tree )

	# Update species names
	species_tree$tip.label = sub( '_', ' ', species_tree$tip.label ) %>%  sub( '_', '', . )

	# Make ultrametric tree
	# By default, root has a depth of 1
	
	species_ultrametric <- ape::chronos( species_tree, lambda=1, model="correlated", quiet=TRUE )
	class( species_ultrametric ) <- "phylo"
	
	# Subsample the species tree to only those for which expression data are available
	species_tree_focal <- species_ultrametric
	tips <- species_tree_focal$tip.label
	
	for( tip in tips ) {
		if ( ! tip %in% focal_species ) {
			species_tree_focal = ape::drop.tip( species_tree_focal, tip )
		}
	}

	# Now grab the phyldog species tree that shows how node 
	# numbering in gene trees corresponds to gene tree nodes
#	phyldog_species_numbered_tree = read.tree(phyldog_species_numbered_tree)
	
	# Parse the species node numbers from tips and clean up tip names
	
	S_tips <- seq(1:length(phyldog_species_numbered_tree$tip.label))
	
	phyldog_species_numbered_tree$tip.label <- 
		phyldog_species_numbered_tree$tip.label %>% sub( '_', ' ', .) %>%  sub( '_', '', . )
	
	# Make sure tip labels are the same
	stopifnot( 
		setequal( phyldog_species_numbered_tree$tip.label, species_tree$tip.label ) 
	)
	
	stopifnot(
		ape::all.equal.phylo( phyldog_species_numbered_tree, species_tree,  use.edge.length=FALSE )
	)
	
	# Get the phyldog species tree node numbers, these are specified by phyldog for
	# tips and internal nodes on the species tree and are different than ape phylo 
	# node numbers
	S_raw <- c( S_tips, phyldog_species_numbered_tree$node.label )
	
	# Get the correspondence of nodes between the trees
	node_correspondence <- all.equal( phyldog_species_numbered_tree, species_tree, use.edge.length=FALSE, index.return=TRUE )
	colnames( node_correspondence ) <- c( "species", "numbered" )
	node_correspondence <- node_correspondence[ order( node_correspondence[,1] ), ]
	
	# Some rows duplicated in all.equal result, remove them
	node_correspondence <- node_correspondence[ ! duplicated( node_correspondence[ ,1 ] ), ]
	
	# Check that the tips correspond correctly
	ntips <- length( species_tree$tip.label )
	stopifnot(
		all( 
			species_tree$tip.label == 
		 phyldog_species_numbered_tree$tip.label[ node_correspondence[1:ntips, 2] ] 
		)
	)
	
	# Order S to correspond to species tree
	S <- S_raw[ node_correspondence[, 2] ]
	
	numbered_tree <-phyldog_species_numbered_tree
	
	numbered_tree$tip.label<-S_tips
	
	node_comparison <- tibble(species=c(phyldog_species_numbered_tree$tip.label, phyldog_species_numbered_tree$node.label), 
	                              numbered =c(numbered_tree$tip.label, numbered_tree$node.label))
	
	#it's important node_comparison is a tibble. dataframes behave strangely.
	
	## Read in gene trees
	
	gene_trees <- lapply(gene_trees_paths, read_gene_trees) 
  names(gene_trees) <- gene_trees_names
	
	# Get vector with values of age of each internal node and names S
	get_age <- function ( node, tree ) {
		distance_matrix = ape::dist.nodes( tree )
		tips = hutan::tip_descendants( tree, node )
		distances = distance_matrix[ node, tips ]
		
		# Make sure the distance from the node to each tip is the same
		stopifnot( near( max( distances ), min( distances ) ) )
		
		return( distances[ 1 ] )
	}
	
	Nnodes<- length( species_ultrametric$tip.label ) + species_ultrametric$Nnode
	node_ages<-sapply( 1:Nnodes, get_age, tree=species_ultrametric )
	
	names( node_ages ) <- S

	get_edge_summary_by_tree = function( edges ) {

		edge_summary_by_tree = edges %>% 
			dplyr::group_by( gene_tree ) %>%
			dplyr::summarise( 
					n_edges = n(), 
					max_length = max( length ),
					sd_length = sd( length ),
					mean_length = mean( length ),
					fraction_default_length = mean( default_length ),
					n_not_default = sum( ! default_length )
			)

		return( edge_summary_by_tree )

	}

	raw_gene_trees_n = length( gene_trees )

	# Get edge summary stats
	edges =  
		parallel::mclapply( 
			gene_trees, 
			summarize_edges_new, 
			default_length_val=default_length_val, 
			mc.cores=cores 
		) %>% 
		dplyr::bind_rows()
	
	gene_tree_hashes = lapply(gene_trees, digest::digest) %>% unlist()
	edge_summary_by_tree = get_edge_summary_by_tree( edges )

	# Filter by specified thresholds
	exclude_hashes = edge_summary_by_tree$gene_tree[
		( edge_summary_by_tree$max_length > edge_length_max_threshold ) |
		( edge_summary_by_tree$fraction_default_length > fraction_default_max_threshold )
	]

	gene_trees = gene_trees[ ! ( gene_tree_hashes %in% exclude_hashes ) ]

	# Recalculate summary statistics on subsampled trees
	edges =  
		parallel::mclapply( 
			gene_trees, 
			summarize_edges_new, 
			default_length_val = default_length_val, 
			mc.cores=cores 
		) %>% 
		dplyr::bind_rows()

	gene_tree_hashes = lapply( gene_trees, digest ) %>% unlist()
	edge_summary_by_tree <- get_edge_summary_by_tree( edges )

	# Get a data frame of summary node statistics for all gene trees
	nodes_raw <- 
	  parallel::mclapply( 
	  	gene_trees, 
	  	agalmar::summarize_nodes, 
	  	default_length_val = default_length_val, 
	  	mc.cores=cores 
	  ) %>% 
	  dplyr::bind_rows()
	

	# Summarize number of species nodes, internal nodes and species nodes that don't match the species tree
	
	node_summary_by_tree<-nodes_raw %>% dplyr::group_by(gene_tree) %>% dplyr::summarize(nodes=n(), internal_nodes=length(node_depth[node_depth!=1]), n_x=length(S[S=="XX"]), xperinternalnode=n_x/internal_nodes)

	exclude_nodes = node_summary_by_tree$gene_tree[node_summary_by_tree$xperinternalnode > max_x_per_node]
	
	gene_trees = gene_trees[ ! ( gene_tree_hashes %in% exclude_nodes ) ]
	
	###########I used this data and resulting plots to determine which threshold to apply.##########################
	
	#threshold=c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9)
	#nodes_retained=sapply(threshold,function(threshold){
	#  return(sum(node_summary_by_tree[node_summary_by_tree$xperinternalnode<threshold,]$internal_nodes))
	#})
	#trees_retained=sapply(threshold,function(threshold){
	#  return(length(node_summary_by_tree[node_summary_by_tree$xperinternalnode<threshold,]$internal_nodes))
	#})
	
	#threshold_data=tibble(
	#  threshold,
	#  nodes_retained,
	#  trees_retained
	#)

	# threshold_data %>% ggplot()+geom_point(aes(threshold, nodes_retained))
	# threshold_data %>% ggplot()+geom_point(aes(threshold, nodes_retained))
	
#Generate dataframe with FC expression data

focal_treatments<-c("Gastrozooid developing","Palpons mature","Gastrozooid mature","Nectophore developing","Pneumatophore","Bract developing","Gonodendron male", "Gonodendron female")

pairs<-combn(focal_treatments,2)

dge_results<-lapply(dge,pairsdge)

#change names in Bargmannia
names(dge_results$`HWI-ST625-51-C02UNACXX-8-BARGMANNIA`[[2]])<-sub("Gaswhimat", "Gasmat",names(dge_results$`HWI-ST625-51-C02UNACXX-8-BARGMANNIA`[[2]])) 
names(dge_results$`HWI-ST625-51-C02UNACXX-8-BARGMANNIA`[[14]])<-sub("Gaswhimat", "Gasmat",names(dge_results$`HWI-ST625-51-C02UNACXX-8-BARGMANNIA`[[14]])) 
names(dge_results$`HWI-ST625-51-C02UNACXX-8-BARGMANNIA`[[15]])<-sub("Gaswhimat", "Gasmat",names(dge_results$`HWI-ST625-51-C02UNACXX-8-BARGMANNIA`[[15]])) 
names(dge_results$`HWI-ST625-51-C02UNACXX-8-BARGMANNIA`[[17]])<-sub("Gaswhimat", "Gasmat",names(dge_results$`HWI-ST625-51-C02UNACXX-8-BARGMANNIA`[[17]]))

dge_results_sorted<-sort_by_pair(dge_results)

tpm_results<-lapply(e,tpmtreat)

names(tpm_results$`HWI-ST625-51-C02UNACXX-8-BARGMANNIA`[[3]])<-sub("Gaswhimat", "Gasmat",names(tpm_results$`HWI-ST625-51-C02UNACXX-8-BARGMANNIA`[[3]]))

#sort tpm results by zooid type and bind rows
tpm_results_sorted<-sort_tpm(tpm_results)

#add the blast hits
Fulldataframe <-dplyr::left_join(nodes_raw,dge_results_sorted[[1]], by="sequence_ids")

for(i in 2:length(dge_results_sorted)){
Fulldataframe <-dplyr::left_join(Fulldataframe,dge_results_sorted[[i]], by="sequence_ids")
}

for(i in 1:length(tpm_results_sorted)){
Fulldataframe <-dplyr::left_join(Fulldataframe,tpm_results_sorted[[i]], by="sequence_ids")
}

Fulldataframe<-Fulldataframe[!is.na(Fulldataframe$sequence_ids),]

#Subset to remove extra data from nodes_raw

Fulldataframe<-Fulldataframe[,c(8,12:ncol(Fulldataframe))]

blast_lookup$sequence_ids<-as.numeric(blast_lookup$sequence_ids)
Fulldataframe<-dplyr::left_join(Fulldataframe,blast_lookup)

#if the treatment is present & was sampled, but the tpm expression value is NA, this should be set to 0 as gene trees will be pruned to remove NA values. Thus affecting the topology of the tree. 

#remove second instance of "sequence_ids"

stopifnot(all(ncol(Fulldataframe)==(length(dge_results_sorted)*6+length(tpm_results_sorted)+2)))


# Parse expression data into @data
	gene_trees_annotated <- parallel::mclapply( 
		gene_trees,
	  add_de_to_nhx,
		combined_de=Fulldataframe, 
		mc.cores=cores
	)

#Create a vector of the number of tips with expression data for each tree
	n_expression_tips <- 
		parallel::mclapply( 
			gene_trees_annotated, 
			function( x ) sum(rowSums(!is.na(x@data[,12:ncol(x@data)])) >0 ), 
			mc.cores=cores 
		) %>%
	  unlist()
	
#Retain only those trees with the minimum number of tips with expression data
	gene_trees_annotated <- gene_trees_annotated[ n_expression_tips >= min_tips ]
	
#Calibrate the gene tree branch lengths to the species tree branch lengths
#Some of these calibrations fail, need to accommodate these failures in 
#later steps
	gene_trees_calibrated <- 
		parallel::mclapply( 
			gene_trees_annotated, 
			function( x ) clone_edge_lengths( x, node_ages=node_ages ), 
			mc.cores=cores
		)
	
#Add the calibrated node ages to the @data
	gene_trees_calibrated <- lapply( gene_trees_calibrated, store_node_age )
	
	# Exclude trees with suspiciously deep roots, which can be indicative of calibration problems
  # Trees that exceed this are replaced with NA to keep indices the same
	gene_trees_calibrated <- 
		parallel::mclapply( 
			gene_trees_calibrated, 
			function( nhx ) {
				if ( class( nhx ) != "treedata" ) {
					return( nhx )
				}
				
				if ( max( nhx@data$node_age ) > max_root_depth ) {
					return( NA )
				}
				else{
					return( nhx )
				}
			}, 
			mc.cores=cores
		)

#Add column with original rownames
 gene_trees_calibrated <- lapply(gene_trees_calibrated,add_node_name)
 
 save.image("Manuscript_prep.RData")
 
#This will now be used by ancestral_trait_recon.R