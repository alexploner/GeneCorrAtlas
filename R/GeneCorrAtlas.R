#' Display the Atlas of Genetic Correlations 
#'
#' Bulik-Sullivan et al. recently (2105) published an atlas of genetic
#' correlations for 50 different traits. This package wraps the more than 1000
#' correlation coefficients involved and offers some convenient functions to
#' extract and display the data
#'
#'
#' @source \url{http://www.nature.com/ng/journal/vaop/ncurrent/extref/ng.3406-S2.csv}
#' @references Bulik-Sullivan B, Finucane HK, Anttila V, Gusev A, Day FR, Loh P-R,
#'    et al. An atlas of genetic correlations across human diseases and traits.
#'    Nat Genet [Internet]. 2015 Sep 28 [cited 2015 Oct 5];advance online publication
#'
#' @docType package
#' @name GeneCorrAtlas
NULL


#' Atlas of genetic correlations
#'
#' A data frame containing the genetic correlations between 50 different traits
#' and diseases presented in Bulik-Sullivan et al. (2015).
#'
#' \itemize{
#'   \item Trait1, name of first phenotype
#'   \item Trait2, name of second phenotype
#'   \item rg, the estimated correlation coefficient
#'   \item se, the standard error of the estimate
#'   \item z, the corresponding Wald statistic
#'   \item p, the corresponding p-value
#' }
#'
#' @format A data frame with 1176 rows and 6 variables. The correlations are
#' given in the order as the appear in the lower triangle of the full correlation
#' matrix, by column. 
#' @source \url{https://raw.githubusercontent.com/bulik/gencor/master/all_rg.csv}
#' @references Bulik-Sullivan B, Finucane HK, Anttila V, Gusev A, Day FR, Loh P-R,
#'    et al. An atlas of genetic correlations across human diseases and traits.
#'    Nat Genet [Internet]. 2015 Sep 28 [cited 2015 Oct 5];advance online publication
#' @name gcatlas
NULL


#' Read Atlas of Genetic Correlations
#'
#' Read the supplementary file from disk or from the web. This is the function
#' originally used to generate the pre-packaged data object \code{\link{gcatlas}}.
#' Only use this function if you want to roll your own.
#'
#' @param x name of the supplementary file (by default, a link to the github repository)
#' @param ... extra arguments to \code{read.csv}, in case of weirdness
#'   
#'
#' @return The data frame described at \code{\link{gcatlas}}
#'
#' @source \url{https://raw.githubusercontent.com/bulik/gencor/master/all_rg.csv}
#' @references Bulik-Sullivan B, Finucane HK, Anttila V, Gusev A, Day FR, Loh P-R,
#'    et al. An atlas of genetic correlations across human diseases and traits.
#'    Nat Genet [Internet]. 2015 Sep 28 [cited 2015 Oct 5];advance online publication
read_gcatlas = function(x, ...)
{
	if (missing(x)) x = "https://raw.githubusercontent.com/bulik/gencor/master/all_rg.csv"
	read.csv(x, ...)
}	

#' Extract trait names from the atlas of genetic correlations
#'
#' A simple helper function that returns a sorted list of traits found
#' in a \code{gcatlas}-type of object
#'
#' @param data the object from which to extract trait names (default: 
#'             the full atlas)
#'
#' @return A character vector of trait names, in the order they appear
#'         in \code{Trait1}/\code{Trait2} in \code{data}. This is \emph{not}
#'         sorted, only almost sorted
#'
#' @seealso \code{\link{gcatlas}}
#'
#' @examples
#' traits()
traits = function(data=gcatlas) 
{
	# DON NOT TRY TO OPTIMIZE THIS BIT; YOU *WILL* FUCK UP THE ORDER
	# OF COLUMNS IN gcmatrix
	unique(c(as.character(data$Trait1), as.character(data$Trait2)))
}	

#' Convert Atlas of Genetic Correlations to a matrix
#'
#' This function extracts any of the variables in the \code{gcatlas} 
#' object and returns it as a square matrix, with traits as rows and columns.
#'
#' @param data The object from which to extract the variable - by default, 
#'             the full \code{gcatlas} object, but a valid subset will 
#'             work, too
#' @param type The variable to extract, one of \code{rg} (default), \code{se}
#'             \code{z} and \code{p}
#' @param diag_value An optional value for the diagonal of the matrix; 
#'                   if missing, the function will use a reasonable 
#'                   default, based on the value of \code{type}
#'
#' @return A square numerical matrix
#'
#' @seealso \code{\link{gcatlas}}
gcmatrix = function(data=gcatlas, type=c("rg", "se", "z", "p"), diag_value)
{
	type = match.arg(type)
	if (missing(diag_value)) {
		diag_value = if (type=="rg") 1 
		             else if (type=="se") 0
		             else NA
	}
	
	traits = traits(data)
	n = length(traits)
	m = matrix(0, nrow=n, ncol=n)
	m[lower.tri(m)] = data[, type, drop=TRUE]
	m = m + t(m)
	diag(m) = diag_value
	colnames(m) = rownames(m) = traits
	m
}	

#' Match trait names to a data set
#'
#' This is a utility function that takes a vector of possibly abbreviated trait
#' names and returns the matching trait names from a data frame with the same
#' format as \code{\link{gcatlas}} (which is the default).
#'
#' @param x a vector of (possibly abbreviated) trait names
#' @param data the data frame against which the trait names are compared
#' @param unique Logical switch indicating whether a unique match is
#'               required; if \code{TRUE}, thorws an error if \code{x}
#'               cannot be matched, or matches more than one trait
#'
#' @seealso \code{\link{gcatlas}}
match_trait_names = function(x, data = gcatlas, unique=FALSE)
{
	refs = traits(data)
	## Trickery: none of the default matching functions does a list
	## of patterns and multiple matches
	ndx = lapply(x, grep, refs)
	len = sapply(ndx, length)
	isna = len == 0
	if (any(isna)) {
		warning("Cannot match ", x[isna], " - excluding mismatch(es)")
		ndx = ndx[!isna]
	}
	ndx = unlist(ndx)
	if (unique) {
		nn = length(ndx)
		if (nn == 0) {
			stop("Cannot match the trait name")
		} else if (nn > 1) {
			stop("Trait name must be a unique match")
		}
	}

	## That's all right then
	refs[ndx]
}	


#' Select traits from the full atlas
#'
#' Subsets a data frame with the same format as \code{\link{gcatlas}} to
#' include only traits from a specified list of trait names. These names can be
#' abbreviated.
#'
#' @param ... a list of trait names (possibly abbreviated)
#' @param data the name of the object to subset; by default the full data set
#' @param drop logical expression indicating whether to return the specified traits
#' (\code{FALSE}, default) or to exclude them from the data (\code{TRUE}).
#'
#' @return An object of the same structure as \code{gcatlas}
#'
#' @seealso \code{\link{gcatlas}}
#'
#' @examples
#' x = sel_trt("ADHD", "Alz", "T2D")
#' x
#' gcmatrix(x)
sel_trt = function(..., data=gcatlas, drop=FALSE)
{
	args = list(...)
	if (length(args) == 0) stop("Need to specify at least one trait name")
	
	full = unlist(sapply(args, match_trait_names, data=data))
	if (drop) {
		full = setdiff(traits(data=data), full)
	}
	
	droplevels(subset(data, (Trait1 %in% full) & (Trait2 %in% full)))
}

#' Convert a correlation matrix to a distance matrix
#'
#' Two simple hepler functions that convert a matrix of correlations to
#' a distance/dissimilarity matrix
#'
#' @param x A square matrix of correlations
#'
#' @return An object of class \code{dist}: \code{cor2dist} returns \eqn{1-r},
#'         \code{cor2dist2} returns \eqn{1-|r|}. 
#'
#' @seealso \code{\link{stats::dist}}, \code{\link{gcheatmap}}
cor2dist = function(x)
{
	as.dist(1 - x)
}
#' @rdname cor2dist
cor2dist_abs = function(x)
{
	as.dist(1 - abs(x))
}
	
#' Heatmap of genetic correlations
#' 
#' Plot a heatmap of genetic correlations
#'
#' @param data An object of the same structure as \code{gcatlas}
#' @param distfun The function to convert a correlation matrix to a 
#'                distance matrix for generating the dendrograms
#' @param symm Logical switch indicating a heatmap of a symmetric matrix; 
#'        passed on to \code{heatmap} and probably best left unchanged
#' @param margins Numerical vector of length two, indicating how much 
#'        space is reserved for the row/column labels. This is the same
#'        argument as for \code{heatmap}, but with a different default
#'        so that the trait names are reasonably readable with normal
#'        device size 
#' @param fix Logical expression indicating whether correlations outside 
#'            of [-1,1] should be set to the limits
#' @param ... Extra arguments passed on to \code{heatmap}
#'
#' @return Invisibly, the same as function \code{heatmap}. This function
#'         is generally invoked to generate a plot.
#'
#' @seealso \code{\link{heatmap}}, \code{\link{cor2dist}}, \code{\link{gcatlas}}
#'
#' @examples
#' ## Default heatmap of the full atlas
#' gcheatmap()	
gcheatmap = function(data=gcatlas, distfun=cor2dist, symm=TRUE, margins=c(9,9), fix=TRUE,...)
{
	rr = gcmatrix(data)
    if (fix) {
        rr[rr>1]  = 1
        rr[rr< (-1)] = -1
    }
    
    heatmap(rr, distfun=distfun, symm=symm, margins=margins, ...)
}

#' Plot all correlations with a specific trait
#'
#' Given the atlas of genetic correlations or a valid subset, plus a trait
#' name, this function plots all correlations of the specified trait and 
#' all other traits in the correlation set. The plotting data (including
#' confidence intervals) is returned invisibly.
#'
#' @param trait The name of the trait for which to plot the correlations;
#'              can be abbreviated, but must be unique.
#' @param data Either the full atlas of genetic correlations (\code{gcatlas},
#'             the default), or a valid subset, from which to draw correlations.
#' @param conf The confidence level for the confidence intervals
#' @param sort Logical switch indicating whether to sort the traits by 
#'             by correlation (default: \code{TRUE})
#' @param rlim Cutoff values for plotting the confidence limits, so that
#'             wide confidence limits do not dominate the plot
#'
#' @return The plot as an ggplot2 object. Note that the plotting data 
#'         can be easily extracted for further processing, see Examples.
#' 
#' @seealso \code{\link{gcatlas}}, \code{\link{sel_trt}}
#'
#' @examples
#' ## Default invocation
#' plotTrait("ADHD")
#' 
#' ## Subset of the full atlas
#' plotTrait("Alz", sel_trt("Obesity"))
#' 
#' ## Extract the plotting data
#' x = plotTrait("Cor")
#' x
#' ## Useful extra attributes
#' attributes(x)
plotTrait = function(trait, data = gcatlas, conf=0.95, sort=TRUE, rlim=c(-0.5, 1.0))
{
	## Check input
	tr = match_trait_names(trait, data=data, unique = TRUE)
	if ((conf <= 0) | (conf >= 1)) stop("confidence levels must be in (0,1)")
	
	## Extract correlations and standard errors
	rr = gcmatrix(data = data)
	se = gcmatrix(type = "se", data = data)
		
	## Cut the data down to the desired column
	ndx = match(tr, rownames(rr))
	nam = rownames(rr)[-ndx]
	rr = rr[-ndx, ndx]
	se = se[-ndx, ndx]
	
	## Do the confidence limits
	quant = abs( qnorm((1-conf)/2) )
	lc = rr - quant*se
	uc = rr + quant*se
	
	## Put it all into a data frame (with extra information)
	dat = data.frame(Trait=nam, R=rr, LCL=lc, UCL=uc)
	rownames(dat) = NULL
	attr(dat, "RefTrait")  = tr
	attr(dat, "ConfLevel") = conf
	
	## The plot
	pl = ggplot(dat, aes(x = if(sort) stats::reorder(Trait, R) else Trait, 
	                     y=R, ymin=pmax(LCL, min(rlim)), ymax=pmin(UCL, max(rlim))), 
	                     environment=environment())
    pl = pl + geom_hline(yintercept=0, size=1, col="grey70")
    pl = pl + geom_pointrange() 
    pl = pl + labs(x="Traits", y="Genetic correlation", title=tr)
    pl = pl + coord_flip()
    pl = pl + cowplot::theme_cowplot()
    pl = pl + theme(axis.text.y = element_text(size=12))
    pl
}

#' Scatterplot of genetic correlations for two traits
#'
#' Given the atlas of genetic correlations or a valid subset, plus two trait
#' names, this function produces a scatterplot of the genetic correlations
#' of one trait against the genetic correlations of the other. 
#'
#' @param xtr,ytr The names of the traits to be displayed on the x- and 
#'                y-axis; can be abbreviated, but must be unique.
#' @param data Either the full atlas of genetic correlations (\code{gcatlas},
#'             the default), or a valid subset, from which to draw correlations.
#' @param co The cutoff under which to consider p-values to be
#'           statistically significant
#'
#' @return The plot as an ggplot2 object. Note that the plotting data 
#'         can be easily extracted for further processing, see Examples.
#' 
#' @seealso \code{\link{gcatlas}}, \code{\link{sel_trt}}
#'
#' @examples
#' ## Default invocation
#' plotPairedTraits("ADHD", "Alz")
#' 
#' ## Change the cutoff 
#' plotPairedTraits("ADHD", "Alz", co=0.01)
#'
#' ## Extract the data
#' x = plotPairedTraits("Cor", "T2D")
#' x$data
#' ## Useful extra attributes
#' attributes(x)
plotPairedTraits = function(xtr, ytr, data = gcatlas, co = 0.05)
{
	## Check input
	xtr = match_trait_names(xtr, data=data, unique = TRUE)
	ytr = match_trait_names(ytr, data=data, unique = TRUE)	
	if ((co <= 0) | (co >= 1)) stop("cutoff must be in (0,1)")
	
	## Extract correlations and p-values
	rr = gcmatrix(data = data)
	pp = gcmatrix(type = "p", data = data)
		
	## Cut the data down to the desired column
	xndx = match(xtr, rownames(rr))
	yndx = match(ytr, rownames(rr))	
	ndx = c(xndx, yndx)
	nam = rownames(rr)[-ndx]	
	rr  = rr[-ndx, ndx]
	pp  = pp[-ndx, ndx]
	ss  = pp <= co
	ss  = 2*ss[,1] + ss[,2]
	signif = factor(ss, levels=0:3, labels=c("None", ytr, xtr, "Both"))
	dat = data.frame(Trait=nam, R=rr, p=pp, Signif=signif)
	rownames(dat) = NULL
	colnames(dat)[2:3] = paste0("R_", 1:2)
	colnames(dat)[4:5] = paste0("p_", 1:2)	
	attr(dat, "RefTrait")  = c(xtr, ytr)
	attr(dat, "Cutoff") = co

	## The plot
	pl = ggplot(dat, aes(x=R_1, y=R_2, label=Trait, shape=Signif, color=Signif, vjust=-0.4))
	pl = pl + geom_hline(yintercept=0, size=1, col="grey70")
	pl = pl + geom_vline(xintercept=0, size=1, col="grey70")
	pl = pl + geom_point() + geom_text(size=3, show.legend=FALSE)
	pl = pl + labs(x=xtr, y=ytr)
	pl = pl + cowplot::theme_cowplot()
	pl = pl + theme(legend.position="top") 
	pl = pl + scale_colour_manual(name="Statistical significance", values=c("grey70", "blue", "red", "purple"), drop=FALSE)
	pl = pl + scale_shape_manual(name="Statistical significance", values=c(20, 16, 17, 15), drop=FALSE)
	pl 
}	
	

