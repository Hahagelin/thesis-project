# -----------------------------------------------------------------------------#
#' Antibody Specific Probabilistic Quotient Normaliation
#' 
#' Normalize the data applying a variation of Probabilistic Quotient method
#' introduced in Dieterle et al.(2006), called AbsPQN.
#' 
#' @param X a \code{\link{matrix}} or a \code{\link{BAf-class}} object
#' @param ... not used 
#' 
#' @return \code{\link{matrix}} (or \code{\link{BAf-class}} object) after the
#'   normalization
#' 
#' @references 
#' Dieterle et al. (2006) Probabilistic quotient normalization.. - Anal.Chem
#' 
#' @seealso
#' \code{\link{apply_per_group}}
#' \code{\link{norm_PQN}}
#' 
#' @author Mun-Gwan Hong <\email{mun-gwan.hong@scilifelab.se}>
#' @rdname abspqn
#' @export
# -----------------------------------------------------------------------------#
# created  : 2017-10-02 by Mun-Gwan
# modified : 
# -----------------------------------------------------------------------------#
setGeneric("abspqn", function(X, ...) standardGeneric("abspqn"));

# -----------------------------------------------------------------------------#
#' @rdname abspqn
#' @export
# -----------------------------------------------------------------------------#
setMethod("abspqn", signature(X = "matrix"), function(X, ...) {
  
  # Stop when any infinite value included
  stopifnot(all(is.finite(X[!is.na(X)])))
  
  # reference per target
  ref <- apply(X, 2, median, na.rm = T)
  # normalized for each target
  quo <- t( t(X)/ref )
  
  # probablistic quotient : normalizing factor for each sample
  pq <- apply(quo, 1, median, na.rm = T)

  # any sample of which the median value is 0
  pq[pq == 0] <- NA_real_
  
  # if NAs were given to samples not to individual data points
  use <- apply(X, 2, function(x) which(is.na(x))) %>% 
    sapply(., function(y) identical(.[[1]], y)) %>% {
      if(all(.)) "complete.obs" else "pairwise.complete.obs"
    }

  # correlation with the 'pq'
  cor_pq <- cbind(pq, X) %>% 
    cor(method= "spearman", use= use) %>% 
    .[-1, 1]    # skip the 'pq'
  
  # Antibody specific PQ
  as_pq <- ((pq - 1) %*% t(cor_pq)) + 1
  
  out <- X / as_pq
  
  attr(out, "pq") <- pq
  attr(out, "cor_pq") <- cor_pq
  
  return(out)
})

# -----------------------------------------------------------------------------#
#' @param by_s,by_b same as \code{by_s} and \code{by_s} in
#'   \code{\link{apply_per_group}}. If any of these are given the \code{pqn} is
#'   applied per group divided by these variables.
#' 
#' @examples
#' data(sba)
#' sba2 <- abspqn(sba[sba@sinfo$cohort != "EMPTY", ])
#' plot_QC_sample_signal_boxplot(sba2)
#' 
#' @rdname abspqn
#' @export
# -----------------------------------------------------------------------------#

setMethod(
  "abspqn", signature(X = "BAf"), 
  function(X,
           by_s = NULL,
           by_b = NULL,
           ...) {
    
    apply_per_group(X, function(EA) {
      EAn <- EA
      NM <- abspqn(sX(EA))
      sX(EAn) <- NM
      ## compute PQN denominator and ASF. Store them in @assy_s and @assy_b
      bbat <- batch(EA, "binder") %>% unique()
      sA(EAn, "sinfo")[["pqn_denominator"]] <-
        sapply(bbat, function(x) attr(NM, "pq")) %>% 
        as.data.frame() %>% 
        `colnames<-`(bbat)
      
      sbat <- batch(EA, "sinfo") %>% unique()
      sA(EAn, "binder")[["corr_pq"]] <- 
        sapply(sbat, function(x) attr(NM, "cor_pq")) %>% 
        as.data.frame() %>% 
        `colnames<-`(sbat)
      return(EAn)
    }, by_s= by_s, by_b= by_b, passBAf= T)
  })

