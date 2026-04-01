#' @title
#' Predictor of Differentiation State for DNAm data
#'
#' @aliases RunEpiDiff
#'
#' @description
#' This function takes as input an Illumina normalized DNAm data matrix
#' with rows annotated to CpGs and returns the probabilistic differentiation
#' state assignments as well as continuous differentiation index (DI) that is
#' theoretically bounded between 0 (pluripotency) and 1 (fully differentiated).
#'
#' @param data.m
#' A normalized DNAm data matrix with rows labeling CpGs (Illumina probes)
#' and columns labeling samples.
#'
#' @param w.v
#' A weight vector to specify how the DI is to be computed. The DI is computed
#' as w_1*p0 + w_2*p1 + w_3*p2, where the probabilities p0, p1 and p2 are
#' the probalities that a sample belongs to the pluripotent, multipotent and
#' differentiated states, respectively. The default weights for 
#' w_1(=0) and w_3(=1) should never be changed. The default value for w_2 is
#' 0.25 which is the optimal choice and should also not be changed.
#'
#' 
#' @return A list containing the following elements
#'
#' @return di
#' A vector containing the DI estimates of the samples for the optimal
#' predictor (ntrees=400). DI-values are bounded between 0 and 1. Typically,
#' values less than 0.25 may indicate pluripotency, values between 0.25 and
#' 0.55 indicate multipotency and higher values indicating increasingly
#' differentiated states.
#' 
#' @return diALL
#' A list of vectors containing the DI estimates of the samples, with each
#' entry in the list labeling a different number of trees (there are 9
#' elements, with number of trees varying from 100 to 500 trees in units of 50)
#'
#' @return class
#' Predicted classes for the optimal predictor (0=pluripotent, 1=multipotent, 2=differentiated)
#'
#' @return classALL
#' Predicted classes for all predictors parameterized by the number of trees.
#'
#' @return prob
#' A matrix containing the estimated probabilities for each sample (row)
#' to belong to the pluripotent, multipotent or differentiated state. Row sums
#' add to 1.
#'
#' @return probALL
#' A list of matrices containing the estimated differentiation state
#' probabilities, one entry for each choice of number of trees.
#' 
#' @references
#' Teschendorff A, Tong H, Guo X.
#' \emph{DNA methylation-based prediction of differentiation state at cell-type resolution} Submitted.
#'
#'
#' @examples
#' data(dataCD34);
#' ediff.o <- doDiffRF(dataCD34.m)
#'
#' @importFrom randomForest predict
#'
#' @export
#' 


DoDiffRF <- function(data.m,w.v=c(0,0.25,1)){
  data("rfTRdiff2"); ### contains candDiffALL2.v and rfTRdiff2.lo
  match(candDiffALL2.v,rownames(data.m)) -> tmp.idx;
  na.idx <- which(is.na(tmp.idx));
  print(paste("There are ",length(na.idx)," missing CpGs out of a total of ",length(candDiffALL2.v)," CpGs",sep=""));
  if(length(na.idx)>0){ ### missing cpgs, so impute and augment matrix
  print("Imputing...");
  aug.m <- matrix(mean(data.m),nrow=length(na.idx),ncol=ncol(data.m));
  rownames(aug.m) <- candDiffALL2.v[na.idx];
  tmp.m <- rbind(data.m,aug.m);
  }
  else if (length(na.idx)==0){
  print("No need to impute");
      tmp.m <- data.m;
  }
  match(candDiffALL2.v,rownames(tmp.m)) -> map.idx;
  predRF.lv <- list();probRF.lv <- list();
  diffIdx.lv <- list();
  for(nti in 1:length(rfTRdiff2.lo)){
    predRF.lv[[nti]] <- as.vector(predict(rfTRdiff2.lo[[nti]],newdata=t(tmp.m[map.idx,])));
    probRF.lv[[nti]] <- predict(rfTRdiff2.lo[[nti]],newdata=t(tmp.m[map.idx,]),type="prob");
    diffIdx.lv[[nti]] <- w.v[1]*probRF.lv[[nti]][,1] + w.v[2]*probRF.lv[[nti]][,2] + w.v[3]*probRF.lv[[nti]][,3];
  }
  
  return(list(di=diffIdx.lv[[7]],diALL=diffIdx.lv,class=predRF.lv[[7]],classALL=predRF.lv,prob=probRF.lv[[7]],probALL=probRF.lv));
}

