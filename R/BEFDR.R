BEFDR <- function(minimal.lm, maximal.lm, FDR.q, mfactor=1){
  #This function remains the same for FS & BE
  compute.Lambda<-function(k, m, Q) {i<-c(1:k)
  return( (1/(k+1)) * sum(qnorm((Q/2) * (i/(m+1-i*(1-Q))))^2))  }

  #This function remains the same for FS & BE
  get.model.size <- function(a.lm) {require(MASS);
    return(extractAIC(a.lm)[1]-1) }

  require(MASS);

  #Scope remains the same for FS & BE
  the.scope <- list(lower = minimal.lm, upper = maximal.lm)
  # Size of maximal.lm multiplied by mfactor
  # as BE procedure is initiated when 1/mfactor of all the genotyped regions are already removed
  # need to correct for all the tests that are actually performed when building
  # the final model
  m <- mfactor*get.model.size(maximal.lm)
  new.model.size <- get.model.size(maximal.lm)
  for (i in 1:m)
  {
    #test i <- 3
    old.model.size <- new.model.size
    Lambda <- compute.Lambda(k = old.model.size-1, m, Q = FDR.q)

    new.model <- stepAIC(maximal.lm, direction="backward", scope=the.scope, k=Lambda, trace=FALSE)
    new.model.size <- get.model.size(new.model);

    if (new.model.size >= old.model.size) break;
  }
  new.lm <- lm(new.model)
  return(new.lm) }
