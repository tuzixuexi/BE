\name{BEFDR}
\alias{BEFDR}
\title{
Backward elimination in multi-locus association analysis
}
\description{
Performes multi-locus backward elimination using adaptive FDR
}
\usage{
BEFDR( minimal.lm, maximal.lm, FDR.q,mfactor=1)
}
\arguments{
  \item{minimal.lm}{
A minimum model returned from lm, stating fixed effects and SNPs that will always be kept
}
 \item{maximal.lm}{
A maximal model returned from lm, stating fixed effects and all SNPs. Note that the number of SNPs can not exceed the number of responses to avoid overfitting. Under overfitting, a chromosome based pre-selction can be performed
}
\item{FDR.q}{
Desired FDR value, usually set to 0.05-0.20.
}
\item{mfactor}{
A factor created to account for the ture model size. When a model selection started on the whole SNP set, set this factor to 1. Under the situation of pre-selection, the number of variants(n2) after pre selection is smaller than the whole set(n1). set m= n1/n2
}
}
\details{
This function perform Backward analysis based on an adaptive FDR. This model selection criterion is flexable to the input model size. Detials can be viewed in Sheng et al.(2015)

}
\value{
Object of class lm
}
\references{

Abramovich, F., Benjamini, Y., & Donoho, D. L. (2006). Special Invited Lecture: Adapting to Unknown Sparsity by Controlling the False Discovery Rate on JSTOR. The Annals of Statistics,Vol34.No.2584-653

Gavrilov, Y., Benjamini, Y., & Sarkar, S. K. (2009). JSTOR: The Annals of Statistics, Vol. 37, No. 2 (Apr., 2009), pp. 619-629. The Annals of Statistics.

Sheng Z, Pettersson ME, Honaker CF, Siegel PB, Carlborg Ö. 2015. Standing genetic variation as a major contributor to adaptation in the Virginia chicken lines selection experiment. Genome Biol 16: 219.

}
\author{
yanjun zan
}
\note{
%%  ~~further notes~~
}


\seealso{
}
\examples{
 data(Chicken) # load required data


 sigMrk_sub<- sigMrk[!(sigMrk \%in\% mrksT)] # this is the fixed effect we want to keep in the minimal model
    # find the id without missing values
    colnames(mrk_geno)<- mrks_info[,1]
    nna<-complete.cases(mrk_geno[,mrks_info[,1] \%in\% c(mrksT,sigMrk_sub)])
    geno_reg<-mrk_geno[nna,c(mrksT,sigMrk_sub)] ## pick out the genotype info for all the included marker

    # construct full model including all markers within the reg
    id.full<- which(colnames(geno_reg) \%in\% c(mrksT,sigMrk_sub))
    geno_add <- paste("as.numeric(geno_reg[,",id.full,"])",sep="",collapse="+")
    test_phe<-phe[nna]
    test_fx1<-fx1[nna]
    fm <- as.formula(paste("test_phe ~ test_fx1 +",geno_add,sep=""))
    reg.lm <- lm(fm, y=TRUE)

    # min model
    id.min<-c(1:ncol(geno_reg))[colnames(geno_reg) \%in\% sigMrk_sub]
    geno_bc.add<-paste("as.numeric(geno_reg[,",id.min,"])",sep="",collapse="+")
    fm.min<- as.formula(paste("test_phe ~ test_fx1 +",geno_bc.add,sep=""))
    min.lm <- lm(fm.min,y=TRUE)

    #Perform Backward-Elimination in the original data at different adaptive FDR thresholds
    fitFDR5 <- BEFDR( minimal.lm = min.lm, maximal.lm = reg.lm, FDR.q = 0.05,mfactor=1)

    ### get out the name of kept snp
    if(sum(grepl("geno_reg",rownames(summary(fitFDR5)$coefficients)))>0){
      terms5<-rownames(summary(fitFDR5)$coefficients)[grep("geno_reg",rownames(summary(fitFDR5)$coefficients))]
      idx5<-unlist(strsplit(terms5,",|[]]"))
      name5<-colnames(geno_reg)[as.numeric(idx5[seq(2,length(idx5),3)])]
      id.in5<- name5[!(name5 \%in\% sigMrk_sub)] # this is the maker pass the FDR
    }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BEFDR}
