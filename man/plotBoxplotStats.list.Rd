%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  plotBoxplotStats.list.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{plotBoxplotStats.list}
\alias{plotBoxplotStats.list}
\alias{plotBoxplotStats.list}


 \title{Plots a (merged or non-merged) list of boxplot.stats() elements}

 \description{
  Plots a (merged or non-merged) list of boxplot.stats() elements.
 }

 \usage{\method{plotBoxplotStats}{list}(stats, ylim=NULL, outline=FALSE, las=2, ...)}

 \arguments{
   \item{stats}{A (merged or non-merged) \code{\link[base]{list}} of
                                  \code{\link[grDevices]{boxplot.stats}} elements.}
   \item{ylim, outline, las, ...}{Arguments passed to \code{\link[graphics]{bxp}}.}
 }

 \value{
   Returns (invisibly) the merged \code{\link[grDevices]{boxplot.stats}} structure.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{methods}
\keyword{internal}
