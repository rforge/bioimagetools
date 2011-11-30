\name{readBMP}
\alias{readBMP}
\title{Read Bitmap files}
\description{
Function to read 2D grey-value BMP files
}
\usage{
read.bmp(file) }
\arguments{
  \item{file}{A character vector of file names or URLs.}
}
\value{
 Returns a matrix with BMP data as integer.
}
\examples{
bi<-readBMP("http://bioimagetools.r-forge.r-project.org/images/bi.bmp")
image(result,col=grey(seq(1,0,length=100)))
}
\author{Volker J. Schmid \email{volkerschmid@users.sourceforge.net}}


