##
## Copyright (c) 2011 Volker Schmid
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * The names of the authors may not be used to endorse or promote
##       products derived from this software without specific prior
##       written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 

filterImage3d<- function(img, filter="var", window, z.scale=1, silent=FALSE) 
{

dims<-dim(img)
N <- prod(dims)
img<-as.vector(img)
filtered<-rep(0,N)

if (filter=="var")
{
    filtered <- .C("varfilter",
                    as.double(img),
                    as.double(filtered),
                    as.double(c(window,z.scale)),
                    as.integer(dims), 
                    as.integer(filtered),
                    as.double(c(0,0)),
                    as.integer(ifelse(silent,1,0)),
                    PACKAGE="bioimagetools")
    minmax <- filtered[[6]]
    filtered <- filtered[[5]]
}
if (filter=="max")
{
    filtered <- .C("maxfilter",
                    as.double(img),
                    as.double(filtered),
                    as.double(c(window,z.scale)),
                    as.integer(dims), 
                    as.integer(filtered),
                    as.double(c(0,0)),
                   as.integer(ifelse(silent,1,0)),
                   PACKAGE="bioimagetools")
    minmax <- filtered[[6]]
    filtered <- filtered[[5]]
}
if (filter=="min")
{
    filtered <- .C("minfilter",
                    as.double(img),
                    as.double(filtered),
                    as.double(c(window,z.scale)),
                    as.integer(dims), 
                    as.integer(filtered),
                    as.double(c(0,0)),
                   as.integer(ifelse(silent,1,0)),
                   PACKAGE="bioimagetools")
    minmax <- filtered[[6]]
    filtered <- filtered[[5]]
}
filtered<-array(filtered,dims)
filtered<-minmax[1]+(minmax[2]-minmax[1])*filtered/65535
return(filtered)
}	

# 2D convolution-based linear filter for images and matrix data

# Copyright (c) 2007 Gregoire Pau, Oleg Sklyar

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.

# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# See the GNU Lesser General Public License for more details.
# LGPL license wording: http://www.gnu.org/licenses/lgpl.html

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

filterImage2d = function(x, filter) {
  validObject(x)
  validObject(filter)

  dx = dim(x)
  df = dim(filter)

  if (any(df%%2==0)) stop("dimensions of 'filter' matrix must be odd")
  if (any(dx[1:2]<df)) stop("dimensions of 'x' must be bigger than 'filter'")

  ## find centres of x and filter
  cx = dx%/%2
  cf = df%/%2

  ## create fft filter matrix
  wf = matrix(0.0, nrow=dx[1], ncol=dx[2])
  wf[(cx[1]-cf[1]):(cx[1]+cf[1]),(cx[2]-cf[2]):(cx[2]+cf[2])] = filter
  wf = fft(wf)

  ## convert to a frame-based 3D array
  dim(x) = c(dx[1:2],prod(dx)/prod(dx[1:2]))

  index1 = c(cx[1]:dx[1],1:(cx[1]-1))
  index2 = c(cx[2]:dx[2],1:(cx[2]-1))
  pdx = prod(dim(x)[1:2])
  y = apply(x, 3, function(xx) {
    dim(xx) = dx[1:2]
    Re(fft(fft(xx)*wf, inverse=TRUE)/pdx)[index1, index2]
  })
  dim(y) = dx

  ## convert it back
  return(y)

}



