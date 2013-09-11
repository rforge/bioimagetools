### R code from vignette source 'EBImage-introduction.Rnw'

###################################################
### code chunk number 1: library
###################################################
library("EBImage")


###################################################
### code chunk number 2: display-hack
###################################################
display = function(...) if (interactive()) EBImage::display(...)


###################################################
### code chunk number 3: readImage1
###################################################
f = system.file("images", "lena.png", package="EBImage")
lena = readImage(f)


###################################################
### code chunk number 4: display
###################################################
display(lena)


###################################################
### code chunk number 5: readImage2
###################################################
lenac = readImage(system.file("images", "lena-color.png", package="EBImage"))
display(lenac)
nuc = readImage(system.file('images', 'nuclei.tif', package='EBImage'))
display(nuc)


###################################################
### code chunk number 6: readImage2h
###################################################
writeImage(nuc, 'nuc.jpeg', quality=85)


###################################################
### code chunk number 7: writeImage
###################################################
writeImage(lena,  'lena.jpeg', quality=85)
writeImage(lenac, 'lenac.jpeg', quality=85)


###################################################
### code chunk number 8: print
###################################################
print(lena)


###################################################
### code chunk number 9: math1
###################################################
lena1 = lena+0.5
lena2 = 3*lena
lena3 = (0.2+lena)^3


###################################################
### code chunk number 10: math1h
###################################################
writeImage(lena1, 'lena1.jpeg', quality=85)
writeImage(lena2, 'lena2.jpeg', quality=85)
writeImage(lena3, 'lena3.jpeg', quality=85)


###################################################
### code chunk number 11: math2
###################################################
lena4 = lena[299:376, 224:301]
lena5 = lena>0.5
lena6 = t(lena)
print(median(lena))


###################################################
### code chunk number 12: math2h
###################################################
writeImage(lena4, 'lena4.jpeg', quality=85)
writeImage(lena5, 'lena5.png')
writeImage(lena6, 'lena6.jpeg', quality=85)


###################################################
### code chunk number 13: combine
###################################################
lenacomb = combine(lena, lena*2, lena*3, lena*4)
display(lenacomb)


###################################################
### code chunk number 14: combineh
###################################################
writeImage(lenacomb, 'lenacomb.jpeg', quality=85)


###################################################
### code chunk number 15: spatial
###################################################
lena7 = rotate(lena, 30)
lena8 = translate(lena, c(40, 70))
lena9 = flip(lena)


###################################################
### code chunk number 16: spatialh
###################################################
writeImage(lena7, 'lena7.jpeg', quality=85)
writeImage(lena8, 'lena8.jpeg', quality=85)
writeImage(lena9, 'lena9.jpeg', quality=85)


###################################################
### code chunk number 17: print
###################################################
print(lenac)


###################################################
### code chunk number 18: colorMode
###################################################
colorMode(lenac) = Grayscale
display(lenac)


###################################################
### code chunk number 19: colorModeh
###################################################
writeImage(lenac, 'lenac.jpeg', quality=85)


###################################################
### code chunk number 20: colorMode2
###################################################
colorMode(lenac) = Color


###################################################
### code chunk number 21: channel
###################################################
lenak = channel(lena, 'rgb')
lenak[236:276, 106:146, 1] = 1
lenak[236:276, 156:196, 2] = 1
lenak[236:276, 206:246, 3] = 1
lenab = rgbImage(red=lena, green=flip(lena), blue=flop(lena))


###################################################
### code chunk number 22: channelh
###################################################
writeImage(lenak, 'lenak.jpeg', quality=85)
writeImage(lenab, 'lenab.jpeg', quality=85)


###################################################
### code chunk number 23: filter
###################################################
flo = makeBrush(21, shape='disc', step=FALSE)^2
flo = flo/sum(flo)
lenaflo = filter2(lenac, flo)

fhi =  matrix(1, nc=3, nr=3)
fhi[2,2] = -8
lenafhi = filter2(lenac, fhi)


###################################################
### code chunk number 24: filterh
###################################################
writeImage(lenaflo, 'lenaflo.jpeg', quality=85)
writeImage(lenafhi, 'lenafhi.jpeg', quality=85)


###################################################
### code chunk number 25: morpho
###################################################
ei = readImage(system.file('images', 'shapes.png', package='EBImage'))
ei = ei[110:512,1:130]
display(ei)

kern = makeBrush(5, shape='diamond')
eierode = erode(ei, kern)
eidilat = dilate(ei, kern)


###################################################
### code chunk number 26: morphoh
###################################################
writeImage(ei, 'ei.png')
writeImage(kern, 'kern.png')
writeImage(eierode, 'eierode.png')
writeImage(eidilat, 'eidilat.png')


###################################################
### code chunk number 27: segmentation
###################################################
eilabel = bwlabel(ei)
cat('Number of objects=', max(eilabel),'\n')

nuct = nuc[,,1]>0.2
nuclabel = bwlabel(nuct)
cat('Number of nuclei=', max(nuclabel),'\n')


###################################################
### code chunk number 28: segmentationh
###################################################
writeImage(eilabel/max(eilabel), 'eilabel.png')
writeImage(nuclabel/max(nuclabel), 'nuclabel.png')


###################################################
### code chunk number 29: segmentation2
###################################################
nuct2 =  thresh(nuc[,,1], w=10, h=10, offset=0.05)
kern = makeBrush(5, shape='disc')
nuct2 = dilate(erode(nuct2, kern), kern)
nuclabel2 = bwlabel(nuct2)
cat('Number of nuclei=', max(nuclabel2),'\n')


###################################################
### code chunk number 30: segmentation2h
###################################################
writeImage(nuclabel2/max(nuclabel2), 'nuclabel2.png')


###################################################
### code chunk number 31: manip
###################################################
nucgray = channel(nuc[,,1], 'rgb')
nuch1 = paintObjects(nuclabel2, nucgray, col='#ff00ff')
nuclabel3 = fillHull(nuclabel2)
nuch2 = paintObjects(nuclabel3, nucgray, col='#ff00ff')


###################################################
### code chunk number 32: maniph
###################################################
writeImage(nuch1, 'nuch1.jpeg', quality=85)
writeImage(nuch2, 'nuch2.jpeg', quality=85)


###################################################
### code chunk number 33: manip2
###################################################
xy = computeFeatures.moment(nuclabel3)[, c("m.cx", "m.cy")]
xy[1:4,]


###################################################
### code chunk number 34: cs1
###################################################
nuc = readImage(system.file('images', 'nuclei.tif', package='EBImage'))
cel = readImage(system.file('images', 'cells.tif', package='EBImage'))
img = rgbImage(green=1.5*cel, blue=nuc)


###################################################
### code chunk number 35: cs1h
###################################################
writeImage(cel, 'cel.jpeg', quality=85)
writeImage(img, 'img.jpeg', quality=85)


###################################################
### code chunk number 36: cs2
###################################################
nmask = thresh(nuc, w=10, h=10, offset=0.05)
nmask = opening(nmask, makeBrush(5, shape='disc'))
nmask = fillHull(nmask)
nmask = bwlabel(nmask)


###################################################
### code chunk number 37: cs2h
###################################################
writeImage(nmask/max(nmask), 'nmask.png')


###################################################
### code chunk number 38: cs3
###################################################
ctmask = opening(cel>0.1, makeBrush(5, shape='disc'))
cmask = propagate(cel, seeds=nmask, mask=ctmask)


###################################################
### code chunk number 39: cs3h
###################################################
writeImage(cmask/max(cmask), 'cmask.png')


###################################################
### code chunk number 40: cs4
###################################################
res = paintObjects(cmask, img, col='#ff00ff')
res = paintObjects(nmask, res, col='#ffff00')


###################################################
### code chunk number 41: cs4h
###################################################
writeImage(res, 'res.jpeg', quality=85)


