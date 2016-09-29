## ----data01--------------------------------------------------------------
library(GeneCorrAtlas)
head(gcatlas)

## ----data02--------------------------------------------------------------
mygc = read_gcatlas()
all.equal(gcatlas, mygc)

## ----subset01------------------------------------------------------------
traits()

## ----subset02------------------------------------------------------------
sel_trt("Alz", "ADHD", "Bip", "Schizo")

## ----subset03------------------------------------------------------------
shortlist = traits()[1:5]
shortlist
sel_trt(shortlist)

## ----subset04------------------------------------------------------------
sel_trt("BMI", "Height")

## ----subset05------------------------------------------------------------
sel_trt("^BMI$", "^Height$")

## ----subset06------------------------------------------------------------
sub1 = sel_trt("BMI", drop=TRUE)
traits(sub1)
sel_trt("Alz", "Bip", "T2D", data=sub1)

## ----matrix01------------------------------------------------------------
gcmatrix()[1:5, 1:5]

## ----matrix02------------------------------------------------------------
gcmatrix(sel_trt("Depr", "Schiz", "Bipo", "ADHD", "Alz"))

## ----matrix03------------------------------------------------------------
gcmatrix(sel_trt("Depr", "Schiz", "Bipo", "ADHD", "Alz"), type="p")

## ----heatmap01, out.width="750px", dpi=300, fig.width=10, fig.height=10----
gcheatmap()

## ----heatmap02, out.width="500px", dpi=300, fig.width=8, fig.height=8----
gcheatmap(sel_trt(c("BMI", "Obesity", "Overweight", "Height", "Birth")))

## ----plotTrait01, out.width="500px", dpi=300, fig.width=8, fig.height=8----
plotTrait("Cor")

## ----plotTrait02, out.width="500px", dpi=300, fig.width=7, fig.height=7----
plotTrait("Cor", sel_trt("BMI", "Height", "Cor"))

## ----plotPairedTrait01, out.width="500px", dpi=300, fig.width=5, fig.height=5----
plotPairedTraits("Depression", "Bipolar")

## ----plotPairedTrait02, out.width="500px", dpi=300, fig.width=5, fig.height=5----
plotPairedTraits("Depression", "Bipolar", co=0.05)

