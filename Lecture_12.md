Lecture 12
================

Inspecting input structure
--------------------------

``` r
library(bio3d)
file.name <- get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

``` r
hiv <- read.pdb(file.name)
```

Prepare separate protein and ligand files
-----------------------------------------

``` r
prot <- trim.pdb(hiv, "protein")
lig <- trim.pdb(hiv, "ligand")

write.pdb(prot, file="1hsg_protein.pdb")
write.pdb(lig, "1hsg_ligand.pdb")
```

Inspect docking results
-----------------------

``` r
res <- read.pdb("all.pdbqt", multi=TRUE)

write.pdb(res, "results.pdb")

ori <- read.pdb("1hsg_ligand.pdbqt")

rmsd(ori, res)
```

    ##  [1]  0.590 11.163 10.531  4.364 11.040  3.682  5.741  3.864  5.442 10.920
    ## [11]  4.318  6.249 11.084  8.929

RMSD, no hydrogen
-----------------

``` r
ori_noh <- trim.pdb(ori, "noh")
res_noh <- trim.pdb(res, "noh")

rmsd(ori_noh, res_noh)
```

    ##  [1]  0.458 11.021 10.374  4.301 10.891  3.717  5.764  3.791  5.498 10.759
    ## [11]  4.224  6.308 10.889  8.776

Alternate RMSD, no hydrogen method
----------------------------------

``` r
ins <- atom.select(ori, "noh")
ins_res <- atom.select(res, "noh")

rmsd(ori$xyz[, ins$xyz], res$xyz[, ins_res$xyz])
```

    ##  [1]  0.458 11.021 10.374  4.301 10.891  3.717  5.764  3.791  5.498 10.759
    ## [11]  4.224  6.308 10.889  8.776
