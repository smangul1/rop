[To be updated for ImReP output]

Based on IgBLAST result one can quantify immune response based on the combinatorial diversity of BCR and TCR loci. We consider reads spanning V(D)J recombinations and estimate the relative proportion of each recombination. Per sample immune diversity is estimated  measured	using	the	alpha diversity (measured by Shannon entropy)	and incorporates total	number	of	V(D)J	combinations	and	their	relative	
proportions. 

Please use `/rop/source/diversity/alphaGeneral.py` to calculate the total number of reads spanning the V(D)J recombinations, number of V(D)J recombinations and alpha diversity measured by Shannon entropy for all the samples in 'raw/VJ_recombination/IGK_VJ/' directory. The results are saved into `alphaIGK_VJ.csv` located in 
`alpha` directory. Make sure to provide the extension of the files. 

```
python <ropDir>/rop/source/diversity/alphaGeneral.py raw/VJ_recombination/IGK_VJ/ alpha IGK_VJ 0 1 IGK_igblast.IGK_recomb
```
