Here we explain hot to extract reads spanning VDJ recombinations

```
python ~/code2/rop/source/iprofile/extract_VDJ.py unmapped_108318A_IGL_igblast.csv $PWD/VJ/ IGL 1e-05
```

If case you have many samples, which are saved in a directory 

```
for f in IGK/* ; do python ~/code2/rop/source/iprofile/extract_VDJ.py $f $PWD/VJ/ IGL 1e-05;done
```
