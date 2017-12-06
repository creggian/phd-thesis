# Roadmap data pipeline preprocessing

## Check file creation and logs for errors
```
qsub roadmapGzToFixedGz
```

## Check stdout log to check integrity test
```
qsub roadmapIntegrityGzToFixedGz
```

## create the .pp files (links) from either the original or fixed files
```
qsub roadmapCreatePrepLinks
```

## MACS2 whole genome peak calling
```
qsub -N FE_4me1_4me3_27ac -v TISSUE_LIST=Fetal_Brain,MARKER_GREP_STRING="Histone_H3K4me1\\\|Histone_H3K4me3\\\|Histone_H3K27ac" roadmapWgMacs
```

## prepare narrowPeak files to be moved into Hadoop cluster
```
qsub roadmapNarrowPeakMetadata
```
