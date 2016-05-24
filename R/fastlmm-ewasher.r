#################################################################
# call FLE by pointing to data frames containing relevant inputs
#################################################################
FLE = function(rSrcDir, data, phen, covar)
{        
    #write each input to file, then call RunFLEFromFiles
    write.table(data, file = "datafileTmp.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(phen, file = "phenfileTmp.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(covar, file = "covarfileTmp.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    RunFLEFromFiles(rSrcDir, "datafileTmp.txt", "phenfileTmp.txt", "covarfileTmp.txt")
}

#################################################################
# call FLE by pointing to the input data files
#################################################################
RunFLEFromFiles = function(rSrcDir, dataFileName, phenotypeFileName, covarFileName = NULL, mapFileName = NULL)
{
    #### these need to be switched to params
    
    #constitutively unmethylated/methylated filter values
    minFilterVal = 0.2; maxFilterVal = 0.8
    
    #values for fastlmm autoselect to try for number of loci in GSM
    #may be advantageous to try a finer, or more enlarged grid, but this default setting seems reasonable.
    asVals = append(append(c(1), seq(100, 900, 100)), seq(1000, 10000, 1000))
    
    # maximum number of pcs to allow.
    maxPC = 10
    
    #threshold on lambda for when to stop adding PCs
    lamdaThresh = 1.0
    
    #### end of params
    
    
    
    source(file.path(rSrcDir, 'utils.r'))
    
    # normalize filenames before we switch the working dir
    rSrcDir = normalizePath(rSrcDir)
    data_file = normalizePath(dataFileName)
    phenotype_file = normalizePath(phenotypeFileName)
    covar_file = covarFileName
    if (length(covar_file) != 0)
    {
        covar_file = normalizePath(covar_file)
    }
    map_file = mapFileName
    if (length(map_file) != 0)
    {
        map_file = normalizePath(map_file)	
    }
    
    # create dirs and set the temp dir as the working dir
    work_dir = getwd()
    results_dir = file.path(work_dir, 'results')
    temp_dir = file.path(work_dir, 'tmp')
    if (!file.exists(results_dir))
    {
        dir.create(results_dir)
    }
    if (!file.exists(temp_dir))
    {
        dir.create(temp_dir)
    }
    setwd(temp_dir)
    
    # get the filename of the executable for this OS
    if (substr(version$platform, 1, 11) == 'x86_64-w64-')
    {
        fastlmmFileName = file.path(rSrcDir, 'fastlmmC.exe')
    } else if (version$os == "linux")
    {
        fastlmmFileName = file.path(rSrcDir, 'fastlmmC')
    } else
    {
        stop('unsupported architecture or operating system')
    }
    
    cat('preparing inputs for FaST-LMM-EWASher...\n')
    CreateFastlmmInputFiles(data_file, phenotype_file, map_file, minFilterVal, maxFilterVal, asVals, rSrcDir)
    if (length(covar_file) == 0)
    {
        covar_file = CreateCovarFile()
    }
    
    cat('running the main loop of FaST-LMM-EWASher...\n')
    TK = RunMainLoop(fastlmmFileName, covar_file, maxPC, lamdaThresh)
    T = TK$T
    
    cat('cleaning up...')
    file.copy('out_linreg.txt', file.path(results_dir, 'out_linreg.txt'), overwrite = TRUE)
    file.copy('qq_linreg.png', file.path(results_dir, 'qq_linreg.png'), overwrite = TRUE)
    file.copy(paste('out_lmm_', T, '.txt', sep = ''),
        file.path(results_dir, 'out_ewasher.txt'), overwrite = TRUE)
    file.copy('qq_plot.png', file.path(results_dir, 'qq_ewasher.png'), overwrite = TRUE)
    file.copy(paste('K_lmm_', T, '.txt', sep = ''),
        file.path(results_dir, 'similarity.txt'), overwrite = TRUE)
    setwd(results_dir)
    
    WriteSummary(T, TK$K)
}

#if user specified mapping file, read it here, otherwise read in the default annotation info
ReadMapFile = function(mapFileName, rSrcDir)
{
	if (length(mapFileName) == 0)
	{
		cat("  loading default map file for Illumina methylation chip\n")
		map = read.delim(file.path(rSrcDir, 'annot450k.tab'), FALSE)
		chrDict = map$V2
		names(chrDict) = map$V1
		posDict = map$V3
		names(posDict) = map$V1
	} else
	{
		cat("  reading in map file:", mapFileName, '\n')
		map = read.delim(mapFileName, FALSE)
		chrDict = map$V1
		names(chrDict) = map$V2
		posDict = map$V3
		names(posDict) = map$V2
	}
  
    return(list(ChrDict = chrDict, PosDict = posDict))
}

# Input data: a matrix of beta values.
# Each row correspond to one probe and each column correspond to one sample.
# The first column lists probe IDs and the first row lists sample IDs. 
CreateFastlmmInputFiles = function(dataFileName, phenotypeFileName, mapFileName, minFilterVal, maxFilterVal, asVals, rSrcDir)
{
	dataWithIDs = read.delim(dataFileName, TRUE)	
	probeIDs = dataWithIDs$ID
	data = dataWithIDs[, 2:ncol(dataWithIDs)]
	sampleIDs = names(data)
	
	phenotypes = read.delim(phenotypeFileName, FALSE)
	names(phenotypes) = c('SampleID', 'Phenotype')
	phenoDict = phenotypes$Phenotype
    names(phenoDict) = gsub('-', '.', phenotypes$SampleID)
	sampleOriginalIDDict = as.character(phenotypes$SampleID)
    names(sampleOriginalIDDict) = names(phenoDict)
    
	# create phenotype file and fam file
    outPheno = t(sapply(sampleIDs,
        function(id) c(1, sampleOriginalIDDict[id], phenoDict[id])))
	write.table(outPheno, file = 'phenotype.txt', quote = FALSE, sep = '\t',
        row.names = FALSE, col.names = FALSE)
    outFam = t(sapply(sampleIDs,
        function(id) c(1, sampleOriginalIDDict[id], 0, 0, 0, phenoDict[id])))
    write.table(outFam, file = 'testmarkers.fam', quote = FALSE, sep = '\t',
        row.names = FALSE, col.names = FALSE)

    # remove rows with probe ID starting with rs or ch
    probeIDPrefixes = substr(probeIDs, 1, 2)
	rowsToKeep = probeIDPrefixes != 'rs' & probeIDPrefixes != 'ch'
	probeIDs = probeIDs[rowsToKeep]
	data = data[rowsToKeep, ]
    
    # remove rows with probe ID not found in map file
	map = ReadMapFile(mapFileName, rSrcDir)
	rowsToKeep = probeIDs %in% names(map$ChrDict)
	numFound = sum(rowsToKeep)
	numMissing = sum(!rowsToKeep)
	probeIDs = probeIDs[rowsToKeep]
	data = data[rowsToKeep, ]
	cat("  num missing map markers =", numMissing, '\n')
	cat("  num found =", numFound, '\n')
	if (numFound == 0)
    {
        stop("no mapping information found for any markers")
	}
    
    # remove rows where more than half of the data is missing
    rowsToKeep = t(apply(data, 1, function(row) sum(!is.na(row)) >= length(row) / 2))
    probeIDs = probeIDs[rowsToKeep]
	data = data[rowsToKeep, ]
	
	# remove rows where the average is out of range
	averages = t(apply(data, 1, function(row) mean(row, na.rm = TRUE)))
    rowsToKeep = averages >= minFilterVal & averages <= maxFilterVal
	probeIDs = probeIDs[rowsToKeep]
	data = data[rowsToKeep, ]
	
	# get dictionary of data file contents
    dataDict = as.data.frame(apply(data, 1, function(row)
    {
        min = min(row, na.rm = TRUE)
        range = max(row, na.rm = TRUE) - min
        average = mean(row, na.rm = TRUE)
        tempRow = replace(row, is.na(row), average)
        return(append(c('A', 'C'), round(2 * (tempRow - min) / range, 3)))
    }))
    names(dataDict) = probeIDs
    rownames(dataDict) = append(c('A', 'C'), sampleIDs)
    
	# get sorted map file contents
    rowsToKeep = names(map$ChrDict) %in% probeIDs
	chr = map$ChrDict[rowsToKeep]
	pos = map$PosDict[rowsToKeep]
	tempMap = data.frame(Chr = chr, ProbeID = names(chr), PosA = pos, PosB = pos)
	sortedMap = suppressWarnings(tempMap[order(as.numeric(as.vector(tempMap$Chr)),
        tempMap$Chr, tempMap$PosA, tempMap$ProbeID), ])
	
    # create map file and data file
    write.table(sortedMap, file = 'testmarkers.map', quote = FALSE, sep = '\t',
        row.names = FALSE, col.names = FALSE)
    outData = t(apply(sortedMap, 1, function(row) as.vector(dataDict[, row[2]])))
    write.table(outData, file = 'testmarkers.dat', quote = FALSE, sep = '\t',
        row.names = TRUE, col.names = FALSE)
    
    # create AS values file
    numLoci = nrow(sortedMap)
    asVals = asVals[asVals <= numLoci]
    cat("  using autoselect values:", asVals, '\n')
    write(asVals, file = 'ASvalues.txt', ncolumns = length(asVals), sep = ' ')
    
    # create principal components file
    cat('computing principal components...\n')
    pcData = outData[, 3:ncol(outData)]
    mode(pcData) = 'numeric'
	dataM = apply(pcData, 1, function(row)
	{
	    average = mean(row, na.rm = TRUE)
	    return(row - average)
	})
	X = svd(dataM) # spectral decomposition
	scores = t(apply(X$u, 1, function(row) row * sqrt(X$d)))
    write.table(scores[, 1:100], file = 'pc.txt', quote = FALSE, sep = '\t',
        row.names = FALSE, col.names = FALSE)
}

CreateCovarFile = function()
{
    covar = read.delim('testmarkers.fam', FALSE)
    covarOut ='covariates.txt'
    write.table(covar[, 1:2], file = covarOut, quote = FALSE, sep = '\t',
        row.names = FALSE, col.names = FALSE)
    return(covarOut)
}

RunMainLoop = function(fastlmmFileName, covarFileName, maxPC, lamdaThresh)
{
	# run linear regression
	cmd = paste(fastlmmFileName, ' -dfile1 testmarkers -pheno phenotype.txt -SetOutputPrecision 3 -linreg -out out_linreg.txt -REML -covar ', covarFileName, ' 2>log.txt', sep = '')
	shell(cmd)
	TryMakeQQPlot('out_linreg.txt','qq_linreg.png')
    
    for (T in 0:maxPC)
    {
		cat("  trying", T, "principal components...\n")
		pcCovarFileName = WritePCCovar(covarFileName, 'pc.txt', T)
		
		cmd = paste(fastlmmFileName, ' -autoSelect ASout -autoSelectSearchValues ASvalues.txt -dfile1sim testmarkers -pheno phenotype.txt -REML -covar ', pcCovarFileName, ' 2>log.txt', sep = '')
		shell(cmd)
		
		numMarkers = nrow(read.delim('ASout.snps.txt', FALSE))
		cat("  num markers selected =", numMarkers, '\n')
		
		outFileName = paste('out_lmm_', T, '.txt', sep = '')
		simOutFileName = paste('K_lmm_', T, '.txt', sep = '')
		cmd = paste(fastlmmFileName, ' -dfile1 testmarkers -dfile1Sim testmarkers -pheno phenotype.txt -excludebygeneticdistance 50000 -REML -SetOutputPrecision 3 -verboseOut -extractSim ASout.snps.txt -out ', outFileName, ' -simOut ', simOutFileName, ' -covar ', pcCovarFileName, ' 2>log.txt', sep = '')
		shell(cmd)
		
		L = EstimateLambda(outFileName)
        cat("  lambda =", L, '\n')
		if (T == maxPC || L < lamdaThresh)
		{
		    TryMakeQQPlot(outFileName, 'qq_plot.png')
			break
		}
    }
    
	return(list(T = T, K = numMarkers))
}

WriteSummary = function(T, K)
{
    numLoci = nrow(read.delim('out_ewasher.txt', TRUE))
	write(paste('number of loci analyzed: ', numLoci,
        '\nnumber of loci selected for the LMM: ', K,
        '\nnumber of principal components used by FaST-LMM-EWASher: ', T, '\n',
        sep = ''), file = 'summary.txt')
}

# to run the demo using source(), first uncomment the following line and set the
# working directory to ../demo
#RunFLEFromFiles('../R', 'input_data.txt', 'input_phenotype.txt', 'input_covariates.txt')
