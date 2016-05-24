#################################################################
# wrapper of FaST-LMM-EWASHer (FLE)
# 
# usage: FLE-script.r data_file phenotype_file [-covar covar_file] [-map map_file]
#################################################################

# get the executing directory
fileArg = commandArgs()[substring(commandArgs(), 1, 7) == '--file=']
thisFile = normalizePath(substring(fileArg, 8))
rSrcDir = substring(thisFile, 1, nchar(thisFile) - 13)

#import the main functions
source(file.path(rSrcDir, 'fastlmm-ewasher.r'))

# parse the cmd args
cmdArgs = commandArgs(trailingOnly = TRUE)
if (!any(length(cmdArgs) == c(1, 2, 4, 6)))
{
    stop('invalid number of args')
}

dataFileName = cmdArgs[1]
phenotypeFileName = cmdArgs[2]
covarFileName = NULL
mapFileName = NULL

if (length(cmdArgs) > 2)
{
    if (cmdArgs[3] == '-covar')
    {
        covarFileName = cmdArgs[4]
    } else if (cmdArgs[3] == '-map')
    {
        mapFileName = cmdArgs[4]
    } else
    {
        stop(paste('invalid third arg "', cmdArgs[3], '" should be "-covar" or "-map"', sep = ''))
    }
    
    if (length(cmdArgs) == 6)
    {
        if (cmdArgs[5] == '-covar')
        {
            covarFileName = cmdArgs[6]
        } else if (cmdArgs[5] == '-map')
        {
            mapFileName = cmdArgs[6]
        } else
        {
            stop(paste('invalid fifth arg "', cmdArgs[5], '" should be "-covar" or "-map"', sep = ''))
        }
    }
}

# call FLE with the cmd args given
RunFLEFromFiles(rSrcDir, dataFileName, phenotypeFileName, covarFileName, mapFileName)
