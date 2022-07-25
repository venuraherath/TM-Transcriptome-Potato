#########################################################RCODE###########################

#loading the libary

library(IsoformSwitchAnalyzeR)


### Import Kallisto example data in R package
KallistoQuant <- importIsoformExpression(
      parentDir = "./",
      addIsofomIdAsColumn = TRUE
)



### Make design matrix
myDesign <- read.table("mydesign.tsv.txt", header = TRUE)


### conditions
cond <- data.frame(condition_1 = c('Mock2hr','Mock5hr'),
                   condition_2 = c('TM2hr','TM5hr'))


### Create switchAnalyzeRlist
aSwitchList <- importRdata(
  isoformCountMatrix   = KallistoQuant$counts,
  isoformRepExpression = KallistoQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "cr.working_models.pm.locus_assign_clean_kallisto.gtf.gz",
  isoformNtFasta       = "cr.working_models.pm.cdna_isoremoved.fasta",
  comparisonsToMake = cond,
  showProgress = FALSE
)

#checking the list

head(aSwitchList$isoformFeatures,2)

head(aSwitchList$exons,2)

head(aSwitchList$ntSequence,2)

#prefiltering 

#relaxed filtering
aSwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 2,
  removeSingleIsoformGenes = TRUE
)


#Testing for Isoform Switches via DEXSeq

# pre-filter

aSwitchListFiltered <- preFilter(aSwitchListFiltered) # preFilter to remove lowly expressed features


# Perform test
aSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = aSwitchListFiltered,
  reduceToSwitchingGenes=TRUE
)

extractSwitchSummary(aSwitchListAnalyzed)

#check whehther ORF are added

addORFfromGTF(aSwitchListAnalyzed) 


#discovery of both known and novel transcripts #optional step

analyzeNovelIsoformORF(aSwitchListAnalyzed, 
                       analysisAllIsoformsWithoutORF=TRUE, # also analyse all those annoatated as without CDS in ref annottaion
                       genomeObject = NULL,
                       
                       ### Advanced argument
                       minORFlength = 100,
                       orfMethod = 'longest.AnnotatedWhenPossible',
                       PTCDistance = 50,
                       startCodons = "ATG",
                       stopCodons = c("TAA", "TAG", "TGA"),
                       showProgress = TRUE,
                       quiet = FALSE
)


#Extracting Nucleotide and Amino Acid Sequences

## This example relies on the example data from the 'Analyzing Open Reading Frames' section above 

aSwitchListAnalyzed <- extractSequence(
  aSwitchListAnalyzed, 
  alsoSplitFastaFile = TRUE, #to facilitate the pfam submission
  pathToOutput = "./",
  writeToFile=TRUE 
)

#importing from external tools
### Load test data (matching the external sequence analysis results)
#data("exampleSwitchListIntermediary")
#exampleSwitchListIntermediary


### Add CPAT analysis #one is enough CPC2 was carried out
#aSwitchListAnalyzed <- analyzeCPAT(
  #switchAnalyzeRlist   = aSwitchListAnalyzed,
  #pathToCPATresultFile = system.file("extdata/cpat_results.txt", package = "IsoformSwitchAnalyzeR"),#needs to be changed
  #codingCutoff         = 0.725, # the coding potential cutoff we suggested for human
  #removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
#)

### Add CPC2 analysis
aSwitchListAnalyzed <- analyzeCPC2(
  switchAnalyzeRlist   = aSwitchListAnalyzed,
  pathToCPC2resultFile = "cpc2_result.txt",#needs to be changed
  removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
)

### Add PFAM analysis
aSwitchListAnalyzed <- analyzePFAM(
  switchAnalyzeRlist   = aSwitchListAnalyzed,
  pathToPFAMresultFile = "pfam_results.txt",#needs to be changed
  showProgress=FALSE
)

### Add SignalP analysis
aSwitchListAnalyzed <- analyzeSignalP(
  switchAnalyzeRlist       = aSwitchListAnalyzed,
  pathToSignalPresultFile  = "SingalP5_output_protein_type.txt" #needs to be changed
)

### Add NetSurfP2 analysis
aSwitchListAnalyzed <- analyzeIUPred2A(
  switchAnalyzeRlist        = aSwitchListAnalyzed,
  pathToIUPred2AresultFile = "IUPRED2_isoformSwitchAnalyzeR_isoform_AA_complete.result",#needs to be changed
  showProgress = FALSE
)

aSwitchListAnalyzed


#Predicting Alternative Splicing

aSwitchListAnalyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = aSwitchListAnalyzed,
  quiet=TRUE
)

### overview of number of intron retentions (IR)
table( aSwitchListAnalyzed$AlternativeSplicingAnalysis$IR )

#Predicting Switch Consequences

# the consequences highlighted in the text above
consequencesOfInterest <- c('intron_retention','coding_potential','NMD_status','domains_identified','ORF_seq_similarity')

aSwitchListAnalyzed <- analyzeSwitchConsequences(
  aSwitchListAnalyzed,
  consequencesToAnalyze = consequencesOfInterest, 
  showProgress=FALSE
)

extractSwitchSummary(aSwitchListAnalyzed, dIFcutoff = 0.1, filterForConsequences = FALSE)

extractSwitchSummary(aSwitchListAnalyzed, dIFcutoff = 0.1, filterForConsequences = TRUE)

#Post Analysis of Isoform Switches with Consequences

#Analysis of Individual Isoform Switching


### Let's reduce the switchAnalyzeRlist to only one condition
aSwitchListAnalyzedSubset2hr <- subsetSwitchAnalyzeRlist(
  aSwitchListAnalyzed, 
  aSwitchListAnalyzed$isoformFeatures$condition_1 == 'Mock2hr' #?????????
)
aSwitchListAnalyzedSubset2hr


### Extract top switching genes (by q-value)
extractTopSwitches(
  aSwitchListAnalyzedSubset2hr, 
  filterForConsequences = TRUE, 
  n = NA, ####change-NA will give all significant
  sortByQvals = TRUE
)


### Extract data.frame with all switching isoforms to study individual genes
switchingIso <- extractTopSwitches( 
  aSwitchListAnalyzedSubset2hr, 
  filterForConsequences = TRUE, 
  n = NA,                  # n=NA: all features are returned
  extractGenes = FALSE,    # when FALSE isoforms are returned
  sortByQvals = TRUE
)

#looking at individual gene
subset(switchingIso, gene_id == 'Soltu.Cru.05_1G002100')

#plotting the gene of interest

switchPlot(aSwitchListAnalyzedSubset2hr, gene = 'Soltu.Cru.05_1G002100')
switchPlot(aSwitchListAnalyzedSubset2hr, gene = 'Soltu.Cru.S125310')
switchPlot(aSwitchListAnalyzedSubset2hr, gene = 'Soltu.Cru.06_1G003100')

#screening we have implemented switchPlotTopSwitches() which will extract the top n genes (controlled by the n argument) with significant switches (as defined by alpha and dIFcutoff) and output a pdf or png version

switchPlotTopSwitches(
  switchAnalyzeRlist = aSwitchListAnalyzed, 
  n = Inf,                                             # Set to Inf for all
  filterForConsequences = TRUE,
  fileType = "pdf",                                   # alternative is "png"
  pathToOutput = "./"
)


#Genome-Wide Analysis of Isoform Switching  

#Global summary statistics #above codes are focusing on one condition based subset

aSwitchListAnalyzed
extractConsequenceSummary(aSwitchListAnalyzed)

extractSplicingSummary(aSwitchListAnalyzed)

#Analysis of splicing/consequence enrichment

extractConsequenceEnrichment(aSwitchListAnalyzed)
extractSplicingEnrichment(aSwitchListAnalyzed)

#Analysis of splicing/consequence enrichment

extractConsequenceEnrichment(aSwitchListAnalyzed)
extractSplicingEnrichment(aSwitchListAnalyzed)

#Comparison of enrichment
extractConsequenceEnrichmentComparison(aSwitchListAnalyzed)
extractSplicingEnrichmentComparison(aSwitchListAnalyzed)

#Analysis of genome-wide changes in isoform usage

extractConsequenceGenomeWide(aSwitchListAnalyzed) 
extractSplicingGenomeWide(aSwitchListAnalyzed)

#individual gene visualization

switchPlotTranscript(aSwitchListAnalyzedSubset, gene = '???')

switchPlotGeneExp (aSwitchListAnalyzedSubset, gene = '???')

switchPlotIsoExp(aSwitchListAnalyzedSubset, gene = '???')

switchPlotIsoUsage(aSwitchListAnalyzedSubset, gene = '???')

#Analyzing Alternative Splicing  - First 4 steps are given above so no need to redo
aSwitchListAnalyzed

extractSplicingSummary(
  aSwitchListAnalyzed,
  asFractionTotal = FALSE,
  plotGenes=FALSE
)

#This function summarizes the uneven usage within each comparison by for each alternative splicing type calculate the fraction of events being gains (as opposed to loss) and perform a statistical analysis of this fraction
splicingEnrichment <- extractSplicingEnrichment(
  aSwitchListAnalyzed,
  splicingToAnalyze='all',
  returnResult=TRUE,
  returnSummary=TRUE
)

#some alternative splicing types are preferentially used in the two cancer types, and from the statistical summary also created by extractSplicingEnrichment() (due to the returnSummary=TRUE argument):
subset(splicingEnrichment, splicingEnrichment$AStype == 'ATSS')

#Enrichment difference although the overall trend in usage of alternative splicing is the same the two cancer types there seem to be some differences
extractSplicingEnrichmentComparison(
  aSwitchListAnalyzed,
  splicingToAnalyze=c('A3','MES','ATSS'), # Those significant in COAD in the previous analysis
  returnResult=FALSE # Preventing the summary statistics to be returned as a data.frame
)

# the analysis of alternative splicing can also be done from the perspective of the individual splice types

extractSplicingGenomeWide(
  aSwitchListAnalyzedSubset,
  featureToExtract = 'all',                 # all isoforms stored in the switchAnalyzeRlist
  splicingToAnalyze = c('???','???','???'), # Splice types significantly enriched in COAD
  plot=TRUE,
  returnResult=FALSE  # Preventing the summary statistics to be returned as a data.frame
)

#Analyzing the Biological Mechanisms Behind Isoform Switching

### Reduce datasize for fast runtime
selectedGenes <- unique(aSwitchListAnalyzedSubset$isoformFeatures$gene_id)[50:100]
aSwitchListAnalyzedSubsetSubset <- subsetSwitchAnalyzeRlist(
  aSwitchListAnalyzedSubset,
  aSwitchListAnalyzedSubset$isoformFeatures$gene_id %in% selectedGenes
)

### analyze the biological mechanisms

bioMechanismeAnalysis <- analyzeSwitchConsequences(
  aSwitchListAnalyzedSubsetSubset, 
  consequencesToAnalyze = c('tss','tts','intron_structure'),
  showProgress = FALSE
)$switchConsequence # only the consequences are interesting here

### subset to those with differences
bioMechanismeAnalysis <- bioMechanismeAnalysis[which(bioMechanismeAnalysis$isoformsDifferent),]

### extract the consequences of interest already stored in the switchAnalyzeRlist
myConsequences <- aSwitchListAnalyzedSubsetSubset$switchConsequence
myConsequences <- myConsequences[which(myConsequences$isoformsDifferent),]
myConsequences$isoPair <- paste(myConsequences$isoformUpregulated, myConsequences$isoformDownregulated) # id for specific iso comparison

### Obtain the mechanisms of the isoform switches with consequences
bioMechanismeAnalysis$isoPair <- paste(bioMechanismeAnalysis$isoformUpregulated, bioMechanismeAnalysis$isoformDownregulated)
bioMechanismeAnalysis <- bioMechanismeAnalysis[which(bioMechanismeAnalysis$isoPair %in% myConsequences$isoPair),]  # id for specific iso comparison

### Create list with the isoPair ids for each consequence
AS   <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'intron_structure')]
aTSS <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'tss'             )]
aTTS <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'tts'             )]

mechList <- list(
  AS=AS,
  aTSS=aTSS,
  aTTS=aTTS
)

### Create Venn diagram
library(VennDiagram)
#> Loading required package: grid
#> Loading required package: futile.logger
myVenn <- venn.diagram(
  x = mechList,
  col='transparent',
  alpha=0.4,
  fill=RColorBrewer::brewer.pal(n=3,name='Dark2'),
  filename=NULL
)

### Plot the venn diagram
grid.newpage() ; grid.draw(myVenn)


### Volcano like plot:
ggplot(data=aSwitchListAnalyzed$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  facet_wrap( ~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw()

### Switch vs Gene changes:
ggplot(data=aSwitchListAnalyzed$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) + 
  facet_wrap(~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='Gene log2 fold change', y='dIF') +
  theme_bw()


#############################End##########################
