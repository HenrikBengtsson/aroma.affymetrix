######################################################################
# Downloads all files needed for testing the Aroma Framework
######################################################################
library("aroma.core");

ar <- AromaRepository(verbose=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ANNOTATION DATA FILES
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

chipType <- "Cytogenetics_Array";   # [x]
#downloadCDF(ar, chipType); # [x]
##downloadACS(ar, chipType, tags="HB20080915");
##downloadUFL(ar, chipType, tags="na26,HB20080916");
downloadUGP(ar, chipType, tags="na28,HB20090519"); # [x]

chipType <- "Mapping10K_Xba142";   # [x]
#downloadCDF(ar, chipType);   # [x]
downloadACS(ar, chipType, tags="HB20080915");
##downloadUFL(ar, chipType, tags="na26,HB20080916");
##downloadUGP(ar, chipType, tags="na26,HB20080916");

chipType <- "Mapping50K_Hind240";   # [x]
downloadACS(ar, chipType, tags="HB20080710");
#downloadCDF(ar, chipType);
downloadUFL(ar, chipType, tags="na31,hg19,HB20110328");
downloadUGP(ar, chipType, tags="na31,hg19,HB20110328");

chipType <- "Mapping50K_Xba240";
downloadACS(ar, chipType, tags="HB20080710");
#downloadCDF(ar, chipType);
downloadUFL(ar, chipType, tags="na31,hg19,HB20110328");
downloadUGP(ar, chipType, tags="na31,hg19,HB20110328");

chipType <- "Mapping250K_Nsp";
downloadACS(ar, chipType, tags="HB20080710");
#downloadCDF(ar, chipType);
downloadUFL(ar, chipType, tags="na31,HB20101007");
downloadUGP(ar, chipType, tags="na31,HB20101007");

chipType <- "Mapping250K_Sty";
downloadACS(ar, chipType, tags="HB20080710");
#downloadCDF(ar, chipType);
downloadUFL(ar, chipType, tags="na31,HB20101007");
downloadUGP(ar, chipType, tags="na31,HB20101007");

chipType <- "GenomeWideSNP_5";
downloadACS(ar, chipType, tags="HB20080710");

chipType <- "GenomeWideSNP_5,Full,r2";  # [x]
#downloadCDF(ar, chipType);  # [x]
##downloadUFL(ar, chipType, tags="na26,HB20080822");
##downloadUGP(ar, chipType, tags="na26,HB20080822");

chipType <- "GenomeWideSNP_6"; # [x]
#downloadCDF(ar, chipType); # [x]
downloadACS(ar, chipType, tags="HB20080710"); # [x]
downloadUFL(ar, chipType, tags="na31,hg19,HB20110328"); # [x]
downloadUGP(ar, chipType, tags="na31,hg19,HB20110328"); # [x]

chipType <- "GenomeWideSNP_6,Full"; # [x]
#downloadCDF(ar, chipType); # [x]
downloadUFL(ar, chipType, tags="na31,hg19,HB20110328"); # [x]
downloadUGP(ar, chipType, tags="na31,hg19,HB20110328"); # [x]

chipType <- "HG-U133_Plus_2"; # [x]
#downloadCDF(ar, chipType);

chipType <- "Hs_PromPR_v02"; # [x]
#downloadCDF(ar, chipType);
downloadACS(ar, chipType);
downloadACM(ar, chipType);
downloadACP(ar, chipType, tags="unique");
downloadACC(ar, chipType, tags="unique");

chipType <- "HuEx-1_0-st-v2"; # [x]
downloadCDF(ar, chipType, tags="coreR3,A20071112,EP"); # [x]
downloadCDF(ar, chipType, tags="fullR3,A20071112,EP"); # [x]

chipType <- "Mm_PromPR_v02"; # [x]
#downloadCDF(ar, chipType); # [x]

chipType <- "MoEx-1_0-st-v1"; # [x]
downloadACS(ar, chipType, tags="HB20100926"); # [x]
downloadCDF(ar, chipType, tags="coreR1,A20080718,MR"); # [x]

chipType <- "MoGene-1_0-st-v1"; # [x]
downloadCDF(ar, chipType, tags="r3"); # [x]
downloadACS(ar, chipType, tags="r3,HB20110325"); # [x]

chipType <- "MOUSEDIVm520650"; # [x]
# downloadCDF(ar, chipType); # [x]
downloadACS(ar, chipType, tags="HB20100530"); # [x]
downloadUFL(ar, chipType, tags="na30,mm9,HB20100603"); # [x]
downloadUGP(ar, chipType, tags="na30,mm9,HB20100603"); # [x]

chipType <- "Test3"; # [x]
# downloadCDF(ar, chipType); # [x]


######################################################################
# HISTORY:
# 2011-09-29
# o Created.
######################################################################
