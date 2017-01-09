#!/bin/bash

RCODE=/home/gosia/R/package_devel/drimseq_package_devel_tests
RWD=/home/gosia/multinomial_project/package_devel/Test_DRIMSeq_1_3_2_beta_binomial
ROUT=$RWD/Rout

mkdir -p $ROUT

R CMD BATCH --no-save --no-restore "--args rwd='$RWD'" $RCODE/Test_DRIMSeq_0.3.3_auto_moderation.R $ROUT/Test_DRIMSeq_0.3.3_auto_moderation_${count_method}.Rout








###
