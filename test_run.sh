#!/bin/bash

RCODE=/home/gosia/R/package_devel/tests
RWD=/home/gosia/multinomial_project/package_devel/Test_DRIMSeq_0.3.3_auto_moderation
ROUT=$RWD/Rout

mkdir $ROUT


count_method='kallistofiltered5'

R32dev2 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}'" $RCODE/Test_DRIMSeq_0.3.3_auto_moderation.R $ROUT/Test_DRIMSeq_0.3.3_auto_moderation_${count_method}.Rout


count_method='htseqprefiltered5'

R32dev2 CMD BATCH --no-save --no-restore "--args rwd='$RWD' count_method='${count_method}'" $RCODE/Test_DRIMSeq_0.3.3_auto_moderation.R $ROUT/Test_DRIMSeq_0.3.3_auto_moderation_${count_method}.Rout






