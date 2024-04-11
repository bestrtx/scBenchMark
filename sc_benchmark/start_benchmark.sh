#!/bin/bash
mkdir output
mkdir output/step1
mkdir output/step2

mkdir output/step1/sctransform
Rscript R_conf.R
Rscript "./scripts/step1/sctransform.R"
python3 evaluation.py './output/step1/nomor_output.txt'

mkdir output/step1/scran
Rscript R_conf.R
Rscript "./scripts/step1/scran.R"
python3 evaluation.py './output/step1/nomor_output.txt'

mkdir output/step1/Dino
Rscript R_conf.R
Rscript "./scripts/step1/Dino.R"
python3 evaluation.py './output/step1/nomor_output.txt'

mkdir output/step1/Log
Rscript R_conf.R
Rscript "./scripts/step1/Log.R"
python3 evaluation.py './output/step1/nomor_output.txt'

mkdir output/step1/LIGER
Rscript R_conf.R
Rscript "./scripts/step1/LIGER.R"
python3 evaluation.py './output/step1/nomor_output.txt'


#插补方法进行测试
mkdir output/step2/deepImpute
Rscript R_conf.R
Rscript "./scripts/step2/pre_to_DeepImpute.R"
deepImpute --cores 50 ./output/step2/deepImpute/DeepImpute_log.csv -o ./output/step2/deepImpute/imputed.csv
Rscript "./scripts/step2/deepImpute.R"
python3 evaluation.py './output/step2/imputation_output.txt'

mkdir output/step2/DrImpute
Rscript R_conf.R
Rscript "./scripts/step2/DrImpute.R"
python3 evaluation.py './output/step2/imputation_output.txt'

mkdir output/step2/Rmagic
Rscript R_conf.R
Rscript "./scripts/step2/Rmagic.R"
python3 evaluation.py './output/step2/imputation_output.txt'

mkdir output/step2/SAVER
Rscript R_conf.R
Rscript "./scripts/step2/SAVER.R"
python3 evaluation.py './output/step2/imputation_output.txt'

mkdir output/step2/ALRA
Rscript R_conf.R
Rscript "./scripts/step2/ALRA.R"
python3 evaluation.py './output/step2/imputation_output.txt'

#此插补方法放在最后测试
mkdir output/step2/scImpute
Rscript R_conf.R
Rscript "./scripts/step2/scImpute.R"
python3 evaluation.py './output/step2/imputation_output.txt'

#此标准化方法放在最后测试
mkdir output/step1/SCnorm
Rscript R_conf.R
Rscript "./scripts/step1/SCnorm.R"
python3 evaluation.py './output/step1/nomor_output.txt'