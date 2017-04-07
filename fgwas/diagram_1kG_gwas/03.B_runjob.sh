
#$ -N build
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e /well/got2d/jason/projects/t2d-integration/fgwas/diagram_1kG_gwas/logs/build.error
#$ -o /well/got2d/jason/projects/t2d-integration/fgwas/diagram_1kG_gwas/logs/build.out
#$ -V

module load R/3.3.1
Rscript	--vanilla /well/got2d/jason/projects/t2d-integration/fgwas/diagram_1kG_gwas/03.A_build_fgwas_annot.R
