conda activate freyja-env
#make dev
#mkdir -p data
#freyja update --buildlocal --outdir data
##rm my_* cwap* resorted.bam

#sed -i 's/NC_045512v2/NC_045512.2/g' SRR30906372_alignment.sam 
##sed -i 's/MN908947.3/NC_045512.2/g' my_vcf_alignment.sam 
#samtools view -Sb SRR30906372_alignment.sam > SRR30906372_alignment.bam
#samtools sort SRR30906372_alignment.bam -o resorted.bam 

freyja variants resorted.bam --variants cwap_variants.tsv --depths cwap_depth.tsv --minq 20
#freyja demix cwap_variants.tsv cwap_depth.tsv --barcodes data/usher_barcodes.csv --output my_output_latest.txt --eps 0.005
#python correct_format.py
conda deactivate
