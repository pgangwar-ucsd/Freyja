if [ $# -eq 1 ]; then
    file_path=$1
    source ~/.bashrc
    conda activate freyja-env
    #make dev

    #sed -i 's/NC_045512v2/NC_045512.2/g' my_vcf_alignment.sam 
    ##sed -i 's/MN908947.3/NC_045512.2/g' my_vcf_alignment.sam 
    #samtools view -Sb my_vcf_alignment.sam > my_vcf_alignment.bam
    #samtools sort my_vcf_alignment.bam -o resorted.bam 

    freyja variants ${file_path}/resorted.bam --variants ${file_path}/cwap_variants.tsv --depths ${file_path}/cwap_depth.tsv --minq 20
    #freyja demix ${file_path}/cwap_variants.tsv ${file_path}/cwap_depth.tsv --barcodes ${file_path}/barcodes.csv --output ${file_path}/freyja_output_latest.txt --eps 0.005
    #python correct_format.py
    conda deactivate
fi