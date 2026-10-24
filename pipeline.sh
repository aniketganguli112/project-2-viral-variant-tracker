#!/bin/bash
read -p "give the accession no. of the samples,separated by space" -a acc_list
read -p "enterbiosample id,in same order of accession numer seperated by space:" -a biosample
echo "acc_no,dates" > dates.csv
read -p "give the path to ref_genome" genome
read -p "enter the name of the virus" name
bwa index "$genome"
mkdir QC_report
for acc in "${acc_list[@]}"
do
	if [[ "$acc" != "" ]]
	then
		echo "here is your fastq file:"
		fasterq-dump -p "$acc" -O samples
		ls samples
		fq1="samples/${acc}_1.fastq"
		fq2="samples/${acc}_2.fastq"
		fastqc "$fq1" "$fq2" --outdir QC_report
		echo "here is your fastq report"|ls QC_report
		read -p "give the name of output file1:" op1
		read -p "give the name of output file2:" op2
		fastp\
		-i "$fq1"\
		-I "$fq2"\
		-o "$op1"\
		-O "$op2"\
		-q 20\
		-u 30\
		--thread 4\
		-h "samples/${acc}.html"\
		-j "samples/${acc}.report.json"\
		-l 50\
		--n_base_limit 5\
		--cut_right\
		--cut_right_mean_quality 20
		bwa mem "$genome" "$op1" "$op2"|samtools sort -o "samples/${acc}_bam_output.bam"
		samtools index "samples/${acc}_bam_output.bam"
		samtools depth -a "samples/${acc}_bam_output.bam" > "samples/${acc}_depth"
		awk '$3<10 {print $1"\t"$2-1"\t"$2}' "samples/${acc}_depth" > "samples/${acc}_low_coverage.bed"
		bcftools mpileup -a AD,DP -d 1000 -q 20 -E --threads 4 -Ou -f "$genome" "samples/${acc}_bam_output.bam"|bcftools call --ploidy 1 -mv -Ob -o "samples/${acc}_vcf_output"
		bcftools view "samples/${acc}_vcf_output" > "samples/${acc}_vcf_op.vcf"
		echo "here is a sample of your vcf file"
		head -n 20 "samples/${acc}_vcf_op.vcf"
		bcftools filter -i "QUAL>20 && INFO/DP>=10" "samples/${acc}_vcf_op.vcf" -Ov -o "samples/${acc}_filtered_variants.vcf"
		bgzip -f "samples/${acc}_filtered_variants.vcf"
		tabix -p vcf "samples/${acc}_filtered_variants.vcf.gz"
		bcftools consensus -f "$genome" -m "samples/${acc}_low_coverage.bed"  "samples/${acc}_filtered_variants.vcf.gz" -o "samples/${acc}_consensus.fasta"
		sed  -i "1s/.*/>${acc}/" "samples/${acc}_consensus.fasta"
		head -n 5 "samples/${acc}_consensus.fasta"
		fi
done
cat samples/*_consensus.fasta > all_consensus.fasta
nextclade dataset get --name "$name" --output-dir dataset
nextclade run --input-dataset dataset --output-tsv nc.tsv all_consensus.fasta
cut -f1,2,3 nc.tsv|head
for bio in "${!acc_list[@]}"
do
        d=$(efetch -db biosample -id "${biosample[$bio]}" -format xml|grep -i "collection_date"|sed 's/.*>//')
        echo  "${acc_list[$bio]},$d" >> dates.csv
done
cat  dates.csv
python3 script2.py
