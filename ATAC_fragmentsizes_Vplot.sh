#!/bin/bash

####SBATCH --job-name Atac  
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

# --mem-per-cpu 6000 

# sbatch --array=0-7 script.sh    # 8 control samples

# Modules required

module load parallel
module load samtools/1.9
module load macs/2.1.0
module load R/3.6.0
module load macs/2.1.0
#module load UCSC-tools/v373 # bedGraphToBigWig
#module load bedops

if true; then # set of modules required to load deepTools:
module load bedtools 
module load BigWig_tools 
module load bx-python
module load deepTools/2.4.2
fi

module load UCSC-tools/v373  # bamtobigwig, liftover..

adir="bam"
#files=(${adir}/Sample_HeLa_interphase ${adir}/Sample_HeLa_mitosis)
files=(`ls ${adir}/*bam`)
file=${files[SLURM_ARRAY_TASK_ID]}
hg38="/path/to/your/hg38/folder/bowtie/2.1.0/hg38.fa"
mkdir -p pdf/


# stats on alignments

if false; then # sbatch --array=0-9 script.sh

cd bam
files=($(ls -d Sample*)) # directories to get basenames
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=`echo "${file//Sample_/}"`
cd ${file}
echo -e "sample\traw_reads_fastq\tnb_align_5max\tbw2_uniq\tbw2_hc_uniq\tfinal_nb_unique_aligned_pairs\tpercentage_unique_alignment" > stats_pipeline.txt
ndir=sorted_files_for_stats
mkdir -p $ndir

samtools sort -@ \$SLURM_CPUS_PER_TASK bw2_${file}.bam ${ndir}/bw2_${file}_sort 2> ${ndir}/bw2_${file}_sort.err
samtools sort -@ \$SLURM_CPUS_PER_TASK bw2_uniq_${file}.bam ${ndir}/bw2_uniq_${file}_sort 2> ${ndir}/bw2_uniq_${file}_sort.err
samtools sort -@ \$SLURM_CPUS_PER_TASK bw2_hc_uniq_${file}.bam ${ndir}/bw2_hc_uniq_${file}_sort 2> ${ndir}/bw2_hc_uniq_${file}_sort.err # old -- to trash

samtools index ${ndir}/bw2_${file}_sort.bam 2>> ${ndir}/bw2_${file}_sort.err
samtools index ${ndir}/bw2_uniq_${file}_sort.bam 2>> ${ndir}/bw2_uniq_${file}_sort.err
samtools index ${ndir}/bw2_hc_uniq_${file}_sort.bam 2>> ${ndir}/bw2_hc_uniq_${file}_sort.err # old -- to trash

tmp=(`cat ../../stats_reads.txt | grep ${NAME}`)
nfastq=${tmp[1]}
bw2=`samtools view ${ndir}/bw2_${file}_sort.bam | wc -l `
bw2_uniq=`samtools view ${ndir}/bw2_uniq_${file}_sort.bam | wc -l `
bw2_hc_uniq=`samtools view bw2_hc_uniq_${file}.bam | wc -l `
final=`samtools view ../${NAME}_full.bam | wc -l `
percent=`echo ${final}/${nfastq} | bc -l`
echo -e "${NAME}\t${nfastq}\t${bw2}\t${bw2_uniq}\t${bw2_hc_uniq}\t${final}\t${percent}" >> stats_pipeline.txt

cd ../..

# TO DO final step: merging all info into one single file 

fi 


# filter for fragments sizes on bam files

if false; then
cd bam
files=($(ls *_full.bam)) 
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=$(basename $file _full.bam)
samtools view -H ${file} > header_${SLURM_ARRAY_TASK_ID}
samtools view ${file} | awk '($9*$9)<19881' | cat header_${SLURM_ARRAY_TASK_ID} - | samtools view -bS - > ${NAME}_0-140.bam # fragments lengths all positive in new pipeline
samtools sort ${NAME}_0-140.bam > tmp_${SLURM_ARRAY_TASK_ID}.bam && mv tmp_${SLURM_ARRAY_TASK_ID}.bam ${NAME}_0-140.bam 
samtools index ${NAME}_0-140.bam 

rm header_${SLURM_ARRAY_TASK_ID}
cd -
fi

# generation of bigwigs from bam files: slurm_bigwigs_deeptools.sh


# convert bam into bed (for downstream analysis)    # TO DO AGAIN ON NEW _full.bam FILE

if false; then
# sbatch --array=0-9 script.sh 
mkdir -p bed

cd bam

files=($(ls *_full.bam))
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=$(basename "${file}" _full.bam)

mkdir -p ../bed/${NAME}
samtools view ${file} | \
awk '$4>0 {if($9>0){print $3 "\t" $4-1 "\t" $4+$9-1 "\t" $9}else{print $3 "\t" $8-1 "\t" $8-$9-1 "\t" $9*-1} }' | \
sortBed -i - > ../bed/${NAME}/${NAME}.bed # 0-based/half-open bed format 
#sortBed -i tmp_${SLURM_ARRAY_TASK_ID} > ../bed/${NAME}/${NAME}.bed # 0-based/half-open bed format  --- DEBUG
cd -
fi


# distribution of fragment sizes: # GIVE A LIST OF BAM FILES AND GET A PAGE/ FILE IN THE FINAL PDF

if false; then # to do on bed files  # sbatch --mem-per-cpu 8000 -q fast -p common,dedicated --array=0-9 script.sh
mkdir -p pdf

cd R/
#files=($(ls ../bed/*/*.bed)) 
#file=${files[SLURM_ARRAY_TASK_ID]}
file="../bed/C2C12_I2/C2C12_I2.bed"  # debug
Rscript fragment_size_distribution.R $file
cd - 

fi




# cut sites and midpoints of bed files 

if false; then
#sbatch --array=0-11 --mem-per-cpu=6000 -c 2 -q fast -p common,dedicated script.sh

cd bam
files=($(ls *_full.bam))
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=$(basename $file _full.bam)
bed="../bed/${NAME}/${NAME}.bed"
mkdir -p ../bed/${NAME}

#cut sites
if false;then
awk '{ print $1 "\t"$2 "\t"$2+1"\t" $4}' $bed > tmp_${SLURM_ARRAY_TASK_ID}
awk '{ print $1 "\t"$3-1 "\t"$3"\t" $4}' $bed | \
cat tmp_${SLURM_ARRAY_TASK_ID} - | sort -k1,1 -k2,2n > ../bed/${NAME}/${NAME}_cut.bed
fi

# midpoints
bed="../bed/${NAME}/${NAME}.bed"
awk '{ if(($4%2)==1){mp=$2+(($4+1)/2); print $1 "\t" mp "\t" mp+1 "\t" $4}else{ mp=$2+($4/2); print $1 "\t" mp "\t" mp+1 "\t" $4} }' $bed  \
| sort -V -k1,1 -k2,2n > ../bed/${NAME}/${NAME}_midpoint.bed 

# midpoints for small fragments 0-120
awk '$4<121' ../bed/${NAME}/${NAME}_midpoint.bed > ../bed/${NAME}/${NAME}_midpoint_0-120.bed
# midpoints for nucleosomes 180-250
awk '$4>139 && $4<251' ../bed/${NAME}/${NAME}_midpoint.bed > ../bed/${NAME}/${NAME}_midpoint_140-250.bed   # old way
awk '$4>179 && $4<251' ../bed/${NAME}/${NAME}_midpoint.bed > ../bed/${NAME}/${NAME}_midpoint_180-250.bed

fi


# merging of bam files 

if false;then 
#sbatch --array=0-1 --mem-per-cpu=6000 -c 2 -q fast -p common,dedicated script.sh
cd bam

files=($(ls C2C12_*2_full.bam)); 
file=${files[SLURM_ARRAY_TASK_ID]}
echo "now file ${file}"
samtools merge -f tmp_${SLURM_ARRAY_TASK_ID} ${file} ${file/2_full/3_full}
samtools sort tmp_${SLURM_ARRAY_TASK_ID} > ${file/2_full/_merged_full}
samtools index ${file/2_full/_merged_full}; #rm tmp_${SLURM_ARRAY_TASK_ID}
fi


# peak calling on bam files (small and full fragments , separate and merged replicates)

if false; then
#sbatch --array=0-21 -q fast -p common,dedicated --mem-per-cpu=4000 script.sh  # 10 samples * 2 types of frag length + 2 merged C2C12 files
mkdir -p peaks/q0.05

cd bam
#files=(M*.bam)  #human
#files=($(ls C*.bam  T*.bam)) # mouse
files=($(ls C*.bam)) # mouse
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=$(basename "${file}" .bam)

#macs2 callpeak -t $file -g hs -n $NAME -q 0.05 --outdir ../peaks/q0.05  #Human
macs2 callpeak -t $file -g mm -n $NAME -q 0.05 --outdir ../peaks/q0.05  #mouse

cd -
fi


# merging of bed files 

if false;then # sbatch -q fast -p common,dedicated --mem-per-cpu=64000 --array=0-1 script.sh

mkdir -p bed/C2C12_I_merged; mkdir -p bed/C2C12_M_merged; 

if [[ $SLURM_ARRAY_TASK_ID == 0 ]]; then
for file in bed/C2C12_I2/*midpoint*.bed; do
NAME=$(basename "${file}")
cat ${file} bed/C2C12_I3/${NAME/I2/I3} > bed/C2C12_I_merged/tmp_${SLURM_ARRAY_TASK_ID} 
sort -V -k1,1 -k2,2n bed/C2C12_I_merged/tmp_${SLURM_ARRAY_TASK_ID} > bed/C2C12_I_merged/${NAME/I2/I_merged}
done
fi
if [[ $SLURM_ARRAY_TASK_ID == 1 ]]; then
for file in bed/C2C12_M2/*midpoint*.bed; do
NAME=$(basename "${file}")
cat ${file} bed/C2C12_M3/${NAME/M2/M3} > bed/C2C12_M_merged/tmp_${SLURM_ARRAY_TASK_ID}
sort -V -k1,1 -k2,2n bed/C2C12_M_merged/tmp_${SLURM_ARRAY_TASK_ID} > bed/C2C12_M_merged/${NAME/M2/M_merged}
done
fi

fi

# formatting ATAC peak files

if false; then
mkdir -p bed; mkdir -p bed/peaks
cd bed/peaks
files=($(ls ../../peaks/*_full_summits.bed)) # more stringent
l=${#files[*]}
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=$(basename $file .bed)
awk '$1!="chrM" {print $1 "\t"$2-1000 "\t"$3+1000}' $file | sort -V -k1,1 -k2,2n > ${NAME}_1kb.bed
cd -
fi

# p300 peaks +/- 1kb

if false; then
mkdir -p analysis
cd analysis/
p300="p300_mm10_encode_opt_idr_hmm_enc.bed"  # file from p300 peaks/chromhmm intersection stored in ~/Epic_gaia/resource/mm10/enhancers/
awk  'NR>1{summit=$2+$10-1; print $1"\t"summit-1000"\t"summit+1001}' $p300 | sort -V -k1,1 -k2,2n > p300_enh_1kb.bed
fi


# V-plot at peaks present in all interphase samples

if false; then # part 1 
# sbatch --array=0-9 script.sh
cd analysis/

mkdir -p tmp

#midpoints

files=($(ls ../bed/*/*_midpoint.bed)) 
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=$(basename $file .bed)

peaks=$(basename $file _midpoint.bed)
if [[ "$peaks" =~ "C2C12_M" ]]; then peaks="C2C12_I2"; fi
if [[ "$peaks" =~ "MCF10a_M" ]]; then peaks="MCF10a_I"; fi
if [[ "$peaks" =~ "MCF7_M" ]]; then peaks="MCF7_I_Vallot"; fi
if [[ "$peaks" =~ "12h_M" ]]; then peaks="Tg2a_Noco_12h_I"; fi
peaks=${peaks}_full

# midpoints
if false;then
# all ATAC peaks 
bedtools intersect -a ${file} -b ../bed/peaks/${peaks}_summits_1kb.bed -wa -wb > tmp_${SLURM_ARRAY_TASK_ID}
awk '{print $0 "\t"$3-$6}' tmp_${SLURM_ARRAY_TASK_ID} > tmp/${NAME}_peak.bed

# p300 all
bedtools intersect -a ${file} -b enhancers.bed -wa -wb > tmp_${SLURM_ARRAY_TASK_ID}
awk '{print $0 "\t"$3-$6}' tmp_${SLURM_ARRAY_TASK_ID} > tmp/${NAME}_p300.bed

# p300 overlapping with ATAC peaks
bedtools intersect -a enhancers.bed -b ../bed/peaks/${peaks}_summits_1kb.bed -wa > tmp_${SLURM_ARRAY_TASK_ID}
bedtools intersect -a ${file} -b tmp_${SLURM_ARRAY_TASK_ID} -wa -wb | \
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$3-$6}' > tmp/${NAME}_peak_p300.bed
fi

# CTCF motif all
bedtools intersect -a ${file} -b Ctcf_all_mm10.bed -wa -wb > tmp_${SLURM_ARRAY_TASK_ID}
awk '{if($8=="+"){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$3-$6"\t+"}else{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"2001-($3-$6)"\t-"}}' tmp_${SLURM_ARRAY_TASK_ID} > tmp/${NAME}_Ctcf.bed

# CTCF overlapping with ATAC peaks
bedtools intersect -a Ctcf_all_mm10.bed -b ../bed/peaks/${peaks}_summits_1kb.bed -wa > tmp_${SLURM_ARRAY_TASK_ID}
bedtools intersect -a ${file} -b tmp_${SLURM_ARRAY_TASK_ID} -wa -wb | \
awk '{if($8=="+"){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$3-$6"\t+"}else{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"2001-($3-$6)"\t-"}}' > tmp/${NAME}_peak_Ctcf.bed

# cutsites
if false; then
files=($(ls ../bed/*/*_cut.bed)) 
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=$(basename $file .bed)

# all ATAC peaks 
bedtools intersect -a ${file} -b ../bed/peaks/${peaks}_summits_1kb.bed -wa -wb > tmp_${SLURM_ARRAY_TASK_ID}
awk '{print $0 "\t"$3-$6}' tmp_${SLURM_ARRAY_TASK_ID} > tmp/${NAME}_peak.bed

# p300 all
bedtools intersect -a ${file} -b enhancers.bed -wa -wb > tmp_${SLURM_ARRAY_TASK_ID}
awk '{print $0 "\t"$3-$6}' tmp_${SLURM_ARRAY_TASK_ID} > tmp/${NAME}_p300.bed

# p300 overlapping with ATAC peaks
bedtools intersect -a enhancers.bed -b ../bed/peaks/${peaks}_summits_1kb.bed -wa > tmp_${SLURM_ARRAY_TASK_ID}
bedtools intersect -a ${file} -b tmp_${SLURM_ARRAY_TASK_ID} -wa -wb | \
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$3-$6}' > tmp/${NAME}_peak_p300.bed
fi
# CTCF motif all
bedtools intersect -a ${file} -b Ctcf_all_mm10.bed -wa -wb > tmp_${SLURM_ARRAY_TASK_ID}
awk '{if($8=="+"){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$3-$6"\t+"}else{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"2001-($3-$6)"\t-"}}' tmp_${SLURM_ARRAY_TASK_ID} > tmp/${NAME}_Ctcf.bed

# CTCF overlapping with ATAC peaks
bedtools intersect -a Ctcf_all_mm10.bed -b ../bed/peaks/${peaks}_summits_1kb.bed -wa > tmp_${SLURM_ARRAY_TASK_ID}
bedtools intersect -a ${file} -b tmp_${SLURM_ARRAY_TASK_ID} -wa -wb | \
awk '{if($8=="+"){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$3-$6"\t+"}else{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"2001-($3-$6)"\t-"}}' > tmp/${NAME}_peak_Ctcf.bed

cd -
fi


if false; then  # part 2
# sbatch -q fast -p common,dedicated --array=0-9 script.sh 
cd R
mkdir -p txt/
files=($(ls ../bam/*_full.bam)) 
l=${#files[*]}

file=${files[SLURM_ARRAY_TASK_ID]}
sample=$(basename $file _full.bam)

#for category in midpoint_peak midpoint_p300 midpoint_peak_p300 midpoint_Ctcf midpoint_peak_Ctcf; do
for category in midpoint_Ctcf midpoint_peak_Ctcf; do
Rscript Vplot.R $sample $category
done

cd -
fi



# metaplots with normalization to total number of reads (and to small fragments for comparison) 

if false; then

mkdir -p pdf/metaplots
cd pdf/metaplots; 
for i in peaks tss ctcf p300-ATAC-overlap p300-ChromHMM-overlap_Owens2019; do mkdir -p $i; done
cd ../../R

if false;then
#sbatch -q fast -p common,dedicated --array=1-4 --mem-per-cpu 24000 -c 2 script.sh
# 1:peaks; 2:tss; 3:ctcf; 4:p300 method p300/ATAC overlap; 5:p300 method Owens2019 p300/ChromHMM
#Rscript metaplots.R $SLURM_ARRAY_TASK_ID C2C12_I2 C2C12_I3 C2C12_M2 C2C12_M3
#Rscript metaplots.R $SLURM_ARRAY_TASK_ID MCF10a_I MCF10a_M
Rscript metaplots.R $SLURM_ARRAY_TASK_ID MCF7_I_Vallot MCF7_M_Vallot
#Rscript metaplots.R $SLURM_ARRAY_TASK_ID Tg2a_Noco_12h_I Tg2a_Noco_12h_M
fi

if true;then # for debugging --- to delete
#sbatch -q fast -p common,dedicated --array=0-3 --mem-per-cpu 24000 -c 2 script.sh
option=4
if [[ "$SLURM_ARRAY_TASK_ID" == 0 ]]; then  
Rscript metaplots.R $option C2C12_I2 C2C12_I3 C2C12_M2 C2C12_M3
fi
if [[ "$SLURM_ARRAY_TASK_ID" == 1 ]]; then  
Rscript metaplots.R $option MCF10a_I MCF10a_M
fi
if [[ "$SLURM_ARRAY_TASK_ID" == 2 ]]; then  
Rscript metaplots.R $option MCF7_I_Vallot MCF7_M_Vallot
fi
if [[ "$SLURM_ARRAY_TASK_ID" == 3 ]]; then  
Rscript metaplots.R $option Tg2a_Noco_12h_I Tg2a_Noco_12h_M
fi
fi

fi



# V-plots/metaplots at CTCF sites


if false;then


ctcf.folder="~/Epic_gaia/beds/"

fi





###########################################################
###########################################################
###########################################################
### OLD 
###########################################################
###########################################################
###########################################################
##########################################################
###########################################################


# Promoter file 

if false; then # not threaded
mdkir -p analysis
cd analysis/
# promoters - UCSC Table browser -- NCBI genes - 1 pb upstream
zcat promoters.bed.gz | awk '$1!="chrM"&&length($1)<6{if($6=="+"){print $1 "\t" $3-999-1 "\t" $3+1001 "\t+"}else{print $1 "\t" $3-1001 "\t" $3+999-1 "\t-"}}' > tmp.bed  
sortBed -i tmp.bed | uniq > promoters.bed # eliminates redundant promoters

cd -
fi

# ChromHMM/ enhancer file 

if false; then # not threaded
mdkir -p analysis
cd analysis/
## wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed.gz
## wget https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
liftOver wgEncodeBroadHmmGm12878HMM.bed.gz hg19ToHg38.over.chain.gz wgEncodeBroadHmmGm12878HMM_hg38.bed.gz unMapped.bed.gz

if false; then
zcat wgEncodeBroadHmmGm12878HMM_hg38.bed.gz |
#awk '{ if(($4%2)==1){mp=$2+(($4+1)/2); print $1 "\t" mp "\t" mp+1 "\t" $4}else{ mp=$2+($4/2); print $1 "\t" mp "\t" mp+1 "\t" $4} }' \
awk '{ if((($3-$2)%2)==1){mp=(($2+$3+1)/2); print $1 "\t" mp-1 "\t" mp "\t" $4}else{ mp=(($2+$3)/2); print $1 "\t" mp-1 "\t" mp "\t" $4} }' \
| sort -k1,1 -k2,2n > wgEncodeBroadHmmGm12878HMM_hg38_midpoint.bed
fi

cd -
fi


# intersection of peaks common to all interphase cancer samples  -- small fragments only

if false; then 
# sbatch --job-name="intersect" script.sh
mkdir -p bed/intersections

peaks=($(ls peaks/*I[12]_0-140_peaks.narrowPeak))
summits=($(ls peaks/*I[12]_0-140_summits.bed))
l=${#summits[*]}

for (( i=0; i<$l; i++ )); do 
NAME=$(basename "${summits[i]}" _summits.bed)
paste ${peaks[i]} <(awk '{print $1 "\t"$2 "\t"$3}' ${summits[i]}) > bed/intersections/${NAME}_peaks_summits.narrowPeak
done
cd bed/intersections
files=(`ls *narrowPeak`)

cp ${files[0]} tmp 
for (( i=1; i<$l; i++ )); do # file with common peaks in all interphase samples
bedtools intersect -a tmp -b ${files[i]} -wa -wb > tmp2 && mv tmp2 tmp
done
mv tmp intersect_all_I.bed

for (( i=0; i<$l; i++ )); do # overlaps of peak summit in each sample with +/-1000bp region around common peaks regions
file=${files[i]}
NAME=$(basename $file _peaks_summits.narrowPeak)
nb1=$((11 + 13*i)); nb2=$((nb1+1)); nb3=$((nb1+2)); 
awk -v col1="$nb1" -v col2="$nb2" -v col3="$nb3" '{print $col1 "\t" $col2 "\t" $col3}' intersect_all_I.bed > intersect_all_I_${NAME}_summits.bed
awk '$1!="chrM"{print $1 "\t" $2-1000 "\t" $3+1000}' intersect_all_I_${NAME}_summits.bed  > intersect_all_I_${NAME}_intervals.bed
done

cd -
fi


# V-plot at peaks present in all interphase samples

if false; then # part 1 
# sbatch --array=0-9 script.sh
cd analysis/

mkdir -p tmp

#midpoints
files=($(ls ../bed/*/*_midpoint.bed)) 
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=$(basename $file .bed)

if [[ "$NAME" =~ "A549_" ]]; then peaks="A549_I1"; fi
if [[ "$NAME" =~ "HCT116_" ]]; then peaks="HCT116_I1"; fi
if [[ "$NAME" =~ "MCF7_" ]]; then peaks="MCF7_I1"; fi

bedtools intersect -a ${file} -b ../bed/intersections/intersect_all_I_${peaks}_intervals.bed -wa -wb > tmp_${SLURM_ARRAY_TASK_ID}
awk '{print $0 "\t"$3-$6}' tmp_${SLURM_ARRAY_TASK_ID} > tmp/${NAME}_peak.bed

#cutsites
files=($(ls ../bed/*/*_cut.bed)) 
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=$(basename $file .bed)

bedtools intersect -a ${file} -b ../bed/intersections/intersect_all_I_${peaks}_intervals.bed -wa -wb > tmp_${SLURM_ARRAY_TASK_ID}
awk '{print $0 "\t"$3-$6}' tmp_${SLURM_ARRAY_TASK_ID} > tmp/${NAME}_peak.bed
rm tmp_${SLURM_ARRAY_TASK_ID} 

cd -
fi

if false; then  # part 2
# sbatch script.sh 
cd R
Rscript Vplot.R

cd -
fi














######################################################################################
######################################################################################
######################################################################################
######################################################################################
########################### OLD stuff
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
############################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################################################
############################################################################################################################################################################



# ATAC profiles /metaplots at interphase peaks (either full fragments size or cut sites)

if false;then # sbatch --array=0-2 script.sh
cd analysis/

if [ "${SLURM_ARRAY_TASK_ID}" == 0 ]; then # ATAC profiles at peaks for A549 - small fragments
computeMatrix reference-point --referencePoint center -R ../bed/intersections/intersect_all_I_A549_I1_0-140_summits.bed \
   -S ../bw/A549_I1_0-140/A549_I1_0-140.bw ../bw/A549_I2_0-140/A549_I2_0-140.bw \
   ../bw/A549_M_Filtered_0-140/A549_M_Filtered_0-140.bw ../bw/A549_M_Not_filtered_0-140/A549_M_Not_filtered_0-140.bw \
   -out Matrix_peaks_A549_0-140 -a 1000 -b 1000

plotProfile --matrixFile Matrix_peaks_A549_0-140 -out Plot_peaks_A549_0-140 --perGroup -T "ATAC small fragments at peaks (A549)" --refPointLabel "0" --regionsLabel " " --samplesLabel I1 I2 M_Filtered M_Not_filtered
fi
if [ "${SLURM_ARRAY_TASK_ID}" == 1 ]; then # ATAC profiles at peaks for HCT116 - small fragments
computeMatrix reference-point --referencePoint center -R ../bed/intersections/intersect_all_I_HCT116_I1_0-140_summits.bed \
   -S ../bw/HCT116_I1_0-140/HCT116_I1_0-140.bw ../bw/HCT116_M1_0-140/HCT116_M1_0-140.bw \
   -out Matrix_peaks_HCT116_0-140 -a 1000 -b 1000

plotProfile --matrixFile Matrix_peaks_HCT116_0-140 -out Plot_peaks_HCT116_0-140 --perGroup -T "ATAC small fragments at peaks (HCT116)" --refPointLabel "0" --regionsLabel " " --samplesLabel I1 M1
fi
if [ "${SLURM_ARRAY_TASK_ID}" == 2 ]; then # ATAC profiles at peaks for MCF7 - small fragments
computeMatrix reference-point --referencePoint center -R ../bed/intersections/intersect_all_I_MCF7_I1_0-140_summits.bed \
   -S ../bw/MCF7_I1_0-140/MCF7_I1_0-140.bw ../bw/MCF7_I2_0-140/MCF7_I2_0-140.bw \
   ../bw/MCF7_M1_0-140/MCF7_M1_0-140.bw ../bw/MCF7_M2_0-140/MCF7_M2_0-140.bw \
   -out Matrix_peaks_MCF7_0-140  -a 1000 -b 1000

plotProfile --matrixFile Matrix_peaks_MCF7_0-140 -out Plot_peaks_MCF7_0-140 --perGroup -T "ATAC small fragments at peaks (MCF7)" --refPointLabel "0" --regionsLabel " " --samplesLabel I1 I2 M1 M2
fi

cd -

fi


# ATAC profiles at promoters, enhancers and peaks that do not overlap with any of these two

if false;then # sbatch --array=0-8 script.sh
cd analysis/

if [ "${SLURM_ARRAY_TASK_ID}" == 0 ]; then # ATAC profiles at promoters for A549
computeMatrix reference-point --referencePoint center -R promoters.bed \
   -S ../bw/A549_I1/A549_I1_0-140.bw ../bw/A549_I2/A549_I2_0-140.bw \
   ../bw/A549_M_Filtered/A549_M_Filtered_0-140.bw ../bw/A549_M_Not_filtered/A549_M_Not_filtered_0-140.bw \
   -out Matrix_promoters_A549 
plotProfile --matrixFile Matrix_promoters_A549 -out Plot_promoters_A549 --perGroup -T "ATAC at promoters (A549)"  -startLabel "-1000" -endLabel "+1000" --samplesLabel I1 I2 M_Filtered M_Not_filtered
fi
if [ "${SLURM_ARRAY_TASK_ID}" == 1 ]; then # ATAC profiles at promoters for HCT116
computeMatrix reference-point --referencePoint center -R promoters.bed \
   -S ../bw/A549_I1/A549_I1_0-140.bw ../bw/A549_I2/A549_I2_0-140.bw \
   ../bw/A549_M_Filtered/A549_M_Filtered_0-140.bw ../bw/A549_M_Not_filtered/A549_M_Not_filtered_0-140.bw \
   -out Matrix_promoters_A549 
plotProfile --matrixFile Matrix_promoters_A549 -out Plot_promoters_A549 --perGroup -T "ATAC at promoters (A549)"  -startLabel "-1000" -endLabel "+1000" --samplesLabel I1 I2 M_Filtered M_Not_filtered
fi

cd -

fi


# debug -- filter for fragments below 40 pb

if false; then
# sbatch --array=0-2 script.sh   # three files to process
cd bam/Sample_A549_I2   # example files
files=($(ls *.bam))
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=$(basename $file .bam)
samtools view -H ${file} > header_${SLURM_ARRAY_TASK_ID}
samtools view ${file} | awk '($9*$9)<1681' | cat header_${SLURM_ARRAY_TASK_ID} - | samtools view -bS - > ${NAME}_0-40.bam # fragments lengths all positive in new pipeline
samtools sort ${NAME}_0-40.bam > tmp_${SLURM_ARRAY_TASK_ID}.bam && mv tmp_${SLURM_ARRAY_TASK_ID}.bam ${NAME}_0-40.bam 
samtools index ${NAME}_0-40.bam 

rm header_${SLURM_ARRAY_TASK_ID}
cd -
fi
if false; then
cd bam
files=($(ls *_full.bam)) 
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=$(basename $file _full.bam)
samtools view -H ${file} > header_${SLURM_ARRAY_TASK_ID}
samtools view ${file} | awk '$9<41' | cat header_${SLURM_ARRAY_TASK_ID} - | samtools view -bS - > ${NAME}_0-40.bam # fragments lengths all positive in new pipeline
samtools sort ${NAME}_0-40.bam > tmp_${SLURM_ARRAY_TASK_ID}.bam && mv tmp_${SLURM_ARRAY_TASK_ID}.bam ${NAME}_0-40.bam 
samtools index ${NAME}_0-40.bam 

rm header_${SLURM_ARRAY_TASK_ID}
cd -
fi


# filter for fragments sizes on bam files instead of bed files

if false; then # old -- from non deduplicated bam files
cd bam
files=($(ls *_full.bam)) 
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=$(basename $file _full.bam)
samtools view -H ${file} > header_${SLURM_ARRAY_TASK_ID}
samtools view ${file} | awk '($9*$9)<10001' | cat header_${SLURM_ARRAY_TASK_ID} - | samtools view -bS - > ${NAME}_0-100.bam 
samtools view ${file} | awk '($9*$9)>10000 && ($9*$9)<19601' | cat header_${SLURM_ARRAY_TASK_ID} - | samtools view -bS - > ${NAME}_101-140.bam 
samtools view ${file} | awk '($9*$9)>19600 && ($9*$9)<40001' | cat header_${SLURM_ARRAY_TASK_ID} - | samtools view -bS - > ${NAME}_141-200.bam 
rm header_${SLURM_ARRAY_TASK_ID}
cd -
fi


# filter for small fragments (0-150bp) 

if false; then
# sbatch --array=0-1 script.sh

samtools view -H ${file}_full.bam > header_${SLURM_ARRAY_TASK_ID}
samtools view ${file}_full.bam | awk '($9*$9)<22500' | cat header_${SLURM_ARRAY_TASK_ID} - | samtools view -bS - > ${file}_0-150.bam 
rm header_${SLURM_ARRAY_TASK_ID}
#mv ${cdir}/${NAME}_0-150.bam $adir
#mv ${cdir}/${NAME}_full.bam $adir  # keep version with all fragments  -- BEFORE LAUNCHING THIS SCRIPT CHANGE FILENAME TO _full.bam
#cd -
fi


# sorting and indexing

if false; then  
# sbatch --array=0-3 script.sh

cd bam

files=($(ls *.bam))
file=${files[SLURM_ARRAY_TASK_ID]}
samtools sort ${file} tmp_${SLURM_ARRAY_TASK_ID} && mv tmp_${SLURM_ARRAY_TASK_ID}.bam ${file}
samtools index ${file}
cd -
fi


# bigwig files

if false; then
# sbatch --array=0-22 script.sh # 46? files here, 2 per/condition (full and 0-150)
cd bam
files=($(ls *.bam))
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=$(basename $file .bam)

mkdir -p ../bed/${NAME}
samtools view ${file} | \
awk '{if($9>0){print $3"\t"$4-1"\t"$4+$9-1}else{print $3"\t"$8-1"\t"$8-$9-1} }' | \
sortBed -i - > ../bed/${NAME}/${NAME}.bed 
#awk '{if($9>0){print $3"\t"$4+1"\t"$4+1+$9-1}else{print $3"\t"$8+1"\t"$8+1-$9-1} }' | \

mkdir -p ../bedGraph/${NAME}
bedtools genomecov -i ../bed/${NAME}/${NAME}.bed -g ~/utils/chrom.sizes.hg38 -dz > ../bedGraph/${NAME}/tmp_${NAME} # -dz uses 0 parity on start & end position - reason for modifs in the awk
awk 'col3=$2+1 {print $1 "\t" $2 "\t" col3 "\t" $3}' ../bedGraph/${NAME}/tmp_${NAME} \
> ../bedGraph/${NAME}/${NAME}.bedGraph 

mkdir -p ../bw/${NAME}
bedGraphToBigWig ../bedGraph/${NAME}/${NAME}.bedGraph ~/utils/chrom.sizes.hg38 ../bw/${NAME}/${NAME}.bw
#rm ../bedGraph/${NAME}/tmp_${NAME} ../bedGraph/${NAME}/${NAME}.bedGraph

cd -
fi


# bigwigs for regions -/+4pb around cut positions/extremities of the fragments

if false; then 
# sbatch --array=0-3 script.sh # 4 files here, 2 per/condition (150+ and 0-150)
cd bam
files=($(ls *.bam))
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=$(basename $file .bam)

samtools view ${file} | \
awk '$3!="chrM"{if($9>0){print $3"\t"$4-1-4"\t"$4+4}else{print $3"\t"$8-1-4"\t"$8+4} }' > tmp_${SLURM_ARRAY_TASK_ID}

samtools view ${file} | \
awk '$3!="chrM"{if($9>0){print $3"\t"$4-1+$9-1-4"\t"$4+$9-1+4}else{print $3"\t"$8-1-$9-1-4"\t"$8-$9-1+4} }' | \
cat tmp_${SLURM_ARRAY_TASK_ID} - | sortBed -i - > ../bed/${NAME}/${NAME}_cut.bed # no merge anymore? 
#cat tmp_$SLURM_ARRAY_TASK_ID - | sortBed -i - | mergeBed -i - > ../bed/${NAME}/${NAME}_cut.bed
rm tmp_$SLURM_ARRAY_TASK_ID

bedtools genomecov -i ../bed/${NAME}/${NAME}_cut.bed -g ~/utils/chrom.sizes.hg38 -dz > ../bedGraph/${NAME}/tmp_${NAME}
# -dz uses 0 parity on start & end position - reason for modifs in the awk
awk 'col3=$2+1 {print $1 "\t" $2 "\t" col3 "\t" $3}' ../bedGraph/${NAME}/tmp_${NAME} \
> ../bedGraph/${NAME}/${NAME}_cut.bedGraph 

bedGraphToBigWig ../bedGraph/${NAME}/${NAME}_cut.bedGraph ~/utils/chrom.sizes.hg38 ../bw/${NAME}/${NAME}_cut.bw

rm ../bed/${NAME}/${NAME}_cut.bed ../bedGraph/${NAME}/tmp_${NAME} ../bedGraph/${NAME}/${NAME}_cut.bedGraph # TO DO AFTER THE FIRST ROUND
cd -
fi

# peak calling for broad peaks on each file

if false; then 
cd bam

for filetype in 0-150.bam full.bam; do
files=($(ls *${filetype}))
file=${files[SLURM_ARRAY_TASK_ID]}
NAME=$(basename $file .bam)
macs2 callpeak -t $file -g mm -n $NAME --broad -q 0.1 --outdir ./peaks/

done
cd -
fi

# merging all peaks -- we don't filter for peaks with 2 replicates in the same conditions cause we check if replication ok

if false; then
cat peaks/*0-150_peaks.broadPeak > peaks/all_merged_peaks_0-150.bed
cat peaks/*full_peaks.broadPeak > peaks/all_merged_peaks_full.bed

sort -k1,1 -k2,2n -k3,3n peaks/all_merged_peaks_0-150.bed | bedtools merge -i - > tmp_${SLURM_ARRAY_TASK_ID} && mv tmp_${SLURM_ARRAY_TASK_ID} peaks/all_merged_peaks_0-150.bed
sort -k1,1 -k2,2n -k3,3n peaks/all_merged_peaks_full.bed | bedtools merge -i - > tmp_${SLURM_ARRAY_TASK_ID} && mv tmp_${SLURM_ARRAY_TASK_ID} peaks/all_merged_peaks_full.bed

cat peaks/all_merged_peaks_0-150.bed | awk '{$3=$3"\t""peak_"NR}1' OFS="\t" > tmp_${SLURM_ARRAY_TASK_ID} && mv tmp_${SLURM_ARRAY_TASK_ID} peaks/all_merged_peaks_0-150.bed
cat peaks/all_merged_peaks_full.bed | awk '{$3=$3"\t""peak_"NR}1' OFS="\t" > tmp_${SLURM_ARRAY_TASK_ID} && mv tmp_${SLURM_ARRAY_TASK_ID} peaks/all_merged_peaks_full.bed

fi

# counting reads at peaks 

if false; then # using bedops
cd bam

# if error try to pre-process with: 
# https://github.com/crazyhottommy/ChIP-seq-analysis/blob/master/part2_Preparing-ChIP-seq-count-table.md

bedtools multicov -f 0.5 -r -bams *0-150.bam -bed ../peaks/all_merged_peaks_0-150.bed > tmp_${SLURM_ARRAY_TASK_ID} && mv tmp_${SLURM_ARRAY_TASK_ID} ../peaks/all_merged_peaks_0-150.bed
bedtools multicov -f 0.5 -r -bams *full.bam -bed ../peaks/all_merged_peaks_full.bed > tmp_${SLURM_ARRAY_TASK_ID} && mv tmp_${SLURM_ARRAY_TASK_ID} ../peaks/all_merged_peaks_full.bed

cd -

# OLD VERSION 
#for filetype in 0-150 150+; do
#files=($(ls *${filetype}.bam))
#file=${files[SLURM_ARRAY_TASK_ID]}
#NAME=$(basename $file .bam)
#
#bam2bed < ${file} | bedmap --echo --count --fraction-map 0.51 peaks/all_merged_peaks_${filetype}.bed - > peaks.txt
#
#done
##awk -F'[\t|]' '{print $1 "\t" $2 "\t" $3}' bw2_sort_hc_uniq_ATAC_A_1.txt > bw2_sort_hc_uniq_ATAC_A_1.bed
##awk -F'[\t|]' '{print $1 "\t" $2 "\t" $3}' bw2_sort_hc_uniq_ATAC_A_2.txt > bw2_sort_hc_uniq_ATAC_A_2.bed
##awk -F'[\t|]' '{print $1 "\t" $2 "\t" $3}' bw2_sort_hc_uniq_ATAC_M_1.txt > bw2_sort_hc_uniq_ATAC_M_1.bed
##awk -F'[\t|]' '{print $1 "\t" $2 "\t" $3}' bw2_sort_hc_uniq_ATAC_M_2.txt > bw2_sort_hc_uniq_ATAC_M_2.bed
#cd -
fi

if false; then # scaling factor based on small fragments 
cd ../data/atac_Teves_data/
#a=`samtools view ../data/atac_Teves_data/bw2_sort_hc_uniq_ATAC_A_1_final.bam | wc -l`
#b=`samtools view ../data/atac_Teves_data/bw2_sort_hc_uniq_ATAC_A_2_final.bam | wc -l`
#c=`samtools view ../data/atac_Teves_data/bw2_sort_hc_uniq_ATAC_M_1_final.bam | wc -l`
#d=`samtools view ../data/atac_Teves_data/bw2_sort_hc_uniq_ATAC_M_2_final.bam | wc -l`

for filetype in 0-150 150+; do
files=($(ls *${filetype}.bam))
depths=(`samtools view -c ${files[0]}` `samtools view -c ${files[1]}` `samtools view -c ${files[2]}` `samtools view -c ${files[3]}`) 
min=`printf "%d\n" "${depth[@]}" | sort -rn | tail -1`
awk -v min="$min" -v a1="${depths[0]}" a2="${depths[1]}" m1="${depths[2]}" m2="${depths[3]}" -F'[\t|]' \
'{print $1 "\t" $2 "\t" $3 "\t" $4*min/a1 "\t" $5*min/a2 "\t" $6*min/m1 "\t" $7*min/m2}' peaks/all_merged_peaks_${filetype}.bed \
> tmp && mv tmp peaks/all_merged_peaks_${filetype}.bed


done

#rm test1.txt test2.txt test3.txt test4.txt bw2_sort_hc_uniq_ATAC_A_1.txt bw2_sort_hc_uniq_ATAC_A_2.txt bw2_sort_hc_uniq_ATAC_M_1.txt bw2_sort_hc_uniq_ATAC_M_2.txt
cd -
fi


# correlation table --> using package R package DiffBind (clearly old)? 



# correlation plots for each replicates (deep tools?) -- at all peaks, promoters and enhancers? 



# PCA at all merged peaks

if false; then 
Rscript pca.R count_table.bed >wig.Rout 2>wig.Rerr
  
fi

# Comparing replicates - @ active promoters

if false; then
cp /path/to/folder/where/ChromHMM/regions/file/is/stored/hg38/tss/hg38_tss_ens_93_hmm_act_prom.tsv tmp.bed # intersection of ChromHMM active promoters and Ensembl data version 93
awk 'FNR>1 {print $1 "\t" $2 "\t" $3 "\t" $4}' tmp.bed | sort -k1,1 -k2,2n -k3,3n > active_prom.bed 
rm tmp.bed
bedtools intersect -a ../data/atac_Teves_data/peaks/bw2_sort_hc_uniq_ATAC_A_1_peaks_peaks.broadPeak -b active_prom.bed -wa -wb > ../data/atac_Teves_data/peaks/intersections/ATAC_A_1_vs_active_prom.bed # intersections ATAC peaks - active promoters
bedtools intersect -a ../data/atac_Teves_data/peaks/bw2_sort_hc_uniq_ATAC_A_2_peaks_peaks.broadPeak -b active_prom.bed -wa -wb > ../data/atac_Teves_data/peaks/intersections/ATAC_A_2_vs_active_prom.bed
bedtools intersect -a ../data/atac_Teves_data/peaks/bw2_sort_hc_uniq_ATAC_M_1_peaks_peaks.broadPeak -b active_prom.bed -wa -wb > ../data/atac_Teves_data/peaks/intersections/ATAC_M_1_vs_active_prom.bed
bedtools intersect -a ../data/atac_Teves_data/peaks/bw2_sort_hc_uniq_ATAC_M_2_peaks_peaks.broadPeak -b active_prom.bed -wa -wb > ../data/atac_Teves_data/peaks/intersections/ATAC_M_2_vs_active_prom.bed
# --> then use Venny to show overlap b. different conditions

Rscript intersect.R ../data/atac_Teves_data/peaks/intersections/ATAC_A_1_vs_active_prom.bed ../data/atac_Teves_data/peaks/intersections/ATAC_A_2_vs_active_prom.bed ../data/atac_Teves_data/peaks/intersections/ATAC_M_1_vs_active_prom.bed ../data/atac_Teves_data/peaks/intersections/ATAC_M_2_vs_active_prom.bed ATAC_samples_overlap_counts_at_active_promoters.txt ATAC_samples_overlap_percent_at_active_promoters.txt

fi

# Comparing replicates - @ enhancers

if false; then
cp /folder/where/enhancer/regions/file/is/stored/hg38/enhancers/p300_hg38_encode_opt_idr_hmm_enc.bed tmp.bed # intersection of ChromHMM enhancers & p300 peaks
awk 'FNR>1 {print $1 "\t" $2 "\t" $3 "\t" $4}' tmp.bed | sort -k1,1 -k2,2n -k3,3n > enhancers.bed 
rm tmp.bed
bedtools intersect -a ../data/atac_Teves_data/peaks/bw2_sort_hc_uniq_ATAC_A_1_peaks_peaks.broadPeak -b enhancers.bed -wa -wb > ../data/atac_Teves_data/peaks/intersections/ATAC_A_1_vs_enhancers.bed # intersections ATAC peaks - active promoters
bedtools intersect -a ../data/atac_Teves_data/peaks/bw2_sort_hc_uniq_ATAC_A_2_peaks_peaks.broadPeak -b enhancers.bed -wa -wb > ../data/atac_Teves_data/peaks/intersections/ATAC_A_2_vs_enhancers.bed
bedtools intersect -a ../data/atac_Teves_data/peaks/bw2_sort_hc_uniq_ATAC_M_1_peaks_peaks.broadPeak -b enhancers.bed -wa -wb > ../data/atac_Teves_data/peaks/intersections/ATAC_M_1_vs_enhancers.bed
bedtools intersect -a ../data/atac_Teves_data/peaks/bw2_sort_hc_uniq_ATAC_M_2_peaks_peaks.broadPeak -b enhancers.bed -wa -wb > ../data/atac_Teves_data/peaks/intersections/ATAC_M_2_vs_enhancers.bed
# --> then use Venny to show overlap b. different conditions

Rscript intersect.R ../data/atac_Teves_data/peaks/intersections/ATAC_A_1_vs_enhancers.bed ../data/atac_Teves_data/peaks/intersections/ATAC_A_2_vs_enhancers.bed ../data/atac_Teves_data/peaks/intersections/ATAC_M_1_vs_enhancers.bed ../data/atac_Teves_data/peaks/intersections/ATAC_M_2_vs_enhancers.bed ATAC_samples_overlap_counts_at_enhancers.txt ATAC_samples_overlap_percent_at_enhancers.txt

fi


# PCA at all active promoters & enhancers

if false; then 

bedtools intersect -a count_table.bed -b active_prom.bed -wa > count_table_at_active_prom_peaks.bed
Rscript pca.R count_table_at_active_prom_peaks.bed >wig.Rout 2>wig.Rerr

bedtools intersect -a count_table.bed -b enhancers.bed -wa > count_table_at_enhancers.bed
Rscript R/pca.R count_table_at_enhancers.bed >wig.Rout 2>wig.Rerr
  
fi

# K-means clustering on the scaled peak value 

if false; then

Rscript R/clustering.R 

fi




# plot profiles of replicates at active promoters

if false; then
awk 'FNR>1 {print $1 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $6}' /folder/where/promoter/regions/file/is/stored/hg38/tss/hg38_tss_ens_93_hmm_act_prom.tsv > active_promoters_genes_coordinates.bed # active genes coordinates 

#if false;then
computeMatrix reference-point -S ../data/atac_Teves_data_Teves_data/wigs/ATAC_interphase_rep1.bw \
                                 ../data/atac_Teves_data_Teves_data/wigs/ATAC_interphase_rep2.bw  \
                                 ../data/atac_Teves_data_Teves_data/wigs/ATAC_mitosis_rep1.bw \
                                 ../data/atac_Teves_data_Teves_data/wigs/ATAC_mitosis_rep2.bw \
                              -R active_promoters_genes_coordinates.bed \
                              --beforeRegionStartLength 1000 \
                              --afterRegionStartLength 1000 \
                              -o profiles_ATAC_active_promoters.mat.gz
#                              --samplesLabel Interphase1 Interphase2 Mitosis1 Mitosis2 \
#fi
# --skipZeros 
plotProfile -m profiles_ATAC_active_promoters.mat.gz \
     --perGroup \
     -out profiles_ATAC_active_promoters.png \
     --colors red pink blue lightblue \
     --plotTitle "ATAC profiles at active promoters"
fi



# genomic distribution of peaks -- using 

if false; then

#awk -F "\t" '{$1="peak_"NR FS$1;$4=$4FS"."}1' > subread.saf
# install featureCounts via http://bioinf.wehi.edu.au/featureCounts/ 
#featureCounts -a subread.saf -F SAF -o counts_subread.txt ../../data/*bam -T 4
echo "test"
fi

####################333 OLD STUFF

if false; then

for file in bw2_sort_hc_uniq_ATAC_A_1 bw2_sort_hc_uniq_ATAC_A_2 bw2_sort_hc_uniq_ATAC_M_1 bw2_sort_hc_uniq_ATAC_M_2 ; do 
bam2bed < ../data/atac_Teves_data/${file}_150+_final.bam > ../data/atac_Teves_data/${file}_150+_final.bed
bedGraphToBigWig ../data/atac_Teves_data/${file}_150+_final.bed ~/chrom.sizes.hg38 ../data/atac_Teves_data/wigs/${file}_150+.bw
#rm ../data/atac_Teves_data/wigs/${file}_150+.bw
done
fi


#depth=($a $b $c $d)
#min=`printf "%d\n" "${depth[@]}" | sort -rn | tail -1`
#awk -v min="$min" -v a1="$a" -F'[\t|]' '{print $1 "\t" $2 "\t" $3 "\t" $4*min/a1}' bw2_sort_hc_uniq_ATAC_A_1.txt > test1.txt
#awk -v min="$min" -v a2="$b" -F'[\t|]' '{print $4*min/a2}' bw2_sort_hc_uniq_ATAC_A_2.txt > test2.txt
#awk -v min="$min" -v m1="$c" -F'[\t|]' '{print $4*min/m1}' bw2_sort_hc_uniq_ATAC_M_1.txt > test3.txt
#awk -v min="$min" -v m2="$d" -F'[\t|]' '{print $4*min/m2}' bw2_sort_hc_uniq_ATAC_M_2.txt > test4.txt
#paste -d'\t' test.txt tes1t2.txt test3.txt test4.txt > count_table.bed








