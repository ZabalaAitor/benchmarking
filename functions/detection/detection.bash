# ------------------------------------------
# Pipeline commands for circular DNA (eccDNA) 
# and circular RNA (circRNA) detection
# ------------------------------------------
#
# The following Nextflow and Python commands run 
# established pipelines and tools to detect eccDNA 
# and circRNA from sequencing data.
#
# For eccDNA:
# - Runs the nf-core/circdna Nextflow pipeline with Apptainer profile.
# - Runs the ecc_finder tool with bwa and minimap2 aligners.
#
# For circRNA:
# - Runs the nf-core/circRNA Nextflow pipeline for both in silico and biological datasets.
#
# Replace placeholders like <<samplesheet.csv>>, <<output_dir>>, 
# <<reference.fasta>>, and <<reference sr.idx>> with your actual file paths.
#
# Adjust CPU and memory resources as needed for your environment.

# eccDNA

nextflow run nf-core/circdna \
	-r 1.1.0 \
	-profile apptainer \
	--input <<samplesheet.csv>> \
	--input_format 'FASTQ' \
	--outdir <<output_dir>> \
	--circle_identifier  ‘circle_map_realign,circle_finder,circexplorer2’ \
	--genome GRCh38 \
	--max_cpus 24 \
	--max_memory 256GB \
	-w <<work_dir>> \
    --save_trimmed \
    -resume 

python run_ecc_finder.py \
  --reference <<reference.fasta>> \
  --reference_idx <<reference sr.idx>> \
  --samplesheet <<samplesheet.csv>> \
  --output_dir <<output_dir>> \
  --aligner bwa
  
python run_ecc_finder.py \
  --reference <<reference.fasta>> \
  --reference_idx <<reference sr.idx>> \
  --samplesheet <<samplesheet.csv>> \
  --output_dir <<output_dir>> \
  --aligner minimap2


# circRNA

## insilico

nextflow run nf-core/circRNA \
	-r  dev \
	-profile apptainer \
	--input <<samplesheet.csv>> \
	--phenotype <<phenotype.csv>> \
	--outdir <<output_dir>>\
	--tool 'ciriquant,circexplorer2,find_circ,circrna_finder,segemehl' \
	--genome GRCh38 \
	--max_cpus 24 \
	--max_memory 256GB \
	--gtf /data/azabala/gtf/UCSC/NCBI/hg38.ncbiRefSeq.gtf  \
	--save_reference false \
	--save_intermediates \
	-resume
	
## biological

nextflow run nf-core/circRNA \
	-r  dev \
	-profile apptainer \
	--input <<samplesheet.csv>> \
	--phenotype <<phenotype.csv>> \
	--outdir <<output_dir>>\
	--tool 'ciriquant,circexplorer2,find_circ,circrna_finder,segemehl' \
	--genome GRCh38 \
	--max_cpus 24 \
	--max_memory 256GB \
	--save_reference false \
	--save_intermediates \
	-resume
	

