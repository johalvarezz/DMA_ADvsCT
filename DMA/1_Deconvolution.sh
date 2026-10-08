# Clone
git clone https://github.com/nloyfer/UXM_deconv.git
cd UXM_deconv/
# Optionally, add to PATH (or link "./uxm" to some directory in $PATH)
export PATH=${PATH}:$PWD
git clone https://github.com/nloyfer/wgbs_tools.git
cd wgbs_tools

# compile
python setup.py

# Load genome reference
wgbstools init_genome hg38 --fasta_path /data/scratch/NAVA1/JAJ/genomic_reference/genome.fa

# Create the conda environment with the following dependencies: python 3+ ( with: pandas version 1.0+, numpy, scipy), samtools
tabix / bgzip, bedtools.

conda create -n uxm python=3.10 -y
conda activate uxm
pip install pandas numpy scipy
conda install -c bioconda samtools bedtools tabix
cd wgbs_tools
python setup.py

# Creamo el bin y lo ponemos en PATH para poder llamar wgstools desde cualquier directorio
mkdir -p ~/bin
ln -s /data/scratch/NAVA1/JAJ/UXM_deconv/wgbs_to/python/wgbs_tools.py ~/bin/wgbstools
export PATH=$HOME/bin:$PATH
source ~/.bashrc
# para uxm
export PATH=/data/scratch/NAVA1/JAJ/UXM_deconv:$PATH
source ~/.bashrc

# Convertimos .bam files con .pat files
wgbstools bam2pat /data/scratch/NAVA1/Results/ADvsCT/bismark/deduplicated/*.bam -o /data/scratch/NAVA1/JAJ/UXM_deconv/pat_files/
Running deconvolution
$ uxm deconv /data/scratch/NAVA1/JAJ/UXM_deconv/pat_files/*pat.gz -o /data/scratch/NAVA1/JAJ/UXM_deconv/results/output_AD.csv --atlas /data/scratch/NAVA1/JAJ/UXM_deconv/supplemental/Atlas.U25.l4.hg38.full.tsv
