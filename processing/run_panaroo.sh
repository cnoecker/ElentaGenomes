#!/bin/bash -l
#SBATCH --time=72:00:00
#SBATCH --ntasks=24
#SBATCH --mem=96g
module load python3

genome_dir=good_gffs/elenta/
results_dir=panaroo_out_fixed/

mkdir -p $results_dir

source activate $HOME/conda_envs/panaroo_v1.5.0

panaroo -i $genome_dir/*.gff3 -o $results_dir --remove-invalid-genes -t 24 -a core --clean-mode moderate

conda deactivate


