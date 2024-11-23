#!/bin/bash
#SBATCH -J sym3_xyk
#SBATCH -p memtest	
#SBATCH -N 1	        
#SBATCH -o stdout.%j    
#SBATCH -e stderr.%j    
#SBATCH --ntasks-per-node=56

module load compilers/intel/oneapi-2023/config

./CTBL_sym3_xyk
