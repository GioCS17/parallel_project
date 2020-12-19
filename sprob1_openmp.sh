#!/bin/bash
#SBATCH -J problem1 # nombre del job
#SBATCH -p investigacion # nombre de la particion 
#SBATCH -N 1 # numero de nodos
#SBATCH --tasks-per-node=16 # numero de tasks por nodo

sizes=(16777216 134217728)
threads=( 1 2 4 8 16)

for i in "${threads[@]}"
do
	for j in "${sizes[@]}"
	do
   		./ct.out $j $i
	done
done

