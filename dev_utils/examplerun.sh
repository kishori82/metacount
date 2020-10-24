#./rpkm -c ~/SAG_project/RPKM_project/data/HYXS/test_data/test_4093112_combined_unique.fasta --r1 ~/SAG_project/RPKM_project/data/HYXS/test_data/test_4093112_artificial_reads.sam -O ~/SAG_project/RPKM_project/data/HYXS/test_data/test_4093113_combined_unique.unannot.gff -f sam-1 --multireads -o ~/SAG_project/RPKM_project/RPKM_out/python/
./metacount  --gff data/lagoon-sample2.unannot.gff \
  --sam /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample2/bwa/lagoon-sample2_batchr1.sam \
  --sam /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample2/bwa/lagoon-sample2.sam \
  --sam /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample2/bwa/lagoon-sample2_a.sam \
  --estimate-type RPKM --out-file x --print-stats --stats-out-file y


#  --sam /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample2/bwa/lagoon-sample2_batchr1.sam \
#  --sam /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample2/bwa/lagoon-sample2_batchr2.sam \
#  --sam /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample2/bwa/lagoon-sample2_batch3.sam \
#  --sam /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample2/bwa/lagoon-sample2_2.sam \
#  --sam /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample2/bwa/lagoon-sample2_4.sam  \
#  --sam /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample2/bwa/lagoon-sample2_batch2.sam \
#  --sam /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample2/bwa/lagoon-sample2.sam \
#  --sam /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample2/bwa/lagoon-sample2_a.sam \
#  --sam /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample2/bwa/lagoon-sample2_1.sam
