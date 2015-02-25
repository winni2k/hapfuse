###########
# Tue Feb 24 12:20:40 GMT 2015
# creating test40 input files
cd test40
for i in chr20_20033136_20177280.STv1.2.13.C100.K100.first2Samp.bin.vcf chr20_20149699_20313076.STv1.2.13.C100.K100.first2Samp.bin.vcf chr20_20281142_20442047.STv1.2.13.C100.K100.first2Samp.bin.vcf; do
    bcftools convert --hapsample test40.$i ../test30/test30.$i
    gunzip test40.$i.hap.gz
done

### using ligateHAPLOTYPES to create truth set
# need to create pseudo "GL" vcf from haps
# create concatenated haps file
perl -ane '@F[5..8] = (0,0,0,0); print join(" ", @F)."\n"' test40.chr20_20*.hap| sort -k 3 -g |uniq > test40.pseudo-concat.hap

bcftools convert --hapsample2vcf <(perl -ane '$F[0] = $F[1]; print join(" ", @F)."\n"' test40.pseudo-concat.hap),test40.chr20_20033136_20177280.STv1.2.13.C100.K100.first2Samp.bin.vcf.sample -Ov -o test40.pseudo-concat.vcf

ligateHAPLOTYPES --vcf test40.pseudo-concat.vcf --chunks test40.chr20_20*.hap --output test40.ligated.hap test40.ligated.sample

# convert to impute 2 hap
~/feng/marchini/scripts/haps2hapLegend.pl  -i test40.ligated.hap -o test40.ligated.i2

# flip second pair of haplotypes
perl -ane '@F = @F[1,0,3,2]; print join(" ", @F)."\n"' test40.ligated.i2.hap >test40.ligated.i2.flip.hap

