#!/bin/bash
#SBATCH -J name change
#SBATCH -t 4:00:00
#SBATCH -c 8

sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/dbEST_CLYHEM.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HISEQ-129-C25BFACXX-8-HYDRACTINIA.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HISEQ-134-C27MTACXX-2-PODOCORYNA.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HISEQ-168-C3DEYACXX-8-ALATINA.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST420-69-D0F2NACXX-3-ECTOPLEURA.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-159-C4MVCACXX-5-AGTTCC.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-159-C4MVCACXX-5-CCGTCC.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-159-C4MVCACXX-5-GTCCGC.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-159-C4MVCACXX-5-TGACCA.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-181-C76K5ACXX-2-ACTTGA.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-181-C76K5ACXX-2-ATCACG.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-181-C76K5ACXX-2-CTTGTA.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-181-C76K5ACXX-2-GGCTAC.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-181-C76K5ACXX-2-TAGCTT.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-51-C02UNACXX-6-FRILLAGALMA.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-51-C02UNACXX-7-NANOMIA.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-51-C02UNACXX-8-BARGMANNIA.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-3-TATC.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-6-ACAGTG.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-6-ATCACG.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-6-CAGATC.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-6-GCCAAT.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-6-TTAGGC.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-7-CGATGT.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-7-CTTGTA.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-7-GATCAG.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-7-GGCTAC.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-7-TAGCTT.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-7-TGACCA.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-8-ACAGTG.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-8-ACTTGA.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-8-ATCACG.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-8-CAGATC.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-8-GCCAAT.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-54-C026EACXX-8-TTAGGC.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-61-C02G3ACxx-6-GGCTAC.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-73-C0JUVACXX-7-AGALMA2.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-73-C0JUVACXX-7-ATCACG.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/HWI-ST625-73-C0JUVACXX-7-TTAGGC.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/JGI_NEMVEC.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/K00162-189-HJTYGBBXX-7-NCAGTG.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/NCBI_HYDMAG.fa
sed -E -i 's/([A-Z]*) ([A-Z]*)/\1_\2/g' /home/cm2477/scratch60/siphonophore2016/OrthoFinder/sequences/SRX231866.fa
