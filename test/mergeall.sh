#!/bin/sh

declare -a arr=(
    toys50_smear0_btag0_gen1_test1
    toys50_smear1_btag0_gen1_test1
    toys50_smear1_btag1_gen1_test1
    toys50_smear1_btag2_gen1_test1
    toys50_smear0_btag0_gen6_test1
    toys50_smear1_btag0_gen6_test1
    toys50_smear1_btag1_gen6_test1
    toys50_smear1_btag2_gen6_test1
    toys50_smear0_btag0_gen2_test2
    toys50_smear1_btag0_gen2_test2
    toys50_smear1_btag1_gen2_test2
    toys50_smear1_btag2_gen2_test2
    toys50_smear0_btag0_gen7_test2
    toys50_smear1_btag0_gen7_test2
    toys50_smear1_btag1_gen7_test2
    toys50_smear1_btag2_gen7_test2
    toys50_smear0_btag0_gen3_test3
    toys50_smear1_btag0_gen3_test3
    toys50_smear1_btag1_gen3_test3
    toys50_smear1_btag2_gen3_test3
    toys50_smear0_btag0_gen8_test3
    toys50_smear1_btag0_gen8_test3
    toys50_smear1_btag1_gen8_test3
    toys50_smear1_btag2_gen8_test3
)

for i in ${arr[@]}
do
  ./merge.sh ${i}
done