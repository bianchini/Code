#!/bin/sh

declare -a arr=(
    #toys100_smear1_btag0_gen0_test0
    #toys100_smear1_btag1_gen0_test0
    #toys100_smear1_btag2_gen0_test0
    #toys100_smear1_btag0_gen5_test0 
    #toys100_smear1_btag1_gen5_test0
    #toys100_smear1_btag2_gen5_test0
    #toys50_smear0_btag0_gen1_test1
    #toys50_smear1_btag0_gen1_test1
    #toys50_smear1_btag1_gen1_test1
    #toys50_smear1_btag2_gen1_test1
    #toys50_smear0_btag0_gen6_test1
    #toys50_smear1_btag0_gen6_test1
    #toys50_smear1_btag1_gen6_test1
    #toys50_smear1_btag2_gen6_test1
    #toys50_smear0_btag0_gen2_test2
    #toys50_smear1_btag0_gen2_test2
    #toys50_smear1_btag1_gen2_test2
    #toys50_smear1_btag2_gen2_test2
    #toys50_smear0_btag0_gen7_test2
    #toys50_smear1_btag0_gen7_test2
    #toys50_smear1_btag1_gen7_test2
    #toys50_smear1_btag2_gen7_test2
    #toys50_smear0_btag0_gen3_test3
    #toys50_smear1_btag0_gen3_test3
    #toys50_smear1_btag1_gen3_test3
    #toys50_smear1_btag2_gen3_test3
    #toys50_smear0_btag0_gen8_test3
    #toys50_smear1_btag0_gen8_test3
    #toys50_smear1_btag1_gen8_test3
    #toys50_smear1_btag2_gen8_test3
    #toys25_smear1_btag2_gen4_test4
    #toys25_smear1_btag2_gen9_test4
    #toys25_smear1_btag1_gen10_test5
    #toys25_smear1_btag1_gen11_test5
    #tH_tL_Hh_smear1_btag1_gen2_test2
    #tH_tL_Hh_smear1_btag1_gen7_test2
    #tH_tL_Hh_smear1_btag1_gen12_test2
    tH_tL_Hh-Ha_smear1_btag1_gen11_test6
    tH_tL_Hh-Hb1_smear1_btag1_gen12_test6
    tH_tL_Hh-Hb2_smear1_btag1_gen7_test6
    tH_tL_Hh-Hs_smear1_btag1_gen2_test6
)

for i in ${arr[@]}
do
  ./merge.sh ${i}
done