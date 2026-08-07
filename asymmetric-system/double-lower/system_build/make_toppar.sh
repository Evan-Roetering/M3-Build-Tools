#!/bin/bash
# Script to generate Martini 3 lipid itps using the lipid-generator.py script
# Last updated by:
#   Evan Roetering
#   Klauda Lab
#   University of Maryland
#   8/7/2026


mkdir gen_itps
cp toppar/lipid-generator.py ./gen_itps
cd gen_itps

# insert from itp generation sheet here

# Example format:
# python3 lipid-generator.py -alname '' -alhead '' -allink '' -altail '' -o .itp
#   alname - Use any 3-4 letter code, must be unique for each lipid
#   alhead - headgroup types options are as follows:
#       P - phosphate [PO4] {Q5 Bead}
#       O - phosphatidyl [PO4] {Q5 Bead w/ -2 charge} (use P for phosphatidyl headgroup)
#       C - choline [NC3] {Q1 Bead} (use PC for phosphatidylcholine)
#       E - ethanolamine [NH3] {Q4p bead} (use PE for phosphatidylethanolamine)
#       G - glycerol [GL0] {P4r bead} (use PG for phosphatidylglycerol)
#       S - serine [CNO] {P6 Bead} (use PS for phosphatidylserine)
#       COH - Cap for diacylglycerol and ceramide [COH] {TN6 bead}
#       PI - Phosphatidylinositol (There is no inositol only headgroup option)
#       P1 - PIP1(3) (Phosphatidylinositol 3-phosphate -3 charge)
#       P2 - PIP2(3,4) (Phosphatidylinositol 3,4-bisphosphate -4 charge)
#       P3 - PIP3(3,4,5) (Phosphatidylinositol 3,4,5-trisphosphate -5 charge)
#       P4 - PIP1(4) (Phosphatidylinositol 4-phosphate -3 charge)
#       P5 - PIP1(5) (Phosphatidylinositol 5-phosphate -3 charge)
#       P6 - PIP2(4,5) (Phosphatidylinositol 4,5-bisphosphate -4 charge)
#       P7 - PIP2(3,5) (Phosphatidylinositol 3,5-bisphosphate -4 charge)
#   alink - Use in format "linker1 linker2" 
#           linker types options are as follows:
#       G - glycerol [GL1/GL2] {SN4a bead}
#       L - Plasmalogen [PL1/PL2] {SN4as bead}
#       E - Ether [GL1/GL2] {SN3a bead}
#       Y - Lysolipid glycerol termination (use for SN2 chain) [LY2] {TN6 bead}
#       A - Sphingosine [OH1/AM2] {SP1/SP2 bead} (Can not be mixed with G/L/E/Y)
#   altail - Use in format "tail1 tail2"    
#            tail types options are as follows:
#       c - short saturated chain segment (-CH2-CH2-) [SC1] 
#       C - long saturated chain segment (-CH2-CH2-CH2-CH2-) [C1]
#       D - long unsaturated chain segment (-CH2-CH=CH-CH2-) [D1]
#       F - long polyunsaturated chain segment (-CH2-CH=CH-CH=) [D1]
#       T - trans unsaturated chain segment (-CH2-CH=CH-CH2-) [D1] (use D for cis unsaturated chain)
#       t - trans short unsaturated chain segment (use for sphingosine) (-CH2-CH=CH-CH2-) [SC1] (use c for cis unsaturated chain)
#       - - nothing (use for lysolipid)
# example usages:
#   DPPC - python3 lipid-generator.py -alname 'DPPC' -alhead 'C P' -allink 'G G' -altail 'cCCC cCCC' -o DPPC.itp
#   POPC - python3 lipid-generator.py -alname 'POPC' -alhead 'C P' -allink 'G G' -altail 'cCCC CDCC' -o POPC.itp
#   PLPC - python3 lipid-generator.py -alname 'PLPC' -alhead 'C P' -allink 'G G' -altail 'cCCC CDDC' -o PLPC.itp
#   PAPC - python3 lipid-generator.py -alname 'PAPC' -alhead 'C P' -allink 'G G' -altail 'cCCC cFFDC' -o PAPC.itp
#   SN1 ether POPC - python3 lipid-generator.py -alname 'PO1C' -alhead 'C P' -allink 'E G' -altail 'cCCC CDCC' -o PO1C.itp
#   SN1 Plasmalogen POPC - python3 lipid-generator.py -alname 'POLC' -alhead 'C P' -allink 'L G' -altail 'cCCC CDCC' -o POLC.itp
#   POPE - python3 lipid-generator.py -alname 'SOPE' -alhead 'E P' -allink 'G G' -altail 'cCCC CDCC' -o POPE.itp
#   POPI - python3 lipid-generator.py -alname 'POPI' -alhead 'PI' -allink 'G G' -altail 'cCCC CDCC' -o POPI.itp
#   POPS - python3 lipid-generator.py -alname 'POPS' -alhead 'S P' -allink 'G G' -altail 'cCCC CDCC' -o POPS.itp
#   PSM - python3 lipid-generator.py -alname 'PSM' -alhead 'C P' -allink 'A1 A2' -altail 'tCCC cCCC' -o PSM.itp
#   PCER - python3 lipid-generator.py -alname 'PCER' -alhead 'COH' -allink 'A1 A2' -altail 'tCCC cCCC' -o PCER.itp

cat ./*.itp > lipids.itp

cd ../
cp ./gen_itps/lipids.itp ./toppar
cd toppar

