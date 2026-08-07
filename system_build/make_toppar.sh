#!/bin/bash

mkdir gen_itps
cp toppar/lipid-generator.py ./gen_itps
cd gen_itps

# insert from itp generation sheet here
python3 lipid-generator.py -alname '' -alhead '' -allink '' -altail '' -o .itp
python3 lipid-generator.py -alname '' -alhead '' -allink '' -altail '' -o .itp
python3 lipid-generator.py -alname 'DPPC' -alhead 'C P' -allink 'G G' -altail 'cCCC cCCC' -o DPPC.itp
python3 lipid-generator.py -alname 'POPC' -alhead 'C P' -allink 'G G' -altail 'cCCC CDCC' -o POPC.itp
python3 lipid-generator.py -alname 'PLPC' -alhead 'C P' -allink 'G G' -altail 'cCCC CDDC' -o PLPC.itp
python3 lipid-generator.py -alname 'PAPC' -alhead 'C P' -allink 'G G' -altail 'cCCC cFFDC' -o PAPC.itp
python3 lipid-generator.py -alname 'DOPC' -alhead 'C P' -allink 'G G' -altail 'CDCC CDCC' -o DOPC.itp
python3 lipid-generator.py -alname 'OAPC' -alhead 'C P' -allink 'G G' -altail 'CDCC cFFDC' -o OAPC.itp

python3 lipid-generator.py -alname 'DP1C' -alhead 'C P' -allink 'E G' -altail 'cCCC cCCC' -o DP1C.itp
python3 lipid-generator.py -alname 'PO1C' -alhead 'C P' -allink 'E G' -altail 'cCCC CDCC' -o PO1C.itp
python3 lipid-generator.py -alname 'PA1C' -alhead 'C P' -allink 'E G' -altail 'CCCC cFFDC' -o PA1C.itp

python3 lipid-generator.py -alname 'PSLC' -alhead 'C P' -allink 'L G' -altail 'cCCC CCCC' -o PSLC.itp
python3 lipid-generator.py -alname 'SALC' -alhead 'C P' -allink 'L G' -altail 'CCCC cFFDC' -o SALC.itp

python3 lipid-generator.py -alname 'SOPE' -alhead 'E P' -allink 'G G' -altail 'CCCC CDCC' -o SOPE.itp
python3 lipid-generator.py -alname 'SLPE' -alhead 'E P' -allink 'G G' -altail 'CCCC CDDC' -o SLPE.itp
python3 lipid-generator.py -alname 'SAPE' -alhead 'E P' -allink 'G G' -altail 'CCCC cFFDC' -o SAPE.itp
python3 lipid-generator.py -alname 'OAPE' -alhead 'E P' -allink 'G G' -altail 'CDCC cFFDC' -o OAPE.itp

python3 lipid-generator.py -alname 'POLE' -alhead 'E P' -allink 'L G' -altail 'cCCC CDCC' -o POLE.itp
python3 lipid-generator.py -alname 'SALE' -alhead 'E P' -allink 'L G' -altail 'CCCC cFFDC' -o SALE.itp
python3 lipid-generator.py -alname 'OALE' -alhead 'E P' -allink 'L G' -altail 'CDCC cFFDC' -o OALE.itp

python3 lipid-generator.py -alname 'SOPI' -alhead 'PI' -allink 'G G' -altail 'CCCC CDCC' -o SOPI.itp
python3 lipid-generator.py -alname 'SAPI' -alhead 'PI' -allink 'G G' -altail 'CCCC cFFDC' -o SAPI.itp
python3 lipid-generator.py -alname 'OAPI' -alhead 'PI' -allink 'G G' -altail 'CDCC cFFDC' -o OAPI.itp
python3 lipid-generator.py -alname 'SOPS' -alhead 'S P' -allink 'G G' -altail 'CCCC CDCC' -o SOPS.itp
python3 lipid-generator.py -alname 'SAPS' -alhead 'S P' -allink 'G G' -altail 'CCCC cFFDC' -o SAPS.itp
python3 lipid-generator.py -alname 'OAPS' -alhead 'S P' -allink 'G G' -altail 'CDCC cFFDC' -o OAPS.itp
python3 lipid-generator.py -alname 'PSM' -alhead 'C P' -allink 'A1 A2' -altail 'tCCC cCCC' -o PSM.itp
python3 lipid-generator.py -alname 'XSM' -alhead 'C P' -allink 'A1 A2' -altail 'tCCC CCCCCC' -o XSM.itp
python3 lipid-generator.py -alname 'NSM' -alhead 'C P' -allink 'A1 A2' -altail 'tCCC CCCDCC' -o NSM.itp
python3 lipid-generator.py -alname 'TLCR' -alhead 'P' -allink 'A1 A2' -altail 'tCCC cC' -o TLCR.itp
python3 lipid-generator.py -alname 'PLCR' -alhead 'P' -allink 'A1 A2' -altail 'tCCC cCCC' -o PLCR.itp
python3 lipid-generator.py -alname 'NLCR' -alhead 'P' -allink 'A1 A2' -altail 'tCCC CCCDCC' -o NLCR.itp
python3 lipid-generator.py -alname 'ZLCR' -alhead 'P' -allink 'A1 A2' -altail 'tCCC CCCCCCC' -o ZLCR.itp

cat ./*.itp > lipids.itp

cd ../
cp ./gen_itps/lipids.itp ./toppar
cd toppar

