curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python get-pip.py
pip install pyfaidx
faidx --transform bed vcfout.fasta > vcfout.bed
