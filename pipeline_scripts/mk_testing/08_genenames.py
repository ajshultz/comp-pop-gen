## A line reader to extract gene abrevitions from bed files with the contents imbeded as ;gene=CONTENT;

#!/usr/bin/env python3

import fileinput
import sys

with fileinput.input('corCor.consensus.bed') as intake:
    for line in intake:
        parts = line.split(';')
        abrvs = [part.rsplit('=')[1] for part in parts if 'gene=' in part]
        abrvs = set(abrvs)
        print(",".join(abrvs), end='\n', file=sys.stdout)
        
## for easiest output, run on command line: 

# python3 08_genenames.py > corCor.genenames.txt
