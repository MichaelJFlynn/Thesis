from rdatkit import datahandlers
import sys

## Usage:
## python getMutatedSeqs.py <input .rdat> <output filename>
rdat = datahandlers.RDATFile()
rdat.load(open(sys.argv[1]))
sections = rdat.constructs.values()[0]

f = open(sys.argv[2], 'w')
for datum in sections.data:
    f.write(datum.annotations['mutation'][0]+ '\t' + '\t'.join(map(str, datum.values)) + '\n') 
f.close()
