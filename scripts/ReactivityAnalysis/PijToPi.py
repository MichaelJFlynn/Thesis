import sys

if len(sys.argv) != 2:
    print "Usage: python PijToPi.py <pairs file>"
    sys.exit()

f = open(sys.argv[1], 'r')

# Get sequence
sequence = next(f).strip()
# Valid .pairs files should follow this formatting
rows = [[float(x) for x in line.split("\t")] for line in f]

# bins to hold probabilities
output = [0] * len(sequence)
for line in rows:
    output[int(line[0])-1] += line[2]
    output[int(line[1])-1] += line[2]

for i in range(len(output)):
    print str(i+1) + "\t" + str(output[i])

    
