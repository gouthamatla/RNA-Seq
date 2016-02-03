'''
Usage: script.py count_matrix.txt min_count min_samples > filtered_out.txt

count matrix:
Geneid    322-11A    322-11B    322-11C    322-12A    322-12B    322-12C    322-1A    322-1B    322-1C    322-2A    322-2B    322-2C    322-3A    322-3B    322-3C
ENSG00000000457    401    423    593    490    474    581    459    404    436    435    512    525    526    444    575
ENSG00000000460    237    258    257    272    248    310    249    263    202    285    291    277    294    265    321
...
...
'''

import sys

infile=sys.argv[1]
min_count=sys.argv[2]
min_samples=sys.argv[3]

with open(infile) as counts:
    print counts.readline().strip()
    for line in counts:
        val = line.split()
        if len([x for x in val[1:] if int(x) > min_count ]) >= min_samples :
            print line.strip()
