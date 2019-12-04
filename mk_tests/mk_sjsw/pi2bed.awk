# converts pi output from VCFtools to BED format
#
# usage:
#
# $ awk -f pi2bed.awk input.pi > output.bed
#

BEGIN {
    # splitting character is the tab
    FS = OFS = "\t"
}

# skip  header
NR <= 1 {
    next;
}

# transform the output
{
    print $1, $2-1, $2, $3 
}
