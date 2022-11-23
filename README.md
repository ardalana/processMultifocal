# processMultifocal
Procedure to list denovo mutations for each tumor in a group of clustered tumors

This script takes in intersected mutation lists created by using a number of variant callers for each and every tumor in one group of clustered tumors (i.e. from one patient). The script then removes positions mutating nonspecifically across different tumors of the group, finds denovo (i.e. non-dbsnp-listed) mutations in the group (dbsnp file supplied), finds non-repeat-region denovo mutations in the group (repeatmasker file supplied), and prints them for each tumor in inclusive and exclusive lists, also separately for indels and snps.
