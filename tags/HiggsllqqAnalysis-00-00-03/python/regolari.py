# author: Valerio Ippolito <valerio.ippolito@cern.ch>
# this script converts a comma separated list of datasets into an one-dataset-per-row file

# import python regular expressions
import re
# define p, a regular expression
p = re.compile(', *') # a comma followed by >=0 spaces

# open input.txt
f = open("input.txt", "r")
# read its content
dataset_list = f.read()
# close it
f.close()

# open it back as an output file
g = open("input.txt", "w")

# substitute, in the previous file content, each comma followed by >=0 spaces with a newline char
new_list = p.sub('\n', dataset_list)

# write the file
g.write(new_list)

# and close it
g.close()
