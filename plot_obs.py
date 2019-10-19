import matplotlib.pyplot as plt


infile = open("output_hot.txt", "r")

x = []
for line in infile:
    line.strip('[]\n')
    x.append(line.split())

###################################
# data = infile.read().split("]") #
# for x in data:                  #
#     x.strip("\n")               #
###################################

