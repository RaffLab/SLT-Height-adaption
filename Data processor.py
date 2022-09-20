import csv
from statistics import mean



FRAMES = 17
DATAPATH = "/Users/thomas/Documents/PhD/Data/Project Laura/Ana2/Highspeed/Center vs outside test/S3rNB_meanRectIntensity unlinked for frame averging.csv"
EXPORTPATH = "S3rNB"

dataArrs = []
for i in range(FRAMES):
    dataArrs.append([])

#read it and split it
csvfile = open(DATAPATH, newline='')
reader = csv.reader(csvfile, delimiter=',')
i = 0
for row in reader:
    formattedRow = []
    for entry in row:
        if entry: formattedRow.append(float(entry))
    dataArrs[i].append(formattedRow)
    i += 1
    if i > 16: i = 0
csvfile.close()

for i in range(len(dataArrs)):
    file = open("/Users/thomas/Documents/PhD/Data/Project Laura/Ana2/Highspeed/Center vs outside test/"+EXPORTPATH+"_Z"+str(i)+".csv",mode="x", newline="")
    writer = csv.writer(file, dialect="excel")
    for j in range(len(dataArrs[i])):
        line = [str(j)]

        
        line.append(mean(dataArrs[i][j]))
        writer.writerow(line)
    
    file.close()




