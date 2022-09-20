#height tool, works for bright embryos
#requires 8bit greyscale tiffs
#requires input of N of slices (cos I'm lazy and cba to go metadata hunting again)


### DEPENDENCIES ###
import os 
import sys
from typing import OrderedDict
import cv2 as cv
import pandas as pd
import tifffile as tf
import numpy as np
import math
from PIL import Image
import csv


### Globals ew ###
G_CURRFRAME = 0
G_CURRZ = 0

### Classy ones ###


#object to represent the individual dots 
class splodgeOfInterest:
    
    #calculate the intensity across the dimensioning rectangle
    def calculateMeanCntIntensity(self,image):
        mask = np.zeros(image.shape,np.uint8)
        cv.drawContours(mask, self.contour,0,255,-1)
        return cv.mean(image, mask)[0]

    def calculateMeanRectIntesity(self,image):
        mask = np.zeros(image.shape, np.uint8)
        box = cv.boxPoints(self.rect)
        box = np.int0(box)
        cv.drawContours(mask, [box],0,255,-1)
        return cv.mean(image, mask)[0]   
    
    #calculate the centre of the contour; if the contour is open we probs don't want it but we're just gonna calculate 
    #the geometric centre ie. ave x ave y instead
    def calculateSplodgeCentroid(self):
        MomentMat = cv.moments(self.contour)
        #for domaining

        if MomentMat['m00'] == 0: #not enclosed contour
            centre = np.mean(self.contour, axis=0)
            return centre[0][0], centre[0][1]
            
        cx = MomentMat['m10']/MomentMat['m00']
        cy = MomentMat['m01']/MomentMat['m00']
        return cx,cy


    def __init__(self,ID,cnt,frame,zPos,image):
    #bit inefficent storing frame here as well but it'll jut help keep data together
        self.splodgeID = ID #id - we do not touch this cos this binds the hash
        self.contour = cnt #contour data
        self.frame = frame #id of frame this splode exists in
        self.excluded = False #bool which can be flicked on and off to hide
        self.rect = cv.minAreaRect(cnt) #rectangle
        self.zPos = zPos
        if cnt.shape[0] >= 6:
            (x,y),(self.elpsWidth,self.elpsLength),self.angle = cv.fitEllipse(cnt) #elipse 
        else:
            self.elpsLength = "N/A"
            self.elpsWidth = "N/A"
            self.angle = "N/A"
        self.rectArea = self.rect[1][0]*self.rect[1][1]
        self.maxDim = max(self.rect[1])
        self.minDim = min(self.rect[1])
        if self.minDim == 0:
            self.ratio = "n/a"
        else:
            self.ratio = self.maxDim/self.minDim


        if self.elpsWidth == 0 or self.elpsLength == "N/A":
            self.elpsRatio = "n/a"
        else:
            self.elpsRatio = self.elpsLength/self.elpsWidth
        #self.meanCntIntensity = self.calculateMeanCntIntensity(image)
        #currently unusued to save time, we should add a paramenter to control how much
        #is calculated
        self.meanRectIntensity = self.calculateMeanRectIntesity(image)
        self.rectX = self.rect[0][0]
        self.rectY = self.rect[0][1]
        self.cntX, self.cntY = self.calculateSplodgeCentroid()
        

    #this is a bit of code just to make this object hashable so we can use it as a key 
    #to a dict, which I may or may not use.... 
    def __hash__(self):
        return hash(self.splodgeID)
    
    def __eq__(self, other):
        return (self.__class__, self.splodgeID) == (other.__class__, other.splodgeID)

    #function to evaluate if the splodge warrents exclusion based on a filter dict
    #filterdict structured as a dictionary keyed by parameter names containing subdicts 
    #keyed by min, max, and equal doing pretty much what they say on the tin 
    def evaluateExclusion(self, filterDict):
        self.excluded = False #base state

        #loop through the parameters
        for paramKey, subFilterDict in filterDict.items():
        #first we can shortcut if there is an equal 
            testAttribute = getattr(self, paramKey)
            if '=' in subFilterDict.keys():
                if testAttribute != subFilterDict['=']:
                    #its not equal and excluded
                    self.excluded = True
                    return True #just incase we ever want to assign with it, plus short cuts the function
                else: continue #can just continue with the loop as equals is most prescriptive
            if '>' in subFilterDict.keys(): #more likely to be greater than so lets get this shortcut out the way faster
                if (testAttribute > subFilterDict['>']) == False:
                    self.excluded = True
                    return True
            if '<' in subFilterDict.keys():
                if (testAttribute < subFilterDict['<']) == False:
                    self.excluded = True
                    return True
        
        #so if we've got here then its defo false but to be safe, I'm just gonna return the excluded
        return self.excluded #note: false is good :D

    #function to handle ids, use to do some funky stuff, guess not needed now
    #should make something to get as an arr though maybe?
    def setTrackID(self,id):
        self.trackID = id

    #function to evaluate the pseudogrid position of a particular splodge
    #returns tuple of rowIndex columnIndex 
    def getPseudoGridPos(self, pseudoGridMetaData):

        #calculate column
        #normalise wrt to the image bounds
        normalisedX = (self.cntX - pseudoGridMetaData['minX'])/pseudoGridMetaData['xDiff']
        xPos = math.floor(normalisedX * pseudoGridMetaData['columns']) #scale to columns then we want to round down
        
        #calculate row
        normalisedY = (self.cntY - pseudoGridMetaData['minY'])/pseudoGridMetaData['yDiff']
        yPos = math.floor(normalisedY * pseudoGridMetaData['rows']) #scale to rows and again round donw

        return yPos,xPos

    #function to set the current parent info to link together splodges
    def setParentData(self, currParentCandidate):
        self.parentCandidate = currParentCandidate
        #{parentSplodge, euDist, tNormEuDist, currSplodge:self, linkDisp}



#obj to contain the individual SoIs 
#contains handling code to access and create them more easily
#note: this an SoI have been ripped and adapted from the SORA tool so may contain some unneccesary stuff
class splodgeContainer:
    def __init__(self): #really basic initating constructor
        self.splodges = {}
        #dict {frame:[splodges]}
        #data then the splodge objects
        self.currSplodgeID = 0
        self.splodgeStats = {
            'minRectIntensity':300,#images are 8bit so this is safe
            'maxRectIntensity':0,
            'minCntIntensity':300,
            'maxCntIntensity':0,
            'areas':[],
            'mincX':2147483647, #intlimit
            'maxcX':0,
            'mincY':2147483647,
            'maxcY':0,
            'minLength':2147483647,
            'maxLength':0,
            'minWidth':2147483647,
            'maxWidth':0,
            'minZ':2147483647,
            'maxZ':-1
        }
        self.filterDict = {}
        #dictionary standardising filters
        #labled by splodge parameters
        #subdicts then indicate min,max and equal parameters

        #adds splodges from a list of contours
    #checks to make sure contours are concave to check we're looking at a
    #circular blob, lets check this is consistent though
    def appendSplodges(self,cntList,frame,zPos, image):
        currData = []
        for cnt in cntList:
            if cv.isContourConvex(cnt): #so something looks weird
                #check this holds 
                 continue
            #so it looks good lets add it to the splodge list
            splodge = splodgeOfInterest(self.currSplodgeID, cnt,frame,zPos,image)
            currData.append(splodge)
            
            #append stats: look at automating this
            if self.splodgeStats['minRectIntensity'] > splodge.meanRectIntensity:
                self.splodgeStats['minRectIntensity'] = splodge.meanRectIntensity
            elif self.splodgeStats['maxRectIntensity'] < splodge.meanRectIntensity:
                self.splodgeStats['maxRectIntensity'] = splodge.meanRectIntensity

            if self.splodgeStats['mincX'] > splodge.cntX:
                self.splodgeStats['mincX'] = splodge.cntX
            elif self.splodgeStats['maxcX'] < splodge.cntX:
                self.splodgeStats['maxcX'] = splodge.cntX
            
            if self.splodgeStats['mincY'] > splodge.cntY:
                self.splodgeStats['mincY'] = splodge.cntY
            elif self.splodgeStats['maxcY'] < splodge.cntY:
                self.splodgeStats['maxcY'] = splodge.cntY

            if self.splodgeStats['minLength'] > splodge.maxDim:
                self.splodgeStats['minLength'] = splodge.maxDim
            elif self.splodgeStats['maxLength'] < splodge.maxDim:
                self.splodgeStats['maxLength'] = splodge.maxDim

            if self.splodgeStats['minWidth'] > splodge.minDim:
                self.splodgeStats['mincWidth'] = splodge.minDim
            elif self.splodgeStats['maxcWidth'] < splodge.minDim:
                self.splodgeStats['maxcWidth'] = splodge.minDim

            if self.splodgeStats['minZ'] > splodge.zPos:
                self.splodgeStats['minZ'] = splodge.zPos
            elif self.splodgeStats['maxZ'] < splodge.zPos:
                self.splodgeStats['maxZ'] = splodge.zPos    
            
            
            #if self.splodgeStats['minCntIntensity'] > splodge.meanCntIntensity:
            #    self.splodgeStats['minCntIntensity'] = splodge.meanCntIntensity
            #elif self.splodgeStats['maxCntIntensity'] < splodge.meanCntIntensity:
            #    self.splodgeStats['maxCntIntensity'] = splodge.meanCntIntensity

            self.splodgeStats['areas'].append(splodge.rectArea)

            #if self.splodgeStats['minRectArea'] > splodge.rectArea:
            #    self.splodgeStats['minRectArea'] = splodge.rectArea
            #elif self.splodgeStats['maxRectArea'] < splodge.rectArea:
            #    self.splodgeStats['maxRectArea'] = splodge.rectArea

            self.currSplodgeID += 1 #increment this to maintain unique IDs
        #so at this point we've got them all nicely packaged into objects in an array; this handles multiple additions to the same section 
        #wtithout fucking up the preivous code assuming one layer of organisation
        if frame not in self.splodges.keys():
            self.splodges[frame] = []
        self.splodges[frame].append(currData)

    #function to output data as csv
    def outputRawData(self, path):
        outputDict = {
            "Splodge_ID":[],
            "Frame_ID":[],
            "Length":[],
            "Width":[]
        }
        for frameDict in self.splodges.values():
            for splodge in frameDict:
                if splodge.excluded: #ignores if data told to be ignored
                    continue
                outputDict['Splodge_ID'].append(splodge.splodgeID)
                outputDict['Frame_ID'].append(splodge.frame)
                outputDict['Length'].append(splodge.maxDim)
                outputDict['Width'].append(splodge.minDim)

        #this will construct the dict
        outputDF = pd.DataFrame(outputDict)
        #this then exports the panda as csv
        outputDF.to_csv(path, index=False)

    #function to output a table of single data values
    #this thus is transposed and aligned with frame indexs as rows and repeats as columns
    def outputDataValue(self, dataVal, path):
        import csv 
        with open(path, 'w',newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for frameNum, frameDict in self.splodges.items():
                rowArr = [frameNum]
                for splodge in frameDict:
                    if splodge.excluded: #ignores if data told to be ignored
                        continue
                    rowArr.append(getattr(splodge, dataVal))
                writer.writerow(rowArr)

    #as a dict this is really easy lol, why didn't I do this first /facepalm
    def getFrameArr(self,frame):
        return(self.splodges[frame])


    #function to generate an array of the rect contours to be plotted onto images
    def getRectContours(self,frame,zRef):
        rectCnts = []
        frameArr = self.getFrameArr(frame)[zRef]
        for splodge in frameArr:
            #convert the rect param into contours to plot
            if splodge.excluded: #ignores if data ignored
                continue
            box = cv.boxPoints(splodge.rect)
            rectCnts.append(np.int0(box))
        return rectCnts

    #function to merge/add new filters onto the filter list
    def mergeFilterDict(self,filterChanges):
        for paramKey, newFilterDict in filterChanges.items(): #loop through the filter changes
            if (paramKey in self.filterDict.keys()) == False: #ie the key isn't currently there
                self.filterDict[paramKey] = {} #initialise it
            for filterType, filterVal in newFilterDict.items(): #loop through the filter parameter stuff
                if filterVal == -1: #remove filter
                    self.filterDict[paramKey].pop(filterType)
                else: #update the filter
                    self.filterDict[paramKey][filterType] = filterVal
            if bool(self.filterDict[paramKey]) == False: #ie the dict is empty
                self.filterDict.pop([paramKey]) #remove the filter from the array to speed up further evals


    #generalised function which goes through the data and selects data based on filters
    #theoretically filterchanges can be optional?
    def excludeData(self, filterChanges):
        #first we need to merge filters 
        self.mergeFilterDict(filterChanges)
        #now we need to go through and reevaluate with regards to the new features
        #let reevaluate everything otherwise its gonna be complex for when we're removing features
        #this is arguably a bit overly complex as we're going through everything again but otherwise
        #its difficult to know whats changed and and what was dependent on what
        #hypothetically each splodge could have an identifier to show which key lost it but I feel going through
        #those comparisons on each splodge will take a similar amount of time at the cost of more memory
        for frameArr in self.splodges.values():
            for zArr in frameArr:
                for splodge in zArr:
                    splodge.evaluateExclusion(self.filterDict)


#class to define the spots refering to the same centriole 
class centrioleLinkage():
    def __init__(self,initialSplodge):
        self.splodges = {} #dict of splodges composing the linkage by zStack
        self.splodges[initialSplodge.zPos] = initialSplodge
        self.recentlyAddedTo = True
        self.prevZ = initialSplodge.zPos
        self.active = True
        self.excluded = False

    def appendSplodge(self, splodge):
        self.splodges[splodge.zPos] = splodge
        self.recentlyAddedTo = True
        self.prevZ = splodge.zPos

    #function to check if the centriole has been recently added to or if its finished, this shortcuts it preventing excessive comparisons or 
    #linkages across large Z-distances. Assumes a centriole will be reliably tracked vertically across its length w/o gaps
    def checkActive(self):
        if self.active == False: return False #inactive so we don't care anymore
        if self.recentlyAddedTo == False: #becomes inactive
            self.active = False
            return False
        else:
            self.recentlyAddedTo = False #flip this for the next time it checks
            return True #return I'm interesting

    #function to get the pseudoGridpositon of a specific splodge by Zpos
    def getPseudoGridPosition(self, targetZPos, pseudoGridMetaData):
        return self.splodges[targetZPos].getPseudoGridPos(pseudoGridMetaData)

    def getPrevCentroidPos(self):
        return (self.splodges[self.prevZ].cntX, self.splodges[self.prevZ].cntY)

    def getEuDistToPrevCentroid(self,comparisonSplodge):
        xDist = self.splodges[self.prevZ].cntX-comparisonSplodge.cntX
        yDist = self.splodges[self.prevZ].cntY-comparisonSplodge.cntY
        return math.hypot(xDist,yDist)

    #function to calculate and remmeber the hihgest intensity stack
    def calculatePrimaryZStack(self):
        currMaxIntensity = 0
        currBestFrame = -1
        for frameRef, splodge in self.splodges.items():
            if splodge.meanRectIntensity > currMaxIntensity:
                currBestFrame = frameRef
                currMaxIntensity = splodge.meanRectIntensity
        
        self.primaryZStack = currBestFrame
    
    def getOutputParameters(self, parameters):
        targetSplodge = self.splodges[self.primaryZStack]
        outputParameters = [str(self.primaryZStack)]
        for parameter in parameters:
            outputParameters.append(str(getattr(targetSplodge, parameter)))
        outputParameters.append(str(len(self.splodges)))
        return outputParameters



#container class for the centriole linkages
class centrioleContationer():
    

    #function to divide the splodges into a dict of lists seperating each by their zPos
    #note: using a predefined ordered dict protects us if theres an intermediary layer where nothing is identified
    #turns out this is a bit irrelevent but we're just gonna leave it in for time as we've accidentally already subordered them
    def splitSplodgesByZStacks(self, splodgeArr, minZ, maxZ):
        splitDict = OrderedDict.fromkeys(range(minZ,maxZ+1)) #+1 as its not inclusive at the top
        for key in splitDict.keys():
            splitDict[key] = [] #precondition

        for zArr in splodgeArr:
            for splodge in zArr:
                if splodge.excluded: continue #exclude excluded splodges
                splitDict[splodge.zPos].append(splodge) #slap it in the right spot
        return splitDict

    #function to generate and populate a pseudogrid off of the the SC parameters for optimised comparison
    #cell size defined by maxLength from sc, generates psuedogrid and metadata returns a dict
    def generatePseudoGrid(self,zPos, centrioleArr,sc):

        #First we gotta calcualte the grid parameters
        xRange = sc.splodgeStats['maxcX'] - sc.splodgeStats['mincX']
        yRange = sc.splodgeStats['maxcY'] - sc.splodgeStats['mincY']

        sideLength = sc.splodgeStats['maxLength']
        if 'rectArea' in sc.filterDict.keys(): #override max length if we've limited area
            if '<' in sc.filterDict['rectArea'].keys():
                sideLength = math.sqrt(sc.filterDict['rectArea']['<']) # should be a reasonable estimate of the max length of toleranced shapes

        nOfColumns = math.ceil(xRange/sideLength)
        nOfRows = math.ceil(yRange/sideLength)

        pseudoGridMetaData = {'minX':sc.splodgeStats['mincX'],'minY':sc.splodgeStats['mincY'],'xDiff':xRange,
            'yDiff':yRange, 'rows':nOfRows,'columns':nOfColumns}

        #now lets build the grid
        pseudoGrid = []
        for i in range(nOfRows+1):
            rowList = [] #initialise the row
            for j in range(nOfColumns+1):
                rowList.append([]) #initalises the empty containers
            pseudoGrid.append(rowList)
        
        #now we populate it
        for centriole in centrioleArr:
            rowIndex, columnIndex = centriole.getPseudoGridPosition(zPos,pseudoGridMetaData)
            pseudoGrid[rowIndex][columnIndex].append(centriole)
        
        return {'metadata':pseudoGridMetaData,'grid':pseudoGrid}


    #function to pickup and grab everything in a pseudogrid cell
    def getSpeciesInPseudoGrid(self, rowIndex,columnIndex,pseudoGrid):
        rtnArr = []
        for item in pseudoGrid[rowIndex][columnIndex]:
            rtnArr.append(item)
        return rtnArr

    #function to return a list of centrioles that maybe the parent of a splodge
    def getNearbyCentrioles(self, queryPos, pseudoGrid):
        
        candidateArr = []
        #so we're gonna check a 9 sqr about the querypos
        rowIndex = queryPos[0] #up/down
        colIndex = queryPos[1] #left/right
        
        topRow = rowIndex-1
        btmRow = rowIndex+1
        leftCol = colIndex-1
        rightCol = colIndex+1

        if topRow >= 0: #can we go above
            candidateArr.extend(self.getSpeciesInPseudoGrid(topRow,colIndex,pseudoGrid['grid']))
            if leftCol >= 0: #can we go left
                candidateArr.extend(self.getSpeciesInPseudoGrid(topRow,leftCol,pseudoGrid['grid']))
            if rightCol <= pseudoGrid['metadata']['columns']: #can we go right
                candidateArr.extend(self.getSpeciesInPseudoGrid(topRow,rightCol,pseudoGrid['grid']))
        
        if leftCol>= 0:
            candidateArr.extend(self.getSpeciesInPseudoGrid(rowIndex,leftCol,pseudoGrid['grid']))
        
        candidateArr.extend(self.getSpeciesInPseudoGrid(rowIndex,colIndex,pseudoGrid['grid']))

        if rightCol<=pseudoGrid['metadata']['columns']: #can we go right
            candidateArr.extend(self.getSpeciesInPseudoGrid(rowIndex,rightCol,pseudoGrid['grid']))

        if btmRow<=pseudoGrid['metadata']['rows']: #can we go down
            candidateArr.extend(self.getSpeciesInPseudoGrid(btmRow,colIndex,pseudoGrid['grid']))
            if leftCol >= 0: #can we go left
                candidateArr.extend(self.getSpeciesInPseudoGrid(btmRow,leftCol,pseudoGrid['grid']))
            if rightCol <= pseudoGrid['metadata']['columns']: #can we go right
                candidateArr.extend(self.getSpeciesInPseudoGrid(btmRow,rightCol,pseudoGrid['grid']))

        return candidateArr

    #function to find the parent centriole from a list of nearby centrioles
    #if we have none we'll say we've found none
    #if we have 1 with a centriod pos within the child then we pass that
    #if we have multiple with the centrioid within, we find the one with the lowest euclidean distance between centroids
    def getParentCentriole(self, orphanCentriole, parentCandidates):
        if parentCandidates == False: return False #cos returning something is truthy, this is falsey so allows handling later


        secondPhaseCandidates = [] #candidates that have a centroid in the child contour
        for centriole in parentCandidates:
            if cv.pointPolygonTest(orphanCentriole.contour, centriole.getPrevCentroidPos(),False) >= 0:
                secondPhaseCandidates.append(centriole)
        
        match len(secondPhaseCandidates):
            case 0:
                return False
            case 1:
                return secondPhaseCandidates[0]


        #now time to find the best of the second phase candidates 
        closestCentriole = None
        closestEuDist = 2147483647 #int limit
        for centriole in secondPhaseCandidates:
            currEuDist = centriole.getEuDistToPrevCentroid(orphanCentriole)
            if currEuDist < closestEuDist:
                closestEuDist = currEuDist
                closestCentriole = centriole
        
        return closestCentriole

    def __init__(self, sc):
        self.centrioles = {} #centrioles kept in a dict organised by frame

        #first we gotta loop through each of the frames
        for frameRef, splodgeArr in sc.splodges.items():
            #while splodgeArr contains multiple different Z-stacks, it is organised such that Z-stacks are in blocks and only go up
            #so we can split it reliably into subsets here to be smort :3
            splodgeDict = self.splitSplodgesByZStacks(splodgeArr,sc.splodgeStats['minZ'], sc.splodgeStats['maxZ'])
            for zRef, subSplodgeArr in splodgeDict.items(): #now we're looping through keys
                if zRef == sc.splodgeStats['minZ']: #for the first case we want to shortcut into just shoving everything into the centriole linkages as we only do vertical 
                    #linkages hence everything is gonna be the start of a centriole linkage, also allows us to inialise the arr
                    self.centrioles[frameRef] = []
                    for splodge in subSplodgeArr:
                        a = centrioleLinkage(splodge)
                        self.centrioles[frameRef].append(a)

                    continue #skip the rest of the stuff, yes we could put it in an else but this just makes it cleaner imo

                #so if we're here we're beyond the inital case, hence we need to start linking

                #first we're gonna check to see if the centrioles are still active, and grab them for comparaison
                centrioleLinakgesForComparison = []
                for centriole in self.centrioles[frameRef]:
                    if centriole.checkActive(): centrioleLinakgesForComparison.append(centriole)

                #now we're gonna create a pseudogrid out of the previous layer of active centrioles to optimise comparisons w/ the splodgeArr
                pseudoGridData = self.generatePseudoGrid(zRef-1, centrioleLinakgesForComparison, sc)

                #now lets look thorugh our new splodges to see what we can link up to centrioles or what makes new centrioles
                for splodge in subSplodgeArr:
                    splodgePosTuple = splodge.getPseudoGridPos(pseudoGridData['metadata'])

                    #so can we find something that aligns 
                    parentCentriole = self.getParentCentriole(splodge,self.getNearbyCentrioles(splodgePosTuple,pseudoGridData))
                    if parentCentriole:#ie there is something
                        parentCentriole.appendSplodge(splodge)
                    else:
                        self.centrioles[frameRef].append(centrioleLinkage(splodge))
        
        self.processPrimaryZStacks()

    #for each centriole go thoruhg and find the brightest 'most representative' zstack
    def processPrimaryZStacks(self):
        for frame, centrioleArr in self.centrioles.items():
            for centriole in centrioleArr:
                centriole.calculatePrimaryZStack()

    #output as Frame/Z/PARAMETERS
    def outputData(self,exportPath, filePrefix, outputParameters):
        with open(os.path.join(exportPath,filePrefix+".csv"),'x') as csvfile:
            writer = csv.writer(csvfile,dialect="excel")
            for frame, centrioleArr in self.centrioles.items():
                for centriole in centrioleArr:
                    row = [str(frame)]
                    row.extend(centriole.getOutputParameters(outputParameters))
                    writer.writerow(row)
                

                



    

### FUNCTIONS N'STUFF ###

#I mean they do what they say on the tin really...
def validateDirectory(dirPath):
    if os.path.isdir(dirPath) == False:
        print("Error: output directory not found" + "\n"
        +"path: "+dirPath)
        sys.exit()

#checks both if it exists and if its got the right file extention
def validateFile(filePath, fileExt):
    if os.path.isfile(filePath) == False:
        print("Error: file cannot be found" + "\n"
        +"path: "+filePath)
        sys.exit()
    
    trialExt = os.path.splitext(filePath)[1]
    if (trialExt == ("."+fileExt)) == False:
        print("Error: unexpected file type, ."+fileExt+" expected; " + trialExt+" recieved.")
        sys.exit()

#extracts the image and packages it into a 2d nested list of [f1[z1,z2,z3....],f2[...],...]
def extractInputImage(pilImage, nOfZStacks):
    imgData = []
    #prime the array
    for i in range(int(pilImage.n_frames/nOfZStacks)): #should check this maths at some point
        imgData.append([])

    #fill the array
    for i in range(0,pilImage.n_frames,nOfZStacks): #this takes each time poijnt
        for j in range(nOfZStacks): #this then takes the Z stacks within each time point
            pilImage.seek(i+j)
            imgData[int(i/17)].append(np.asarray(pilImage))
    
    return imgData

#function to process individual frames to identify contours around centrioles
#if blurkernal is 0, we don't do it 
#note: blurKernal must be a posiitve odd integer
def generateContours(frameMat, blurKernal,tophatKernal):
    newFrame = []
    if blurKernal == 0: #ie we don't have blur to do
        newFrame = frameMat
    else:
        if blurKernal % 2 != 1: #cos otherwise it'll break
            print("Error: blurKernal must be Odd")
        else:
            newFrame = cv.GaussianBlur(frameMat,(blurKernal,blurKernal),0) #blurr
    
    #cv.imshow("blurred frame",newFrame)
    #cv.waitKey(0)

    tophattedFrame =[]
    if not("tophatKernal" in locals()): #does the parameter exist?
        tophattedFrame = newFrame
    else: #if it does exist and we have hats to do
        tophattedFrame = cv.morphologyEx(newFrame, cv.MORPH_TOPHAT, tophatKernal)


    #threshold
    ret, thrFrame = cv.threshold(tophattedFrame,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
    #cv.imshow("thrFrame", thrFrame)
    #cv.waitKey(0)

    #find the contours
    contours, hierachy = cv.findContours(thrFrame, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)
    return contours
        

#function to produce a tif with the current modifications
#consider swapping this for frames and doing it at that level
def generatePreviewTif(inputImg,sc):
    newImg = []
    for frameRef in range(len(inputImg)):
        newFrame = []
        for zRef in range(Z_STACKS):
            currFrame = inputImg[frameRef][zRef].copy() #create a new object
            currRectCnts = sc.getRectContours(frameRef,zRef) #get the contours
            cv.drawContours(currFrame, currRectCnts, -1,(255,255,255),1) # draw contours
            newFrame.append(currFrame) #draw the frame
        newImg.append(newFrame)
    return newImg

#callbackfunction to change the displayed frame
def changeDisplayedFrame(value):
    #print("displaying frame: "+str(value))
    global G_CURRFRAME
    G_CURRFRAME = value
    cv.imshow(currWindow,previewImg[G_CURRFRAME][G_CURRZ])

#callbackfunction to change the displayed Z
def changeDisplayedZ(value):
    #print("displaying frame: "+str(value))
    global G_CURRZ
    G_CURRZ = value
    cv.imshow(currWindow,previewImg[G_CURRFRAME][G_CURRZ])

#callbackfunction to tell the splodecontainer alter the displayed frame
#so lets currently looking at doing this on intensity from max to min
#and lets go for contour intensity as I'd rather use that and remove the other
def excludeMinAveIntensity(value):
    #first we want to convert value which is between 0 and 100 into 
    #a real number based on the min/max
    minIntensity = sc.splodgeStats['minRectIntensity']
    maxIntensity = sc.splodgeStats['maxRectIntensity']
    diffIntensity = maxIntensity - minIntensity
    scaledVal = ((value/100)*diffIntensity)+minIntensity
    sc.excludeData({'meanRectIntensity':{'>':scaledVal}})
    global previewImg
    previewImg = generatePreviewTif(img, sc)
    changeDisplayedFrame(cv.getTrackbarPos('frame','previewWindow'))

def excludeMaxAveIntensity(value):
    #first we want to convert value which is between 0 and 100 into 
    #a real number based on the min/max
    minIntensity = sc.splodgeStats['minRectIntensity']
    maxIntensity = sc.splodgeStats['maxRectIntensity']
    diffIntensity = maxIntensity - minIntensity
    scaledVal = ((value/100)*diffIntensity)+minIntensity
    sc.excludeData({'meanRectIntensity':{'<':scaledVal}})
    global previewImg
    previewImg = generatePreviewTif(img, sc)
    changeDisplayedFrame(cv.getTrackbarPos('frame','previewWindow'))


#computed as percentile because range to great
def excludeMinArea(value):
    scaledVal = np.percentile(sc.splodgeStats['areas'],value)
    sc.excludeData({'rectArea':{'>':scaledVal}})
    global previewImg
    previewImg = generatePreviewTif(img, sc)
    changeDisplayedFrame(cv.getTrackbarPos('frame','previewWindow'))

#computed as percentile because range to great
def excludeMaxArea(value):
    scaledVal = np.percentile(sc.splodgeStats['areas'],value)
    sc.excludeData({'rectArea':{'<':scaledVal}})
    global previewImg
    previewImg = generatePreviewTif(img, sc)
    changeDisplayedFrame(cv.getTrackbarPos('frame','previewWindow'))


### SETTINGS ###
INPUT_PATH = "/Users/thomas/Documents/PhD/Data/Project Laura/Ana2/Highspeed/Center vs outside test/C2-S4rNbGreyStack.tif" #file of interest
OUTPUT_FOLDER = "/Users/thomas/Documents/PhD/Data/Project Laura/Ana2/Highspeed/Center vs outside test/Individual Z tracking" #dir to export into
OUTPUT_PREFIX = "AndorS4rNb_2" #prefix to put on output files

PRODUCE_OUTPUT_TIF = False #tif summarising identified species

Z_STACKS = 17 #n of Z stacks in the image (Note: this should be the 1/n number in FIJI)

BLUR = 11 #size of gaussian blur px - must be odd
STRUCTURING_FACTOR = 11 #size of structuring factor px diameter (I believe???)

EXPORT_PARAMETERS = ["meanRectIntensity",'cntX','cntY'] #parameters to epxort

MINIMUM_CENTRIOLE_FILTER = 1

### CODE ###

#validation of parameters
validateDirectory(OUTPUT_FOLDER)
validateFile(INPUT_PATH, 'tif')

pilImage = Image.open(INPUT_PATH)
img = extractInputImage(pilImage, Z_STACKS)
sc = splodgeContainer()

tophatKernal = cv.getStructuringElement(cv.MORPH_RECT,(STRUCTURING_FACTOR,STRUCTURING_FACTOR))
nOfFrames = len(img)
for frameRef in range(nOfFrames):
    for zRef in range(Z_STACKS):
        print("Current Analysis > FRAME: "+str(frameRef)+" Z: "+str(zRef))
        contours = generateContours(img[frameRef][zRef],BLUR,tophatKernal)
        sc.appendSplodges(contours,frameRef,zRef, img[frameRef][zRef])
        
previewImg = generatePreviewTif(img, sc)
currWindow = 'previewWindow'
#preview data
cv.namedWindow('previewWindow')
cv.createTrackbar('frame','previewWindow',0,nOfFrames-1,changeDisplayedFrame)
cv.createTrackbar('z','previewWindow',0,Z_STACKS,changeDisplayedZ)
cv.createTrackbar('minAveIntensity','previewWindow',0,100,excludeMinAveIntensity)
cv.createTrackbar('maxAveIntensity','previewWindow',0,100,excludeMaxAveIntensity)
cv.createTrackbar('minArea','previewWindow',0,100,excludeMinArea)
cv.createTrackbar('maxArea','previewWindow',0,100,excludeMaxArea)

cv.imshow('previewWindow', previewImg[0][0])
cv.waitKey(0)
cv.destroyWindow('previewWindow')

#right so at this point it should be thresholded and shizzle, so now what we need to do is the compression to centrioles
cc = centrioleContationer(sc)

cc.outputData(OUTPUT_FOLDER,OUTPUT_PREFIX,EXPORT_PARAMETERS)
