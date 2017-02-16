from ij import IJ
from ij import *
from ij import ImagePlus
from ij import ImageStack
from ij.measure import *
from ij.measure import ResultsTable
from ij.plugin import *
from ij.process import *
from ij.plugin import ChannelSplitter
from ij.measure import Measurements
from ij.plugin import ImageCalculator
from ij.plugin.frame import RoiManager
from ij.process import ImageProcessor
from ij.process import ImageStatistics
from ij.gui import Roi
from ij.gui import PointRoi
from ij.gui import OvalRoi
import math 
from java.awt import * 
from java.awt import Font
import itertools 
import os
import glob


start_number = 0	#Starting position for numbering of cells in report

noiseC1 = 50   #values from find maxima
noiseC2 = 50    ##values from find maxima
pixelsize = 0.092  ### insert pixel size here in µm
#pixelsize = imp.getCalibration().pixelWidth  ### Alternative: if image is calibrated, this can be activated instead
7
noisec1_1 = '"noise=' + str(noiseC1) + ' output=[Point Selection]"'   # generates the parameters for the find maxima dialog in Channel 1
noisec2_1 = '"noise=' + str(noiseC2) + ' output=[Point Selection]"'   # generates the parameters for the find maxima dialog in Channel 2

automatic_save_results = True

savepath = "L:/KRK/KRF"   # storage of the data , use "/" instead of "\" as path separator

### Do MaxEntropy thresholding prior to Maximum Detection? 
### Im proves MaxEntropy detection
####Attention: Use not in combinaton with intensity measurements, since this will eliminate cellular structures below the threshold
### Use in conjuction with Pearson might generate artifacts or division by 0 errors (should be catched by the pearsons calculation, but use at your own risk!)
DoThreshold = False


#### Do pearsons correlation between the 2 channels per cell or per spot?

DoPearson_per_Cell = False
DoPearson_per_Spot = False

Measure_Spot_Intensity = True

Intensity_Roi_Size = 5

### Should only the spot in Channel2 that is closest to a spot in Channel1 be measured and displayed?
Nearest_Neighbour_Only = True


## Should galleries be generated?
Generate_Gallery_Area = False
Generate_Gallery_Length = False #True
Generate_Gallery_SpotDistance = False #True
Generate_Gallery_Spotnumber_C1 = False #True
Generate_Gallery_Spotnumber_C2 = False #True




### Should Galleries be labelled by the ROI ID?
Label_Gallery = False

## Should the cell outlines be drawn?

DrawOutline = False


OutlineMaxima = False




# how big is the desired ROI around the cell?
roisize = 75

#### Initialize filtering criteria for gallery generation

## Filter by area, in um^2
Filter_Area_Min = 2
Filter_Area_Max = 10

## filter by cell length in um
Filter_Length_Min = 0
Filter_Length_Max = 10

## filter by spots under / over a certain distance
Filter_Distance_Min = 0
Filter_Distance_Max = 10

## filter by number of spots in Channel 1
C1_Filter_Spots_Min = 0
C1_Filter_Spots_Max = 5
#C1_Filter_Spots = 3

## filter by number of spots in Channel 2
C2_Filter_Spots_Min = 0
C2_Filter_Spots_Max = 5
##C2_Filter_Spots = 3




def ThresholdMaxEntropy(imp0):
    """Thresholds image and returns thresholded image, merge code still quite clumsy but functional"""
    imp0 = IJ.getImage()
    impthres = imp0.duplicate()
    imp01 = ImagePlus("Channel1", ChannelSplitter.getChannel(imp0, 1))
    imp02 = ImagePlus("Channel2", ChannelSplitter.getChannel(imp0, 2))
    imp001 = imp01.duplicate()
    imp002 = imp02.duplicate()
    IJ.setAutoThreshold(imp001, "MaxEntropy dark")
    IJ.run(imp001, "Convert to Mask", "");
    IJ.run(imp001, "Divide...", "value=255");
    IJ.setAutoThreshold(imp002, "MaxEntropy dark")
    IJ.run(imp002, "Convert to Mask", "");
    IJ.run(imp002, "Divide...", "value=255");
    ic = ImageCalculator()
    imp0001 = ic.run("Multiply create", imp01, imp001)
    ic2 = ImageCalculator()
    imp0002 = ic2.run("Multiply create", imp02, imp002)
    imp0001.copy()
    impthres.setC(1)
    impthres.paste()
    imp0002.copy()
    impthres.setC(2)
    impthres.paste()
    imp01.close()
    imp02.close()
    imp001.close()
    imp002.close()
    imp0001.close()
    imp0002.close()
    return impthres




#dist_min = 1 
#dist_max = 10 
 
def dist(p0, p1): 
    """ Calculates the  distance between two xy coordinates, each
    each coordinated supplied by a tupel"""
    return math.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2) 

#x = [] 
#y = [] 
 
# pearson code from from http://stackoverflow.com/questions/3949226/calculating-pearson-correlation-and-significance-in-python 
 
def average(x): 
    """calculates the averade of a list x"""
    assert len(x) > 0 
    return float(sum(x)) / len(x) 
 
def pearson_def(x, y):
    """Calculate pearson correlation between list x and y, lists need to have the same length"""
    assert len(x) == len(y) 
    n = len(x) 
    assert n > 0 
    avg_x = average(x) 
    avg_y = average(y) 
    diffprod = 0 
    xdiff2 = 0 
    ydiff2 = 0 
    for idx in range(n): 
        xdiff = x[idx] - avg_x 
        ydiff = y[idx] - avg_y 
        diffprod += xdiff * ydiff 
        xdiff2 += xdiff * xdiff 
        ydiff2 += ydiff * ydiff 
    if (xdiff2 * ydiff2) != 0:   # checks if one of the pixels is 0  intensity, avoids div0 error
        return round((diffprod / math.sqrt(xdiff2 * ydiff2) ),3)
    else:
        return "n/a" #returns "n/a" if one pixel value is 0
        

# function to get the pixel values in a ROI as a list 
# gets the bounding box and the corresponding mask for polygon ROI, then retrieves the intensity values of all pixels within the ROI as a list 
# Code adapted from Java implementation in Burger and Burge "Digital image processing" 
 
def pixel_values(roi): 
    """extract each pixel value from a polygon roi into a list
    does not work for square rois!"""
    pixel = [] 
    mask =  roi.getMask()    # polygon rois are defined by a mask
    box = roi.getBounds() 
    boxLeft = box.x 
    boxRight = boxLeft + box.width 
    boxTop = box.y 
    boxBottom = boxTop + box.height 
    for v in range (boxTop, boxBottom): 
        for u in range (boxLeft, boxRight):       
            if mask.getPixel(u - boxLeft, v - boxTop) > 0: 
                pixel.append(imp.getProcessor().getPixel(u,v)) 
    return pixel 
 
def pixel_values_rect(roi): 
    """extract each pixel value from a square roi into a list"""
    pixel = [] 
    mask =  roi.getMask() 
    box = roi.getBounds() 
    boxLeft = box.x 
    boxRight = boxLeft + box.width 
    boxTop = box.y 
    boxBottom = boxTop + box.height 
    for v in range (boxTop, boxBottom): 
        for u in range (boxLeft, boxRight):
            #if mask.getPixel(u - boxLeft, v - boxTop) > 0: 
            pixel.append(imp.getProcessor().getPixel(u,v)) 
    return pixel 
 


def GenerateGallery(imp, CoordinateList):  
    """ Generates a gallery from a list of coordinates:
    input values are:
    a coordinate list (x,y as list of tupels)
    an image from which the gallery is generated (imp),
    Image must (for now) be a 8 or 16 bit composite
    a ROI size for the gallery (roisize)"""
    if CoordinateList != []:    
        CellNumber = len(CoordinateList)
        channelnumber = imp.getNChannels()
        slicenumber = imp.getNSlices()
        bitdepth = imp.getBitDepth()
        imp2 = IJ.createHyperStack("Gallery",roisize, roisize, channelnumber, slicenumber, CellNumber, bitdepth)
        if bitdepth != 24:
            imp2.copyLuts(imp)
        timer = 0
        for cells in CoordinateList:
            timer = timer + 1
            roi2 = Roi(cells[0]-(roisize/2), cells[1]-roisize/2, roisize,roisize)
            imp.setC(1)
            imp.setRoi(roi2)
            imp.copy()
            imp2.setT(timer)
            imp2.setC(1)
            imp2.paste()
            if channelnumber > 1:
                imp.setC(2)
                imp.copy()
                imp2.setC(2)
                imp2.paste()
                if channelnumber > 2:
                   imp.setC(3)
                   imp.copy()
                   imp2.setC(3)
                   imp2.paste()
        imp2.show()
        if Label_Gallery:
            timer2 = 0
            for cells in CoordinateList:
                
                timer2 = timer2 + 1
                imp2.setT(timer2)
                imp2.setC(1)
                ip = imp2.getProcessor()
                Labelfont = Font("Arial", Font.PLAIN, 12)
                ip.setFont(Labelfont)
                ip.setColor(65000)
                ROINumber = str(cells[2])
                ip.drawString(ROINumber[:4], 5, 14)
                imp.updateAndDraw()
            
        return imp2




######################################################
### End of function definitions ######################
######################################################


################################################################
###  Get original image, then either threshold it first using 
###  MaxEntropy autothreshold or work directly on the original #
### All work is done none-destrcutively on a copy of the original


imp_orig = IJ.getImage()
IJ.run(imp_orig, "Select None", "");
if DoThreshold:
    imp = ThresholdMaxEntropy (imp_orig)
else:
    imp = imp_orig.duplicate()
imp.show()
IJ.run(imp, "Select None", "");





###########################################################

## Initialize lists for gallery generation
Filtered_Cells_Area = []   # Cells filtered by cell area
Filtered_Cells_Length = []   # cells filtered by cell length
Filtered_Cells_SpotDistance = []   # Cell filtered by minimal / maximal distance between spots
C1_Filtered_Cells_Spots = []  # Number of spots in Channel 1
C2_Filtered_Cells_Spots = [] # Number of spots in Channel 2
##  End of list initialisation

## set pixel calibration based on the entry on top
cal = imp.getCalibration()
cal.pixelHeight = pixelsize
cal.pixelWidth = pixelsize
pixelWidth = imp.getCalibration().pixelWidth
print pixelWidth


coordinates = [] 

### setting up the results tables ####

ort = ResultsTable() 
ort.setPrecision(0) 
#print ort.getCounter 
ort.setHeading(0, "Cell") 
ort.setHeading(1, "Point_C1") 
ort.setHeading(2, "Point_C2") 
ort.setHeading(3, "Distance in um") 
ort.setHeading(4, "Distance in pixel") 
ort.setHeading(5, "Area in um") 
ort.setHeading(6, "Feret in um") 

#pixelWidth = imp.getCalibration().pixelWidth 

 
ort2 = ResultsTable() 
ort2.setPrecision(3) 
ort2.setHeading(0, "Cell") 
ort2.setHeading(1, "Channel") 
ort2.setHeading(2, "Point #") 
ort2.setHeading(3, "Coordinates") 
ort2.setHeading(4, "Distance from center in um") 
if DoPearson_per_Spot:
    ort2.setHeading(5, "Pearson of spot") 

threshold = ImageProcessor.NO_THRESHOLD
 
allPearson = [] 

cell = start_number   # used as a counter for each measured cell

######### End of initial declarations
### Start of the real program
### Each ROI in the ROI manager is actiavted, measured and  maxima are identified. Diverse Measurements are made in addition


for roi in RoiManager.getInstance().getRoisAsArray(): 
    HasFilteredSpots = False    # will become True in case some of the spots correspond the the filter criteria defined above
    roiname=RoiManager.getInstance().getName(cell)   # Find the name of the currently active ROI
    item1counter = 0 
    item2counter = 0 
    pixelC1 = [] 
    pixelC2 = [] 
    cell = cell+1 
    imp.setC(1)   
    imp.setRoi(roi) 
    if DoPearson_per_Cell:
        pixelC1 = pixel_values(roi) 
    stats = imp.getStatistics(Measurements.MEAN | Measurements.AREA | Measurements.FERET | Measurements.CENTROID | Measurements.ELLIPSE)    ## Measurements using the ImageJ Stats module
    IJ.run(imp, "Find Maxima...", noisec1_1)    # Finds the maxima in Channel 1, gets them as array of point ROIs
    points_C1 = []
    proi =  IJ.getImage().getRoi() 
    if proi.getClass() == PointRoi: 
        #### For each point ROI identified by FindMaxima, the coordinates are added to a list
        px = proi.getXCoordinates() 
        py = proi.getYCoordinates() 
        bounds = proi.getBounds() 
        #points_C1 = [] 
        for i in range(proi.getNCoordinates()): 
           points_C1.append((bounds.x + px[i], bounds.y + py[i])) 
           #points.append((px[i], py[i])) 
        #print points_C1 
    ### reading out the measured statistics
    means = stats.mean       
    area = stats.area 
    ellipse_length = stats.major
    ellipse_width = stats.minor
    feret = roi.getFeretsDiameter() 
    boundingbox = roi.getBounds() 
    centerx = round(stats.xCentroid/pixelWidth , 3)  
    centery = round(stats.yCentroid/pixelWidth , 3)    
    center = (centerx,centery)  
    center_ROInumber = (centerx,centery,roiname)
    print cell, center
    print cell, center_ROInumber

    ###reading out the measured statistics
    IJ.run(imp, "Select None", "");
    imp.setC(2)
    imp.setRoi(roi)
    if DoPearson_per_Cell:
        pixelC2 = pixel_values(roi)
        pearson = (pearson_def(pixelC1, pixelC2)) # calculate pearson of the cell ROI in channel 1 and 2
        allPearson.append(pearson) 
    IJ.run(imp, "Find Maxima...", noisec2_1)
    proi2 =  IJ.getImage().getRoi()
    points_C2 = []
    if proi2.getClass() == PointRoi:
        px2 = proi2.getXCoordinates()
        py2 = proi2.getYCoordinates()
        bounds2 = proi2.getBounds()
        #points_C2 = []
        for i in range(proi2.getNCoordinates()):
           points_C2.append((bounds2.x + px2[i], bounds2.y + py2[i]))
           #points.append((px[i], py[i]))
        #print points_C2
    if points_C1 == [] or points_C2 == []:
        ort.incrementCounter()
        ort.addValue("Cell", cell)
        if points_C1 == []:
            print 'No spots found in Channel1 in Cell' + str(cell)
        if points_C2 == []:
            print 'No spots found in Channel2 in Cell' + str(cell)
        ort.addValue("Point_C1", 'Invalid')
        ort.addValue("Distance in um", 'Invalid')
        ort.addValue("Distance in pixel",'Invalid')
        ort.addValue("Area in um", area)
        ort.addValue("Feret in um", feret)
        ort.addValue("Fittet Elipse - Length", str(round(ellipse_length, 3)))
        ort.addValue("Fittet Elipse - Width", str(round(ellipse_width, 3)))        
        ort.addValue("Center", str(center))
        ort.addValue("Center - Point_C1 in um", 'Invalid')
        ort.addValue("Center - Point_C2 in um", 'Invalid')
        ort.addValue("Center - Point_C1", 'Invalid')
        ort.addValue("Center - Point_C2", 'Invalid')
        if DoPearson_per_Cell:
            ort.addValue("Pearson per cell", str(pearson))
        if DoPearson_per_Spot:
            ort.addValue("Pearson_Spot_C1", 'Invalid')
            ort.addValue("Pearson_Spot_C2", 'Invalid') 
        ort.addValue("Point_C2", 'Invalid')  
        ort.addValue("Point_C1", 'Invalid')
    

## point-based measurements and results tables



    if points_C1 != [] and points_C2 != []:
        for item1 in points_C1:
            pixels_Spots1_C1 = []
            pixels_Spots1_C2 = []
            imp.setC(1)
            roi1 = Roi(item1[0]-2, item1[1]-2,5,5)
            if DoPearson_per_Spot:
                imp.setRoi(roi1)
                pixels_Spots1_C1 = pixel_values_rect(roi1)
                imp.setC(2)
                imp.setRoi(roi1)
                pixels_Spots1_C2 = pixel_values_rect(roi1)
                pearsons_spot1 = (pearson_def(pixels_Spots1_C1, pixels_Spots1_C2))
            #Distance from center
            distance_center_spot1 = round(dist(item1, center),3)
            
            if Nearest_Neighbour_Only:
            #Create new list with distances to points in C2 and find index to shortest one
                f = lambda pC2: dist(item1, pC2)
                points_C2_min = map(f, points_C2)
                index_C2 = points_C2_min.index(min(points_C2_min))
                #item2 is set to contain coordinate to nearest point in C2
                item2 = points_C2[index_C2]
                pixels_Spots2_C1 = []
                pixels_Spots2_C2 = []
                if DoPearson_per_Spot:
                        imp.setC(1)
                        roi2 = Roi(item2[0]-3, item2[1]-3, 5,5)
                        imp.setRoi(roi2)
                        pixels_Spots2_C1 = pixel_values_rect(roi2)
                        imp.setC(2)
                        imp.setRoi(roi2)
                        pixels_Spots2_C2 = pixel_values_rect(roi2)
                        pearsons_spot2 = (pearson_def(pixels_Spots2_C2, pixels_Spots2_C1))
                distance_center_spot2 = round(dist(item2, center),3)
                distance = dist(item1, item2)
                distance_in_um = (distance * pixelWidth)
                result = 'Cell '+ str(cell), item1, item2, distance
                ort.incrementCounter()
                ort.addValue("Cell", cell)
                ort.addValue("Point_C1", str(item1))
                ort.addValue("Point_C2", str(item2))
                ort.addValue("Distance in um", str(round((distance*pixelWidth),3)))
                ort.addValue("Distance in pixel", distance)
                ort.addValue("Area in um", area)
                ort.addValue("Feret in um", feret)
                ort.addValue("Fittet Elipse - Length", str(round(ellipse_length, 3)))
                ort.addValue("Fittet Elipse - Width", str(round(ellipse_width, 3)))                       
                ort.addValue("Center", str(center))
                ort.addValue("Center - Point_C1 in um", str(round((distance_center_spot1*pixelWidth),3)))
                ort.addValue("Center - Point_C2 in um", str(round((distance_center_spot2*pixelWidth),3)))
                ort.addValue("Center - Point_C1", str(distance_center_spot1))
                ort.addValue("Center - Point_C2", str(distance_center_spot2))
                if DoPearson_per_Cell:
                    ort.addValue("Pearson per cell", str(pearson))
                if DoPearson_per_Spot:
                    ort.addValue("Pearson_Spot_C1", str(pearsons_spot1))
                    ort.addValue("Pearson_Spot_C2", str(pearsons_spot2))      
                if (distance_in_um > Filter_Distance_Min) and (distance_in_um < Filter_Distance_Max):
                         HasFilteredSpots = True
                
                
            else:    
    
                for item2 in points_C2:
                    pixels_Spots2_C1 = []
                    pixels_Spots2_C2 = []
                    if DoPearson_per_Spot:
                        imp.setC(1)
                        roi2 = Roi(item2[0]-2, item2[1]-2, 5,5)
                        imp.setRoi(roi2)
                        pixels_Spots2_C1 = pixel_values_rect(roi2)
                        imp.setC(2)
                        imp.setRoi(roi2)
                        pixels_Spots2_C2 = pixel_values_rect(roi2)
                        pearsons_spot2 = (pearson_def(pixels_Spots2_C2, pixels_Spots2_C1))
                    distance_center_spot2 = round(dist(item2, center),3)
                    distance = dist(item1, item2)
                    distance_in_um = (distance * pixelWidth)
                    result = 'Cell '+ str(cell), item1, item2, distance
                    ort.incrementCounter()
                    ort.addValue("Cell", cell)
                    ort.addValue("Point_C1", str(item1))
                    ort.addValue("Point_C2", str(item2))
                    ort.addValue("Distance in um", str(round((distance*pixelWidth),3)))
                    ort.addValue("Distance in pixel", distance)
                    ort.addValue("Area in um", area)
                    ort.addValue("Feret in um", feret)
                    ort.addValue("Fittet Elipse - Length", str(round(ellipse_length, 3)))
                    ort.addValue("Fittet Elipse - Width", str(round(ellipse_width, 3)))       
                    ort.addValue("Center", str(center))
                    ort.addValue("Center - Point_C1 in um", str(round((distance_center_spot1*pixelWidth),3)))
                    ort.addValue("Center - Point_C2 in um", str(round((distance_center_spot2*pixelWidth),3)))
                    ort.addValue("Center - Point_C1", str(distance_center_spot1))
                    ort.addValue("Center - Point_C2", str(distance_center_spot2))
                    if DoPearson_per_Cell:
                        ort.addValue("Pearson per cell", str(pearson))
                    if DoPearson_per_Spot:
                        ort.addValue("Pearson_Spot_C1", str(pearsons_spot1))
                        ort.addValue("Pearson_Spot_C2", str(pearsons_spot2))      
                    if (distance_in_um > Filter_Distance_Min) and (distance_in_um < Filter_Distance_Max):
                         HasFilteredSpots = True
                
                
    if points_C1 != []:
        for item1 in points_C1:
            item1counter = item1counter + 1
            pixels_Spots1_C1 = []
            pixels_Spots1_C2 = []
            imp.setC(1)
            roi1 = Roi(item1[0]-2, item1[1]-2,5,5)
            Intensity_ROI = OvalRoi(item1[0]-(Intensity_Roi_Size/2), item1[1]-(Intensity_Roi_Size/2), Intensity_Roi_Size,Intensity_Roi_Size)
            if DoPearson_per_Spot:
                imp.setRoi(roi1)
                pixels_Spots1_C1 = pixel_values_rect(roi1)
                imp.setC(2)
                imp.setRoi(roi1)
                pixels_Spots1_C2 = pixel_values_rect(roi1)
                pearsons_spot1 = (pearson_def(pixels_Spots1_C1, pixels_Spots1_C2))
            if Measure_Spot_Intensity:
                imp.setC(1)
                imp.setRoi(Intensity_ROI)
                stats = imp.getStatistics(Measurements.MEAN | Measurements.AREA | Measurements.FERET | Measurements.CENTROID)
                spot_mean1 =stats.mean
            
            #Distance from center
            distance_center_spot1 = round(dist(item1, center),3)
            ort2.incrementCounter()
            ort2.addValue("Cell", cell)
            ort2.addValue("Length", ellipse_length)
            ort2.addValue("Channel", '1')
            ort2.addValue("Point #", str(item1counter))
            ort2.addValue("Coordinates", (str(item1[0])+ ";" + str(item1[1])))
            ort2.addValue("Distance from center in um", str(round((distance_center_spot1*pixelWidth),3)))
            if DoPearson_per_Spot:
                ort2.addValue("Pearson of spot", str(pearsons_spot1))
            if Measure_Spot_Intensity:
                ort2.addValue("Mean Intensity of spot", str(spot_mean1))
    else:
        ort2.incrementCounter()
        ort2.addValue("Cell", cell)
        ort2.addValue("Length", ellipse_length)
        ort2.addValue("Channel", '1')
        ort2.addValue("Point #", 'No maximum found')
        ort2.addValue("Coordinates", 'No maximum found')
        ort2.addValue("Distance from center in um", 'No maximum found')
        if DoPearson_per_Spot:
            ort2.addValue("Pearson of spot", 'No maximum found')
        if Measure_Spot_Intensity:
            ort2.addValue("Mean Intensity of spot", '')


    if points_C2 != []:        
        for item2 in points_C2:
            item2counter = item2counter + 1
            pixels_Spots2_C1 = []
            pixels_Spots2_C2 = []
            Intensity_ROI = OvalRoi(item2[0]-(Intensity_Roi_Size/2), item2[1]-(Intensity_Roi_Size/2), Intensity_Roi_Size,Intensity_Roi_Size)
            if DoPearson_per_Spot:
                imp.setC(1)
                roi2 = Roi(item2[0]-2, item2[1]-2,5,5)
                imp.setRoi(roi2)
                pixels_Spots2_C1 = pixel_values_rect(roi2)
                imp.setC(2)
                imp.setRoi(roi2)
                pixels_Spots2_C2 = pixel_values_rect(roi2)
                pearsons_spot2 = (pearson_def(pixels_Spots2_C2, pixels_Spots2_C1))
            #Distance from center
            if Measure_Spot_Intensity:
                imp.setC(2)
                imp.setRoi(Intensity_ROI)
                stats = imp.getStatistics(Measurements.MEAN | Measurements.AREA | Measurements.FERET | Measurements.CENTROID)
                spot_mean2 =stats.mean
            
            distance_center_spot2 = round(dist(item2, center),3)
            ort2.incrementCounter()
            ort2.addValue("Cell", cell)
            ort2.addValue("Length", ellipse_length)
            ort2.addValue("Channel", '2')
            ort2.addValue("Point #", str(item2counter))
            ort2.addValue("Coordinates", (str(item2[0])+ ";" + str(item2[1])))
            ort2.addValue("Distance from center in um", str(round((distance_center_spot2*pixelWidth),3)))
            if DoPearson_per_Spot:
                ort2.addValue("Pearson of spot", str(pearsons_spot2))
            if Measure_Spot_Intensity:
                ort2.addValue("Mean Intensity of spot", str(spot_mean2))
    else:
        ort2.incrementCounter()
        ort2.addValue("Cell", cell)
        ort2.addValue("Length", ellipse_length)
        ort2.addValue("Channel", '2')
        ort2.addValue("Point #", 'No maximum found')
        ort2.addValue("Coordinates", 'No maximum found')
        ort2.addValue("Distance from center in um", 'No maximum found')
        if DoPearson_per_Spot:
            ort2.addValue("Pearson of spot", 'No maximum found')
        if Measure_Spot_Intensity:
            ort2.addValue("Mean Intensity of spot", '')
    if DrawOutline:
        imp.setC(1)
        imp.setRoi(roi)
        IJ.run("Draw")
        imp.setC(2)
        imp.setRoi(roi)
        IJ.run("Draw")


     
    if OutlineMaxima:
        if points_C1 != []:        
            for item1 in points_C1:
                roi_outline_C1 = Roi(item1[0]-2, item1[1]-2,5,5)
                imageproc = imp.getProcessor()
                imageproc.setRoi(roi_outline_C1)
                imageproc.setValue(65535)
                imageproc.draw(roi_outline_C1)
        if points_C2 != []:        
            for item2 in points_C2:
                roi_outline_C2 = Roi(item2[0]-2, item2[1]-2,5,5)
                imp.setC(2)
                imageproc = imp.getProcessor()
                imageproc.setRoi(roi_outline_C2)
                imageproc.setValue(65535)
                imageproc.draw(roi_outline_C2)                
   
    
    if (area > Filter_Area_Min) and (area < Filter_Area_Max):
        Filtered_Cells_Area.append(center_ROInumber)
        
    if (area > Filter_Length_Min) and (area < Filter_Length_Max):
        Filtered_Cells_Length.append(center_ROInumber)
     
    if (len(points_C1) > C1_Filter_Spots_Min) and (len(points_C1)< C1_Filter_Spots_Max):
        C1_Filtered_Cells_Spots.append(center_ROInumber)
    if (len(points_C2) > C2_Filter_Spots_Min) and (len(points_C2)< C2_Filter_Spots_Max):
        C2_Filtered_Cells_Spots.append(center_ROInumber)
        
    if HasFilteredSpots:
        Filtered_Cells_SpotDistance.append(center_ROInumber)
        

 

## Generates a gallery for cells filtered by area
if Generate_Gallery_Area:
    print Filtered_Cells_Area
    imp2 = GenerateGallery(imp, Filtered_Cells_Area)
    if imp2 != None:
        imp2.setTitle("Cells filtered by cell area")

if Generate_Gallery_Length:
    imp3 = GenerateGallery(imp, Filtered_Cells_Length)
    if imp3 != None:    
        imp3.setTitle("Cells filtered by cell length")

if Generate_Gallery_SpotDistance:
    imp4 = GenerateGallery(imp,Filtered_Cells_SpotDistance)
    if imp4 != None:
        imp4.setTitle("Cells filtered by Spot Distance")
    
if Generate_Gallery_Spotnumber_C1:
    imp5 = GenerateGallery(imp,C1_Filtered_Cells_Spots)
    if imp5 != None:
        imp5.setTitle("Cells filtered by number of spots in Channel 1")

if Generate_Gallery_Spotnumber_C2:
    imp6 = GenerateGallery(imp,C2_Filtered_Cells_Spots)
    if imp6 != None:
        imp6.setTitle("Cells filtered by number of spots in Channel 2")

#print "Cells filtered by Area: " + str(Filtered_Cells_Area)

#print "Cells filtered by length: " + str(Filtered_Cells_Length)

#print "Cells filtered by number of spots in Channel 1: " + str(C1_Filtered_Cells_Spots)

#print "Cells filtered by number of spots in Channel 2: " + str(C2_Filtered_Cells_Spots)

#print "Cells filtered by distance: " + str(Filtered_Cells_SpotDistance)

#print allPearson

#print average(allPearson)
ort.show("Distance map")
ort2.show("Point measurements")

"""
Saving the Results tables:
Result tables will be saved to the directory specified in "savepath"


"""
dataname = imp_orig.getShortTitle()

if automatic_save_results:
    # Gather filenames of the savedirectory
    filename_ort = dataname+"_Distance_Map_001.csv"
    savename_ort = savepath+"/"+filename_ort # Generate complete savepath
    print savename_ort
    ort.saveAs(savename_ort) # save

        # Gather filenames of the savedirectory
    # Gather filenames of the savedirectory
    filename_ort2 = dataname+"_PointMeasurements_001.csv"
    savename_ort2 = savepath+"/"+filename_ort2 # Generate complete savepath
    print savename_ort2
    ort2.saveAs(savename_ort2) # save


IJ.run(imp, "Select None", "");
imp.changes = False


    

