#!/usr/bin/env python

#TODO - We need to check waht type of image we have (IR, VNIR, Other)
#TODO - Atmospheric parameters (http://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=MRO-M-CRISM-3-RDR-TARGETED-V1.0)
#TODO - Get non-implemented coded

'''
CRISM dataset info:

The S detector is the VNIR detector for TRR3 data products.  With the S detector
 the following are supported:
 
 ['R770', 'RBR', 'BD530', 'SH600', 'BD640', 'BD860', 'BD920', 'RPEAK1', 'BDI1000VIS','R440', 'IRR1']
 
 The L detector is the IR detector for TRR3 data products.  With the L detector
  the following are supported:
  
  ['BD11000IR', 'IRA', 'OLINDEX', 'LCPINDEX', 'HCPINDEX', 'VAR', 'ISLOPE1', 'BD1435', 'BD1500', 'ICER1', 'BD1750', 'BD1900', 'BDI2000', 'BD2100', 'BD2210', 'BD2290', 'Ds2300', 'SINDEX', 'ICER2', 'BDCARB','BD3000', 'BD3100', 'BD3200', 'BD3400', 'CINDEX', 'BD1270O2', 'BD1400H2O', 'BD2000CO2', 'BD2350', 'BD2600', 'IRR2', 'R2700', 'BD2700', 'IRR3']

The J detector is unsupported currently.
'''

#Internal imports
import sys
import argparse
import string
import os

#Module Imports
import csas_algorithms as algorithms

#External imports
try:
    from osgeo import gdal
    from osgeo.gdalconst import *
    version_num = int(gdal.VersionInfo('VERSION_NUM'))
    if version_num <1800 :
        print 'ERROR: Python bindings of GDAL version 1.8.0 or later required'
        raise
    else:
        pass
except ImportError:
    print "GDAL and the GDAL python bindings must be installed."
    raise

try:
    import numpy
except ImportError:
    print "NumPY must be installed."
    raise

try:
    import scipy
    import scipy.stats as stats
except ImportError:
    print "Some functionality will not work without scipy installed."

class GdalIO(object):
    
    def __init__(self, inputdataset):
        """Initialize Class Variables"""
        self.inputds = inputdataset
    
    def load(self):
        """Method to open any GDAL supported raster dataset"""

        #Open the dataset read only using GDAL
        dataset = gdal.Open(self.inputds, gdal.GA_ReadOnly)
        return dataset
    
    def info(self,dataset):
        xsize = dataset.RasterXSize
        ysize = dataset.RasterYSize
        bands = dataset.RasterCount
        metadata = dataset.GetMetadata()
        return xsize, ysize, bands, metadata
    
    def create_output(self,xsize,ysize,name,bands=1,dtype=gdal.GDT_Float32):

        """Method to create an output of the same type, size, projection, and transformation as the input dataset"""
        driver = gdal.GetDriverByName('GTIFF')
        outdataset = driver.Create(name,xsize,ysize,bands,dtype)
        #outdataset = driver.Create(outputname, xsize, ysize)
        
        return outdataset 

def getarray(band):
    #Get NDV, GetBand, MaskNDV
    NDV = band.GetNoDataValue()
    if NDV != None:
        array = band.ReadAsArray().astype(numpy.float64)
        array = numpy.ma.masked_where(array == NDV, array, copy=False)
        array = numpy.ma.masked_where(array == 65535, array, copy=False)
    else:
        array = band.ReadAsArray().astype(numpy.float64)
        array = numpy.ma.masked_where(array == 65535, array, copy=False)
    return array

def parseddtype(dtype):
    _dtypelookup = {'Byte': gdal.GDT_Byte,
                   'Int16' : gdal.GDT_Int16,
                   'UInt16' : gdal.GDT_UInt16,
                   'Int32' : gdal.GDT_Int32,
                   'UInt32' : gdal.GDT_UInt32,
                   'Float32' : gdal.GDT_Float32 
                   }
    dtype = _dtypelookup[dtype]
    return dtype

def getbandnumbers(wavelengths, *args):
    '''
    This parses the wavelenth list,finds the mean wavelength closest to the 
    provided wavelength, and returns the index of that value.  One (1) is added
    to the index to grab the correct band.
    
    Parameters
    ----------
    wavelengths: A list of wavelengths, 0 based indexing
    *args: A variable number of input wavelengths to map to bands
    
    Returns
    -------
    bands: A variable length list of bands.  These are in the same order they are
    provided in.  Beware that altering the order will cause unexpected results.
    
    '''
    bands = []
    for x in args:
        bands.append(min(range(len(wavelengths)), key=lambda i: abs(wavelengths[i]-x)))
    return bands

def get_wavelengths(newlbl,sw):
        sensor_id = None ; slk = None
        for line in newlbl:
            if "MRO:SENSOR_ID" in line:
                sensor_id =  line.split("\"")[1]
            if "SPACECRAFT_CLOCK_START_COUNT" in line:
                slk = int(line.split("\"")[1].split("/")[1].split(".")[0])
                
        #Select only those SW files that match the sensor    
        sw = [cdr6 for cdr6 in sw if sensor_id in cdr6.split("_")[4]]
        #Select only those SW files that are newest in version
        maxversion = 0
        for cdr6 in sw:
            maxversion = int(cdr6.split("_")[5].split(".")[0])   
        sw = [cdr6 for cdr6 in sw if str(maxversion) in cdr6.split("_")[5].split(".")[0]]
        
        #Select the newest CDR6 with a clock time <= image slk
        wavelength_tab = []
        for cdr6 in sw:
            cdr6_slk = int(cdr6.split("_")[2])
            if cdr6_slk < slk:
                wavelength_tab.append(cdr6)
        newlbl.close()
        
        wavelength_tab = wavelength_tab[-1] #Get the most recent version
        return wv2band('SW/' + wavelength_tab, sensor_id),sensor_id
    
def wv2band(wvlength,sensor_id): 
    '''This function opens the wavelength tab file and parses the band
    lookups.  A list is returned where the wavelength index is equal to 
    the band number -1.  This is because python lists are 0 based and GDAL's
    band reading is 1 based.
    '''
    wv = open(wvlength)
    wv2bnd_list = []
    for line in wv:
        linelst = string.split(line, ',')
        wv2bnd_list.append(float(linelst[1]))
    wv.close()
    if sensor_id == "L":
        wv2bnd_list.reverse()
        
    wv2bnd_out = []
    for x in wv2bnd_list:
        if x != 65535.0:
            wv2bnd_out.append(x)
    print wv2bnd_list
    return wv2bnd_list

def metadata2band(metadata):
    #The metadata for each IMG contains a dict with the band mappings, we just need to parse it.
    wv2band = []
    for key, value in metadata.iteritems():
        wv2band.append(float(value))
    wv2band.sort(key=int)
    return wv2band

def scale(array, scalemin, scalemax):
    '''
    This function is used to scale the data between a user defined range [c,d].  By default this range is between 1 and 255.  This maintains 0 as a special no data value should the user wish to set it.
    
    Scale from (a,b) to c,d)
        y = mx+b
        x = ((x-a)(d-c)/(b-a))+c
        
    Returns a scaled ndarray
    '''
    arraymin = array.mean() - 4 * array.std()
    arraymax = array.mean() + 4 * array.std()
    array = (((scalemax - scalemin)*(array-arraymin)) / (arraymax-arraymin)) + scalemin
    return array

def linear(array,scalemin,scalemax):
    '''
    This function exists to perform a linear stretch on the data to 2std from the mean.
    '''
    scalemax = array.mean() + (2* array.std())
    scalemin = array.mean() - (2* array.std())
    array = (array - array.min()) * ((scalemax - scalemin)/ (array.max()-array.min())) + scalemin 
    return array

def parseargs():
    '''
    This function utilizes argparse to parse user input arguments, store them to
    variables, and handle default values.
    
    Returns
    --------
    Instance of the parser, accessed vis args.keyword, where keyword is the dest variable.
    '''
#Setup the arg parser
    parser = argparse.ArgumentParser(description='CRISM Spectral Analysis Tool')
    parser.add_argument('input_data', action='store', help='The input CRISM data to be processed.')
    parser.add_argument('outputname', action='store', nargs='?',default=None, help='The output file name. Defaults are provided if this is not specified.')
    parser.add_argument('--wavelength_tab', action='store',nargs=1,default=None, help='The wavelength lookup table for the input_cube.')
    
    #Helper Parameters
    helpergroup = parser.add_argument_group('Helper Functionality')
    
    helpergroup.add_argument('--scale', '-s', action='store_true', dest='scale', default=False, help='Linearly remap the data to between two values.  The default is 0 - 255.')
    helpergroup.add_argument('--range', action='store', dest='scalerange',nargs=2, default=[0,255], help='Linearly remap the data to between two values.  The default is 0 - 255.')
    helpergroup.add_argument('--NDV', action='store', dest='newNDV',nargs=1, help='Override the intrinsic no data value or create a new no data value.')
    helpergroup.add_argument('-ot', action='store', dest='dtype',default=None, help='A GDAL output format. (Byte, Int16, Float32 are likely candidates.' )
    
    #Algorithms
    #Surface Parameters
    surfacegroup = parser.add_argument_group('Surface Parameters')
    
    surfacegroup.add_argument('--rockdust1','--R770', action='store_true', dest='rockdust1', default=False, help='Parameter: 0.77micron band, Rationale: rock/dust.')
    surfacegroup.add_argument('--rockdust2','--RBR', action='store_true', dest='rockdust2', default=False, help='Parameter: red/blue ratio, Rationale: rock.dust')
    surfacegroup.add_argument('--BD530', action='store_true', dest='bd530', default=False, help='Parameter: 0.53 micron band depth')
    surfacegroup.add_argument('--SH600', action='store_true', dest='sh600', default=False, help='Parameter: 0.60 micron shoulder height')    
    surfacegroup.add_argument('--BD640', action='store_true', dest='bd640', default=False, help='Parameter: 0.64 micron band depth')     
    surfacegroup.add_argument('--BD860', action='store_true', dest='bd860', default=False, help='Parameter: select ferric minerals')   
    surfacegroup.add_argument('--BD920', action='store_true', dest='bd920', default=False, help='Parameter: select ferric minerals')  
    surfacegroup.add_argument('--IRAlbedo', '--IRA', action='store_true', dest='ira', default=False, help='Parameter: 1.3micron reflectance, Rationale: IR albedo.') 
    surfacegroup.add_argument('--OlIndex', '--olivine_index', action='store_true', dest='olivine_index', default=False, help='Parameter:Olivine Index, Rationale: Olivine will be strongly positive; based on fayalite.')    
    surfacegroup.add_argument('--OlIndex2', '--olivine_index2', action='store_true', dest='olivine_index2', default=False, help='Parameter:Olivine Index, Rationale: Olivine will be strongly positive; based on fayalite. TRDRv3 only.')  
    surfacegroup.add_argument('--LCP', '--lcp_index', action='store_true', dest='lcp_index', default=False, help='Parameter: pyroxene index, Rationale:pyroxene will be strongly positive; favors LCP')
    surfacegroup.add_argument('--HCP', '--hcp_index', action='store_true', dest='hcp_index', default=False, help='Parameter: spectral variance.')
    surfacegroup.add_argument('--ISLOPE1', '--ferric1', action='store_true', dest='islope1', default=False, help='Parameter: -1*Spectral Slope_1, Rationale: Ferric coating on dark rocks.')
    surfacegroup.add_argument('--BD1435', '--co2_surface_ice', action='store_true', dest='bd1435', default=False, help='Parameter: 1.435 micron band depth, Rationale: CO2 surface ice.') 
    surfacegroup.add_argument('--BD1500', '--h20_surface_ice', action='store_true', dest='bd1500', default=False, help='Parameter: 1.5 micron band depth, Rationale: H2O surface ice.') 
    surfacegroup.add_argument('--ICER1', '--h2o_co2_ice', action='store_true', dest='icer1', default=False, help='Parameter:1.5micron and 1.43micron band ratio,Rationale: CO2, H2O Ice Mixtures.')
    surfacegroup.add_argument('--BD1750', '--gypsum', action='store_true', dest='bd1750', default=False, help='Parameter: 1.75 micron band depth, Rationale: Gypsum.') 
    surfacegroup.add_argument('--BD1900',action='store_true', dest='bd1900', default=False, help='Parameter: 1.9 micron band depth, Rationale: H2O, chemically bound or absorbed.')
    surfacegroup.add_argument('--BD2100',action='store_true', dest='bd2100', default=False, help='Parameter: 2.1 micron band depth, Rationale: monohydrated minerals.')   
    surfacegroup.add_argument('--BD2210',action='store_true', dest='bd2210', default=False, help='Parameter: 2.21 micron band depth, Rationale: Al-OH minerals.') 
    surfacegroup.add_argument('--BD2290',action='store_true', dest='bd2290', default=False, help='Parameter: 2.9 micron band depth, Rationale: Mg,Fe-OH minerals (at 2.3); also CO2 ice(at 2.292  microns).')     
    surfacegroup.add_argument('--D2300',action='store_true', dest='d2300', default=False, help='Parameter: 2.3 micron drop, Rationale: hydrated minerals; particularly clays.') 
    surfacegroup.add_argument('--SINDEX',action='store_true', dest='sindex', default=False, help='Parameter: Convexity at 2.29 microns  due to absorptions at 1.9/2.1 microns and 2.4 microns, Rationale: hydrated minerals; particularly sulfates.')
    surfacegroup.add_argument('--ICER2', action='store_true', dest='icer2', default=False, help='Parameter:gauge 2.7micron band,Rationale: CO2 ice will be >> 1; H2O ice and soil will be ~ 1.')    
    surfacegroup.add_argument('--BDCARB', action='store_true', dest='bdcarb', default=False, help='Parameter:overtone band depth,Rationale: carbonate overtones.')    
    surfacegroup.add_argument('--BD3000', action='store_true', dest='bd3000', default=False, help='Parameter:3 micron band depth,Rationale: H2O, chemically bound or adsorbed.') 
    surfacegroup.add_argument('--BD3100', action='store_true', dest='bd3100', default=False, help='Parameter:3.1 micron band depth,Rationale: H2Oice.')
    surfacegroup.add_argument('--BD3200', action='store_true', dest='bd3200', default=False, help='Parameter:3.2 micron band depth,Rationale: CO2ice.')
    surfacegroup.add_argument('--BD3400', action='store_true', dest='bd3400', default=False, help='Parameter:3.4 micron band depth,Rationale: Carbonates;organics.')
    surfacegroup.add_argument('--CINDEX', action='store_true', dest='cindex', default=False, help='Parameter:gauge 3.9 micron band,Rationale: carbonates.')
    surfacegroup.add_argument('--RPeak1', action='store_true', dest='rpeak1', default=False, help='Parameter: reflectance peak 1')
    surfacegroup.add_argument('--BDI1000VIS', action='store_true', dest='bdi1000VIS', default=False, help='Parameter: 1 micron integrated band depth; VIS wavelengths')
    
    #The following will raise not implemented errors.
    surfacegroup.add_argument('--VAR', '--spectral_variance', action='store_true', dest='var', default=False, help='Parameter: pyroxene index, Rationale: pyroxene will be strongly positive; favors HCP.')
    surfacegroup.add_argument('--BDI1000IR', action='store_true', dest='bdi1000IR', default=False, help='Parameter: 1 micron integrated band depth; IR wavelengths')
    surfacegroup.add_argument('--BDI2000', action='store_true', dest='bdi2000', default=False, help='Parameter: 2 micron integrated band depth, Rationale: pyroxene abundance and particle size')

    return (parser.parse_args()) 

    

def main():
    args = parseargs()
    name = args.outputname
    
    #Get a list of available CDR6 lookup tables
    sw = os.listdir("SW")
    sw = [x for x in sw if not ".LBL" in x.split("_")[5]]
    
    if args.input_data.split('.')[1] == 'img':
        print "\nYou provided the .img file.  We will read from the .lbl file\n"
        args.input_data = args.input_data.split('.')[0] + '.lbl'

    try:
        input_data = string.split(args.input_data, '.')
        input_data = input_data[0] + '.lbl'
        dataset = GdalIO(input_data)
        raster = dataset.load()
    except:
        #Here we intercept the input .IMG, alter the .LBL so that GDAL can read it and point at the .LBL.
        input_data = string.split(args.input_data, '.')
        input_data = input_data[0] + '.lbl'
        newlbl = open(input_data.split('.')[0] + '_fixed.lbl', 'w')
        for line in open(input_data,'r'):
            #print line
            if "OBJECT          = FILE" in line:
                line = "/* OBJECT = FILE */\n"

                print "\nA new GDAL readable .lbl file has been created with _fixed appended to the file name.  If you are running this multiple times on the same image, we are overwriting the label each time.\n"
            if "LINES" in line:
                lines = int(line.split("=")[1])
            if "LINE_SAMPLES" in line:
                samples = int(line.split("=")[1])
            newlbl.write(line)
        newlbl.close()
        dataset = GdalIO(input_data.split('.')[0] + '_fixed.lbl') #We have to read off disk because GDAL wants a string and newlbl is type file
        raster = dataset.load()
    print lines, samples
    basename =  os.path.basename(input_data).split('.')[0]
    dirname =  os.path.dirname(input_data)
    hdr = open(dirname + '/' + basename +'.hdr', 'w')
    hdr_out = ('''ENVI
samples = {}
lines   = {}
bands   = 438
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bil
byte order = 0''').format(samples,lines)
    hdr.write(hdr_out)
    print "ENVI .hdr file created to allow GDAL to read the input image.\n"
    hdr.close()
    
    dataset = GdalIO(dirname + '/' + basename +'.img') #We have to aim GDAL at the .img.
    raster = dataset.load()
    xsize, ysize, bands, metadata = dataset.info(raster)

    #We are assuming that the NDV is homogeneous over the entire cube
    NDV = raster.GetRasterBand(1).GetNoDataValue()

    try:
        wavelengths,sensor_id = get_wavelengths(open(input_dataset, 'r'),sw)
    except:
        wavelengths,sensor_id = get_wavelengths(open(input_data.split('.')[0] + '_fixed.lbl','r'),sw)
    
    #Check to make sure that the sensor supports the algorithm.
    if sensor_id == 'S':
        print "Input image uses the IR sensor.  Supported algorithms with this dataset are:"
        print "['R770', 'RBR', 'BD530', 'SH600', 'BD640', 'BD860', 'BD920', 'RPEAK1', 'BDI1000VIS','R440', 'IRR1']"
    
    elif sensor_id =='L':
        print "Input image uses the VNIR sensor.  Supported algorithms with this dataset are:"
        print "['BD11000IR', 'IRA', 'OLINDEX', 'LCPINDEX', 'HCPINDEX', 'VAR', 'ISLOPE1', 'BD1435', 'BD1500', 'ICER1'," 
        print " 'BD1750', 'BD1900', 'BDI2000', 'BD2100', 'BD2210', 'BD2290', 'D2300', 'SINDEX', 'ICER2', 'BDCARB',"
        print " 'BD3000', 'BD3100', 'BD3200', 'BD3400', 'CINDEX', 'BD1270O2', 'BD1400H2O', 'BD2000CO2', 'BD2350', 'BD2600',"
        print " 'IRR2', 'R2700', 'BD2700', 'IRR3']\n"
        
    #Start parsing the args
    if args.rockdust1 == True:
        array_out = algorithms.rockdust1(raster, wavelengths)
        
        if args.outputname == None:
            name = 'R770.tif'
        
    if args.rockdust2 == True:
        array_out = algorithms.rockdust2(raster, wavelengths)
        
        if args.outputname == None:
            name = 'R770divR440.tif'
            
    if args.bd530 == True:
        array_out = algorithms.bd530(raster, wavelengths)
        if args.outputname == None:
            name = 'BD530.tif'
                    
    if args.sh600 == True:
        array_out = algorithms.sh600(raster, wavelengths)
        if args.outputname == None:
            name = 'SH600.tif'
            
    if args.bd640 == True:
        array_out = algorithms.bd640(raster, wavelengths)
        if args.outputname == None:
            name = 'BD640.tif'               

    if args.bd860 == True:
        array_out = algorithms.bd860(raster, wavelengths)
        if args.outputname == None:
            name = 'BD860.tif' 

    if args.bd920 == True:
        array_out = algorithms.bd920(raster, wavelengths)
        if args.outputname == None:
            name = 'BD920.tif' 
    
    if args.ira == True:
        array_out = algorithms.ira(raster, wavelengths)
        if args.outputname == None:
            name='IRAlbedo.tif'
            
    if args.olivine_index == True:
        array_out = algorithms.olivine_index(raster, wavelengths)
        if args.outputname == None:
            name = 'OlivineIndex.tif'

    if args.olivine_index2 == True:
        array_out = algorithms.olivine_index2(raster, wavelengths)
        if args.outputname == None:
            name = 'OlivineIndex2.tif'
            
    if args.lcp_index == True:
        array_out = algorithms.hcp_index(raster, wavelengths)
        if args.outputname == None:
            name = 'LCPIndex.tif'
        
    if args.hcp_index == True:
        array_out = algorithms.hcp_index(raster, wavelengths)
        if args.outputname == None:
            name = 'HCPIndex.tif'

    if args.var == True:
        array_out = algorithms.var(raster, wavelengths)
        if args.outputname == None:
            name = 'SpectralVariance.tif'

    if args.rpeak1 == True:
        array_out = algorithms.rpeak1(raster, wavelengths)
        if args.outputname == None:
            name = 'RPeak1.tif'

    if args.bdi1000VIS == True:
        array_out = algorithms.bdi1000VIS(raster, wavelengths)
        if args.outputname == None:
            name = 'BDI1000VIS.tif'        

    if args.bdi1000IR == True:
        array_out = algorithms.bdi1000IR(raster, wavelengths)
        if args.outputname == None:
            name = 'BDI1000IR.tif' 

    if args.islope1 == True:
        array_out = algorithms.islope1(raster, wavelengths)
        if args.outputname == None:
            name = 'ISlope1.tif' 

    if args.bd1435 == True:
        array_out = algorithms.bd1435(raster, wavelengths)
        if args.outputname == None:
            name = 'BD1435.tif' 

    if args.bd1500 == True:
        array_out = algorithms.bd1500(raster, wavelengths)
        if args.outputname == None:
            name = 'BD1500.tif'
    
    if args.icer1 == True:
        array_out = algorithms.icer1(raster, wavelengths)
        if args.outputname == None:
            name = 'ICER1.tif' 
            
    if args.bd1750 == True:
        array_out = algorithms.bd1750(raster, wavelengths)
        if args.outputname == None:
            name = 'BD1750.tif' 
            
    if args.bd1900 == True:
        array_out = algorithms.bd1900(raster, wavelengths)
        if args.outputname == None:
            name = 'BD1900.tif' 

    if args.bdi2000 == True:
        array_out = algorithms.bdi2000(raster, wavelengths)
        if args.outputname == None:
            name = 'BDI2000.tif'   
            
    if args.bd2100 == True:
        array_out = algorithms.bd2100(raster, wavelengths)
        if args.outputname == None:
            name = 'BD2100.tif'     

    if args.bd2210 == True:
        array_out = algorithms.bd2210(raster, wavelengths)
        if args.outputname == None:
            name = 'BD2210.tif'       
            
    if args.bd2290 == True:
        array_out = algorithms.bd2290(raster, wavelengths)
        if args.outputname == None:
            name = 'BD2290.tif'   

    if args.d2300 == True:
        array_out = algorithms.d2300(raster, wavelengths)
        if args.outputname == None:
            name = 'D2300.tif'   

    if args.sindex == True:
        array_out = algorithms.sindex(raster, wavelengths)
        if args.outputname == None:
            name = 'SIndex.tif'    
  
    if args.icer2 == True:
        array_out = algorithms.icer2(raster, wavelengths)
        if args.outputname == None:
            name = 'ICER2.tif'   
            
    if args.bdcarb == True:
        array_out = algorithms.bdcarb(raster, wavelengths)
        if args.outputname == None:
            name = 'BDCARB.tif' 
            
    if args.bd3000 == True:
        array_out = algorithms.bd3000(raster, wavelengths)
        if args.outputname == None:
            name = 'BD3000.tif'     
            
    if args.bd3100 == True:
        array_out = algorithms.bd3100(raster, wavelengths)
        if args.outputname == None:
            name = 'BD3100.tif'  
    
    if args.bd3200 == True:
        array_out = algorithms.bd3200(raster, wavelengths)
        if args.outputname == None:
            name = 'BD3200.tif'  
            
    if args.bd3400 == True:
        array_out = algorithms.bd3400(raster, wavelengths)
        if args.outputname == None:
            name = 'BD3400.tif'  
            
    if args.cindex == True:
        array_out = algorithms.cindex(raster, wavelengths)
        if args.outputname == None:
            name = 'CIndex.tif'        

    #Scale the data if necessary
    if args.scale ==True:
        scalemin = float(args.scalerange[0])
        scalemax = float(args.scalerange[1])
        array_out = scale(array_out, scalemin, scalemax) #Scale the data
        array_out = linear(array_out, scalemin, scalemax) #Stretch the data
    else:
        pass

    #Check to see if the array_out has been masked and if so, fill the mask
    if isinstance(array_out,numpy.ma.core.MaskedArray) == True:
        if args.newNDV == None:
            array_out = numpy.ma.filled(array_out, NDV)
        else:
            array_out = numpy.ma.filled(array_out, float(args.newNDV[0]))
    else:
        pass
    
    #Create the output
    if args.dtype == None:
        dtype = gdal.GetDataTypeName(raster.GetRasterBand(1).DataType)
    else:
        dtype = args.dtype
    output = dataset.create_output(xsize, ysize, name, 1, parseddtype(dtype))    
    output.GetRasterBand(1).WriteArray(array_out)
    
    if NDV or args.newNDV:
        if args.newNDV == None:
            output.GetRasterBand(1).SetNoDataValue(NDV)
        else:
            output.GetRasterBand(1).SetNoDataValue(float(args.newNDV[0]))
            
if __name__ == '__main__':
    gdal.SetCacheMax(805306368) #768MB
    gdal.UseExceptions() #This allows try / except statements with gdal.Open to fail in a pythonic way.
    #gdal.SetConfigOption('CPL_DEBUG', 'ON')
    
    main()
