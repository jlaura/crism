#Internal imports
import sys
import argparse
import string

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
        return xsize, ysize, bands
    
    def create_output(self,xsize,ysize,name,bands=1,dtype=gdal.GDT_Float32):

        """Method to create an output of the same type, size, projection, and transformation as the input dataset"""
        driver = gdal.GetDriverByName('GTIFF')
        outdataset = driver.Create(name,xsize,ysize,bands,dtype)
        #outdataset = driver.Create(outputname, xsize, ysize)
        
        return outdataset    
    
def ferric_coating(raster):
    '''
    Name: ISPLOE1
    Parameter: -1 * Spectral Slope_1
    Formulation: (R1815-R2530)/(R2530-R1815)
    Rationale: ferric coating on dark rock
    Band(s): 321, 213
    '''
            
    array1815 = getarray(raster.GetRasterBand(321))
    array2530 = getarray(raster.GetRasterBand(213))

    #Algorithm
    array_out = (array1815-array2530)/(array2530-array1815)
    
    return array_out

def carbonates(raster):
    '''
    Name: CINDEX
    Parameter: gauge 3.9micron band
    Formulation: (R3750 + (R3750-R3630)/((R3750-R3630)*(R3950-R3750)))/R3950-1
    Rationale: Carbonates.
    Band(s):29, 47, 1, !!!1 is empty - remapped to 5 to get output!!!!!
    '''
    
    R3750 = getarray(raster.GetRasterBand(29))
    R3630 = getarray(raster.GetRasterBand(47))
    R3950 = getarray(raster.GetRasterBand(5))

    #Algorithm
    array_out = (R3750 + (R3750-R3630)/((R3750-R3630)*(R3950-R3750)))/R3950-1
    return array_out
    
def hcp_index(raster):
        '''
        Name: HCPINDEX
        Parameter: pyroxene index
        Formulation: ((R1470-R1050)/(R1470+R1050))*((R1470-R1050)/(R1470+R2067))
        Rationale: pyroxene will be strongly positive;favors HCP
        Band(s): 374, 438 !!!RESET TO 436, the last band with valid data!!!, 253
        '''
            
        array1470 = getarray(raster.GetRasterBand(374))
        array1050 = getarray(raster.GetRasterBand(436))
        array2067 = getarray(raster.GetRasterBand(253))
        
        #Algorithm
        array_out = ((array1470-array1050)/(array1470+array1050))*((array1470-array1050)/(array1470+array2067))
        
        return array_out
    
def h2o1(raster):
    
    '''
    Name: BD3000
    Parameter: 3micron band depth
    Formulation: 1 - (R3000/(R2530*(R2530/R2210)))
    Rationale: h2o.
    Band(s): 142, 213, 262
    '''
        
    R3000 = getarray(raster.GetRasterBand(142))
    R2530 = getarray(raster.GetRasterBand(213))
    R2210 = getarray(raster.GetRasterBand(262))

    #Algorithm
    array_out = 1 - (R3000/(R2530*(R2530/R2210)))
    
    return array_out    
    
def h2o_co2_ice(raster):
    
    '''
    Name: ICER1
    Parameter: 1.5micron and 1.43micron band ration
    Formulation: R1510 / R1430
    Rationale: CO2, H2O Ice Mixtures.
    Band(s): 368, 380
    '''
        
    array1510 = getarray(raster.GetRasterBand(368))
    array1430 = getarray(raster.GetRasterBand(380))

    #Algorithm
    array_out = array1510/array1430
    
    return array_out

def h2o_co2_ice2(raster):
    
    '''
    Name: ICER2
    Parameter: gauge 2.7 micron band
    Formulation: R2530 / R2600
    Rationale: CO2 ice will be >> 1; H2O ice and soil will be ~ 1.
    Band(s): 213, 203
    '''
        
    array2530 = getarray(raster.GetRasterBand(213))
    array2600 = getarray(raster.GetRasterBand(203))

    #Algorithm
    array_out = array2530/array2600
    
    return array_out

def ira(raster):
    '''
    Name:IRA
    Parameter: 1.3micron reflectance
    Formulation: R1330
    Rationale: IR albedo
    Band(s): 395
    '''
        
    array1330 = getarray(raster.GetRasterBand(395))
    
    return array1330

def lcp_index(raster):
        '''
        Name: LCPINDEX
        Parameter: pyroxene index
        Formulation: ((R1330-R1050)/(R1330+R1050))*((R1330-R1815)/(R1330+R1815))
        Rationale: pyroxene will be strongly positive;favors LCP
        Band(s): 395, 438 !!!RESET TO 436, the last band with valid data!!!, 321, 
        '''
            
        array1330 = getarray(raster.GetRasterBand(395))
        array1050 = getarray(raster.GetRasterBand(436))
        array1815 = getarray(raster.GetRasterBand(321))
        
        #Algorithm
        array_out = ((array1330-array1050)/(array1330+array1050))*((array1330-array1815)/(array1330+array1815))
        
        return array_out
    
def olivine_index(raster):
    '''
    Name: OLINDEX
    Parameter: Olivine Index
    Formulation: (R1695/(0.1*R1050 + 0.1*R1210 + 0.4*R1330 + 0.4*R1470))-1
    Rationale: Olivine will be strongly positive; based on fayalite
    Band(s): 340, 438 !!!438 remapped to 436 for valid values!!!, 414, 395, 374
    '''
    
    NDV = raster.GetRasterBand(1).GetNoDataValue()
        
    array1695 = getarray(raster.GetRasterBand(340))
    array1050 = getarray(raster.GetRasterBand(436))
    array1210 = getarray(raster.GetRasterBand(414))
    array1330 = getarray(raster.GetRasterBand(395))
    array1470 = getarray(raster.GetRasterBand(374))
    
    #Algorithm
    array_out = (array1695/((0.1*array1050)+(0.1*array1210)+(0.4*array1330)+(0.4*array1470)))-1
    
    return array_out

    
def rockdust1(raster):
        '''
        Name: R770
        Parameter: 0.77micron reflectance
        Formulation: R770
        Rationale: rock/dust
        Band(s): 248
        '''
            
        array770 = getarray(raster.GetRasterBand(247+1))
        
        return array770


def rockdust2(raster):
    '''
    Name: RBR
    Parameter: Red/Blue Ratio
    Formulation: R770 / R440
    Rationale: rock/dust
    Bands: 248, 196
    '''
        
    band770 = getarray(raster.GetRasterBand(247+1))
    band440 = getarray(raster.GetRasterBand(194+1))
    
    #Algorithm
    array_out = band770 / band440
    
    return array_out

    
def getarray(band):
    #Get NDV, GetBand, MaskNDV
    NDV = band.GetNoDataValue()
    array = band.ReadAsArray().astype(numpy.float64)
    array = numpy.ma.masked_where(array == NDV, array, copy=False)
    return array

def wv2band(wvlength): 
    wv = open(wvlength)
    wv2bnd = {}
    for line in wv:
        linelst = string.split(line, ',')
        wv2bnd[float(linelst[1].rstrip())] = int(linelst[0])
    wv.close
    return wv2bnd

def scale(array, scalemin, scalemax):
    '''
    This function is used to scale the data between a user defined range [c,d].  By default this range is between 1 and 255.  This maintains 0 as a special no data value should the user wish to set it.
    
    Scale from (a,b) to c,d)
        y = mx+b
        x = ((x-a)(d-c)/(b-a))+c
        
    Returns a scaled ndarray
    '''
    arraymin = numpy.amin(array)
    arraymax = numpy.amax(array)
    array -= arraymin
    array *= (scalemax-scalemin)
    array *= 1.0/(arraymax-arraymin)
    array += scalemin
    
    return array

def parseargs():
#Setup the arg parser
    parser = argparse.ArgumentParser(description='CRISM Spectral Analysis Tool')
    parser.add_argument('input_cube', action='store', help='The input CRISM data cube to be processed.')
    parser.add_argument('outputname', action='store', nargs='?',default=None, help='The output file name. Defaults are provided if this is not specified.')
    parser.add_argument('wavelength_tab', action='store', help='The wavelength lookup table for the input_cube.  Values are currently hard coded.')
    
    #Helper Parameters
    helpergroup = parser.add_argument_group('Helper Functionality')
    
    helpergroup.add_argument('--scale', '-s', action='store_true', dest='scale', default=False, help='Scale the data.  Typically from 1-255. NDV is set to 0.')
    
    #Surface Parameters
    surfacegroup = parser.add_argument_group('Surface Parameters')
    
    surfacegroup.add_argument('--rockdust1','--R770', action='store_true', dest='rockdust1', default=False, help='Parameter: 0.77micron band, Rationale: rock/dust.')
    surfacegroup.add_argument('--rockdust2', action='store_true', dest='rockdust2', default=False, help='Parameter: red/blue ratio, Rationale: rock.dust')
    surfacegroup.add_argument('--IRAlbedo', '--IRA', action='store_true', dest='ira', default=False, help='Parameter: 1.3micron reflectance, Rationale: IR albedo.')
    surfacegroup.add_argument('--OlIndex', '--olivine_index', action='store_true', dest='olivine_index', default=False, help='Parameter:Olivine Index, Rationale: Olivine will be strongly positive; based on fayalite.')
    surfacegroup.add_argument('--LCP', '--lcp_index', action='store_true', dest='lcp_index', default=False, help='Parameter: pyroxene index, Rationale:pyroxene will be strongly positive; favors LCP')
    surfacegroup.add_argument('--HCP', '--hcp_index', action='store_true', dest='hcp_index', default=False, help='Parameter: pyroxene index, Rationale: pyroxene will be strongly positive; favors HCP.')
    surfacegroup.add_argument('--ISLOPE1', '--ferric1', action='store_true', dest='ferric_coating', default=False, help='Parameter: -1*Spectral Slope_1, Rationale: Ferric coating on dark rocks.')
    surfacegroup.add_argument('--ICER1', '--h2o_co2_ice', action='store_true', dest='h2o_co2_ice', default=False, help='Parameter:1.5micron and 1.43micron band ratio,Rationale: CO2, H2O Ice Mixtures.')
    surfacegroup.add_argument('--ICER2', '--h2o_co2_ice2', action='store_true', dest='h2o_co2_ice2', default=False, help='Parameter:guge 2.7micron band,Rationale: CO2 ice will be >> 1; H2O ice and soil will be ~ 1.')
    surfacegroup.add_argument('--BD3000', '--h2o1', action='store_true', dest='h2o1', default=False, help='Parameter: 3micron band depth, Rationale: h2o.')
    surfacegroup.add_argument('--CINDEX', '--carbonates', action='store_true', dest='carbonates', default=False, help='Parameter: gauge 3.9micron band, Rationale:Carbonates.')
    return (parser.parse_args())

def main():
    args = parseargs()
    #Generate a wavelength to band number lookup table
    wavelengths = wv2band(args.wavelength_tab)
    
    #Open the input dataset & get info
    dataset = GdalIO(args.input_cube)
    raster = dataset.load()
    xsize, ysize, bands = dataset.info(raster)
    NDV = raster.GetRasterBand(1).GetNoDataValue()

    if args.rockdust1 == True:
        array_out = rockdust1(raster)
        
        if args.outputname == None:
            name = 'R770.tif'
        
    if args.rockdust2 == True:
        array_out = rockdust2(raster)
        
        if args.outputname == None:
            name = 'R770divR440.tif'
        
    if args.ira == True:
        array_out = ira(raster)
        
        if args.outputname == None:
            name='IRAlbedo.tif'
    
    if args.olivine_index == True:
        array_out = olivine_index(raster)
        
        if args.outputname == None:
            name = 'OlivineIndex.tif'
        
    if args.lcp_index == True:
        array_out = hcp_index(raster)
        if args.outputname == None:
            name = 'LCPIndex.tif'
        
    if args.hcp_index == True:
        array_out = hcp_index(raster)
        if args.outputname == None:
            name = 'HCPIndex.tif'
            
    if args.ferric_coating == True:
        array_out = ferric_coating(raster)
        if args.outputname == None:
            name = 'FerricCoating.tif'
            
    if args.h2o_co2_ice == True:
        array_out = h2o_co2_ice(raster)
        if args.outputname == None:
            name = 'h2o_co2_ice.tif'

    if args.h2o_co2_ice2 == True:
        array_out = h2o_co2_ice2(raster)
        if args.outputname == None:
            name = 'h2o_co2_ice2.tif'
            
    if args.h2o1 == True:
        array_out = h2o1(raster)
        if args.outputname == None:
            name = 'h2o1.tif'
            
    if args.carbonates == True:
        array_out = carbonates(raster)
        if args.outputname == None:
            name = 'Carbonates.tif'
            
    #Scale the data if necessary
    if args.scale ==True:
        array_out = scale(array_out, 1, 255)
        array_out = numpy.ma.filled(array_out, 0)
        output = dataset.create_output(xsize, ysize, name, 1, gdal.GDT_Byte)
    else:
        array_out = numpy.ma.filled(array_out, NDV)
        output = dataset.create_output(xsize, ysize, name)

    output.GetRasterBand(1).WriteArray(array_out)

    #Setup the output NDV
    if args.scale == True:
        output.GetRasterBand(1).SetNoDataValue(0)
    else:
        output.GetRasterBand(1).SetNoDataValue(NDV)
    
if __name__ == '__main__':
    main()
