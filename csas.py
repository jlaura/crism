#!/usr/bin/env python

#Internal imports
import sys
import argparse
import string

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

def wv2band(wvlength): 
    '''This function opens the wavelength tab file and parses the band
    lookups.  A list is returned where the wavelength index is equal to 
    the band number -1.  This is because python lists are 0 based and GDAL's
    band reading is 1 based.
    '''
    wv = open(wvlength)
    wv2bnd = []
    for line in wv:
        linelst = string.split(line, ',')
        wv2bnd.append(float(linelst[2]))
    wv.close
    return wv2bnd

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
    print array.mean(), array.max(), array.min(), array.std()
    array = (((scalemax - scalemin)*(array-arraymin)) / (arraymax-arraymin)) + scalemin
    return array

def linear(array,scalemin,scalemax):
    '''
    This function exists to perform a linear stretch on the data to 2std from the mean.
    '''
    scalemax = array.mean() + (2* array.std())
    scalemin = array.mean() - (2* array.std())
    array = (array - array.min()) * ((scalemax - scalemin)/ (array.max()-array.min())) + scalemin 
    print array.mean(), array.max(), array.min(), array.std()
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
    name = args.outputname
    
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
                line = "/* OBJECT = FILE */"

                print "\nA new GDAL readable .lbl file has been created with _fixed appended to the file name.  If you are running this multiple times on the same image, we are overwriting the label each time.\n"
            newlbl.write(line)
        newlbl.close()

        dataset = GdalIO(input_data.split('.')[0] + '_fixed.lbl') #We have to read off disk because GDAL wants a string and newlbl is type file
        raster = dataset.load()
        
    xsize, ysize, bands, metadata = dataset.info(raster)

    #We are assuming that the NDV is homogeneous over the entire cube
    NDV = raster.GetRasterBand(1).GetNoDataValue()
    
    #Generate a wavelength to band number lookup table via either the provided lookup or the image metadata.
    if args.wavelength_tab == None:
        #Attempt to parse the wv.tab file that is shipped with an MRDR
        wavelength_lookup = args.input_data.replace('mrrif', 'mrrwv').split('.')
        wavelengths = wv2band(wavelength_lookup[0]+'.tab')
    print wavelengths
    #Start parsing the args
    if args.rockdust1 == True:
        array_out = algorithms.rockdust1(raster, wavelengths)
        
        if args.outputname == None:
            name = 'R770.tif'
        
    if args.rockdust2 == True:
        array_out = algorithms.rockdust2(raster, wavelengths)
        
        if args.outputname == None:
            name = 'R770divR440.tif'
        
    if args.ira == True:
        array_out = algorithms.ira(raster, wavelengths)
        
        if args.outputname == None:
            name='IRAlbedo.tif'
    
    if args.olivine_index == True:
        array_out = algorithms.olivine_index(raster, wavelengths)
        
        if args.outputname == None:
            name = 'OlivineIndex.tif'
        
    if args.lcp_index == True:
        array_out = algorithms.hcp_index(raster, wavelengths)
        if args.outputname == None:
            name = 'LCPIndex.tif'
        
    if args.hcp_index == True:
        array_out = algorithms.hcp_index(raster, wavelengths)
        if args.outputname == None:
            name = 'HCPIndex.tif'
            
    if args.ferric_coating == True:
        array_out = algorithms.ferric_coating(raster, wavelengths)
        if args.outputname == None:
            name = 'FerricCoating.tif'
            
    if args.h2o_co2_ice == True:
        array_out = algorithms.h2o_co2_ice(raster, wavelengths)
        if args.outputname == None:
            name = 'h2o_co2_ice.tif'

    if args.h2o_co2_ice2 == True:
        array_out = algorithms.h2o_co2_ice2(raster, wavelengths)
        if args.outputname == None:
            name = 'h2o_co2_ice2.tif'
            
    if args.h2o1 == True:
        array_out = algorithms.h2o1(raster, wavelengths)
        if args.outputname == None:
            name = 'h2o1.tif'
            
    if args.carbonates == True:
        array_out = algorithms.carbonates(raster, wavelengths)
        if args.outputname == None:
            name = 'Carbonates.tif'
    
    
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
    
    if args.newNDV == None:
        output.GetRasterBand(1).SetNoDataValue(NDV)
    else:
        output.GetRasterBand(1).SetNoDataValue(float(args.newNDV[0]))
            
if __name__ == '__main__':
    gdal.SetCacheMax(805306368) #768MB
    gdal.UseExceptions() #This allows try / except statements with gdal.Open to fail in a pythonic way.
    #gdal.SetConfigOption('CPL_DEBUG', 'ON')
    
    main()
