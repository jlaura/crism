#!/usr/bin/env python

import csas
from osgeo import gdal
def rockdust1(raster, wavelengths):
    '''
    Name: R770
    Parameter: 0.77micron reflectance
    Formulation: R770
    Rationale: rock/dust
    '''
    bands = csas.getbandnumbers(wavelengths, 770)     
    print raster, raster.GetRasterBand(247).ReadRaster(0,0,320,450,1,gdal.GDT_Float64)
    array770 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    
    return array770


def rockdust2(raster, wavelengths):
    '''
    Name: RBR
    Parameter: Red/Blue Ratio
    Formulation: R770 / R440
    Rationale: rock/dust
    '''
    bands = csas.getbandnumbers(wavelengths, 440,770)  
    band440 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band770 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    #Algorithm
    array_out = band770 / band440
    
    return array_out

def bd530(raster, wavelengths):
    '''
    NAME: BD530
    PARAMETER: 0.53 micron band depth
    FORMULATION *: 1 - (R530/(a*R709+b*R440))
    RATIONALE: Crystalline ferric minerals
    '''
    bands = csas.getbandnumbers(wavelengths, 440, 530, 709)
    band440 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band530 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band709 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    wv440 = wavelengths[bands[0]]
    wv530 = wavelengths[bands[1]]
    wv709 = wavelengths[bands[2]]
    #Algorithm
    a = (wv530 - wv440) / (wv709 - wv440)
    b = 1.0 - a
    array_out = 1.0 - (band530/((a*band709)+(b*band440)))
    return array_out

def sh600(raster, wavelengths):
    '''
    NAME: SH600
    PARAMETER: 0.60 micron shoulder height
    FORMULATION *: R600/(a*R530+b*R709)
    RATIONALE: select ferric minerals
    '''
    bands = csas.getbandnumbers(wavelengths, 533, 600, 710)
    band533 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band600 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band710 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    wv533 = wavelengths[bands[0]]
    wv600 = wavelengths[bands[1]]
    wv710 = wavelengths[bands[2]]
    #Algorithm
    a = (wv600 - wv533) / (wv710 - wv533)
    b = 1.0 - a
    array_out = 1.0 - (((b * band533)+(a*band710))/band600)
    return array_out

def bd640(raster, wavelengths):
    '''
    NAME: BD640
    PARAMETER: 0.64 micron band depth
    FORMULATION *: 1 - (R648/(a*R600+b*R709))
    RATIONALE: select ferric minerals, especially maghemite
    '''
    bands = csas.getbandnumbers(wavelengths, 600, 648, 709)
    band600 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band648 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band709 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    wv600 = wavelengths[bands[0]]
    wv648 = wavelengths[bands[1]]
    wv709 = wavelengths[bands[2]]
    #Algorithm
    a = (wv648 - wv600) / (wv709 - wv600)
    b = 1.0 - a
    array_out = 1.0 - (band648/((b*band600)+(a*band709)))
    return array_out

def bd860(raster, wavelengths):
    '''
    NAME: BD860
    PARAMETER: 0.86 micron band depth
    FORMULATION *: 1 - (R860/(a*R800+b*R984))
    RATIONALE: select ferric minerals ('hematite band')
    '''
    bands = csas.getbandnumbers(wavelengths, 800, 860, 984)
    band800 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band860 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band984 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    wv800 = wavelengths[bands[0]]
    wv860 = wavelengths[bands[1]]
    wv984 = wavelengths[bands[2]]
    #Algorithm
    a = (wv860 - wv800) / (wv984 - wv800)
    b = 1.0 - a
    array_out = 1.0 - (band860/((b*band800)+(a*band984)))
    return array_out

def bd920(raster, wavelengths):
    '''
    NAME: BD920
    PARAMETER: 0.92 micron band depth
    FORMULATION *: 1 - ( R920 / (a*R800+b*R984) )
    RATIONALE: select ferric minerals ('Pseudo BDI1000 VIS')
    '''
    bands = csas.getbandnumbers(wavelengths, 800, 920, 984)
    band800 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band920 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band984 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    wv800 = wavelengths[bands[0]]
    wv920 = wavelengths[bands[1]]
    wv984 = wavelengths[bands[2]]
    #Algorithm
    a = (wv920 - wv800) / (wv984 - wv800)
    b = 1.0 - a
    array_out = 1.0 - (band920/((b*band800)+(a*band984)))
    return array_out

def rpeak1(raster, wavelengths):
    import numpy as np
    import sys
    import time
    import multiprocessing as mp
    starttime = time.clock()
    '''
    NAME: RPEAK1
    PARAMETER: reflectance peak 1
    FORMULATION *: wavelength where 1st derivative=0 of 5th order
      polynomial fit to R600, R648, R680, R710, R740, R770, R800, R830
    RATIONALE: Fe mineralogy
    '''
    bands = csas.getbandnumbers(wavelengths, 442,533,600,710,740,775,800,833,860,893,925)
    band442 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band533 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band600 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    band710 = csas.getarray(raster.GetRasterBand(bands[3]+1))
    band740 = csas.getarray(raster.GetRasterBand(bands[4]+1))
    band775 = csas.getarray(raster.GetRasterBand(bands[5]+1))
    band800 = csas.getarray(raster.GetRasterBand(bands[6]+1))
    band833 = csas.getarray(raster.GetRasterBand(bands[7]+1))
    band860 = csas.getarray(raster.GetRasterBand(bands[8]+1))
    band893 = csas.getarray(raster.GetRasterBand(bands[9]+1))
    band925 = csas.getarray(raster.GetRasterBand(bands[10]+1))
    
    wavelength_index = np.array([442,533,600,710,740,775,800,833,860,893,925])
    wavelength_vector = np.arange(len(wavelength_index), dtype=np.float64)
    for index,band in enumerate(bands):
        wavelength_vector[index] = wavelengths[band]
    num_model_points = 1001
    poly_degree = 4.0 #4th order polynomial fit
    model_wv_vector = np.arange(num_model_points, dtype=np.float64) / (num_model_points - 1.0)*((wavelength_vector.max() - wavelength_vector.min())+wavelength_vector.min())

    #Create a multi dimensional data cube from the selected bands 
    rpeak_cube = np.dstack((band442,band533,band600,band710,band740,band775, band800, band833, band860, band893, band925))
    
    #Create output arrays.
    rpeak_wavelength = np.ones(shape=(band442.shape),dtype=np.float32)
    rpeak_value = np.ones(shape=(band442.shape),dtype=np.float32)
    
    #Now we need to iterate over each pixel in the band depth dimension (11)
    for x in range(0,rpeak_cube.shape[1]):
        sys.stdout.write("Processed column %i of %i \r" %(x,rpeak_cube.shape[1]))
        sys.stdout.flush()
        
        for y in range(0,rpeak_cube.shape[0]):
            spec_vec = rpeak_cube[y][x]
            rpeak_params = np.polyfit(wavelength_vector,spec_vec,poly_degree)
            #model_spec = np.polyval(model_wv_vector, rpeak_params)
            model_spec = np.zeros(num_model_points)
            for m in range(0,int(poly_degree)+1):
                model_spec = model_spec + rpeak_params[m] * (model_wv_vector**float(m))
            model_spec_max = model_spec.max()
            model_spec_max_index = np.argmax(model_spec)            
            model_spec_max_wv = model_wv_vector[model_spec_max_index]
            rpeak_wavelength[y][x] =  model_spec_max_wv
            print rpeak_wavelength[y][x]
            rpeak_value[y][x] = model_spec_max
    array_out = rpeak_wavelength / 1000.0
    
    stoptime = time.clock()
    print "Total time to perform Rpeak1 %s" %(stoptime-starttime)
    return array_out
    
def bdi1000VIS(raster, wavelengths):
    import numpy as np
    from scipy import integrate
    '''
    NAME: BDI1000VIS
    PARAMETER: 1 micron integrated band depth; VIS wavelengths
    FORMULATION *: divide R830, R860, R890, R915 by RPEAK1 then
      integrate over (1 -  normalized radiances)
    RATIONALE: crystalline Fe+2 or Fe+3 minerals
    '''
    bands = csas.getbandnumbers(wavelengths, 833, 860, 892, 925, 951, 984, 1023)
    band833= csas.getarray(raster.GetRasterBand(bands[0]+1))
    band860 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band892 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    band925 = csas.getarray(raster.GetRasterBand(bands[3]+1))
    band951 = csas.getarray(raster.GetRasterBand(bands[4]+1))
    band984 = csas.getarray(raster.GetRasterBand(bands[5]+1))
    band1023 = csas.getarray(raster.GetRasterBand(bands[6]+1))
    
    wavelength_index = np.array([833, 860, 892, 925, 951, 984, 1023])
    wavelength_vector = np.arange(len(wavelength_index), dtype=np.float64)
    for index,band in enumerate(bands):
        wavelength_vector[index] = wavelengths[band]    
    
    wavelength_vector_um = wavelength_vector / 1000.0
    bdivis_value = np.zeros(shape=(band984.shape),dtype=np.float64)
    bdi1000_cube = np.dstack((band883,band860,band892,band925,band951,band984,band1023))
    rpeak_value_cube = np.zeros(shape=(band984.shape),dtype=np.float64)
    print "Computing the rpeak value for each cell."
    for x in range(0,rpeak_cube.shape[1]):
            sys.stdout.write("Processed column %i of %i \r" %(x,rpeak_cube.shape[1]))
            sys.stdout.flush()
            for y in range(0,rpeak_cube.shape[0]):
                spec_vec = rpeak_cube[x][y]
                rpeak_params = np.polyfit(wavelength_vector,spec_vec,poly_degree)
                #model_spec = np.polyval(model_wv_vector, rpeak_params)
                model_spec = np.zeros(num_model_points)
                for m in range(0,int(poly_degree)+1):
                    model_spec = model_spec + rpeak_params[m] * (model_wv_vector**float(m))
                model_spec_max = model_spec.max()
                model_spec_max_index = np.argmax(model_spec)            
                model_spec_max_wv = model_wv_vector[model_spec_max_index]
                rpeak_wavelength[x][y] = model_spec_max_wv
                rpeak_value_cube[x][y] = model_spec_max 
    print "Finished computing rpeak values. Now computing Integrated Band Depth."
    bdi1000_normalized_cube = bdi1000_cube / rpeak_value_cube
    for x in range(0,rpeak_cube.shape[1]+1):
                sys.stdout.write("Processed column %i of %i \r" %(x,rpeak_cube.shape[1]))
                sys.stdout.flush()
                for y in range(0,rpeak_cube.shape[0]+1):
                    spec_vec = bdi1000_normalized_cube[x][y]  
                    check_vec = bdi1000_cube[x][y]
                    bdivis_value[x][y] = scipy.integrate.newton_cotes(wavelength_vecor_um, (1.0-spec_vec))
                    
    return bdvis_value
    #raise NotImplementedError

def bdi1000IR(raster, wavelengths):
    '''
     NAME: BDI1000IR
     PARAMETER: 1 micron integrated band depth; IR wavelengths
     FORMULATION *: divide R1030, R1050, R1080, R1150
       by linear fit from peak R  between 1.3 - 1.87 microns to R2530
       extrapolated backwards, then integrate over (1 -  normalized
       radiances)
     RATIONALE: crystalline Fe+2 minerals; corrected for overlying
       aerosol induced slope
    '''

    raise NotImplementedError

def ira(raster, wavelengths):
    '''
    NAME: IRA
    PARAMETER: 1.3 micron reflectance
    FORMULATION *: R1330
    RATIONALE: IR albedo
    '''
    bands = csas.getbandnumbers(wavelengths, 1330)    
    array1330 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    
    return array1330

def olivine_index(raster, wavelengths):
    '''
    NAME: OLINDEX (prior to TRDR version 3)
    PARAMETER: olivine index
    FORMULATION *: (R1695 / (0.1*R1080 + 0.1*R1210 + 0.4*R1330 +
      0.4*R1470)) - 1
    RATIONALE: olivine will be strongly +; based on fayalite
    '''
    bands = csas.getbandnumbers(wavelengths, 1080,1210, 1470, 1695)  
    band1080 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band1210 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band1470 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    band1695 = csas.getarray(raster.GetRasterBand(bands[3]+1))
    #Algorithm
    array_out = (band1695 / (0.1*band1080 + 0.1*band1210 + 0.4*band1330 +
      0.4*band1470)) - 1
    
    return array_out

def olivine_index2(raster, wavelengths):
    '''
    NAME: OLINDEX2 (beginning with TRDR version 3)
    PARAMETER: olivine index with less sensitivity to illumination
    FORMULATION *: (((RC1054 ? R1054)/RC1054) * 0.1)
      + (((RC1211 ? R1211)/(RC1211) * 0.1)
      + (((RC1329 ? R1329)/RC1329) * 0.4)
      + (((RC1474 ? R1474)/RC1474) * 0.4)
    RATIONALE: olivine will be strongly positive
    '''
    print "Olivine Index 2 (for TRDR v3) has not been implemented in CAT yet."
    exit()

def hcp_index(raster, wavelengths):
    '''
    NAME: HCPXINDEX
    PARAMETER: pyroxene index
    FORMULATION *: 100 * ((R1470 - R1080)/(R1470 + R1080)) * ((R1470 - R2067)/(R1470+R2067))
    RATIONALE: pyroxene is strongly +; favors high-Ca pyroxene
    Algorithm differs from published - coded as per CAT
    '''
    bands = csas.getbandnumbers(wavelengths, 1080,1470, 2067)  
    band1080 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band1470 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band2067 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    
    #Algorithm
    array_out = 100 * ((band1470 - band1080)/(band1470 + band1080)) * ((band1470 - band2067)/(band1470+band2067))
    
    return array_out        

def lcp_index(raster, wavelengths):
    '''
    NAME: LCPINDEX
    PARAMETER: pyroxene index
    FORMULATION *: 100 * ((R1330 - R1080)/(R1330 + R1080)) * ((R1330 - R1815)/(R1330+R1815))
    RATIONALE: pyroxene is strongly +; favors low-Ca pyroxene
    Algorithm differs from published - coded as per CAT
    '''
    bands = csas.getbandnumbers(wavelengths, 1080,1330,1815)  
    band1080 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band1330 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band1815 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    
    #Algorithm
    array_out = 100 * ((band1330 - band1080)/(band1330 + band1080)) * ((band1330 - band1815)/(band1330+band1815))
    
    return array_out  

def var(raster, wavelengths):
    '''
     NAME: VAR
    PARAMETER: spectral variance
    FORMULATION *: find variance from a line fit from 1 - 2.3 micron
      by summing in quadrature over the intervening wavelengths
    RATIONALE: Ol & Px will have high values; Type 2 areas will have
      low values
    '''
    raise NotImplementedError

def islope1(raster, wavelengths):
    '''
    NAME: ISLOPE1
    PARAMETER: -1 * spectral slope1
    FORMULATION *: (R1815-R2530) / (2530-1815)
    RATIONALE: ferric coating on dark rock
    '''
    bands = csas.getbandnumbers(wavelengths, 1815,2530)  
    band1080 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band1330 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    wv1815 = wavelengths[bands[0]]
    wv2530 = wavelengths[bands[1]]
    
    #Algorithm
    array_out = (band1815-band2530)/(wv2530-wv1815)    
    return array_out  

def bd1435(raster, wavelengths):
    '''
    NAME: BD1435
    PARAMETER: 1.435 micron band depth
    FORMULATION *: 1 - ( R1430 / (a*R1370+b*R1470) )
    RATIONALE: CO2 surface ice
    '''
    bands = csas.getbandnumbers(wavelengths, 1370, 1430,1470)
    band1370 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band1430 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band1470 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    wv1370 = wavelengths[bands[0]]
    wv1430 = wavelengths[bands[1]]
    wv1470 = wavelengths[bands[2]]
    #Algorithm
    a = (wv1430 - wv1370) / (wv1470 - wv1370)
    b = 1.0 - a
    array_out = 1.0 - (band1430/((b*band1370)+(a*band1470)))
    return array_out

def bd1500(raster, wavelengths):
    '''
    NAME: BD1500
    PARAMETER: 1.5 micron band depth
    FORMULATION *: 1.0 - ((R1558 + R1505)/(R1808 + R1367))
    RATIONALE: H2O surface ice
    Algorithm differs from published - coded as per CAT (reduced instrument noise)
    '''
    bands = csas.getbandnumbers(wavelengths, 1367, 1505, 1558,1808)
    band1367 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band1505 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band1558 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    band1808 = csas.getarray(raster.GetRasterBand(bands[3]+1))
    #Algorithm
    array_out = 1.0 - ((band1558 + band1505)/(band1808 + band1367))
    return array_out

def icer1(raster, wavelengths):
    
    '''
    NAME: ICER1
    PARAMETER: 1.5 micron and 1.43 micron band ratio
    FORMULATION *: R1510 / R1430
    RATIONALE: CO2, H20 ice mixtures
    '''
    bands = csas.getbandnumbers(wavelengths, 1430, 1510)
    band1430 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band1510 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    #Algorithm
    array_out = band1430 / band1510
    return array_out

def bd1750(raster, wavelengths):
    '''
    NAME: BD1750
    PARAMETER: 1.75 micron band depth
    FORMULATION *: 1 - ( R1750 / (a*R1660+b*R1815) )
    RATIONALE: gypsum
    '''
    bands = csas.getbandnumbers(wavelengths, 1557, 1750, 1815)
    band1557 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band1750 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band1815 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    wv1557 = wavelengths[bands[0]]
    wv1750 = wavelengths[bands[1]]
    wv1815 = wavelengths[bands[2]]
    #Algorithm
    a = (wv1750 - wv1557) / (wv1815 - wv1557)
    b = 1.0 - a
    array_out = 1.0 - (band1750/((b*band1557)+(a*band1815)))
    return array_out

def bd1900(raster, wavelengths):
    '''
    NAME: BD1900
    PARAMETER: 1.9 micron band depth
    FORMULATION *: 1.0 - ((R1972 + R1927)/(R2006 + R1874))
    RATIONALE: H2O, chemically bound or adsorbed
    Algorithm differs from published - coded as per CAT (reduced instrument noise)
    '''
    bands = csas.getbandnumbers(wavelengths, 1874, 1927, 1973, 2006)
    band1874= csas.getarray(raster.GetRasterBand(bands[0]+1))
    band1927 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band1973 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    band2006 = csas.getarray(raster.GetRasterBand(bands[3]+1))
    #Algorithm
    array_out = 1.0 - ((band1972 + band1927)/(band2006 + band1874))
    return array_out

def bdi2000(raster, wavelengths):
    '''
    NAME: BDI2000
    PARAMETER: 2 micron integrated band depth
    FORMULATION *: divide R1660, R1815, R2140, R2210, R2250, R2290,
      R2330, R2350, R2390, R2430, R2460 by linear fit from peak R
      between 1.3 - 1.87 microns to R2530, then integrate over
     (1 -  normalized radiances)
    RATIONALE: pyroxene abundance and particle size
             
    '''
    raise NotImplementedError

def bd2100(raster, wavelengths):
    '''
    NAME: BD2100
    PARAMETER: 2.1 micron band depth
    FORMULATION *: 1 - ( ((R2120+R2140)*0.5) / (a*R1930+b*R2250) )
    RATIONALE: monohydrated minerals
    '''
    bands = csas.getbandnumbers(wavelengths, 1930, 2120, 2140, 2250)
    band1930 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band2120 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band2140 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    band2250 = csas.getarray(raster.GetRasterBand(bands[3]+1))
    wv1930 = wavelengths[bands[0]]
    wv2120 = wavelengths[bands[1]]
    wv2140 = wavelengths[bands[2]]
    wv2250 = wavelengths[bands[2]]
    #Algorithm
    a = (((wv2120 + wv2140) / 2) - wv1930) / (wv2250 - wv1930)
    b = 1.0 - a
    array_out = 1.0 - (((band2120 + band2140)*0.5)/((b*band1930)+(a*band2250)))
    return array_out

def bd2210(raster, wavelengths):
    '''
    NAME: BD2210
    PARAMETER: 2.21 micron band depth
    FORMULATION *: 1 - ( R2210 / (a*R2140+b*R2250) )
    RATIONALE: Al-OH minerals: monohydrated minerals
    '''
    bands = csas.getbandnumbers(wavelengths, 2140, 2210, 2250)
    band2140 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band2210 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band2250 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    wv2140 = wavelengths[bands[0]]
    wv2210 = wavelengths[bands[1]]
    wv2250 = wavelengths[bands[2]]
    #Algorithm
    a = (wv2210 - wv2140) / (wv2250 - wv2140)
    b = 1.0 - a
    array_out = 1.0 - ((band2210)/((b*band2140)+(a*band2250)))
    return array_out

def bd2290(raster, wavelengths):
    '''
    NAME: BD2290
    PARAMETER: 2.29 micron band depth
    FORMULATION *: 1 - ( R2290 / (a*R2250+b*R2350) )
    RATIONALE: Mg,Fe-OH minerals (at 2.3); also CO2 ice
      (at 2.292  microns)
    '''
    bands = csas.getbandnumbers(wavelengths, 2250, 2290, 2350)
    band2250 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band2290 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band2350 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    wv2250 = wavelengths[bands[0]]
    wv2290 = wavelengths[bands[1]]
    wv2350 = wavelengths[bands[2]]
    #Algorithm
    a = (wv2290 - wv2250) / (wv2350 - wv2250)
    b = 1.0 - a
    array_out = 1.0 - ((band2290)/((b*band2250)+(a*band2350)))
    return array_out

def d2300(raster, wavelengths):
    '''
     NAME: D2300
    PARAMETER: 2.3 micron drop
    FORMULATION *: 1 - ( (CR2290+CR2320+CR2330) /
      (CR2140+CR2170+CR2210) ) (CR values are observed R values
      divided by values fit along the slope as determined between 1.8
      and 2.53 microns - essentially continuum corrected))
    RATIONALE: hydrated minerals; particularly clays
    '''
    bands = csas.getbandnumbers(wavelengths, 1815, 2120, 2170, 2210, 2290, 2320,2330, 2530 )
    band1815 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band2120 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band2170 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    band2210 = csas.getarray(raster.GetRasterBand(bands[3]+1))
    band2290 = csas.getarray(raster.GetRasterBand(bands[4]+1))
    band2320 = csas.getarray(raster.GetRasterBand(bands[5]+1))
    band2330 = csas.getarray(raster.GetRasterBand(bands[6]+1))
    band2530 = csas.getarray(raster.GetRasterBand(bands[7]+1))
    wv1815 = wavelengths[bands[0]]
    wv2110 = wavelengths[bands[1]]
    wv2170 = wavelengths[bands[2]]
    wv2210 = wavelengths[bands[3]]
    wv2290 = wavelengths[bands[4]]
    wv2320 = wavelengths[bands[5]]
    wv2330 = wavelengths[bands[6]]
    wv2530 = wavelengths[bands[7]]
    #Algorithm
    #Continuum removal phase
    slope = (band2530-band1815) / (wv2530 - wv1815)
    cr2290 = band1815 + slope * (wv2290 - wv1815)
    cr2320 = band1815 + slope * (wv2320 - wv1815)
    cr2330 = band1815 + slope * (wv2330 - wv1815)
    cr2120 = band1815 + slope * (wv2120 - wv1815)
    cr2170 = band1815 + slope * (wv2170 - wv1815)
    cr2210 = band1815 + slope * (wv2210 - wv1815)
    #Computation phase
    array_out = 1.0 - (((band2290/cr2290)+(band2320/cr2320)+(band2330/cr2330))/((band2120/cr2120)+(band2170/cr2170)+(band2210/cr2210)))
    return array_out

def sindex(raster, wavelengths):
    '''
  NAME: SINDEX
    PARAMETER: Convexity at 2.29 microns  due to absorptions at
      1.9/2.1 microns and 2.4 microns
    FORMULATION *: 1 - (R2100 + R2400) / (2 * R2290) CR
      values are observed R values divided by values fit along the
      slope as determined between 1.8 - 2.53 microns (essentially
      continuum corrected))
    RATIONALE: hydrated minerals; particularly sulfates
    '''
    bands = csas.getbandnumbers(wavelengths, 2100, 2400, 2290)
    band2100 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band2400 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band2290 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    #Algorithm
    array_out = 1.0 - ((band2100 + band2400)/(2*band2290))
    return array_out

def icer2(raster, wavelengths):
    '''
    NAME: ICER2
    PARAMETER: gauge 2.7 micron band
    FORMULATION *: R2530 / R2600
    RATIONALE: CO2 ice will be >>1, H2O ice and soil will be about 1
    '''
    bands = csas.getbandnumbers(wavelengths, 2530, 2600)
    band2530 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band2600 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    #Algorithm
    array_out = band2530 / band2600
    return array_out    

def bdcarb(raster, wavelengths):
    from math import sqrt
    '''
    NAME: BDCARB
    PARAMETER: overtone band depth
    FORMULATION *: 1 - ( sqrt [ ( R2330 / (a*R2230+b*R2390) ) *
      ( R2530/(c*R2390+d*R2600) ) ] )
    RATIONALE: carbonate overtones        
    '''
    bands = csas.getbandnumbers(wavelengths, 2230, 2330, 2390, 2530, 2600)
    band2230 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band2330 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band2390 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    band2530 = csas.getarray(raster.GetRasterBand(bands[3]+1))
    band2600 = csas.getarray(raster.GetRasterBand(bands[4]+1))
    wv2230 = wavelengths[bands[0]]
    wv2330 = wavelengths[bands[1]]
    wv2390 = wavelengths[bands[2]]
    wv2530 = wavelengths[bands[3]]
    wv2600 = wavelengths[bands[4]]
    #Algorithm
    a = (((wv2330 + wv2120)*.5) - wv2230)/ (wv2390 - wv2230)
    b = 1.0 - a
    c = (((wv2530 + wv2120)*.5) - wv2390)/ (wv2600 - wv2390)
    d = 1.0 - c
    array_out = 1 - (sqrt(band2330 / ((b*band2230)+(a*band2390)))* (band2530 / ((d*band2230)+(c*band2600))))
    return array_out    

def bd3000(raster, wavelengths):

    '''
    NAME: BD3000
    PARAMETER: 3 micron band depth
    FORMULATION *: 1 - ( R3000 / (R2530*(R2530/R2210)) )
    RATIONALE: H2O, chemically bound or adsorbed    
    '''
    bands = csas.getbandnumbers(wavelengths, 2210, 2530, 3000)
    band2210 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band2530 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band3000 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    #Algorithm
    array_out = 1 - (band3000 / (band2530 * (band2530 / band2210)))
    return array_out    

def bd3100(raster, wavelengths):
    '''
    NAME: BD3100
    PARAMETER: 3.1 micron band depth
    FORMULATION *: 1 - ( R3120 / (a*R3000+b*R3250) )
    RATIONALE: H2O ice      
    '''
    bands = csas.getbandnumbers(wavelengths, 3000,3120, 3250)
    band3000 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band3120 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band3250 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    wv3000 = wavelengths[bands[0]]
    wv3120 = wavelengths[bands[1]]
    wv3250 = wavelengths[bands[2]]
    #Algorithm
    a = (wv3120 - wv3000)/ (wv3250 - wv3000)
    b = 1.0 - a
    array_out = 1.0 - (band3120/((b*band3000)+(a*band3250)))
    return array_out    

def bd3200(raster, wavelengths):
    '''
    NAME: BD3200
    PARAMETER: 3.2 micron band depth
    FORMULATION *: 1 - ( R3320 / (a*R3250+b*R3390) )
    RATIONALE: CO2 ice          
    '''
    bands = csas.getbandnumbers(wavelengths, 3250, 3320, 3390)
    band3250 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band3320 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band3390 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    wv3250 = wavelengths[bands[0]]
    wv3320 = wavelengths[bands[1]]
    wv3390 = wavelengths[bands[2]]
    #Algorithm
    a = (wv3320 - wv3250)/ (wv3390 - wv3250)
    b = 1.0 - a
    array_out = 1.0 - (band3320/((b*band3250)+(a*band3390)))
    return array_out 

def bd3400(raster, wavelengths):
    '''
    NAME: BD3400
    PARAMETER: 3.4 micron band depth
    FORMULATION *: 1 - ( (a*R3390+b*R3500) / (c*R3250+d*R3630) )
    RATIONALE: carbonates; organics         
    '''
    bands = csas.getbandnumbers(wavelengths, 3250, 3390, 3500, 3630)
    band3250 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band3390 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band3500 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    band3630 = csas.getarray(raster.GetRasterBand(bands[3]+1))
    wv3250 = wavelengths[bands[0]]
    wv3390 = wavelengths[bands[1]]
    wv3500 = wavelengths[bands[2]]
    wv3630 = wavelengths[bands[3]]
    #Algorithm
    c = (((wv3390+wv3500)*0.5)- wv3250)/ (wv3630 - wv3250)
    d = 1.0 - c
    array_out = 1.0 - (((band3390 + band3500)*0.5)/((d*band3250)+(c*band3630)))
    return array_out 

def cindex(raster, wavelengths):
    '''
    NAME: CINDEX
    PARAMETER: gauge 3.9 micron band
    FORMULATION *: ( R3750 + (R3750-R3630) / (3750-3630) *
      (3920-3750) ) / R3920 - 1
    RATIONALE: carbonates 
    Algorithm differs from published - coded as per CAT
    '''
    bands = csas.getbandnumbers(wavelengths, 3630, 3750, 3950)
    band3630 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band3750 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    band3950 = csas.getarray(raster.GetRasterBand(bands[2]+1))
    wv3630 = wavelengths[bands[0]] / 1000 #Need in microns
    wv3750 = wavelengths[bands[1]]
    wv3950 = wavelengths[bands[2]]
    #Algorithm
    array_out = ((band3750+((band3750-band3630)/((wv3750-wv3630)*(wv3920-wv3750))))) / band3920 - 1
    return array_out 