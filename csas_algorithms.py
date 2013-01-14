#!/usr/bin/env python

import csas
def ferric_coating(raster, wavelengths):
    '''
    Name: ISPLOE1
    Parameter: -1 * Spectral Slope_1
    Formulation: (R1815-R2530)/(R2530-R1815)
    Rationale: ferric coating on dark rock
    Band(s): 1815, 2530
    '''

    bands = csas.getbandnumbers(wavelengths, 1815, 2530)
    print csas.getarray(raster.GetRasterBand(bands[0]+1))

    array1815 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    array2530 = csas.getarray(raster.GetRasterBand(bands[1]+1))

    #Algorithm
    array_out = (array1815-array2530)/(array2530-array1815)
    
    return array_out

def carbonates(raster, wavelengths):
    '''
    Name: CINDEX
    Parameter: gauge 3.9micron band
    Formulation: (R3750 + (R3750-R3630)/((R3750-R3630)*(R3950-R3750)))/R3950-1
    Rationale: Carbonates.
    Band(s):29, 47, 1, !!!1 is empty - remapped to 5 to get output!!!!!
    '''
    
    R3750 = csas.getarray(raster.GetRasterBand(29))
    R3630 = csas.getarray(raster.GetRasterBand(47))
    R3950 = csas.getarray(raster.GetRasterBand(5))

    #Algorithm
    array_out = (R3750 + (R3750-R3630)/((R3750-R3630)*(R3950-R3750)))/R3950-1
    return array_out
    
def hcp_index(raster, wavelengths):
        '''
        Name: HCPINDEX
        Parameter: pyroxene index
        Formulation: ((R1470-R1050)/(R1470+R1050))*((R1470-R1050)/(R1470+R2067))
        Rationale: pyroxene will be strongly positive;favors HCP
        Band(s): 374, 438 !!!RESET TO 436, the last band with valid data!!!, 253
        '''
            
        array1470 = csas.getarray(raster.GetRasterBand(374))
        array1050 = csas.getarray(raster.GetRasterBand(436))
        array2067 = csas.getarray(raster.GetRasterBand(253))
        
        #Algorithm
        array_out = ((array1470-array1050)/(array1470+array1050))*((array1470-array1050)/(array1470+array2067))
        
        return array_out
    
def h2o1(raster, wavelengths):
    
    '''
    Name: BD3000
    Parameter: 3micron band depth
    Formulation: 1 - (R3000/(R2530*(R2530/R2210)))
    Rationale: h2o.
    Band(s): 142, 213, 262
    '''
        
    R3000 = csas.getarray(raster.GetRasterBand(142))
    R2530 = csas.getarray(raster.GetRasterBand(213))
    R2210 = csas.getarray(raster.GetRasterBand(262))

    #Algorithm
    array_out = 1 - (R3000/(R2530*(R2530/R2210)))
    
    return array_out    
    
def h2o_co2_ice(raster, wavelengths):
    
    '''
    Name: ICER1
    Parameter: 1.5micron and 1.43micron band ration
    Formulation: R1510 / R1430
    Rationale: CO2, H2O Ice Mixtures.
    Band(s): 368, 380
    '''
        
    array1510 = csas.getarray(raster.GetRasterBand(368))
    array1430 = csas.getarray(raster.GetRasterBand(380))

    #Algorithm
    array_out = array1510/array1430
    
    return array_out

def h2o_co2_ice2(raster, wavelengths):
    
    '''
    Name: ICER2
    Parameter: gauge 2.7 micron band
    Formulation: R2530 / R2600
    Rationale: CO2 ice will be >> 1; H2O ice and soil will be ~ 1.
    Band(s): 213, 203
    '''
        
    array2530 = csas.getarray(raster.GetRasterBand(213))
    array2600 = csas.getarray(raster.GetRasterBand(203))

    #Algorithm
    array_out = array2530/array2600
    
    return array_out

def ira(raster, wavelengths):
    '''
    Name:IRA
    Parameter: 1.3micron reflectance
    Formulation: R1330
    Rationale: IR albedo
    Band(s): 395
    '''
        
    array1330 = csas.getarray(raster.GetRasterBand(395))
    
    return array1330

def lcp_index(raster, wavelengths):
        '''
        Name: LCPINDEX
        Parameter: pyroxene index
        Formulation: ((R1330-R1050)/(R1330+R1050))*((R1330-R1815)/(R1330+R1815))
        Rationale: pyroxene will be strongly positive;favors LCP
        Band(s): 395, 438 !!!RESET TO 436, the last band with valid data!!!, 321, 
        '''
            
        array1330 = csas.getarray(raster.GetRasterBand(395))
        array1050 = csas.getarray(raster.GetRasterBand(436))
        array1815 = csas.getarray(raster.GetRasterBand(321))
        
        #Algorithm
        array_out = ((array1330-array1050)/(array1330+array1050))*((array1330-array1815)/(array1330+array1815))
        
        return array_out
    
def olivine_index(raster, wavelengths):
    '''
    Name: OLINDEX
    Parameter: Olivine Index
    Formulation: (R1695/(0.1*R1050 + 0.1*R1210 + 0.4*R1330 + 0.4*R1470))-1
    Rationale: Olivine will be strongly positive; based on fayalite
    Band(s): 340, 438 !!!438 remapped to 436 for valid values!!!, 414, 395, 374
    '''
    
    NDV = raster.GetRasterBand(1).GetNoDataValue()
        
    array1695 = csas.getarray(raster.GetRasterBand(340))
    array1050 = csas.getarray(raster.GetRasterBand(436))
    array1210 = csas.getarray(raster.GetRasterBand(414))
    array1330 = csas.getarray(raster.GetRasterBand(395))
    array1470 = csas.getarray(raster.GetRasterBand(374))
    
    #Algorithm
    array_out = (array1695/((0.1*array1050)+(0.1*array1210)+(0.4*array1330)+(0.4*array1470)))-1
    
    return array_out

    
def rockdust1(raster, wavelengths):
        '''
        Name: R770
        Parameter: 0.77micron reflectance
        Formulation: R770
        Rationale: rock/dust
        Band(s): 248
        '''
        bands = csas.getbandnumbers(wavelengths, 770)    
        array770 = csas.getarray(raster.GetRasterBand(bands[0]+1))
        
        return array770


def rockdust2(raster, wavelengths):
    '''
    Name: RBR
    Parameter: Red/Blue Ratio
    Formulation: R770 / R440
    Rationale: rock/dust
    Bands: 248, 196
    '''
    bands = csas.getbandnumbers(wavelengths, 440,770)  
    print bands
    band440 = csas.getarray(raster.GetRasterBand(bands[0]+1))
    band770 = csas.getarray(raster.GetRasterBand(bands[1]+1))
    
    
    #Algorithm
    array_out = band770 / band440
    
    return array_out