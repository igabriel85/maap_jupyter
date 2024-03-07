# -*- coding: utf-8 -*-
"""
Property of the European Space Agency (ESA-ESRIN) - contact: clement.albinet@esa.int / nuno.miranda@esa.int
Developed for Python 3.5.1 and GDAL 2.0.2
Date  -  Version  -  Author(s)  -  List of changes
17/05/2016 - V1.0 - Clement Albinet - First version of the code.
01/02/2017 - V1.1 - Clement Albinet - Projectors now work with multi-bands files.
14/04/2017 - V1.2 - Clement Albinet - Update of GrdToSlrProj to work with GDAL instead of numpy (~40 times faster).
"""
########## ########## ORCHESTRATOR OF THE BIOMASS ALGORITHM TEST BED ########## ########## 

from osgeo import gdal
import numpy as np

##########################################################################################
def SlrToGrdProj(slrFile, grdFile, azimuthFile, rangeFile):
    'Projection of an image from Slant Range geometry to Ground Projected geometry'
    # Open original image in slant range geometry:
    slr_image_driver = gdal.Open(slrFile, 0)
    slr_image = slr_image_driver.ReadAsArray(0,0,1,1)
    
    # Open Azimuth coordinates file:
    Azimuth_driver = gdal.Open(azimuthFile, 0)
    Azimuth = Azimuth_driver.ReadAsArray()
    
    # Open Range coordinates file:
    Range_driver = gdal.Open(rangeFile, 0)
    Range = Range_driver.ReadAsArray()
    
    # Mask of the data inside the GRD projected image:
    mask = np.logical_and(Range!=55537, Azimuth!=55537)
    
    # Create an empty image of NaN:
    grd_image = np.full((Range_driver.RasterYSize, Range_driver.RasterXSize), np.NaN, dtype=slr_image.dtype)
    
    # Create the image in the ground projected geometry:
    outdriver = gdal.GetDriverByName('GTiff')
    grd_image_driver = outdriver.Create(grdFile, Range_driver.RasterXSize, Range_driver.RasterYSize, slr_image_driver.RasterCount, slr_image_driver.GetRasterBand(1).DataType)
    grd_image_driver.SetGeoTransform(Range_driver.GetGeoTransform())
    grd_image_driver.SetProjection(Range_driver.GetProjection())
    
    for band in range(slr_image_driver.RasterCount):
        # Read original image in slant range geometry:
        slr_image = slr_image_driver.GetRasterBand(band+1).ReadAsArray()
        
        # Project the image in the ground projected geometry:
        grd_image[mask] = slr_image[Azimuth[mask], Range[mask]]
        
        # Save the corresponding band of the image in the ground projected geometry:
        grd_image_driver.GetRasterBand(band+1).WriteArray(grd_image)
    
    # Close data sets:
    slr_image_driver = None
    Azimuth_driver = None
    Range_driver = None
    grd_image_driver = None

##########################################################################################
def GrdToSlrProj(grdFile, slrFile, azimuthFile, rangeFile, originalImageFile):
    'Projection of an image from Ground Projected geometry to Slant Range geometry'
    'Remark: Computation can need a long time'
    from osgeo import gdal
    from gdalconst import GA_ReadOnly
    
    # Open original image in slant range geometry:
    grd_image_driver = gdal.Open(grdFile, GA_ReadOnly)
    grd_image = grd_image_driver.ReadAsArray(0,0,1,1)
    
    # Open Azimuth coordinates file:
    Azimuth_driver = gdal.Open(azimuthFile, GA_ReadOnly)
    Azimuth = Azimuth_driver.ReadAsArray()
    
    # Open Range coordinates file:
    Range_driver = gdal.Open(rangeFile, GA_ReadOnly)
    Range = Range_driver.ReadAsArray()
    
    # Open original image file:
    original_driver = gdal.Open(originalImageFile, GA_ReadOnly)
    
    # Mask of the data inside the GRD projected image:
    mask = np.logical_and(Range!=55537, Azimuth!=55537)
    
    # Create the image in the slant range geometry:
    outdriver = gdal.GetDriverByName('GTiff')
    slr_image_driver = outdriver.Create(slrFile, original_driver.RasterXSize, original_driver.RasterYSize, grd_image_driver.RasterCount, grd_image_driver.GetRasterBand(1).DataType)
    
    for band in range(grd_image_driver.RasterCount):
        # Read original image in slant range geometry:
        grd_image = grd_image_driver.GetRasterBand(band+1).ReadAsArray()
        
        # Create an empty image of NaN:
        slr_image = np.full((original_driver.RasterYSize, original_driver.RasterXSize), np.NaN, dtype=grd_image.dtype)
        
        # Project the image in the slant range geometry:
        slr_image[Azimuth[mask], Range[mask]] = grd_image[mask]
        
        # Save the corresponding band of the image in the ground projected geometry:
        slrBand = slr_image_driver.GetRasterBand(band+1)
        slrBand.WriteArray(slr_image)
        slrBand.FlushCache()
        
        # Get NaN values to interpolate
        maskNaN = ~np.isnan(slr_image)
        
        # Create and fill a dataset with NaN values previously obtained:
        mask_driver = gdal.GetDriverByName('MEM').Create('', original_driver.RasterXSize, original_driver.RasterYSize, grd_image_driver.RasterCount, grd_image_driver.GetRasterBand(1).DataType)
        maskBand = mask_driver.GetRasterBand(1)
        maskBand.WriteArray(maskNaN)
        maskBand.FlushCache()
        
        # Interpolate missing values:
        gdal.FillNodata(targetBand = slrBand, maskBand = maskBand, maxSearchDist = 5, smoothingIterations = 0)
        
        # Close temporary NaN mask dataset:
        mask_driver = None
    
    # Close data sets:
    grd_image_driver = None
    Azimuth_driver = None
    Range_driver = None
    original_driver = None
    slr_image_driver = None


##########################################################################################
if (__name__ == '__main__'):
    
#    folder = r'C:\Users\Clement Albinet\Desktop\2 Biomass Algo\Test Bed input_data\test_Biosar1'
#    slrFile = folder + r'\i07biosar0105x1_ch1_t01_slc.dat.final_dB.tiff'
#    grdFile = folder + r'\i07biosar0105x1_ch1_t01_slc.dat.final_dB_GTC.tiff'
#    rangeFile = folder + r'\range_slc07biosar0105x1_t01_int.tiff'
#    azimuthFile = folder + r'\azimuth_slc07biosar0105x1_t01_int.tiff'
    
#    folder = r'C:\Users\Clement Albinet\Desktop\2 Biomass Algo\Test Bed input_data\test_Biosar2'
#    slrFile = folder + r'\i08biosar0103x1_ch1_t01_slc_dB.tiff'
#    grdFile = folder + r'\i08biosar0103x1_ch1_t01_slc_dB_GTC.tiff'
#    rangeFile = folder + r'\range_slc08biosar0103x1.tiff'
#    azimuthFile = folder + r'\azimuth_slc08biosar0103x1.tiff'
    
#    folder = r'C:\Users\Clement Albinet\Desktop\2 Biomass Algo\Test Bed input_data\test_Biosar1'
#    slrFile = folder + r'\i07biosar0412x1_ch1_t02_slc_dB.tiff'
#    grdFile = folder + r'\i07biosar0412x1_ch1_t02_slc_dB_GTC.tiff'
#    rangeFile = folder + r'\range_slc07biosar0412x1_t02_int.tiff'
#    azimuthFile = folder + r'\azimuth_slc07biosar0412x1_t02_int.tiff'
    
    # Test processing of azimuth and range maps (tropisar) :
#    folder = r'C:\Users\Clement Albinet\Desktop\2 Biomass Algo\Test Bed input_data\test_Tropisar'
#    slrFile = folder + r'\tropi0402_Pcons_Hh_slc_dB.tiff'
#    grdFile = folder + r'\tropi0402_Pcons_Hh_slc_dB_GTC.tiff'
#    rangeFile = folder + r'\rg.tiff'
#    azimuthFile = folder + r'\az.tiff'
    
    # Test processing of azimuth and range maps (biosar3) :
#    folder = r'C:\Users\Clement Albinet\Desktop\2 Biomass Algo\Test Bed input_data\test_Biosar3'
#    slrFile = folder + r'\tropisar_01_HH_dB.tiff'
##    slrFile = folder + r'\tropisar_01_HH.tiff'
#    grdFile = folder + r'\tropisar_01_HH_dB_GTC.tiff'
##    grdFile = folder + r'\tropisar_01_HH_GTC.tiff'
#    rangeFile = folder + r'\rg.tiff'
#    azimuthFile = folder + r'\az.tiff'
##    slrFile2 = folder + r'\tropisar_01_HH_dB_GTC_RGI.tiff'
##    slrFile2 = folder + r'\tropisar_01_HH_GTC_RGI.tiff'
    
    # Test processing of azimuth and range maps (afrisar) :
    folder = r'C:\Users\Clement Albinet\Desktop\2 Biomass Algo\Test Bed input_data\test_Afrisar'
    slrFile = folder + r'\20150705B-10_sar_UHF50MHzHAM_Hh_rad_dB.tiff'
    grdFile = folder + r'\20150705B-10_sar_UHF50MHzHAM_Hh_rad_dB_GRD.tiff'
    rangeFile = folder + r'\rg_linear.tiff'
    azimuthFile = folder + r'\az_linear.tiff'
    
    
    SlrToGrdProj(slrFile, grdFile, azimuthFile, rangeFile)
    
#    slrFile2 = folder + r'\i07biosar0105x1_ch1_t01_slc.dat.final_dB_GTC_RGI.tiff'
#    slrFile2 = folder + r'\i08biosar0103x1_ch1_t01_slc_dB_GTC_RGI.tiff'
    
#    GrdToSlrProj(grdFile, slrFile2, azimuthFile, rangeFile, slrFile)
    
    print(' - Fin -')

