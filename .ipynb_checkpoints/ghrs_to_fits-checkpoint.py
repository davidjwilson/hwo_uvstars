#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
from astropy.table import Table
from astropy.io import ascii
from scipy.interpolate import interp1d

import astropy.units as u
import astropy.constants as const
from datetime import date

__author__ = 'David Wilson'
__version__ = 0.1
__date__ = '20230209'



"""
Takes the standard GHRS files from MAST, where wavelength, flux and error are in separate files, and combines them into on file. If multiple sub exposures are present then they are coadded into one spectrum. The fits file (for now) uses the header from the c0f (wavelength) file with some keys updated.

The calibrated files used are:
c0f = wavelength
c1f = flux
c2f = error
cqf = data quality flag 
"""

def coadd_flux(f_array, e_array, scale_correct=True):
    """
    Returns a variance-weighted coadd with standard error of the weighted mean (variance weights, scale corrected).
    f_array and e_arrays are collections of flux and error arrays, which should have the same lenth and wavelength scale
    """
    weights = 1 / (e_array**2)
    flux = np.average(f_array, axis =0, weights = weights)
    var = 1 / np.sum(weights, axis=0)
    rcs = np.sum((((flux - f_array)**2) * weights), axis=0) / (len(f_array)-1) #reduced chi-squared
    if scale_correct:
        error = (var * rcs)**0.5
    else:
        error = var**2
    return flux,error

def make_plot(rootname, data):
    fig, ax = plt.subplots(num=rootname)
    ax.step(data['WAVELENGTH'], data['FLUX'], where='mid', label='FLUX')
    ax.step(data['WAVELENGTH'], data['ERROR'], where='mid', alpha=0.5, label='ERROR')
    ax.set_xlabel('Wavelength (\AA)')
    ax.set_ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)')
    ax.set_title(rootname)
    ax.legend()
    fig.tight_layout()


def make_ghrs_file(inpath, outpath, plot=True):
    """
    inpath = where the files from mast are stored
    outpath = where you would like the output files to go.
    """
    
    c0f_list = glob.glob('{}*c0f.fits'.format(inpath))
    roots = []
    for c0f in c0f_list:
        roots.append(fits.getheader(c0f, 0)['ROOTNAME'])
    print('inpath contains c0f (wavelength) files for the following datasets:', roots)
    
    for c0f in c0f_list:
        hdr = fits.getheader(c0f, 0)
        rootname = hdr['ROOTNAME']
        print('working with dataset {}'.format(rootname))
        wavelength_arrays = fits.getdata(c0f, 0)
        c1f = '{}{}/{}_c1f.fits'.format(path, star, rootname.lower())
        try: 
            flux_arrays = fits.getdata(c1f, 0) #test if this works
        except: 
                print('WARNING: no c1f (flux) file found for dataset {}, skipping'.format(rootname))
                continue
        c2f = '{}{}/{}_c2f.fits'.format(path, star, rootname.lower())
        try:
            error_arrays = fits.getdata(c2f, 0)
        except: 
            print('WARNING: no c2f (error) file found for dataset {}, skipping'.format(rootname))
            continue
        cqf = '{}{}/{}_cqf.fits'.format(path, star, rootname.lower())
        try:
            dq_arrays = fits.getdata(cqf, 0)
        except: 
            print('WARNING: no cqf (data quality) file found for dataset {}, skipping'.format(rootname))
            continue
    
        if len(np.shape(wavelength_arrays)) == 1: #check if there are multiple subexposures that need to be combined
            wavelength, flux, error, dq = wavelength_arrays, flux_arrays, error_arrays, dq_arrays
        else:
            print('this dataset has subexposures')
            wavelength = wavelength_arrays[0] #each wavelength array in a subexposure is the same as far as I can tell
            flux, error = coadd_flux(flux_arrays, error_arrays)
            dq = [(np.sum(np.unique(dq_arrays[:,i]))) for i in range(len(dq_arrays[0]))]
        
        filename = '{}_{}_{}_{}.fits'.format(hdr['TARGNAME'].lower(), hdr['INSTRUME'].lower(), hdr['GRATING'].lower(), hdr['ROOTNAME'].lower())
      
        data = Table((wavelength*u.AA, flux*u.erg/u.s/u.cm**2/u.AA, error*u.erg/u.s/u.cm**2/u.AA, dq), names = ['WAVELENGTH', 'FLUX', 'ERROR', 'DQ'])
        
        if plot:
            make_plot(rootname, data)
        
        data_ext = fits.table_to_hdu(data)
        data_ext.header.set('EXTNAME', 'SPECTRUM')
        hdr_new = hdr
        hdr_new.set('FILENAME', filename)
        hdr_new.set('FILETYPE', 'SPECTRUM')
        hdr_new.set('FITSDATE', date.today().strftime('%Y-%m-%d'))
        primary_hdu = fits.PrimaryHDU(header=hdr_new)
        
        hdul = fits.HDUList([primary_hdu, data_ext])
        hdul.writeto('{}/{}'.format(outpath, filename), overwrite=True)
        print('spectrum saved as {}'.format(filename))
        
    if plot:
        plt.show()
        
#test
path = '/media/david/2tb_ext_hd/hddata/ghrs_test/' #have a couple of stars in path
stars = os.listdir(path)
for star in stars:
    inpath = '{}{}/'.format(path, star)
    outpath = 'ghrs_test_fits'
    
    make_ghrs_file(inpath, outpath, plot=True)
    
    
