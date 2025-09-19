import os
import numpy as np
import pandas as pd
import math
import matplotlib ; import matplotlib.pyplot as plt ; import matplotlib as mpl ; from matplotlib import colors ; import matplotlib.patches as mpatches
import matplotlib.ticker as ticker ; from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import warnings
from astropy.io.fits.verify import VerifyWarning ; from astropy.io import fits ; from astropy.time import Time ; from astropy.io import ascii ; from astropy.wcs import WCS ; from astropy.stats import SigmaClip
import photutils as phot ; from photutils.background import Background2D, MMMBackground
import casatasks
import shutil
from pathlib import Path
import glob
import re
warnings.simplefilter("ignore",VerifyWarning)

def duplicate_file(src_file, new_name):
    #Creation of a duplicate of a given file with a new name in the same directory.
    #Parameters:
    #    src_file (str): Path to the source file.
    #    new_name (str): New name for the duplicate file.
    #Returns:
    #    str: Path to the duplicated file.

    src_dir = os.path.dirname(src_file) ; dest_file = os.path.join(src_dir, new_name) ; shutil.copy(src_file, dest_file)
    return dest_file

def duplicate_directory(src_dir, dest_dir):
    #Creation of a duplicate of a given directory. If the destination directory exists, it is removed before copying.
    #Parameters:
    #    src_dir (str): Path to the source directory.
    #    dest_dir (str): Path to the destination directory.
    #Returns:
    #    str: Path to the duplicated directory.
    
    if os.path.exists(dest_dir):
        shutil.rmtree(dest_dir)
    shutil.copytree(src_dir, dest_dir)
    return dest_dir

def contours(l1,b):
    #Calculates contour levels for an image.
    #Parameters:
    #    l1 (float): Base level.
    #    b (int): Power exponent.
    #Returns:
    #    float: Computed contour level.
    
    return l1*np.power(2,b)

def map_properties(file):
    #Extracts properties from a FITS map file.
    #Parameters:
    #    file (str): Path to the FITS file.
    #Returns:
    #    tuple: Map data, pixel scale, increment, central pixel scale, peak flux, peak index, observation date in MJD, beam major/minor axes, beam position angle, and total flux from the map.
    
    hdul_map = fits.open(file) ; hdr_map = hdul_map[0].header
    map_data = hdul_map[0].data ; map_data = np.reshape(map_data,(hdr_map['NAXIS1'],hdr_map['NAXIS1']))
    #Map Properties
    delta = hdr_map['CDELT2'] ; increment = delta*3600000.0  ; c_pix = 1/increment
    f_peak = np.max(map_data) ; peak_id = np.unravel_index(np.argmax(map_data,axis=None),map_data.shape)
    MJD = Time(hdr_map['DATE-OBS'],format='iso').mjd
    #Beam Size
    BMaj,BMin,BPA = hdr_map['BMAJ']*3600,hdr_map['BMIN']*3600,hdr_map['BPA']
    #Total flux (Sum pixels/Beam area)
    cdelt1,cdelt2 = abs(hdr_map['CDELT1'])*3600,abs(hdr_map['CDELT2'])*3600
    pix_area = cdelt1*cdelt2 ; beam_area = (math.pi*BMaj*BMin)/(4*math.log(2))
    f_total = np.sum(map_data)*(pix_area/beam_area)
    return map_data,delta,increment,c_pix,f_peak,peak_id,MJD,BMaj,BMin,BPA,f_total

def detection_limit(treshold,file,map_data):
    #Computes the detection limit based on background noise in the map.
    #Parameters:
    #    threshold (float): Sigma level for detection limit.
    #    file (str): Path to the FITS file.
    #    map_data (numpy.ndarray): Image data from the FITS file.
    #Returns:
    #    float: Computed detection limit in Jy/beam.
    
    mmm_bkg = MMMBackground() ; sigma_clip = SigmaClip(sigma_lower=3,sigma_upper=5)
    bkg1 = Background2D(map_data,(64,64),filter_size=(3,3),sigma_clip=sigma_clip,bkg_estimator=mmm_bkg)
    bkg_err_median = bkg1.background_rms_median
    det_limit = treshold*bkg_err_median
    print(f'Detection Limit of: {treshold}sigma = {det_limit:.4f} Jy/beam')
    return det_limit

def map_fitting(folder,file,f_peak,peak_id,BMaj,BMin,BPA,det_limit):
    #Iteratively fits Gaussian components in the map using CASA until the peak flux is below the detection limit.
    #In each iterarion the output of imfit is directly stored into the directory without a proper return value.
    #Parameters:
    #    folder (str): Directory where output files will be stored.
    #    file (str): Path to the FITS file.
    #    f_peak (float): Peak flux of the map.
    #    peak_id (tuple): Coordinates of the peak flux.
    #    BMaj (float): Major axis of the beam.
    #    BMin (float): Minor axis of the beam.
    #    BPA (float): Beam position angle.
    #    det_limit (float): Detection limit for stopping the fitting process.
    #Returns:
    #    count (int): Number of fitted components.
    #    term (str): Termination status, 'Final_T' if last iteration did not converge, else 'Final'.

    
    count = 0 ; print('Iter.\tResidual.\tx\ty')
    while f_peak >= det_limit:
        #Initial Estimates
        mode = 'w+' if count == 0 else 'a'
        with open(f'{folder}/init{count}.txt',mode) as outfile:
            outfile.write(f'{f_peak},  {peak_id[1]},  {peak_id[0]},  {BMaj} arcsec,  {BMin} arcsec,  {BPA} deg \r\n')
        outfile.close()
        #Fit
        casatasks.imfit(file,residual=f'{folder}/res{count+1}.IMAP',model=f'{folder}/model{count+1}.IMAP',estimates=f'{folder}/init{count}.txt',logfile=f'{folder}/log{count+1}.txt',append=False,newestimates=f'{folder}/iter{count+1}.txt',summary=f'{folder}/sum{count+1}.txt')
        casatasks.exportfits(f'{folder}/res{count+1}.IMAP',fitsimage=f'{folder}/res{count+1}.fits',overwrite=True)
        #Read the residual to find the maximum value and its indices
        residual = f'{folder}/res{count+1}.fits'
        hdul_res = fits.open(residual) ; hdr_res = hdul_res[0].header
        res_data = hdul_res[0].data ; res_data = np.reshape(res_data,(hdr_res['NAXIS1'],hdr_res['NAXIS1']))
        f_peak = np.max(res_data) ; peak_id = np.unravel_index(np.argmax(res_data,axis=None),res_data.shape)
        #Creation of the new estimates
        count += 1 ; print(f'{count}\t{f_peak:.3f}\t{peak_id[1]}\t{peak_id[0]}')
        if f_peak >= det_limit:
            duplicate_file(f'{folder}/iter{count}.txt',f'init{count}.txt')
    print('Fit Done')
    #Checking if the fit was succesfull
    if f_peak <= 0:
        count -= 1 ; term = 'Final_T'
    else:
        term = 'Final'
    return count,term

def map_cleaning(folder,file,count,det_limit,BMaj,BMin,term):
    #Performs post-fitting cleaning by filtering out bad components from the fitted VLBI model.
    #The cleaning is done based on the detection limit, eccentricity, and size of the components.
    #The approved final fit is directly stored into the directory without a proper return value.
    #Parameters:
    #   folder (str): Path to the directory where output files are stored.
    #   file (str): Path to the FITS file.
    #   count (int): Number of fitted components from the previous step.
    #   det_limit (float): Detection threshold below which components are discarded.
    #   BMin (float): Minimum acceptable size for the minor axis of a component.
    #   term (str): Label for final output files.
    #Returns:
    #   len(Itera) (int): Number of fitted final components.
    #   clean_state (str): Status of the cleaning process:
    #                      'AutoClean' if te cleaning was performed successfully. 'CleanToCheck' if the cleaning fit did not converge. 'NoClean' if no cleaning was necessary.

    Summary = (ascii.read(f'{folder}/sum{count}.txt',header_start=1)).to_pandas()
    Itera = (ascii.read(f'{folder}/iter{count}.txt', names=['fpeak','pix_x','pix_y','pix_a','pix_b','pix_p'])).to_pandas()
    b = Itera['pix_b'].apply(lambda x: float(x.split()[0])) ; a = Itera['pix_a'].apply(lambda x: float(x.split()[0]))
    beam_ratio = BMaj/BMin
    gaussian_ratio = a/b
    ecc_factor = 7
    beam_area = np.pi*BMaj*BMin/(4*np.log(2))
    gaussian_area = np.pi*a*b/(4 * np.log(2))
    if (Summary['Peak']<det_limit).any():
        print('There are components below Detection Limit')
    if (gaussian_ratio > ecc_factor*beam_ratio).any():
        print('There are components with eccentricity greater than 3 times that of the beam')
    if (gaussian_area < 0.9*beam_area).any():
        print('There are components smaller than the beam')
    if (Summary['Peak']<det_limit).any() or (gaussian_ratio > ecc_factor*beam_ratio).any() or (gaussian_area < 0.9*beam_area).any():
        with open(f'{folder}/init_{term}.txt','w+') as outfile_final:
            for j in Summary.index:
                if Summary['Peak'][j] >= det_limit and gaussian_ratio[j] <= ecc_factor*beam_ratio and gaussian_area[j] >= 0.9*beam_area:
                    row = ',  '.join(Itera.iloc[j].astype(str).values)
                    outfile_final.write(row+'\n')
        casatasks.imfit(file,residual=f'{folder}/{term}_res.IMAP',model=f'{folder}/{term}_model.IMAP',estimates=f'{folder}/init_{term}.txt',logfile=f'{folder}/{term}_log.txt',append=False,newestimates=f'{folder}/{term}_iter.txt',summary=f'{folder}/{term}_sum.txt')
        casatasks.exportfits(f'{folder}/{term}_model.IMAP',fitsimage=f'{folder}/{term}_model.fits',overwrite=True)
        casatasks.exportfits(f'{folder}/{term}_res.IMAP',fitsimage=f'{folder}/{term}_res.fits',overwrite=True)
        print('Map cleaning done')
        #Read the residual to ensure the convergence of the cleaning
        residual = f'{folder}/{term}_res.fits'
        hdul_res = fits.open(residual) ; hdr_res = hdul_res[0].header
        res_data = hdul_res[0].data ; res_data = np.reshape(res_data,(hdr_res['NAXIS1'],hdr_res['NAXIS1'])) ; f_peak = np.max(res_data)
        if f_peak <= 0:
            print('The fitting of the cleaned map did not converge. Delivering map with all components as final')
            for file_type in ['model', 'res']:
                duplicate_directory(f'{folder}/{file_type}{count}.IMAP', f'{folder}/{term}_{file_type}.IMAP')
            for file_type in ['iter', 'log', 'sum']:
                duplicate_file(f'{folder}/{file_type}{count}.txt',f'{term}_{file_type}.txt')
            casatasks.exportfits(f'{folder}/{term}_model.IMAP',fitsimage=f'{folder}/{term}_model.fits',overwrite=True)
            casatasks.exportfits(f'{folder}/{term}_res.IMAP',fitsimage=f'{folder}/{term}_res.fits',overwrite=True)
            clean_state = 'CleanToCheck'
        else: 
            clean_state = 'AutoClean'
    else:
        print('There are no bad components')
        for file_type in ['model', 'res']:
            duplicate_directory(f'{folder}/{file_type}{count}.IMAP', f'{folder}/{term}_{file_type}.IMAP')
        for file_type in ['iter', 'log', 'sum']:
            duplicate_file(f'{folder}/{file_type}{count}.txt',f'{term}_{file_type}.txt')
        casatasks.exportfits(f'{folder}/{term}_model.IMAP',fitsimage=f'{folder}/{term}_model.fits',overwrite=True)
        casatasks.exportfits(f'{folder}/{term}_res.IMAP',fitsimage=f'{folder}/{term}_res.fits',overwrite=True)
        clean_state = 'NoClean'
    return len(Itera),clean_state

def TicksEstimator(file):
    #Estimates tick positions and their corresponding RA and Dec labels for a FITS image.
    #Parameters:
    #   file (str): Path to the FITS file.
    #Returns:
    #   xticks/yticks (numpy.ndarray): X-axis/Y-axis pixel coordinates for ticks.
    #   RA_short/Dec_short (list of str): Formatted Right Ascension/Declination labels.

    hdul_map = fits.open(file) ; hdr_map = hdul_map[0].header ; wcs = WCS(hdr_map,naxis=2)
    xticks,yticks = np.linspace(0,hdr_map['NAXIS1'],6),np.linspace(0,hdr_map['NAXIS2'],6)
    
    RAs = [wcs.pixel_to_world(x,y).ra.to_string(unit='hourangle', precision=5) for x, y in zip(xticks,yticks)]
    Decs = [wcs.pixel_to_world(x,y).dec.to_string(precision=5) for x, y in zip(xticks,yticks)]
    var_x = next(i for i, val in enumerate(zip(*RAs)) if any(s != RAs[0][i] for s in val)) ; var_y = next(i for i, val in enumerate(zip(*Decs)) if any(s != Decs[0][i] for s in val))
    dec_x = RAs[0].find('.') if '.' in RAs[0] else None ; dec_y = Decs[0].find('.') if '.' in Decs[0] else None
    m_x = RAs[0].find('m') if 'm' in RAs[0] else None ; m_y = Decs[0].find('m') if 'm' in Decs[0] else None
    h_x = RAs[0].find('h') if 'h' in RAs[0] else None ; d_y = Decs[0].find('d') if 'd' in Decs[0] else None
    posx = (0 if var_x < h_x else h_x+1 if var_x < m_x else m_x+1 if var_x < dec_x else dec_x)
    posy = (0 if var_y < d_y else d_y+1 if var_y < m_y else m_y+1 if var_y < dec_y else dec_y)
    RA_short = [RA if i == 0 else RA[posx:] for i, RA in enumerate(RAs)]
    Dec_short = [Dec if i == len(Decs)-1 else Dec[posy:] for i, Dec in enumerate(Decs)]
    return xticks,yticks,RA_short,Dec_short

def RMS_estimations(folder,file,term):
    #Calculates the RMS values for the map, model, and residual FITS files. Also computes the residual-to-map and model-to-map ratios of RMS.
    #Parameters:
    #    folder (str): Directory where the output files are located.
    #    file (str): Path to the FITS file.
    #    term (str): Termination status.
    #Returns:
    #    map_RMS,model_RMS,res_RMS,rate_min,rate_max  (float): RMS value of the map data ,model data, residual data, and ratios.

    hdul_map = fits.open(file) ; hdr_map = hdul_map[0].header ; map_data = hdul_map[0].data
    map_data = np.reshape(map_data,(hdr_map['NAXIS1'],hdr_map['NAXIS1'])) ; map_RMS = np.sqrt(np.mean(map_data**2))

    model = f'{folder}/{term}_model.fits'
    hdul_model = fits.open(model) ; hdr_model = hdul_model[0].header ; model_data = hdul_model[0].data
    model_data = np.reshape(model_data,(hdr_model['NAXIS1'],hdr_model['NAXIS1'])) ; model_RMS = np.sqrt(np.mean(model_data**2))

    residual = f'{folder}/{term}_res.fits'
    hdul_res = fits.open(residual) ; hdr_res = hdul_res[0].header ; res_data = hdul_res[0].data
    res_data = np.reshape(res_data,(hdr_res['NAXIS1'],hdr_res['NAXIS1'])) ; res_RMS = np.sqrt(np.mean(res_data**2))

    rate_min = res_RMS/map_RMS ; rate_max = model_RMS/map_RMS
    return map_RMS,model_RMS,res_RMS,rate_min,rate_max


def ploting(folder,file,term,det_limit,increment,f_peak,BMaj,BMin,BPA):
    #Plots the map, model, residuals, and fitted components in a multi-panel figure. The plot includes flux contours, beam size, and component positions overlaid on the map.
    #Parameters:
    #    folder (str): Directory where the output files are located.
    #    file (str): Path to the FITS file.
    #    term (str): Termination status.
    #    det_limit (float): The detection limit for the contours and components.
    #    increment (float): Pixel size in arcseconds for scaling the beam and components.
    #    f_peak (float): The peak flux value of the map for colormap scaling.
    #    BMaj (float): Major axis of the beam in arcseconds.
    #    BMin (float): Minor axis of the beam in arcseconds.
    #    BPA (float): Position angle of the beam in degrees.
    #Returns:
    #    None: Saves a plot image with the map, model, residuals, and fitted components.

    model = f'{folder}/{term}_model.fits' ; residual = f'{folder}/{term}_res.fits' ; fitted = f'{folder}/{term}_iter.txt' ;summmary = f'{folder}/{term}_sum.txt'
    
    file_name = Path(file).stem ; levels = contours(det_limit,np.arange(0,10,1))
    #Map
    hdul_map = fits.open(file) ; hdr_map = hdul_map[0].header ; wcs = WCS(hdr_map,naxis=2) ; date = hdr_map['DATE-OBS']
    map_data = hdul_map[0].data ; map_data = np.reshape(map_data,(hdr_map['NAXIS1'],hdr_map['NAXIS1']))
    #Model
    hdul_model = fits.open(model) ; hdr_model = hdul_model[0].header
    model_data = hdul_model[0].data ; model_data = np.reshape(model_data,(hdr_model['NAXIS1'],hdr_model['NAXIS1']))
    #Residual
    hdul_res = fits.open(residual) ; hdr_res = hdul_res[0].header
    res_data = hdul_res[0].data ;res_data = np.reshape(res_data,(hdr_res['NAXIS1'],hdr_res['NAXIS1']))
    #Components
    Fitted = (ascii.read(fitted,names=['fpeak','pix_x','pix_y','pix_a','pix_b','pix_p'])).to_pandas()
    b = Fitted['pix_b'].apply(lambda x: float(x.split()[0])) ; a = Fitted['pix_a'].apply(lambda x: float(x.split()[0]))
    Fitted['eccen'] = np.sqrt(1-(b**2/a**2))
    Summary = (ascii.read(summmary,header_start=1)).to_pandas()
    properties = pd.DataFrame({'Integrated':Summary['I'],'IntErr':Summary['Ierr'],'Peak':Summary['Peak'],'PeakErr':Summary['PeakErr'],\
                           'RAJ2000':Summary['RAJ2000'].apply(lambda x:x+360 if x<0 else x),'DecJ2000':Summary['DecJ2000'],\
                           'Pos_X':Fitted['pix_x'],'Pos_Y':Fitted['pix_y'],'MajAx':Summary['ConMaj'],'MajAxErr':Summary['ConMajErr'],'MinAx':Summary['ConMin'],'MinAxErr':Summary['ConMinErr'],\
                           'Eccen':Fitted['eccen'],'PA':Summary['ConPA'],'ConPAErr':Summary['ConPAErr'],'Freq':Summary['Freq']})
    units = ['Jy','Jy','Jy/beam','Jy/beam','deg','arcsec','deg','arcsec','pixel','pixel','arcsec','arcsec','arcsec','arcsec','no_units','deg','deg','GHz'] ; headers = pd.DataFrame([units])
    with open(f'{folder}/{file_name}_Comp_Summary.txt','w') as editor:
        headers.to_csv(editor,index=False,header=False,sep='\t') ; properties.to_csv(editor,index=False,sep='\t') 
    
    Text_xy = [int((map_data.shape[1]*10)/100),int((map_data.shape[0]*85)/100)]
    Beam_xy = [int((map_data.shape[1]*80)/100),int((map_data.shape[0]*20)/100)]
    BeamText_xy = [int((map_data.shape[1]*65)/100),int((map_data.shape[0]*12)/100)]
    colormap = 'CMRmap_r'
    dynamic_colormap = matplotlib.colormaps['tab20'] ; num_colors = dynamic_colormap.N
    
    #Figures Settings and Plotting
    plt.rc('font',size=15) ; matplotlib.rcParams['font.family'] = 'Times New Roman' ;matplotlib.rcParams['text.usetex'] = True
    L1 = det_limit/10 ; L2 = f_peak/2 ; params = {'legend.fontsize':15,'legend.handlelength':1} ; plt.rcParams.update(params)
    xticks,yticks,xlabels,ylabels = TicksEstimator(file)
    fig, axes = plt.subplots(1,4,sharey=True,gridspec_kw={'width_ratios':[1,1,1,1]})
    fig.subplots_adjust(wspace=0.00) , fig.set_size_inches([32,8])

    Data_to_Plot = [map_data.data,model_data.data,res_data.data,0] ; Text_to_Plot = [f'{date}\nObserved','Model','Residual',f'Components\nDetection Limit: {det_limit:.4f} Jy/beam']
    for w, axis in enumerate(axes,start=1):
        axis.set_xlabel('R.A. (J2000)')
        if w > 1:
            axis.set_xticks(xticks[1:]) ; axis.set_xticklabels(xlabels[1:],fontsize=13)
        else:
            axis.set_xticks(xticks); axis.set_xticklabels(xlabels,fontsize=13)
        axis.tick_params(which='major',top=True,right=True,direction='in',length=15,pad=8)
        axis.tick_params(which='minor',top=True,right=True,direction='in',length=5,pad=8)
        axis.minorticks_on() ; axis.grid(color='gray',linestyle='--',linewidth=0.8,alpha=0.4)
        
        beam = mpatches.Ellipse(xy=Beam_xy,width=(BMin*1000)/increment,height=(BMaj*1000)/increment,angle=BPA,edgecolor='black',fc='green',lw=1,alpha=0.8)
        axis.add_patch(beam)
        axis.text(BeamText_xy[0],BeamText_xy[1],f'{BMin*1000:.2f} mas x {BMaj*1000:.2f} mas',fontsize=13)
        if w != 3:
            axis.contour(map_data if w != 2 else model_data,levels=levels,colors='black',alpha=0.8,linewidths=1.00)
        axis.text(Text_xy[0],Text_xy[1],Text_to_Plot[w-1])
        
        if w <= len(Data_to_Plot)-1:
            im = axis.imshow(Data_to_Plot[w-1],origin='lower',cmap=colormap,norm=colors.LogNorm(vmin=L1, vmax=L2)) ; im.set_clim(L1,L2)
        if w == 1:
            axis.set_ylabel('Dec (J2000)')
        if w == 4:
            for i in range(len(Fitted['pix_x'])):
                color = dynamic_colormap(i % num_colors)
                knot = mpatches.Ellipse(xy=[Fitted['pix_x'][i],Fitted['pix_y'][i]],width=Summary['ConMin'][i]*1000/increment,height=Summary['ConMaj'][i]*1000/increment, angle=Summary['ConPA'][i],edgecolor=color,fc='None',lw=1.5, linestyle='--')
                axis.add_patch(knot) ; axis.scatter(Fitted['pix_x'][i],Fitted['pix_y'][i],marker='+',s=20,color=color,label='C'+str(i+1))
            axis.legend()
    for axis in axes:
        axis.set_yticks(yticks) ; axis.set_yticklabels(ylabels,fontsize=13)
    #Colorbar
    cbar_ax = fig.add_axes([0.905,0.13,0.012,0.735])
    cb = fig.colorbar(mpl.cm.ScalarMappable(cmap=colormap,norm=mpl.colors.Normalize(vmin=L1,vmax=L2)),cax=cbar_ax)
    cb.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    cb.ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    cb.ax.tick_params(direction='out') 
    cb.set_label('Flux (Jy/beam)',labelpad=10)
    plt.savefig(f'{folder}/{file_name}_Fit.png',bbox_inches='tight',dpi=300)
    plt.show()