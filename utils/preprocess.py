# Generic Python Modules
import numpy as np
import ants
from scipy import ndimage as ndi

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

class Preprocess:
    def __init__(self):
        pass

    def run(self, data, kernsz=3, qa=False, out_fpath=None):
        self._img_sort(data)
        self._gauss_smooth(data, kernsz)
        if qa:
            ants.from_numpy( data.imdata[:,:,:,0,:], spacing=data.params['spacing'], origin=data.params['origin'], direction=data.params['direction'] ).image_write(out_fpath+'/ase_preproc.nii.gz')
            if  data.params['nechoes']==2:
                ants.from_numpy( data.imdata[:,:,:,1,:], spacing=data.params['spacing'], origin=data.params['origin'], direction=data.params['direction'] ).image_write(out_fpath+'/ase_preproc_e2.nii.gz')
    def _img_sort(self, data):
        '''
        Sort ASE data by increasing tau (displacement of 𝜋 pulse) shift value given by tauvec.csv
        '''
        tauind_sort = np.argsort(data.params['tauvec'])
        tauvec_sort = np.array(data.params['tauvec'])[tauind_sort]
        data.imdata = data.imdata[:,:,:,:,tauind_sort]
        data.taudata = tauvec_sort
    
    def _gauss_smooth(self, data, kernsz):
        '''
        Smooth image data by using a Gaussian smoothing filter with defined kernel size (def=3)
        '''
        kernrad = int((kernsz-1)/2)
        data.imdata = ndi.gaussian_filter(data.imdata, sigma=(1,1,1,0,0), radius=kernrad, order=0)

# =========================================================================================================
class PreprocessANTs(Preprocess):
    def __init__(self):
        pass
    
    def run(self, data, kernsz=3, qa=False, out_fpath=None):
        self._img_sort(data)
        self._gauss_smooth(data, kernsz)
        self._motion_corr(data, out_fpath)
        if qa:
            ants.from_numpy( data.imdata[:,:,:,0,:], spacing=data.params['spacing'], origin=data.params['origin'], direction=data.params['direction'] ).image_write(out_fpath+'/ase_preproc.nii.gz')
            if  data.params['nechoes']==2:
                ants.from_numpy( data.imdata[:,:,:,1,:], spacing=data.params['spacing'], origin=data.params['origin'], direction=data.params['direction'] ).image_write(out_fpath+'/ase_preproc_e2.nii.gz')
            
            # self._motcor_vid(data.imdata, out_fpath)

    def _motion_corr(self, data, out_fpath):
        '''
        Motion correction using ANTsPy Rigid transformation to reference volume (def=0)
        Taken from ANTsPy github with adjustments to registration parameters
        (https://github.com/ANTsX/ANTsPy/blob/master/tutorials/motionCorrectionExample.ipynb)
        '''
        refvol = int( round( data.imdata.shape[-1] / 2 ) ) # Takes middle frame
        ref_ants = ants.from_numpy(data.imdata[:,:,:,0,refvol], spacing=data.params['spacing'][:-1], 
                                   origin=data.params['origin'][:-1], direction=data.params['direction'][:-1] )
        for ii in range ( data.params['nechoes'] ): # Re-iterate motion-correction over nEchoes - AKS
            motcor_list = list()
            for jj in range( data.imdata.shape[4] ):
                if jj != refvol:
                    temp_ants = ants.from_numpy(data.imdata[:,:,:,ii-1,jj], spacing=data.params['spacing'][:-1],
                                                 origin=data.params['origin'][:-1], direction=data.params['direction'][:-1])
                    areg = ants.registration( fixed=ref_ants, 
                                             moving=temp_ants, 
                                             type_of_transform="Rigid", 
                                             aff_random_sampling_rate=0.9,
                                             aff_shrink_factors=(1),
                                             aff_smoothing_sigmas=(0),
                                             aff_iterations=(100),
                                             initial_transform="identity"
                                             )
                    motcor_list.append( areg[ 'warpedmovout' ] )
                elif jj==refvol:
                    motcor_list.append( ref_ants )
            tempmc_ants = ants.make_image( data.imdata.shape[:-1], spacing=data.params['spacing'], origin=data.params['origin'], direction=data.params['direction'] )
            tempmc_ants = ants.list_to_ndimage(tempmc_ants, motcor_list)
            data.imdata[:,:,:,ii,:] = tempmc_ants.numpy()
    
    def _motcor_vid(self, data, out_fpath):
        c_map = cm.get_cmap('jet').copy()
        c_map.set_under(color='black')

        nn = int( round(data.imdata.shape[2] / 3) )
        s_max = round( np.max( data.imdata[:,:,:,0,0], axis=(0,1,2) ) * 0.9 , -3)

        fig, axes=plt.subplots()
        frames=[]
        for ii in range( int( data.imdata.shape[4]) ):
            axS = ndi.rotate( data.imdata[:,:,nn,0,ii], 90 )
            frames.append([axes.imshow(axS, cmap=c_map, vmin=0.01, vmax=s_max, aspect='equal')])
            if ii==0:
                axe=axes.imshow(axS, cmap=c_map, vmin=0.01, vmax=s_max, aspect='equal')
        axes.set_axis_off()
        fig.set_facecolor('black')

        cbar=fig.colorbar(axe,orientation="vertical", ax = axes, format='%.0e')
        cbar.set_label("Signal(a.u.)", color='white', size=12)
        cbar.ax.yaxis.set_tick_params(color='white')
        cbar.ax.tick_params(labelsize=12)
        cbar.outline.set_edgecolor('white')
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
        ani = animation.ArtistAnimation(fig, frames, interval=200, blit=False,
                                    repeat_delay=1000)
        ani.save(out_fpath+'/ase_motcorr_qa.gif')