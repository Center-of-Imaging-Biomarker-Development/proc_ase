from tqdm import trange
import time
import ants
from math import ceil
import numpy as np
from scipy import ndimage as ndi
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import random

from model.asefit import AseFit

class VoxelFit():
    def __init__(self, data):
        self._dims = tuple( int(i) for i in data.params['dimensions'][:3])
        self._r2prime_map = np.zeros(self._dims)
        self._rvcbv_map = np.zeros(self._dims)
        self._roef_map = np.zeros(self._dims)
        self._rsquared_map = np.zeros(self._dims)
        self._rand_vox = None
        self._plot_vox = None

    # ======================== Properties ======================== #
    @property
    def r2prime_map(self):
        return self._r2prime_map
    
    @r2prime_map.setter
    def r2prime_map(self,new_r2prime_map):
        self._r2prime_map = new_r2prime_map

    @property
    def rvcbv_map(self):
        return self._rvcbv_map
    
    @rvcbv_map.setter
    def rvcbv_map(self,new_rvcbv_map):
        self._rvcbv_map = new_rvcbv_map

    @property
    def roef_map(self):
        return self._roef_map
    
    @roef_map.setter
    def roef_map(self,new_roef_map):
        self._roef_map = new_roef_map

    @property
    def rsquared_map(self):
        return self._rsquared_map
    
    @rsquared_map.setter
    def rsquared_map(self,new_rsquared_map):
        self._rsquared_map = new_rsquared_map

    # ======================== Methods ======================== #

    def run(self, data, out_fpath):
        # ----- Prep ----- #
        ase_mask = self._mask_ase(data) # Create mask to remove unnecessary voxels (reduce proc. time)
        x,y,z = np.nonzero(ase_mask)
        maxlen = len(x)
        plot_vox = None # Choose random voxel for report
        rand_n = random.randrange(round(maxlen / 3), round(2 * maxlen / 3), 1)
        
        # ----- Voxelwise Model Fitting ----- #
        for n in trange(maxlen-1):
            i,j,k = int(x[n]), int(y[n]), int(z[n])
            sigvec = data.imdata[i,j,k,:,:]
            b0 = data.b0data[i,j,k]    
            voxfit = AseFit(data, sigvec, b0)
            voxfit.yh_1c_full(data)
            self._r2prime_map[i,j,k], self._rvcbv_map[i,j,k], self._roef_map[i,j,k], self._rsquared_map[i,j,k] = \
                voxfit.r2prime, voxfit.vcbv, voxfit.oef, voxfit.rsquared

            if n == rand_n:
                self._rand_vox = (int(x[n]), int(y[n]), int(z[n]))
                self._plot_vox = voxfit
        
        # ----- Save results ----- #
        self._plot_maps(out_fpath)
        self._plot_randomfit(data, out_fpath)
        self._save_maps(data, out_fpath)

    def _mask_ase(self, data):
        ase_ants = ants.from_numpy( data.imdata[:,:,:,0,0], 
                                   spacing=data.params['spacing'][:-1], 
                                   origin=data.params['origin'][:-1],
                                   direction=data.params['direction'][:-1, :-1])
        ase_mask = ants.get_mask(ase_ants)
        ase_mask = ase_mask.numpy()
        return ase_mask
    
    def _plot_maps(self, out_fpath):
        '''
        Function to plot the voxel-wise maps for QA/report
        
        :param self: Description
        :param out_fpath: Description
        '''
        turbo_map = cm.get_cmap('turbo').copy()
        turbo_map.set_under(color='black')

        image_list = [self.r2prime_map, self.rvcbv_map*100, self.roef_map*100, self.rsquared_map]
        image_title = [u"R$_2$'", u"vCBV", u"OEF", u"R$^2$"]
        image_range = [[0.001, 15], [0.001, 20], [0.001, 100], [0.001, 1]]
        image_cbarax = [[0,5,10,15], [0,5,10,15,20], [0,25,50,75,100], [0,0.25,0.50,0.75,1.00]]
        cbar_units = [u"s$^{-1}$", u"%", u"%", "unitless"]
        z_list = np.linspace(0, self._dims[2]-1, num=4)
        z_list = np.round(z_list)

        fig, axarr = plt.subplots(nrows=4, ncols=4, figsize=(6,6), constrained_layout=True)
        for i in range(4):
            i_map = image_list[i]
            i_range = image_range[i]

            for j in range(4):
                z = int(z_list[j])
                im = axarr[i,j].imshow(ndi.rotate(i_map[:,:,z], 90),
                                       cmap = turbo_map,
                                       interpolation='None',
                                       vmin=i_range[0],
                                       vmax=i_range[1]
                                       )
                if i == 0:
                    axarr[i,j].set_title(u"z={}".format(z), color='white')
                if j == 3:
                    divider = make_axes_locatable(axarr[i,j])
                    cax = divider.append_axes('right', size='5%', pad=0.05)
                    cbar = fig.colorbar(im, cax=cax, orientation='vertical', ticks=image_cbarax[i])
                    cbar.ax.yaxis.set_tick_params(color='white')
                    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
                    cbar.set_label(cbar_units[i], color='white')
                if  j == 0:
                    axarr[i,j].set_ylabel(image_title[i], rotation=0, color='white', size=12)
                else:
                    axarr[i,j].axis('off')
        fig.set_facecolor('black')

        # ----- Save figure ----- #
        fig.savefig(out_fpath + '/voxelmaps_qa.png', dpi=300)

    def _plot_randomfit(self, data, out_fpath):
        tau0_z1 = data.imdata[:, :, self._rand_vox[2], 0, 0]

        fig, axarr = plt.subplots(nrows=1, ncols=2, figsize=(6,3), constrained_layout=True)
        im=axarr[0].imshow( ndi.rotate(tau0_z1,90), cmap='gray', interpolation='None' )
        axarr[0].plot( self._rand_vox[0], self._dims[1]-self._rand_vox[1], '.', markersize=8, color='red' )
        axarr[0].title.set_text(u"Voxel {},{},{}".format(self._rand_vox[0], self._rand_vox[1], self._rand_vox[2]))
        axarr[0].axis('off')

        self._plot_vox.plot_fit(axarr[1], data, u"Voxel Model Fitting")

        # ----- Save figure ----- #
        fig.savefig(out_fpath+'/voxelfit_qa.png', dpi=300)
    
    def _save_maps(self, data, out_fpath):
        map_list = ['R2prime', 'rvCBV', 'rOEF', 'Rsquared']
        data_list = [self._r2prime_map, self._rvcbv_map, self._roef_map, self._rsquared_map]
        for i in range(len(map_list)):
            temp_ants = ants.from_numpy( data_list[i], spacing=data.params['spacing'][:-1], origin=data.params['origin'][:-1],
                                         direction=data.params['direction'][:-1, :-1] )
            temp_ants.image_write(out_fpath+'/'+map_list[i]+'.nii.gz')