# Generic Python Modules
import os
import yaml
import ants
import numpy as np

class AseData:
    '''
    Base class to initialize ASE data object
    '''
    def __init__(self, img_fpath, img_seq, hct, rm_tau, out_fpath=None, e2_fpath=None, b0_fpath=None, yml_fpath=None):
        self._params=None
        
        self._taudata=None
        self._tedata=None
        self._dtedata=None

        self._imdata=None
        self._b0data=None

        self._load_yaml(img_fpath, img_seq, hct, rm_tau, e2_fpath, yml_fpath)
        os.makedirs(out_fpath+'/'+self._params['mrID'], exist_ok=True)
    # ======================== Methods ======================== #
    def _findid(self,img_fpath):
        a = img_fpath.rfind('/')
        if 'e1' in img_fpath:
            mrid = img_fpath[a+1:-7]
        else:
            mrid = img_fpath[a+1:-4]
        return mrid

    def _load_yaml(self, img_fpath, img_seq, hct, rm_tau, e2_fpath, yml_fpath):
        if img_seq == 'ASE_HARMON':
            yaml_fpath = os.path.join(os.path.dirname(__file__), "..", "config", "ASE_HARMON.yaml")
        elif img_seq == 'Other':
            yaml_fpath = yml_fpath
        
        with open(yaml_fpath, 'r') as f:
            self._params = yaml.load(f, Loader=yaml.SafeLoader)
        
        self._params['img_seq'] = img_seq
        self._params['macrohct'] = hct
        self._params['microhct'] = hct * 0.85 # macro-to-microvascular hct ratio (~0.85)
        self._params['mrID'] = self._findid(img_fpath)
        self._params['ase_fpath'] = img_fpath
        self._params['rm_tau'] = np.array(rm_tau)
        
        if rm_tau != None:
            self._params['rm_tau'] = self._params['rm_tau']/1000 # convert ms to s
        if e2_fpath is None:
            self._params['nechoes'] = 1

    def set_fitdata(self):
        '''
        Set tau and dTE data from params for fitting functions
        '''
        self._tedata = np.ones( int(self._params['dimensions'][3]), dtype=np.float32 ) * self._params['te1']
        if self._params['nechoes']==2:
            self._taudata = np.append(self._taudata, self._taudata)
            self._tedata = np.append( self._tedata, np.ones( int(self._params['dimensions'][3]), dtype=np.float32 ) * self._params['te2'] )
        
        self._taudata = self._taudata / 1000 # convert ms to s
        self._tedata = self._tedata / 1000 # convert ms to s
        self._dtedata = self._tedata - np.min( self._tedata )

    # ======================== Properties ======================== #
    @property
    def taudata(self):
        return self._taudata
    
    @taudata.setter
    def taudata(self, new_taudata):
        self._taudata = new_taudata

    @property
    def tedata(self):
        return self._tedata
    
    @property
    def dtedata(self):
        return self._dtedata
        
    @property
    def params(self):
        return self._params

    @params.setter
    def params(self, new_params):
        self._params = new_params

    @property
    def imdata(self):
        return self._imdata
    
    @imdata.setter
    def imdata(self, new_imdata):
        self._imdata = new_imdata
    
    @property
    def b0data(self):
        return self._b0data
    
    @b0data.setter
    def b0data(self, new_b0data):
        self._b0data = new_b0data

# =========================================================================================================
class AseDataNifti1(AseData):
    '''
    Subclass to load NIFTI image data
    '''
    def __init__(self, img_fpath, img_seq, hct, rm_tau, out_fpath=None, e2_fpath=None, b0_fpath=None, yml_fpath=None):
        super().__init__(img_fpath, img_seq, hct, rm_tau, out_fpath, e2_fpath, b0_fpath, yml_fpath)
        self._load_img(img_fpath, out_fpath, e2_fpath, b0_fpath, img_seq)

    def _load_img(self, img_fpath, out_fpath, e2_fpath, b0_fpath, img_seq):
        self._params['img_format'] = 'NIFTI'

        # Image Header
        img_hdr = ants.image_header_info(img_fpath)
        self._params['dimensions'] = img_hdr['dimensions']
        self._params['spacing'] = img_hdr['spacing']
        self._params['origin'] = img_hdr['origin']
        self._params['direction'] = img_hdr['direction']
        
        # Image numpy
        img_ants = ants.image_read(img_fpath)
        (nX, nY, nZ, nDyn) = img_ants.shape
        self._imdata = np.zeros( ( nX, nY, nZ, self._params['nechoes'], nDyn ), dtype=np.float32 )
        self._imdata[:,:,:,0,:] = img_ants.numpy()
        if self._params['nechoes'] == 2:
            self._imdata[:,:,:,1,:] = ants.image_read(e2_fpath).numpy()
        
        if b0_fpath == None:
            self._b0data = np.ones_like( self._imdata[:,:,:,0,0] ) * self._params['B0']
        else:
            self._b0data = ants.image_read( b0_fpath ).numpy()