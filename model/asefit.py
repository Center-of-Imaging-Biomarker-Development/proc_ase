# Generic Python Modules
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

class AseFit():
    def __init__(self, data, sigvec, b0, taucut=0.010): # initialize variables to be used for fitting
        self._taucut=taucut
        self._b0=b0
        self._taudata = data.taudata
        self._dtedata = data.dtedata
        self._sigvec = sigvec
        self._sigvec_std = None
        self._clean_sigvec(data)

        self._longtau_fit={}
        self._shorttau_fit={}

        self._r2prime=None
        self._r2star=None
        self._vcbv=None
        self._oef=None
        self._rsquared=None

    # ======================== Properties ======================== #
    @property
    def r2prime(self):
        return self._r2prime
    
    @r2prime.setter
    def r2prime(self, new_r2prime):
        self._r2prime = new_r2prime

    @property
    def r2star(self):
        return self._r2star
    
    @r2star.setter
    def r2star(self, new_r2star):
        self._r2star = new_r2star

    @property
    def vcbv(self):
        return self._vcbv
    
    @vcbv.setter
    def vcbv(self, new_vcbv):
        self._vcbv = new_vcbv

    @property
    def oef(self):
        return self._oef
    
    @oef.setter
    def oef(self, new_oef):
        self._oef = new_oef

    @property
    def rsquared(self):
        return self._rsquared

    @rsquared.setter
    def rsquared(self, new_rsquared):
        self._rsquared = new_rsquared


    # ======================== Methods ======================== #
    def longtau_fx(self, X, Cl, R2star, R2prime):
        '''
        Long tau function
        '''
        tau= X[:,0]
        dTE = X[:,1]
        return Cl - R2star*dTE - R2prime*(2*tau)

    def shorttau_fx(self, X, Cs, R2primeSqDivLambda):
        '''
        Short tau function
        '''
        tau= X[:,0]
        dTE = X[:,1]
        return Cs - 0.3 * R2primeSqDivLambda * ( (2*tau)**2 )

    def yh_1c_full(self, data):
        '''
        Yablonskiy-Haacke 1-compartment ASE model - Full tau regime
        '''
        
        longtauinds = np.argwhere(self._taudata >= self._taucut).flatten()
        shorttauinds = np.argwhere(self._taudata < self._taucut).flatten()

        self._longtau_fit['sigdata'] = np.log( self._sigvec[longtauinds] )
        self._shorttau_fit['sigdata'] = np.log( self._sigvec[shorttauinds] )
        
        try:
            # ----- Long Tau Regime Fit ------ #
            self._longtau_fit['te1_ind'] = np.argwhere(self._dtedata[longtauinds]==0).flatten()
            self._longtau_fit['te2_ind'] = np.argwhere(self._dtedata[longtauinds]!=0).flatten()

            self._longtau_fit['xdata'] = np.empty( [len(longtauinds), 2] )
            self._longtau_fit['xdata'][:,0] = self._taudata[ longtauinds ]
            self._longtau_fit['xdata'][:,1] = self._dtedata[ longtauinds ]

            self._longtau_fit['p0'] = [ self._longtau_fit['sigdata'][0]*1.02, 50, 20 ] # Define guesses for parameters - [Cl, R2star, R2prime]
            self._longtau_fit['lb'] = [ 0, 10, 1 ] # Define lower bounds - [Cl, R2star, R2prime]
            self._longtau_fit['ub'] = [ self._longtau_fit['sigdata'][0]*1.15, 200, 150 ] # Define upper bounds - [Cl, R2star, R2prime]

            self._longtau_fit['popt'], self._longtau_fit['pcov'] = curve_fit( self.longtau_fx, self._longtau_fit['xdata'], self._longtau_fit['sigdata'],
                                                                            p0=self._longtau_fit['p0'], bounds=(self._longtau_fit['lb'], self._longtau_fit['ub']) )
            self._longtau_fit['sigma_ab'] = np.sqrt(np.diagonal(self._longtau_fit['pcov']))

            # ----- Short Tau Regime Fit  ----- #
            self._shorttau_fit['te1_ind'] = np.argwhere(self._dtedata[shorttauinds]==0).flatten()
            self._shorttau_fit['te2_ind'] = np.argwhere(self._dtedata[shorttauinds]!=0).flatten()
            shorttauinds=shorttauinds[self._shorttau_fit['te1_ind']]
            self._shorttau_fit['sigdata'] = self._shorttau_fit['sigdata'][self._shorttau_fit['te1_ind']]

            self._shorttau_fit['xdata'] = np.empty( [len(shorttauinds), 2] )
            self._shorttau_fit['xdata'][:,0] = self._taudata[ shorttauinds ]
            self._shorttau_fit['xdata'][:,1] = self._dtedata[ shorttauinds ]

            self._shorttau_fit['p0'] = [ self._shorttau_fit['sigdata'][0]*1.02, 8000 ] # Define guesses for parameters - [Cl, R2prime^2/lambda]
            self._shorttau_fit['lb'] = [ 0, 6.4 ] # Define lower bounds - [Cl, R2prime^2/lambda]
            self._shorttau_fit['ub'] = [ self._shorttau_fit['sigdata'][0]*1.15, 2e6 ] # Define upper bounds - [Cl, R2prime^2/lambda]

            self._shorttau_fit['popt'], self._shorttau_fit['pcov'] = curve_fit( self.shorttau_fx, self._shorttau_fit['xdata'], self._shorttau_fit['sigdata'],
                                                                              p0 = self._shorttau_fit['p0'], bounds=[self._shorttau_fit['lb'], self._shorttau_fit['ub']] )
            self._shorttau_fit['sigma_ab'] = np.sqrt( np.diagonal(self._shorttau_fit['pcov']) )

            # ----- Parameter Quantification ----- #
            # R2prime
            self.r2prime = self._longtau_fit['popt'][2]
            self.r2star = self._longtau_fit['popt'][1]

            # vCBV
            lnSl_extrap = self._longtau_fit['popt'][0]
            lnSs0 = self._shorttau_fit['popt'][0]
            self.vcbv = lnSl_extrap - lnSs0

            # OEF
            self.oef = self.r2prime / ( data.params['gamma'] * self.vcbv * (4/3) * np.pi * data.params['dchi0'] * data.params['microhct'] * self._b0 )

            # Rsquared
            residuals = np.zeros( len(self._sigvec) )
            residuals[shorttauinds] = self._shorttau_fit['sigdata'] - self.shorttau_fx(self._shorttau_fit['xdata'], *self._shorttau_fit['popt'])
            residuals[longtauinds] =  self._longtau_fit['sigdata'] - self.longtau_fx( self._longtau_fit['xdata'], * self._longtau_fit['popt'])
            ss_res = np.sum( residuals**2 )
            ss_tot = np.sum( ( np.log(self._sigvec) - np.mean( np.log(self._sigvec) ) ) **2 )
            self.rsquared = 1 - ( ss_res / ss_tot )

        except:
            # ----- Return zeros if fit doesn't work ----- #
            self.r2prime = 0
            self.r2star = 0
            self.vcbv = 0
            self.oef = 0
            self.rsquared = 0

    def plot_fit(self, ax, data, roi_name, save_fig=False, out_fpath=None):
        # ----- Plot short tau fit ----- #
        ax.errorbar(2*self._shorttau_fit['xdata'][:,0], self._shorttau_fit['sigdata'][self._shorttau_fit['te1_ind']], marker='o', color='blue', label='Short Tau (obs.)', ls='none')
        ax.plot(2*self._shorttau_fit['xdata'][:,0], self.shorttau_fx(self._shorttau_fit['xdata'], *self._shorttau_fit['popt']), 'b-', label='Short Tau (fit.)')
        bound_upper = self.shorttau_fx(self._shorttau_fit['xdata'], *(self._shorttau_fit['popt'] + self._shorttau_fit['sigma_ab']))
        bound_lower = self.shorttau_fx(self._shorttau_fit['xdata'], *(self._shorttau_fit['popt'] - self._shorttau_fit['sigma_ab']))
        ax.fill_between(2*self._shorttau_fit['xdata'][:,0], bound_lower, bound_upper, color='blue', alpha=0.15)

        # ----- Plot Long tau fit ----- #
        ax.errorbar(2*self._longtau_fit['xdata'][self._longtau_fit['te1_ind'],0], self._longtau_fit['sigdata'][self._longtau_fit['te1_ind']], marker='o', color='red', label='Long Tau (obs.)', ls='none')
        tau_extrap = np.arange(0, np.max(self._taudata), 0.001)
        te_extrap = np.full(tau_extrap.shape, self._dtedata[self._longtau_fit['te1_ind'][0]])
        x_extrap = np.zeros((len(tau_extrap),2))
        x_extrap[:,0] = tau_extrap
        x_extrap[:,1] = te_extrap
        ax.plot(2*x_extrap[:,0], self.longtau_fx(x_extrap, *self._longtau_fit['popt']), 'r-', label='Long Tau (fit.)')
        bound_upper = self.longtau_fx(self._longtau_fit['xdata'][self._longtau_fit['te1_ind']], *(self._longtau_fit['popt'] + self._longtau_fit['sigma_ab']))
        bound_lower = self.longtau_fx(self._longtau_fit['xdata'][self._longtau_fit['te1_ind']], *(self._longtau_fit['popt'] - self._longtau_fit['sigma_ab']))
        ax.fill_between(2*self._longtau_fit['xdata'][self._longtau_fit['te1_ind'],0], bound_lower, bound_upper, color='red', alpha=0.15)

        if data.params['nechoes']==2:
            ax.errorbar(2*self._longtau_fit['xdata'][self._longtau_fit['te2_ind'],0], self._longtau_fit['sigdata'][self._longtau_fit['te2_ind']], marker='o', color='green', label='Long Tau-e2 (obs.)', ls='none')
            bound_upper = self.longtau_fx(self._longtau_fit['xdata'][self._longtau_fit['te2_ind']], *(self._longtau_fit['popt'] + self._longtau_fit['sigma_ab']))
            bound_lower = self.longtau_fx(self._longtau_fit['xdata'][self._longtau_fit['te2_ind']], *(self._longtau_fit['popt'] - self._longtau_fit['sigma_ab']))
            ax.fill_between(2*self._longtau_fit['xdata'][self._longtau_fit['te2_ind'],0], bound_lower, bound_upper, color='green', alpha=0.15)

        # ----- Misc ----- #
        ax.set_xlabel(r'2$\tau$ (s)')
        ax.set_ylim(np.log(np.min(self._sigvec)*0.95), (self._longtau_fit['popt'][0]*1.01))
        ax.set_ylabel('ln(Signal) (au)')
        ax.legend(loc='upper right', ncols=2, fontsize=6)

        textstr = '\n'.join((
            '$R_2\' =%.2f$' % (self.r2prime, ),
            'vCBV=%.2f' % (self.vcbv*100, ),
            'OEF=%.2f' % (self.oef*100, ),
            '$R^2$=%.2f' % (self.rsquared, )))

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        # place a text box in bottom left in axes coords
        ax.text(0.05, 0.35, textstr, transform=ax.transAxes, fontsize=8,
                verticalalignment='top', bbox=props)
        ax.set_title(roi_name)
        
        if save_fig:
            plt.tight_layout()
            fig_fname = out_fpath + '/' + roi_name + '_yh_1c_full_fit.pdf'
            plt.savefig(fig_fname, dpi=300)

    def _clean_sigvec(self,data):
        '''
        Remove indices where signal(au) = 0 and averages tau signal data where more than 1 dynamic collected
        '''
        ind_rm = np.argwhere( (self._sigvec.flatten() == 0) | (np.isin(self._taudata, data.params['rm_tau'])) ).flatten()
        taudata_rm = np.delete(self._taudata, ind_rm)
        dtedata_rm = np.delete(self._dtedata, ind_rm)
        sigvec_rm = np.delete(self._sigvec.flatten(), ind_rm)

        if (  data.params['nechoes']==1 & len(np.unique(taudata_rm)) != len(taudata_rm)): # Check if there are any repeat tau shifts
            taudata_m = np.unique(taudata_rm)
            dtedata_m = np.zeros( taudata_m.shape )
            sigvec_m = np.zeros( taudata_m.shape )
            sigvec_std = np.zeros( taudata_m.shape )

            for i in range( len(taudata_m) ):
                temp_ind =  np.argwhere( taudata_rm == taudata_m[i] ).flatten()
                dtedata_m[i] = np.mean( dtedata_rm[temp_ind] )
                sigvec_m[i] = np.mean(sigvec_rm[temp_ind])
                sigvec_std[i] = np.std(sigvec_rm[temp_ind])
            
            self._taudata = taudata_m
            self._dtedata = dtedata_m
            self._sigvec = sigvec_m
            self._sigvec_std = sigvec_std
        
        else:
            self._taudata = taudata_rm
            self._dtedata = dtedata_rm
            self._sigvec = sigvec_rm


    