from fpdf import FPDF
from datetime import datetime
from PIL import Image
import os
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import ndimage
import numpy as np
import pandas as pd

def _gen_input(img_fpath, img_seq, hct, out_fpath, qa, e2_fpath, b0_fpath, yml_fpath, rm_tau, kernsz, motcor):
    in_str = "python proc_ase.py --img_fpath {} --img_seq {} --out_fpath {} --hct {} --kern_sz {}".format(img_fpath, img_seq, out_fpath, hct, kernsz)
    
    if qa != False:
        in_str = in_str + " --qa"
    if motcor != False:
        in_str = in_str + " --motcor"
    if e2_fpath != None:
        in_str = in_str + " --e2_fpath {}".format(e2_fpath)
    if b0_fpath != None:
        in_str = in_str + " --b0_fpath {}".format(b0_fpath)
    if yml_fpath != None:
        in_str = in_str + " --yml_fpath {}".format(yml_fpath)
    if rm_tau != None:
        in_str = in_str + " --rm_tau {}".format(rm_tau)

    return in_str

class report(FPDF):
    def header(self):
        self.image(os.path.join(os.path.dirname(__file__), "..","etc", "Donahue_logo.png"), 165, 8, 30)
        self.set_font('Arial', 'B', 14)
        self.cell(80)
        self.text(10, 15, 'Report ASE v0.0.1')
        self.line(10,18,160,18)
        self.ln(20)

    def body(self, params, in_str, voxel_maps, voxel_fit):
        font_size = 9
        self.set_font('Arial', '', font_size)
        now = datetime.now()
        strtime = now.strftime("%m-%d-%Y")

        # Basic information 
        self.text(10, 25, 'Date of report: {}'.format(strtime))
        self.text(10, 30, 'Scan ID: {}'.format(params['mrID']))
        self.text(10, 35, 'ASE Sequence: {}'.format(params['img_seq']))
        self.text(10, 40, 'Hematocrit: {}'.format(params['macrohct']))
        self.text(10, 45, 'Command line input:')
        self.ln(17.5)

        # Create code block for input from command line
        self.set_font('Courier', '', 8)
        self.set_text_color(0, 0, 0)
        self.set_fill_color(240, 240, 240)
        self.multi_cell(0, 5, in_str, 0, 'L', 1)
        self.ln(5)
        y = self.get_y()

        # Example voxel fit
        size_of_image = 150
        x=((self.w-size_of_image)/2)
        self.set_font('Arial', '', 9)
        self.set_fill_color(255, 255, 255)
        self.text(10,y,u'Example voxel fit:')
        self.ln(1)
        self.image(voxel_fit, x=((self.w-size_of_image)/2), w=size_of_image)
        self.ln(5)
        y = self.get_y()

        # Voxel-wise map results
        size_of_image = 120
        self.text(10,y, 'Voxel-wise maps:')
        self.ln(1)
        self.image(voxel_maps,x=((self.w-size_of_image)/2), w=size_of_image)
     
def generate_outputs(outpath, name, data, in_str, voxel_maps, voxel_fit):
    pdf = report()
    pdf.set_title(name)
    pdf.add_page()

    pdf.body(data.params, in_str, voxel_maps, voxel_fit)
    pdf.output('{}/report_{}.pdf'.format(outpath, name), 'F')    