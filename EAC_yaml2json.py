import json
import yaml
import os
from hwo_sci_eng.utils.read_yaml import *
import astropy.units as u
from synphot import SpectralElement
from specutils import Spectrum1D
import numpy as np
import pandas as pd

class EAC:
    def __init__(self, eac_num=1):
        if eac_num == 1:
            self.eac = eac1()
        elif eac_num == 2:
            self.eac = eac2()
        elif eac_num == 3:
            self.eac = eac3()
        else:
            raise ValueError('EAC number must be 1, 2, or 3')
        
        # Read in the CI, HRI, and UVI yaml files
        self.ci = ci()
        self.hri = hri()
        self.uvi = uvi()

        # Read in the template EAC JSON file to update
        # From Armen Tokadjian (A.T.; NASA JPL)
        with open('./data/EAC1_EXOSIMS_Template.json', 'r') as f:
            self.eac_json = json.load(f)


        self.optical_path = self.ci['opticalpath']['full_path']

        # Get path to this directory
        self.path_dir = os.path.dirname(os.path.realpath(__file__))


    def read_refl(self, refl_path):
        with open(os.getenv('SCI_ENG_DIR') + '/obs_config/' + refl_path[2:], 'r') as f:
            refl_file = yaml.safe_load(f)
        return refl_file
    
    def read_det(self, det_path):
        """Reads in the detector yaml file and saves it as a csv file while returning the relative path

        Parameters
        ----------
        det_path : str
            Path to the detector yaml in the hwo_sci_eng directory

        Returns
        -------
        str
            Relative path to the csv file
        """
        with open(os.getenv('SCI_ENG_DIR') + '/obs_config/' + det_path[2:], 'r') as f:
            det_file = yaml.safe_load(f)
        # Need to save this as a csv file
        wl = det_file['wavelength'] # Should be in nm
        qe = det_file['QE']
        det_file = pd.DataFrame({'lambda': wl, 'QE': qe})
        csv_path = './data'+det_path[2:-4]+'csv'
        det_file.to_csv(csv_path, index=False)

        return csv_path[1:]

    
    def optics(self, start_i = 0, end_i = -1, exc_modules = ['unitless']):
        
        optical_path = self.optical_path[start_i:end_i]

        # Remove excluded modules
        optical_path = [module for module in optical_path if module not in exc_modules]

        # First two elements are the Primary and Secondary mirrors, which are housed in eac
        eac_reflectivity = [self.read_refl(self.eac[module]['reflectivity']) for module in optical_path[:2] if 'reflectivity' in self.eac[module]]

        # Loop through each optical element, and if there is a reflectance file, read it
        ci_reflectivity = [self.read_refl(self.ci[module]['reflectivity']) for module in optical_path[2:] if 'reflectivity' in self.ci[module]]

        # Add in the transmission of the filters
        ci_transmission = [self.read_refl(self.ci[module]['transmission']) for module in optical_path[2:] if 'transmission' in self.ci[module]]

        # Concatenate the two lists
        optical_path_refl = eac_reflectivity + ci_reflectivity + ci_transmission[1:] # The first element is from beam splitter, which is already included in the reflectivity

        # For each optical element, make sure they are sorted by wavelength
        for i, refl in enumerate(optical_path_refl):
            optical_path_refl[i]['wavelength'], optical_path_refl[i]['reflectivity'] = zip(*sorted(zip(refl['wavelength'], refl['reflectivity']), key=lambda x: x[0]))

        # Create Spectrum1D objects for each optical element
        optical_path_refl = [SpectralElement.from_spectrum1d(Spectrum1D(flux=refl['reflectivity']*u.dimensionless_unscaled, spectral_axis=refl['wavelength']*u.nm)) for refl in optical_path_refl]

        # Multiply the reflectance of each optical element to get the total reflectance
        full_optical_path = optical_path_refl[0]
        for refl in optical_path_refl[1:]:
            full_optical_path = full_optical_path * refl

        self.full_optical_path = full_optical_path

    def update_imager(self, imager_key='imaging_EAC1_BroadbandVisible_500'):
        """Converts information about the Broad Band Imager to JSON format
        """

        # Grab the instrument with the assigned name
        imager = [instrument for instrument in self.eac_json['scienceInstruments'] if imager_key in instrument['name']]

        # Make sure there is only one instrument with the given name
        assert len(imager) == 1, 'Instrument name not found in scienceInstruments or duplicate instruments found'

        imager = imager[0]

        # Optical path at 500 nm
        optics_500 = self.full_optical_path(500*u.nm)
        imager['optics'] = optics_500.value

        # Pixel size/pitch (in meters)
        imager['pixelSize'] = (self.ci['Visible_Channels']['Detectors']['Broadband_Imager']['pixel_pitch'][0] * u.micron).to(u.m).value

        # Pixel number (square detector)
        imager['pixelNumber'] = self.ci['Visible_Channels']['Detectors']['Broadband_Imager']['detector_format'][0][0]

        # Dark Current (units of e-/pix/s)
        imager['idark'] = self.ci['Visible_Channels']['Detectors']['Broadband_Imager']['DC'][0]

        # Read noise (units of e-/pix)
        imager['sread'] = self.ci['Visible_Channels']['Detectors']['Broadband_Imager']['RN'][0]

        # Following three entries are not in yaml, assuming values from
        # https://github.com/spacetelescope/hwo-tools/blob/main/coron_imaging/exosims_wrapper/exosims_sample_parameters.yml 
        # Clock Induced Charge (units of e-/pix)
        # XXX: This is assumed, do not have a value from the yaml
        imager['CIC'] = 1.3e-3

        # Photon counting efficiency
        # XXX: Assuming 100% efficiency
        imager['PCeff'] = 1

        # Excess noise factor
        # XXX: Assuming 1
        imager['ENF'] = 1

        # Quantum efficiency
        imager['QE'] = self.path_dir + self.read_det(self.ci['Visible_Channels']['Detectors']['Broadband_Imager']['QE'])

        # Exposure time (in seconds)
        imager['texp'] = 1000

        # Pixel/plate scale (comes in mas, need to be in as)
        imager['pixelScale'] = (self.ci['Visible_Channels']['Performance']['Broadband_Imager']['plate_scale'][0] * u.mas).to(u.arcsec).value

        # Half Field of View (in arcseconds)
        # Full FoV = 2 * arctan(npix*dpix/2f) -- dpix and f are in meters, output will be in radians
        focal_length = imager['pixelSize'] / (2 * np.tan(imager['pixelScale']*u.arcsec/2)).value # in meters

        imager['FoV'] = (np.arctan(imager['pixelNumber'] * imager['pixelSize'] \
                                / 2 / focal_length) * u.rad).to(u.arcsec).value
        

        print(imager)

        return
    
    def update_ifs(self, ifs_key="spectro_EAC1_IFS_910", wb='IR'):

        ifs = [instrument for instrument in self.eac_json['scienceInstruments'] if ifs_key in instrument['name']]

        # Make sure there is only one instrument with the given name
        assert len(ifs) == 1, 'Instrument name not found in scienceInstruments or duplicate instruments found'

        ifs = ifs[0]

        # Optical path at 1000 nm
        optics_1000 = self.full_optical_path(1000*u.nm)
        ifs['optics'] = optics_1000.value

        # Pixel size/pitch (in meters)
        ifs['pixelSize'] = (self.ci['NIR_Channels']['Detectors']['IFS']['pixel_pitch'][0] * u.micron).to(u.m).value

        # Pixel number (square detector)
        ifs['pixelNumber'] = self.ci['NIR_Channels']['Detectors']['IFS']['detector_format'][0][0]

        # Dark Current (units of e-/pix/s)
        ifs['idark'] = self.ci['NIR_Channels']['Detectors']['IFS']['DC'][0]

        # Read noise (units of e-/pix)
        ifs['sread'] = self.ci['NIR_Channels']['Detectors']['IFS']['RN'][0]

        # Quantum efficiency
        ifs['QE'] = self.path_dir + self.read_det(self.ci['NIR_Channels']['Detectors']['IFS']['QE'])

        # Exposure time (in seconds)
        ifs['texp'] = 1000

        # Spectral resolution
        if wb == 'IR':
            ifs['Rs'] = 70 
        elif wb == 'VIS':
            ifs['Rs'] = 140

        return

    
    def update_starlightSuppressionSystems(self, system_key='CEC2a'):

        # Grab the instrument with the assigned name
        system = [instrument for instrument in self.eac_json['starlightSuppressionSystems'] if system_key in instrument['name']]

        # Make sure there is only one instrument with the given name
        assert len(system) == 1, 'Instrument name not found in starlightSuppressionSystems or duplicate instruments found'

        system = system[0]

        # Plate scale (in milliarcsec)
        system['core_platescale_units'] = 'mas'
        system['core_platescale'] = self.ci['Visible_Channels']['Performance']['Broadband_Imager']['plate_scale'][0]

        # Core area
        # Compute photometric aperture area pi*phot_aper**2 where phot_aper = cirum_diam/inscr_diam * root(2)/2 * lambda/D
        # In units of lambda/D
        core_area = np.pi/2 * (self.eac['PM']['circumscribing_diameter'][0] / self.eac['PM']['inscribing_diameter'][0])**2
        system['core_area'] = core_area

        # Files for core_throughput and core_mean_intensity and occ_trans
        system['core_thruput'] = self.path_dir + '/data/core_thruput.fits'
        system['core_mean_intensity'] = self.path_dir + '/data/core_mean_intensity.fits'
        system['occ_trans'] = self.path_dir + '/data/occ_trans.fits'

        print(system)
        return

# Main
if __name__ == '__main__':
    eac = EAC(1)
    eac.optics(0, -1, ['unitless', 'filters', 'Detector', 'Apodizer', 'Focal_Plane_Mask', 'Lyot_Stop'])
    eac.update_imager()
    eac.update_starlightSuppressionSystems()
    eac.update_ifs()