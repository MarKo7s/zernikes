from mathutils import cart2pol
import numpy as np
from math import factorial
import numexpr as ne

try:
    from .utils import (
        default_coefficients_dict,
        dict_to_zernike_coefficients,
        load_coefficients_json,
        save_coefficients_json,
        zernike_coefficients_to_dict,
    )
except ImportError:
    from utils import (
        default_coefficients_dict,
        dict_to_zernike_coefficients,
        load_coefficients_json,
        save_coefficients_json,
        zernike_coefficients_to_dict,
    )

class Zernikes:

    def __init__(self, max_zernike_radial_number = 6, mask_size_x_in_pixels = 1024, mask_size_y_in_pixels = 1024, aperture_diameter_in_m = 7.5e-3,
                 SLM_pixel_size_in_um = 9.2e-6, coefs_path = None):

        self.max_zernike_radial_number = max_zernike_radial_number

        # mask midpoint
        half_width = mask_size_x_in_pixels/2.0
        half_height = mask_size_y_in_pixels/2.0

        # Makes a polar grid with r = 1 at r_norm (in pixels). This properly fits the Zernike to the required aperture
        x,y = np.meshgrid(np.linspace(0, mask_size_x_in_pixels - 1, mask_size_x_in_pixels, dtype = np.float32),\
                        np.linspace(0, mask_size_y_in_pixels - 1, mask_size_y_in_pixels, dtype = np.float32))
        rho, phi = cart2pol(x - half_width, y - half_height)

        self.N = rho.shape[0] # Defines the array size on which we calculate the Zernikes (same as mask_size_x_in_pixels & mask_size_y_in_pixels if square)

        r_norm = aperture_diameter_in_m/SLM_pixel_size_in_um/2
        rho = rho/r_norm  

        self.zernike_basis_array = self.generate_zernike_basis(rho, phi, max_zernike_radial_number, mask_size_x_in_pixels, mask_size_y_in_pixels)
        
        if coefs_path is not None:
            self.load_coefs(coefs_path)
        else:
            self.clean()

        self.modify_zernike_coefs()
        self.make_zernike_fields()
        
        

    @staticmethod
    def _radial_coefficients(n, m_abs):
        coeffs = []
        for k in range((n - m_abs) // 2 + 1):
            c = ((-1)**k) * factorial(n - k) / (
                factorial(k) * factorial((n + m_abs) // 2 - k) * factorial((n - m_abs) // 2 - k)
            )
            coeffs.append((np.float32(c), n - 2 * k))
        return coeffs

    def _build_rho_powers(self, rho, max_radial_idx):
        max_power = 2 * max_radial_idx
        rho_powers = [np.ones_like(rho, dtype=np.float32)]
        for power in range(1, max_power + 1):
            rho_powers.append(rho ** power)
        return rho_powers

    def _compute_radial_nm(self, rho, rho_powers, n, m_abs, radial_coeffs_cache):
        r_nm = np.zeros((self.N, self.N), dtype=np.float32)
        key = (n, m_abs)
        if key not in radial_coeffs_cache:
            radial_coeffs_cache[key] = self._radial_coefficients(n, m_abs)
        for coeff, power in radial_coeffs_cache[key]:
            if power == 0:
                r_nm += coeff
            else:
                r_nm += coeff * rho_powers[power]
        r_nm = np.where(rho > 1, 0, r_nm)
        r_nm = np.where(rho == 1, 1, r_nm)
        return r_nm

    def _calculate_modelab_angular_convention(self, r_nm, phi, n, m, odd=True):
        # Conversion factors from standard basis (wiki, Born & Wolf, Wolfram) to Modelab basis.
        if odd:
            if n % 2 == 0:
                z_nm = (-1)**(m // 2 + 1) * r_nm * np.sin(m * phi)
            else:
                z_nm = (-1)**((m + 1) // 2) * r_nm * np.cos(m * phi)
        else:
            if (n % 2 == 0) and (m == 0):
                z_nm = r_nm * np.cos(m * phi) + (-1)**(n // 2 + 1)
            elif (n % 2 == 0) and (m != 0):
                z_nm = (-1)**(m // 2) * r_nm * np.cos(m * phi)
            else:
                z_nm = (-1)**((m + 1) // 2) * r_nm * np.sin(m * phi)
        return z_nm

    def generate_zernike_basis(self, rho, phi, max_radial_idx, mask_size_x_in_pixels, mask_size_y_in_pixels):
        self.zernike_count = (max_radial_idx + 1) ** 2
        zernike_basis_array = np.zeros((self.zernike_count, self.N, self.N), dtype=np.float32)

        rho_powers = self._build_rho_powers(rho, max_radial_idx)
        radial_coeffs_cache = {}
        radial_cache = {}

        for n in range(max_radial_idx + 1):
            for m in range(-n, n + 1, 2):
                m_abs = abs(m)
                radial_key = (n, m_abs)
                if radial_key not in radial_cache:
                    radial_cache[radial_key] = self._compute_radial_nm(
                        rho, rho_powers, n, m_abs, radial_coeffs_cache
                    )
                z_nm = self._calculate_modelab_angular_convention(
                    radial_cache[radial_key], phi, n, m_abs, odd=(m < 0)
                )
                idx = (n * (n + 2) + m) // 2
                zernike_basis_array[idx] = z_nm

        return zernike_basis_array
    
    def _apply_coefficients(self, coefficients):
        self.h_mask_center_x = coefficients["h_mask_center_x"]
        self.h_mask_center_y = coefficients["h_mask_center_y"]
        self.v_mask_center_x = coefficients["v_mask_center_x"]
        self.v_mask_center_y = coefficients["v_mask_center_y"]

        self.piston_h = coefficients["piston_h"]
        self.tilt_x_h = coefficients["tilt_x_h"]
        self.tilt_y_h = coefficients["tilt_y_h"]
        self.piston_v = coefficients["piston_v"]
        self.tilt_x_v = coefficients["tilt_x_v"]
        self.tilt_y_v = coefficients["tilt_y_v"]

        self.w_h = coefficients["w_h"]
        self.w_v = coefficients["w_v"]
        self.n = coefficients["n"]
        self.m = coefficients["m"]
        self.zern_idxs = coefficients["zern_idxs"]

    def clean(self):
        self._apply_coefficients(
            dict_to_zernike_coefficients(
                default_coefficients_dict(self.max_zernike_radial_number),
                max_zernike_radial_number=self.max_zernike_radial_number,
            )
        )

    def load_coefs(self, path):
        data = load_coefficients_json(path)
        coefficients = dict_to_zernike_coefficients(
            data,
            max_zernike_radial_number=self.max_zernike_radial_number,
        )
        self._apply_coefficients(coefficients)

    def save_coefs(self, path):
        data = zernike_coefficients_to_dict(self)
        save_coefficients_json(path, data)

    def write_backup_of_zernike_coefs(self):
        # Write a copy of piston coefs
        self.b_piston_h = self.piston_h
        self.b_piston_v = self.piston_v

        # Write a copy of tilt coefs
        self.b_tilt_x_h = self.tilt_x_h
        self.b_tilt_y_h = self.tilt_y_h

        self.b_tilt_x_v = self.tilt_x_v
        self.b_tilt_y_v = self.tilt_y_v

        # Write a copy of the rest of Zernikes. WATCH it!!! This time we need a real copy of array to avoid modifing the original using the view
        self.b_w_h = np.copy(self.w_h)
        self.b_w_v = np.copy(self.w_v)

    def read_backup_of_zernike_coefs(self):
        # Read a copy of piston coefs
        self.piston_h = self.b_piston_h
        self.piston_v = self.b_piston_v

        # Read a copy of tilt coefs
        self.tilt_x_h = self.b_tilt_x_h
        self.tilt_y_h = self.b_tilt_y_h

        self.tilt_x_v = self.b_tilt_x_v
        self.tilt_y_v = self.b_tilt_y_v

        # Read a copy of the rest of Zernikes
        self.w_h = np.copy(self.b_w_h)
        self.w_v = np.copy(self.b_w_v)

    def modify_zernike_coefs(self, apply_tilt_offset = False, apply_defocus_offset = False, scan_tilt = False, scan_defocus = False, **kwargs):

        # NOTE: I'm doing the calc with numexpr module as it gives 10x speed up. Now it is faster then the SLM so no need to push it further for now
        self.write_backup_of_zernike_coefs() # Save the original Zernike coefs (originals loaded after the Zernike field is calculated)

        # ================ Making sure all the required offsets for tilt and defocus are applied to Zernikes =========================================
        if apply_tilt_offset or apply_defocus_offset:
            if 'MASK_corrections' in kwargs:
                MASK_corrections = kwargs.get('MASK_corrections')

                if apply_tilt_offset:
                    self.tilt_x_h += MASK_corrections['tilt_x_H']
                    self.tilt_y_h += MASK_corrections['tilt_y_H']
                    
                    self.tilt_x_v += MASK_corrections['tilt_x_V']
                    self.tilt_y_v += MASK_corrections['tilt_y_V'] 

                if apply_defocus_offset:       
                    self.w_h[1] += MASK_corrections['defocus_H']
                    self.w_v[1] += MASK_corrections['defocus_V']      
            else:
                raise Exception('"apply_tilt_offset = True" or "apply_defocus_offset = True" requires "MASK_corrections" dictionary.')
        # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        # ================ Let's scan using tilt_x and tilt_y =========================================
        if scan_tilt:
            if all (key in kwargs for key in ('tilt_x','tilt_y')):
                tilt_x = kwargs.get('tilt_x')
                tilt_y = kwargs.get('tilt_y')

                self.tilt_x_h += tilt_x
                self.tilt_y_h += tilt_y
                    
                self.tilt_x_v += tilt_x
                self.tilt_y_v += tilt_y                
            else:
                raise Exception('scan_tilt = "True" requires "tilt_x", "tilt_y" keyword arguments')                
        # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        # ================ Let's scan defocus =========================================      
        if scan_defocus:
            if 'defocus' in kwargs:   
                defocus = kwargs.get('defocus')

                self.w_h[1] += defocus
                self.w_v[1] += defocus             
            else:
                raise Exception('scan_defocus = "True" requires "defocus" keyword arguments') 
        # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
    def make_zernike_fields(self):       
        # set the total Zernike phase array in H,V pol to 0
        zern_phase_h = np.zeros((self.N,self.N), dtype = np.float32)
        zern_phase_v = np.zeros((self.N,self.N), dtype = np.float32)

        
        # calculate the phase coming from Zernikes and modelab coefficients (piston, tilt_x and tilt_y are not part of this, that's why there is  (-3) in the loop)
        for i in range(self.zernike_count - 3):
            zern_phase_h += ne.evaluate('w_h_i*zern_i', {'w_h_i': self.w_h[i], 'zern_i': self.zernike_basis_array[self.zern_idxs[i],:,:]})
            zern_phase_v += ne.evaluate('w_v_i*zern_i', {'w_v_i': self.w_v[i], 'zern_i': self.zernike_basis_array[self.zern_idxs[i],:,:]})

        # add tilts (not piston as I never used it so far)    
        zern_phase_h += ne.evaluate('tilt_x_h*zern_1 + tilt_y_h*zern_2', {'tilt_x_h': self.tilt_x_h, 'tilt_y_h': self.tilt_y_h,\
                                                                        'zern_1': self.zernike_basis_array[1,:,:], 'zern_2': self.zernike_basis_array[2,:,:]})
        zern_phase_v += ne.evaluate('tilt_x_v*zern_1 + tilt_y_v*zern_2', {'tilt_x_v': self.tilt_x_v, 'tilt_y_v': self.tilt_y_v,\
                                                                        'zern_1': self.zernike_basis_array[1,:,:], 'zern_2': self.zernike_basis_array[2,:,:]})
        
        # Make the fields from zern_phase
        self.field_h = ne.evaluate('exp(1j*pi*zern_phase_h)', {'zern_phase_h': zern_phase_h, 'pi': np.pi}).astype(np.complex64)
        self.field_v = ne.evaluate('exp(1j*pi*zern_phase_v)', {'zern_phase_v': zern_phase_v, 'pi': np.pi}).astype(np.complex64)

        self.read_backup_of_zernike_coefs()  # Revert Zernnike coefs back to their original state