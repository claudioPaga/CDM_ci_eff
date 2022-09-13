pro setup_cdm_variables

;;;
;;; CP, 8 Sept 2020
;;; Summary - Writes the file to store trap parameters and properties
;;;           values to be used in the CDM processing
;;;           Trap densities are for EOL estimated CTI
;;;
;;; Output - distort_cdm_properties.sav file
;;;
;;; History -
;;; 24 May 2021 - Trap species upped to 4. Trap parameters updated to
;;;               values provided by Tom Buggey (21/05/2021)
;;;               To better match results and future applications the
;;;               trap parameters are now expressed using physical
;;;               units instead of the convenient units such as
;;;               readout time.
;;;               Readout times units are now seconds instead of units
;;;               of readout periods.
;;;               Added additional parameters used for internal
;;;               transformation of units
;;;               According to Tom the same trap species should be
;;;               present in image/store/serial readout.
;;;
;;; 15 June 2021 - Trap densities from Tom are relative to entire pixel.               
;;;                I want to scale the density to the value measured
;;;                from trap-pumping test at the OU, that reveals a
;;;                turn-up at around ~20K e-, supposedly the edge of
;;;                the SBC. At that flux the density is of approx
;;;                0.005 traps/pixel (see plot density.png)
;;;                Normalise Tom's densities to that value, as
;;;                my code implicitly only accounts for traps in the
;;;                SBC
;;;
;;; 17 Nov 2021 - Set temperature to T=-100C (173K) to test possible
;;;               warmer regimes.  
;;;
;;; 11 March 2022 - Re-evaluated the trap_density_sbc factor.
;;;                 During meeting with OU David noted that the trap
;;;                 density values were relative to the signal flux
;;;                 used in the test, not to the entire pixel. The
;;;                 signal during the test is 10k ele (Tom email, 10
;;;                 March 2022). Therefore, the scaling to the SBC is
;;;                 actually (20k/10k)^(1-charge_volume_coeff)~1.125
;;;                 and not 0.25 (the fraction of traps in a pixel in
;;;                 the SBC).
;;;
;;; 13 September 2022 - Add flags and effective densities to model the effect of CI in updated CDM treatment
  
;;;
  
image_section_lines = 3791
store_section_lines = 719
serial_columns = 4510
overscan_cols = 10
readout_nodes = 2
;ccd_mode_binning can be anything but should be any of 1/6/24
ccd_mode_binning = 6
;Charge injection info
charge_injection_flag = 1
;Convergence threshold (ev) of iterative CTI correction algorithm
iter_th = 1
;Makimum number of attempts for iterative CTI correction algorithm
iter_max = 10
ccd_t_kelvin = 173.0
fwc = 20000. ;;; full well capacity of the supplementary buried channel
vth = 674345.06 * sqrt(ccd_t_kelvin) ;;; Thermal velocity of thermal electrons at the operating CCD temperature.
pixel_volume = 1.0E-10

charge_injection_period_frame = 1000 ; Lines gap between CI lines at the top of image after exposure.
charge_injection_time_gap = 2.32  ; This is time in [s] between the parallel CI readout ends (CI line in serial register) and the start of the image readout. This number might be update at the EDU processing/timing is clarified. The EDU time might be longer than the 0.116s George has in his timing plot, as that might only be the time required to transfer pixels in and not actually the processing.
ima_expo_time = 2.5 ; 
charge_injection_block_lines = 1 ; Number of consecutive lines injected with charge
charge_injection_electrons = 10000 ; Number of electrons per pixel of injected charge.
;trap release timescales are in units of the readout periods, for
;simplicity assumed as 1.
readout_image_time = 2.7E-5 ;[s]
readout_store_time = 2.7E-5 ;[s]
readout_serial_time = 1.0/1.4E+5 ; Serial freq 140Khz
trap_species_parallel = 4
trap_species_serial = 4
trap_species_image_density = [0.08, 0.6, 0.0072, 0.6] ; Densities from Tom's email.
;trap_density_sbc = 0.005 ; Density @20k e- (SBC turnup) from
;                                           trap-pumping OU test
;trap_density_sbc = 0.25 ; Fraction of traps in the SBC accoring to OU estimate (Dave report during telecon)
;UPDATE after OU discussion, March 2022.
;Values in Tom's table are measured at a signal of 10K e-. I had assumed those values were relative to the entire pixel, in reality they are relative to a volume
;The volume seen by a signal of 10K e- compared to the SBC of 20K e- is V10/V20 = (10k/20k)^(1-charge_volume_coeff_image)

charge_volume_coeff_image = 0.83 ; Beta derived from PL fit of CCD370_Meanfit.txt (email 6 Oct 2020)
;charge_volume_coeff_image = 0.2
charge_volume_coeff_store = charge_volume_coeff_image
charge_volume_coeff_serial = charge_volume_coeff_image

trap_density_sbc_scaling = (20000./10000.)^(1.-charge_volume_coeff_image) ; With charge_volume_coeff_image = 0.83 this is approx 1.125

trap_species_image_density = trap_species_image_density * trap_density_sbc_scaling
                                ; Scaling to density fraction in SBC.
trap_species_store_density = trap_species_image_density
trap_species_serial_density = trap_species_image_density

; CI effective densities flags. These are set to 1 in case the density of traps must be updated due to the effect of CI
trap_ci_eff_image_flag = [0, 1, 1, 0] ;t0, t3 are unaffected by CI because they are too fast
trap_ci_eff_store_flag = [1, 1, 1, 1]
trap_ci_eff_serial_flag = [1, 1, 1, 1]

; CI effective densities fixed values flags. These are set to 1 if the effective density can be expressed as a single fixed value as a function of DETY
trap_ci_eff_image_fixed = [0, 1, 1, 0] ;
trap_ci_eff_store_fixed = [1, 1, 1, 1]
trap_ci_eff_serial_fixed = [1, 1, 1, 1]

; CI effective densities fractions. Used to specify the fractional effective density
trap_ci_eff_image_fraction = [0.1, 0.2, 0.3, 0.4] 
trap_ci_eff_store_fraction = [0.4, 0.1, 0.1, 0.1]
trap_ci_eff_serial_fraction = [0.4, 0.3, 0.1, 0.1]

; CI effective densities function. These are set to 1 if the effective density is expressed as a function of DETY using exponentials.
; The exponential function parameters are stored in a look-up table.
trap_ci_eff_image_function = [0, 1, 1, 0] ;t0, t3 are unaffected by CI because they are too fast
trap_ci_eff_store_function = [1, 1, 1, 1]
trap_ci_eff_serial_function = [1, 1, 1, 1]

; CI effective densities DETX flag. Indicate if the effective filling caused by CI need a special treatment as a function of DETY in the serial register.
trap_ci_eff_serial_detx = [0, 1, 0, 0]


trap_energy_level_image = [0.235, 0.456, 0.345, 0.165] ;[eV]
trap_energy_level_store = trap_energy_level_image
trap_energy_level_serial = trap_energy_level_image

;;; Transform energy levels into release timescales using the formula:
;;; tau = exp(Et/kT)/sigma*vth*Nc
;;; kT in eV, sigma*vth*Nc ~2.8E+11 
release_image_time = exp(trap_energy_level_image/(ccd_t_kelvin/11604.45))/2.8E+11 ; [s]
release_store_time = release_image_time
release_serial_time = release_image_time

charge_volume_coeff_image = 0.83 ; Beta derived from PL fit of CCD370_Meanfit.txt (email 6 Oct 2020)
;charge_volume_coeff_image = 0.2
charge_volume_coeff_store = charge_volume_coeff_image
charge_volume_coeff_serial = charge_volume_coeff_image

capture_cross_section_image_cm2 = [1.35E-15, 3.7E-15, 5E-16, 6.1E-15] ;[cm2]

capture_cross_section_image = capture_cross_section_image_cm2 * fwc^(1.0-charge_volume_coeff_image) * vth/ (2 * pixel_volume)
capture_cross_section_store = capture_cross_section_image
capture_cross_section_serial = capture_cross_section_image
;capture_cross_section_image = [0.01, 0.01] ; This values corresponds to a Pc = 0.99 at Fe55, Pc = 0.55 @ 0.6 keV
;capture_cross_section_store = [0.01, 0.01]
;capture_cross_section_serial = [0.01, 0.01]
;capture_cross_section_image = [0.005, 0.005] ; This values corresponds to a Pc = 0.9 at Fe55, Pc = 0.25 @ 0.6 keV
;capture_cross_section_store = [0.005, 0.005]
;capture_cross_section_serial = [0.005, 0.005]

flag_split = 1
flag_unbinned = 1 ; If set, the input file has unbinned PHAS values (likely used George code to spread charge to generate FITS file)
flag_event_rec_ll = 0 ; If set, apply energy LL to pixels in damaged event in output file, using lower limit thresholds for total energy and split pixels enegy
detection_ll_eV = 50 ; Threshold of total PICTI energy of the damaged event
detection_ll_split_eV = 20 ; Split pixel threshold.
flag_retrapping = 1
flag_iteration = 0
threshold = 0.1
short_bernoulli = 0
long_bernoulli = 0
flag_binomial = 0
flag_same_pixel_release = 1 ;When set, code will calculate the charge released in the same pixel of caputre and in the following pixel.

SAVE, /VARIABLES, FILENAME = 'distort_cdm_properties.sav'

end

;list_energy_grid[*,event_index]*1000./3.65
;tab[event_index].RAWY
;index = where(xray_unbinned_18cols gt 0., n)
;dims = size(xray_unbinned_18cols,/dim)
;xi = index mod dims[0]
;yi = index / dims[0]
;xi
;yi
;xray_unbinned_18cols[xi[0], tab[event_index].RAWY*6-6:tab[event_index].RAWY*6+11]*1000./3.65
;xray_unbinned_18cols_pcti_image[xi[0], tab[event_index].RAWY*6-6:tab[event_index].RAWY*6+11]
;xray_binned_18cols_pcti_image[xi[0],tab[event_index].RAWY-1]
;xray_binned_18cols_pcti_image[xi[0],tab[event_index].RAWY]   
;xray_binned_18cols_pcti_image[xi[0],tab[event_index].RAWY+1]   
;xray_binned_18cols_store_input[*,store_section_lines]
;xray_binned_18cols_store_input[*,store_section_lines+1]
;xray_binned_18cols_store_input[*,store_section_lines+2]
