;; $Id$
;;
;*************************************************************************
;+
; NAME:
;       mas_dot__define
;
;
; PURPOSE:
;       
;       Reconstructs HARDI MRI voxel data using the Diffusion
;       Orientation Transformation method outlined in reference [1].
;
; REFERENCES: 
;
;       1. Ozarslan, Evren et al. “Resolution of complex tissue
;          microarchitecture using the diffusion orientation transform
;          (DOT).” NeuroImage 31.3 (2006): 1086-1103.
;
;
; CATEGORY:
;
;       DIFFUSION TOOLS
;
; INPUTS:
;
;
;
; KEYWORDS PARAMETERS:
;
;        L_MAX:        The L-order of the spherical harmonic basis. Note
;                      that this choice is highly dependent on the
;                      quality of the data (SNR, B-value). The default
;                      value is 4.
;
;        B_VALUES:     An array containing the B-values for the
;                      data. The number of elements in this array should
;                      be equal to the number of diffusion weighted
;                      scan. This array is used to distinguish the
;                      "high" b-value scans from the "low" by using the
;                      mean value as the pivotal value. They are also
;                      used to compute the ADC. Expected units are
;                      s/mm^2
;
;        GRADIENTS:    A 2xn (spherical - [theta,phi]) or 3xn (cartesian
;                      - [x,y,z]) array of unit gradient directions. See
;                      DEGREES and SPHERICAL.
;
;        DEGREES:      Set this keyword if the GRADIENTS are in
;                      spherical coordinates and the angles are in
;                      degrees.
;
;        SPHERICAL:    Set this keyword if the GRADIENTS are in
;                      spherical coordinates. Radians are assumed
;                      otherwise.
;
;        DEBUG:        Set this keyword to have extra debugging messages
;                      sent to the IDL console
;
;        RADIUS:       The desired diffusion radius
;
;        DIFF_TIME:    The diffusion time from the pulse sequence
;                      parameters. 
;
; OUTPUTS:
;
;
;
; METHODS:
;
;     setProperty: Sets properties
;
;        *** NOT IMPLEMENTED ***
;
;     getProperty: Gets properties
;
;        dot->getProperty, l_max=lmax
;
;     computeWeights: compute the intermediate weights from the
;                     data. These weights can be saved or used to
;                     compute probabilities.
;
;     getDisplacementProbability: Gets the probability
;
;        pr = dot->getDisplacementProbability(data=data, direction=d, $
;                                            /spherical, radius=r)
;       Or, more efficiently:
;
;        w = dot->computeWeights(data=data, radius=r)
;        pr = dot->getDisplacementProbability(weights=w, direction=d1, $
;                                            /spherical, radius=r)
;        pr = dot->getDisplacementProbability(weights=w, direction=d2, $
;                                            /spherical, radius=r)
;        pr = dot->getDisplacementProbability(weights=w, direction=d3, $
;                                            /spherical, radius=r)

;
; PROCEDURE:
;
;     The user is responsible for preparing the gradients and b-value
;     arrays as indicated in the KEYWORDS section. Once prepared, a
;     typical usage is, given data = an array of diffusion weighted data:
;
;     dot = obj_new('mas_dot', b_values=bvals, gradients=grad, $
;                   l_max=4, diff_time=td, radius=12.0)
;     Pr = dot->getDisplacementProbability(data, direction=dir)
;     obj->destroy, dot
;
;     See also the example procedure at the end of this file.
;
; EXAMPLE:
;
;     See the example procedure at the end of this file.
;
; MODIFICATION HISTORY:
;
;     Written by: Bill Triplett, Feb 2009.
;-

function mas_dot::init, $
   debug=debug, $
   diff_time=diff_time, $
   b_values=b_values, $
   bval_thr=bval_thr, $
   gradients=gradients, $
   radius=radius, $
   l_max=l_max, $
   spherical=spherical, $
   degrees=degrees

    if (keyword_set(debug)) then begin
       self.debug =  byte(debug)
    endif
    
    if (keyword_set(diff_time)) then begin
       self.diff_time = diff_time
    endif

    if (keyword_set(bval_thr)) then begin
        self.b_thr = float(bval_thr)
    endif else begin
        self.b_thr = 400.
    endelse
    
    if (keyword_set(radius)) then begin
       self.diff_radius = radius
    endif else begin
       self.diff_radius = 10.0
    endelse

    if (keyword_set(l_max)) then begin
       self.l_max = l_max
    endif else begin
       self.l_max = 4
    endelse

    self.can_use_weights = 1

    if (n_elements(b_values) ne 0) then begin
       self->prepareBValues, b_values
    endif else begin

    endelse

    if (n_elements(gradients) ne 0) then begin
       
       if (keyword_set(spherical)) then begin
          n_gradients = n_elements(gradients)/2
       endif else begin
          n_gradients = n_elements(gradients)/3
       endelse
       
       if (self.debug) then print, "mas_dot::init: n_gradients = ", n_gradients

       self->prepareGradients, gradients, $
                               spherical=spherical, $
                               degrees=degrees
    endif

    return, 1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;
;     mas_dot::cleanup
;
; PURPOSE:
;
;     Destructor
;
;-

pro mas_dot::cleanup

    if (self.debug) then print, "mas_dot::cleanup: cleaning up"
    ptr_free, self.gradients_sph
    ptr_free, [self.b_values, self.b_high, self.b_low]
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;
;     mas_dot::prepareGradients
;
; PURPOSE:
;
;     Private method that prepares the gradients for computation
;
; ARGUMENTS:
;
;     GRADIENTS - the array of gradients
;
; KEYWORDS:
;     
;     SPHERICAL:  set if the gradients are in spherical coordinates
;     DEGREES:    set if the spherical angles are in degrees
; 
; MODIFICATION HISTORY:
;
;-

pro mas_dot::prepareGradients, gradients, $
                               spherical=spherical, $
                               degrees=degrees

    our_gradients = gradients

    if (keyword_set(spherical)) then begin
       if (self.debug) then print, "mas_dot::prepareGradients: Spherical"
       sz_grad = size(our_gradients, /dimension)
       if (sz_grad[1] eq 2 or sz_grad[1] eq 3) then begin
          if (self.debug) then print, "mas_dot::prepareGradients: Transposing"
          our_gradients = transpose(our_gradients)
       endif

       if (keyword_set(degrees)) then begin
          if (self.debug) then print, "mas_dot::prepareGradients: Degrees"
          our_gradients *= !DTOR
       endif

       ptr_free, self.gradients_sph
       self.gradients_sph = ptr_new(our_gradients, /no_copy)

    endif else begin

       if (self.debug) then print, "mas_dot::prepareGradients: Cartesian"
       our_gradients = (cv_coord(from_rect=gradients, /to_sphere))[0:1,*]
       our_gradients[1,*] = !DPI/2. - our_gradients[1,*]
       our_gradients = our_gradients[[1,0],*]
       ptr_free, self.gradients_sph
       self.gradients_sph = ptr_new(our_gradients, /no_copy)

    endelse

    n_of_lms = ((self.l_max)/2+1)^2
    n_of_is  = n_elements(*self.b_high)

    l_arr=[0,2,2,2,4,4,4,4,4,replicate(6,7),replicate(8,9),replicate(10,11)]
    m_arr=[0,0,1,2,0,1,2,3,4,indgen(7),indgen(9),indgen(11)]

    sph_harm_wts = dcomplexarr(n_of_is, n_of_lms)

    for lm = 0, n_of_lms-1 do begin
       sph_harm_wts[*,lm] = conj(spher_harm((*self.gradients_sph)[0,*self.b_high], $
                                            (*self.gradients_sph)[1,*self.b_high], $
                                            l_arr[lm], m_arr[lm]))
    endfor
    
    ptr_free, self.spher_harm_wts
    self.spher_harm_wts = ptr_new(sph_harm_wts)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;
;     mas_dot::prepareBValues
;
; PURPOSE:
;
;     Private method that prepares the b-value array for computation.
;
; ARGUMENTS:
;
;     BVALUES - the array of b-values
;
; MODIFICATION HISTORY:
;
;-
pro mas_dot::prepareBValues, bvalues

    
    b_thr = self.b_thr * 1e-6 ;; mean(bvalues)
    
    high_b_idx = where(bvalues gt b_thr, complement=low_b_idx, count)

    if (count eq 0) then begin
       message, 'There must be a distinct low/high b-value set'
       return
    endif else begin
       if (self.debug) then print, "mas_dot::prepareBValues: setting b-value information"
       ptr_free, [ self.b_high, self.b_low, self.b_values ]
       self.b_high = ptr_new(high_b_idx)
       self.b_low  = ptr_new(low_b_idx)
       self.b_values = ptr_new(bvalues); * 1e-6)
    endelse
  
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;
;     mas_dot::computeWeights
;
; PURPOSE:
;
;     Computes the p_lm weights from the directional ADC. Evaluating
;     I_lu (eqn. 28) from DOT paper. Note that the A_i and B_i
;     computations are being computing using Horners method of
;     evaluating polynomials.
;
; ARGUMENTS:
;
;     DATA:    The voxel data
;
;     RADIUS:  The desired diffusion radius. If not specified, default
;              value is used either from initialization or object
;              default.
;
; MODIFICATION HISTORY:
;
;-

function mas_dot::computeWeights, data, radius=radius

    if (n_elements(*self.b_low) ne 1) then begin
       low_b_data = (moment(data[*self.b_low]))[0]
    endif else begin
       low_b_data = (moment([data[*self.b_low],data[*self.b_low]]))[0]
    endelse

    if (not keyword_set(radius)) then begin
       radius = self.diff_radius
    endif else begin
       self.diff_radius = radius
    endelse

    high_b_data = reform(data[*self.b_high])/low_b_data

    adc = -alog(high_b_data)/(*self.b_values)[*self.b_high]
    
    ;; negative ADC will produce NANs later in computation. We decided
    ;; that negative ADC, which occurred when the S0 data was greater
    ;; than the diffusion-weighted data, was due to noise. Our remedy
    ;; was to set the ADC in those cases to a number close to zero.
    neg_adc = where(adc le 0, n_neg_adc)
    if (n_neg_adc gt 0) then begin
       adc[neg_adc] = 1e-4
    endif

    bi   = double( (SQRT(self.diff_time)/radius) * SQRT(adc) ) ;; beta inverse
    bi2  = bi*bi

    expo = EXP(-1.D /(4.D * bi2))/(4.D * !dpi * adc * self.diff_time)^(3.D / 2.D)
    errf = ERF( 1.D / bi / 2.D )/(4.D * !dpi * radius^3)

    weights = dcomplexarr((self.l_max/2+1)^2)
    for ll = 0, self.l_max, 2 do begin
       
       case ll of 

          0: begin
             A_l = 1D
             B_l = 0D
          end

          2: begin
             A_l = -(1D + 6D*bi2)
             B_l = 3D
          end

          4: begin
             A_l = 1D + bi2*(20D + 210D*bi2)
             B_l = 15D/2. * (1D - 14D*bi2)
          end

          6: begin
             A_l = -(1D + bi2*(42D + bi2*(1575D/2. + bi2*10395D)))
             B_l = 105D/8. * (1D + bi2*(-36D + bi2*396D))
          end

          8: begin
             A_l = 1 + bi2*(72 + bi2*(10395D/4. + bi2*(45045D + bi2*675675D)))
             B_l = 316D/16D * (1D + bi2*(-66D + bi2*(1716D - bi2*17160D)))
          end

       endcase

       Il_u = A_l * expo + B_l * errf

       for l = 0,ll do begin

          tmp = (-1.)^(ll/2) * 2.0 * total( Il_u * (*self.spher_harm_wts)[*,(ll/2)^2 + l] )
          weights[(ll/2)^2 + l] = tmp

       endfor

    endfor

    return, weights

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;
;     mas_dot::getDisplacementProbability
;
; PURPOSE:
;
;     Compute the displacement probability in a particular
;     direction. Note that for this DOT reconstrutor, it is possible
;     to pass an array of directions (3,N)/(2,N) and get a fast
;     computation for every direction. This behavior is equivalent to
;     using the MAS_DOT_BULK_RECONSTRUCTOR object. This behavior is
;     a side effect of the computation process.
;
; KEYWORDS:
;
;     WEIGHTS:   The p_lm weights computed from Eqn 28. [1] For
;                multiple computations on a single voxel, it is more
;                efficient to compute this first and use it to compute
;                probabilities. 
;
;     DATA:      The array of data for a particular voxel
;
;     RADIUS:    The desired diffusion radius
;
;     DIRECTION: the [x,y,z] direction or [theta,phi] spherical
;                direction.
; 
;     SPHERICAL: Set this keyword if DIRECTION is in spherical
;                coordinates
;
;     DEGREES:   Set this keyword if the spherical angles are in
;                degrees. The default is radians
;
; MODIFICATION HISTORY:
;;
;-

function mas_dot::getDisplacementProbability, $
   weights=weights, $
   data=data, $
   radius=radius, $
   direction=direction, $
   spherical=spherical, $
   degrees=degrees

    if (not keyword_set(spherical)) then begin
       our_dir = (cv_coord(from_rect=direction, /to_sphere))[0:1,*]
       our_dir[1,*] = !DPI/2.0 - our_dir[1,*]
    endif else if keyword_set(degrees) then begin
       our_dir = reform(direction[0:1,*] * !DTOR)
    endif else begin
       our_dir = direction
    endelse
       
    n_dirs = n_elements(our_dir)/2

    if (not keyword_set(weights)) then begin
       if (n_elements(data) ne 0) then begin
          weights = self->computeWeights(data, radius=radius)
       endif else begin
          message, 'mas_dot::getDisplacementProbability: weights or data needed'
          return, -1
       endelse
    endif 
          
    N_lm = ((self.l_max)/2+1)^2
    n=0
    spharm = dcomplexarr(N_lm, n_dirs)
    for l = 0, self.l_max, 2 do begin
       for m = 0, l do begin
          spharm[n++,*] = ((m eq 0) ? 1D : 2D) * spher_harm(our_dir[1,*], $
                                                            our_dir[0,*], $
                                                            l, m)
       endfor
    endfor

    sum = 0D

    for ll = 0, n-1 do begin
       sum += real_part(weights[ll] * spharm[ll,*])
    endfor

    return, reform(sum)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;       mas_dot::getProperty
;
;
; PURPOSE:
;       
;       Receive information about the parameters of the object
;
; KEYWORDS PARAMETERS:
;
;       Each keyword can be set to receive its current value

pro mas_dot::getProperty, $
   debug=debug, $
   diff_radius=diff_radius, $
   can_use_weights=can_use_weights, $
   diff_time=diff_time, $
   l_max=l_max

    l_max = self.l_max
    debug = self.debug
    diff_radius = self.diff_radius
    diff_time = self.diff_time
    can_use_weights = self.can_use_weights

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;       mas_dot__define
;
; PURPOSE: 
;
;       The MAS_DOT object's structure definition.
;-

pro mas_dot__define

    struct = { MAS_DOT, $
               can_use_weights: 0B, $
               debug: 0B, $
               diff_radius: 10D, $
               diff_time: 0D, $
               M_matrix: ptr_new(), $
               l_max: 4L, $
               b_thr: float(0), $
               b_low: ptr_new(), $               
               b_high: ptr_new(), $
               b_values:ptr_new(), $
               spher_harm_wts: ptr_new(), $
               gradients_sph: ptr_new() }
end

pro mas_dot_bulk_reconstructor__define

   struct = { MAS_DOT_BULK_RECONSTRUCTOR, inherits MAS_DOT, $
              reco_vertices: ptr_new() }

end

function mas_dot_bulk_reconstructor::makeBootstrap, $
   data=data
   
   ;; Not supported.
   return, 0
end
    
pro mas_dot_bulk_reconstructor::cleanup

    if (self.debug) then print, "mas_dot_bulk_reconstructor::cleanup"
    ptr_free, self.reco_vertices
    self->mas_dot::cleanup

end

pro mas_dot_bulk_reconstructor::setReconstructionDirections, $
   directions, $
   spherical=spherical, $
   degrees=degrees

    if (not keyword_set(spherical)) then begin
       our_dir = (cv_coord(from_rect=directions, /to_sphere))[0:1,*]
       our_dir[1,*] = !DPI/2.0 - our_dir[1,*]
    endif else if keyword_set(degrees) then begin
       our_dir = reform(directions[0:1,*] * !DTOR)
    endif else begin
       our_dir = directions
    endelse

    n_dirs = n_elements(our_dir)/2

    ptr_free, self.reco_vertices

    N_lm = ((self.l_max)/2+1)^2
    n=0
    spharm = dcomplexarr(N_lm, n_dirs)
    for l = 0, self.l_max, 2 do begin
       for m = 0, l do begin
          spharm[n++,*] = ((m eq 0) ? 1D : 2D) * spher_harm(our_dir[1,*], $
                                                            our_dir[0,*], $
                                                            l, m)
       endfor
    endfor

    ptr_free, self.reco_vertices
    self.reco_vertices = ptr_new(spharm, /no_copy)

end

function mas_dot_bulk_reconstructor::getDisplacementProbability, $
   data=data, $
   weights=weights, $
   radius=radius

    if (not keyword_set(weights)) then begin
       weights = self->mas_dot::computeWeights(data, radius=radius)
    endif
    
    sum = 0D

    for ll = 0, n_elements(weights)-1 do begin
       sum += real_part(weights[ll] * (*self.reco_vertices)[ll,*])
    endfor

    return, reform(sum)

end

pro mas_dot_example

    ;; Theta angles (polar, from z-axis = 0)
    theta = [ 0.00000, 63.4349, 63.4513, 63.4581, 63.4581, 63.4513, $
              21.1526, 21.1266, 21.1443, 21.1443, 21.1266, 37.3656, $
              37.3559, 37.3625, 37.3559, 37.3656, 42.2852, 42.2832, $
              42.2797, 42.2797, 42.2832, 58.8617, 58.8800, 58.8731, $
              58.8843, 58.8695, 58.8695, 58.8843, 58.8731, 58.8800, $
              58.8617, 79.1655, 79.1612, 79.1621, 79.1612, 79.1655, $
              81.0275, 81.0263, 81.0253, 81.0288, 81.0260, 81.0260, $
              81.0288, 81.0253, 81.0263, 81.0275, 0.00000 ]

    ;; Phi angles
    phi = [ 180.000, 0.00000, 72.0309, 144.001, 215.999, 287.969, $
            0.00000, 72.0676, 144.019, 215.981, 287.932, 36.0204, $
            108.047, 180.000, 251.953, 323.980, 0.00000, 71.9958, $
            143.948, 216.052, 288.004, 23.6294, 48.3609, 95.6307, $
            120.367, 167.587, 192.413, 239.633, 264.369, 311.639, $
            336.371, 35.9716, 107.974, 180.000, 252.026, 324.028, $
            12.3904, 59.5933, 84.3644, 131.595, 156.367, 203.633, $
            228.405, 275.636, 300.407, 347.610, 0.00000 ]

    grad = [ [theta], [phi] ]

    ;; b_values
    bvals = [ replicate(800.0, 46), 0.] * 1e-6

    ;; diffusion time
    dt = 17.2 * 1e-3

    ;; some diffusion data
    data = [ 94971.8, 96982.2, 68006.5, 82847.6, 60716.6, 74985.8, $
             83197.1, 104232., 96661.7, 84150.5, 84722.3, 94969.2, $
             99566.2, 83487.4, 72094.5, 87053.2, 92669.6, 95042.6, $
             103382., 84024.1, 84522.4, 86291.2, 66306.5, 81414.1, $
             101041., 89636.8, 76524.8, 66982.3, 66651.2, 86889.2, $
             92382.5, 70348.6, 80515.4, 82981.2, 73442.0, 99426.0, $
             79545.3, 71345.7, 74289.4, 85795.4, 101870., 67104.5, $
             59460.6, 76771.4, 77757.9, 85816.7, 150168. ]
    ;; l_max
    l_max = 6

    ;; dot object
    dot = obj_new('mas_dot', l_max=l_max, b_values=bvals, $
                  radius=12.0, diff_time=dt, $
                  gradients=grad, /spherical, /degrees)

    ;; get some probabilities
    print, "--> Using Data"
    print, "PR[1,0,0]: ", dot->getDisplacementProbability(data=data, $
                                                          direction=[1.,0,0])
    print, "PR[0,1,0]: ", dot->getDisplacementProbability(data=data, $
                                                          direction=[0,1.,0])
    print, "PR[0,0,1]: ", dot->getDisplacementProbability(data=data, $
                                                          direction=[0,0,1.])
    ;; exploit the ease and computational effectiveness of using weights
    print, "--> Using Weights"
    w = dot->computeWeights(data)
    print, "PR[1,0,0]: ", dot->getDisplacementProbability(weights=w, $
                                                          direction=[1.,0,0])
    print, "PR[0,1,0]: ", dot->getDisplacementProbability(weights=w, $
                                                          direction=[0,1.,0])
    print, "PR[0,0,1]: ", dot->getDisplacementProbability(weights=w, $
                                                          direction=[0,0,1.])
    
    obj_destroy, dot
    
    ;; special bulk reconstructor object
    dot = obj_new('mas_dot_bulk_reconstructor', b_values=bvals, $
                  gradients=grad, /spherical, /degrees, $
                  diff_time=dt, radius=12.0, l_max=l_max)

    ;; set the directions
    mesh_obj, 4, vert, poly, fltarr(50,50)+1
    dot->setReconstructionDirections, vert

    ;; get the bulk probabilites
    Pr  = dot->getDisplacementProbability(data=data)

    ;; minmax to enhance shape
    Pr  = (Pr - min(Pr))/(max(Pr)-min(Pr))

    ;; distort sphere to create direction shape
    vert[0,*] *= Pr
    vert[1,*] *= Pr
    vert[2,*] *= Pr

    poly = obj_new('idlgrpolygon', vert, polygon=poly, color=[255,0,0])

    xobjview, poly, background=[0,0,0], /block

    obj_destroy, poly
    obj_destroy, dot

end
