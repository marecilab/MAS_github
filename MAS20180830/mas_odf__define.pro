;; $Id$
;;
;*************************************************************************
;+
; NAME:
;       mas_odf__define
;
;
; PURPOSE:
;       
;       Reconstructs HARDI MRI voxel data using the Q-Ball Orientation
;       Distribution Function method outlined in reference [1].
;
; REFERENCES: 
;
;  1. Christopher P. Hess, Pratik Mukherjee. "Q-ball reconstruction of
;     multimodal fiber orientations using the spherical harmonic basis."
;     Magnetic Resonance in Medicine 56.1 (2006): 104-117.
;
;  2. Need to find this reference... do not use this method
;
;  3. M. Descoteaux, E. Angelino, S. Fitzgibbons, R. Deriche  Regularized,
;     Fast and Robust Analytical Q-Ball Imaging,  Magnetic Resonance in 
;     Medicine,  Volume 58, Issue 3. Pages 497-510.
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
;                      mean value as the pivotal value. B-values are not
;                      used in the calculation.
;
;        BVAL_THR:     the value that separates "low" b-values from "high"
;                      b-values
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
;        USE_TIKHONOV: Set this keyword to use TIKHONOV regularization
;                      as outlined in reference [1].
;
;       USE_LAP_SHARP: Set this keyword to use Laplacian sharpening
;                      as outlined in reference [2].
;
;        USE_LAP_BERT: Set this keyword to use Laplace-Beltrami
;                      sharpening as outlined in reference [3].
;
; OUTPUTS:
;
;
;
; METHODS:
;
;     setProperty: Sets properties
;
;        odf->setProperty, l_max=4
;
;     getProperty: Gets properties
;
;        odf->getProperty, l_max=lmax
;
;     getDisplacementProbability:
;
;        odf->getDisplacementProbability(data=data, direction=d, $
;                                        /spherical)
;
;
; PROCEDURE:
;
;     The user is responsible for preparing the gradients and b-value
;     arrays as indicated in the KEYWORDS section. Once prepared, a
;     typical usage is, given data = an array of diffusion weighted data:
;
;     odf = obj_new('mas_odf', b_values=bvals, gradients=grad, l_max=4)
;     Pr = odf->getDisplacementProbability(data, direction=dir)
;     obj->destroy, odf
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

function mas_odf::init, $
    use_lap_bert=use_lap_bert, $
    use_lap_sharp=use_lap_sharp, $
    use_tikhonov=use_tikhonov, $
    l_max=l_max, $
    debug=debug, $
    b_values=b_values, $
    bval_thr=bval_thr, $
    gradients=gradients, $
    spherical=spherical, $
    degrees=degrees
  
    ;; default = 4
    if (not keyword_set(l_max)) then begin
       self.l_max = 4
    endif else begin
       ;; only even numbers are necessary
       self.l_max = (l_max mod 2 eq 0) ? l_max : l_max+1
    endelse
    
    if (keyword_set(bval_thr)) then begin
        self.b_thr = float(bval_thr)
    endif else begin
        self.b_thr = 400.
    endelse
    
    sharp_methods = 0
    if (keyword_set(use_lap_bert)) then begin
       self.use_lap_bert = 1B
       sharp_methods++
    endif

    if (keyword_set(use_tikhonov)) then begin
       self.use_tikhonov = 1B
       sharp_methods++
    endif

    if (keyword_set(use_lap_sharp)) then begin
       self.use_lap_sharp = 1B
       sharp_methods++
    endif
    
    if (sharp_methods gt 1) then begin
       message, 'only one regularization method may be specified.'
       return, 0
    endif

    if (keyword_set(debug)) then begin
       self.debug =  byte(debug)
    endif
    
    self.can_use_weights = 0

    if (n_elements(gradients) ne 0) then begin
       
       if (keyword_set(spherical)) then begin
          n_gradients = n_elements(gradients)/2
       endif else begin
          n_gradients = n_elements(gradients)/3
       endelse
       
       if (self.debug) then print, "mas_odf::init: n_gradients = ", n_gradients

       self->prepareGradients, gradients, $
                               spherical=spherical, $
                               degrees=degrees
    endif

    if (n_elements(b_values) ne 0) then begin
       self->prepareBValues, b_values
    endif

    self->prepareZQMatrix

  return, 1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;
;     mas_odf::cleanup
;
; PURPOSE:
;
;     Destroys pointer data and cleans up object.
;
; MODIFICATION HISTORY:
;;
;-
pro mas_odf::cleanup

    if (self.debug) then print, 'mas_odf::cleanup'
    ptr_free, [ self.b_high, self.b_low, self.b_values ]
    ptr_free, self.gradients_sph
    ptr_free, [ self.A_matrix, self.Z_Q_dagger_matrix, self.P_lm_matrix ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;
;     mas_odf::computeReconstructionZU
;
; PURPOSE:
;
;     Private method to compute the ZU matrix identified in Ref [1].
;
;-

function mas_odf::computeReconstructionZU, directions

    n_directions = n_elements(directions)/2

    if (self.debug) then print, "mas_odf::computeReconstructionZU: n_directions: ", n_directions 

    N = 1L * (self.l_max+1)*(self.l_max+2)/2

    Z_U = dcomplexarr(N, n_directions)
    ptr = 0

    for el = 0, self.l_max, 2 do begin
        for em = -el, el do begin

           tmp = spher_harm(directions[1,*], directions[0,*], el, em)
           Z_U[ptr,*] = tmp
           
           ptr++

        endfor
    endfor

    return, Z_U

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;
;     mas_odf::getDisplacementProbability
;
; PURPOSE:
;
;     Compute the displacement probability in a particular
;     direction. Note that for this ODF reconstrutor, it is possible
;     to pass an array of directions (3,N)/(2,N) and get a fast
;     computation for every direction. This behavior is equivalent to
;     using the MAS_ODF_BULK_RECONSTRUCTOR object. This behavior is
;     a side effect of the computation process.
;
; KEYWORDS:
;
;     DATA:  The array of data for a particular voxel
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

function mas_odf::getDisplacementProbability, $
   data=data, $
   direction=direction, $
   spherical=spherical, $
   degrees=degrees

    if (not ptr_valid(self.A_matrix)) then begin
       message, 'mas_odf::getDisplacementProbability: '+ $
                'Cannot compute probabilities without first '+$
                'supplying gradient directions and b values.'
    endif

    ;; make it our own 
    our_dir = direction

    if (not keyword_set(spherical)) then begin
       our_dir = cv_coord(from_rect=our_dir, /to_sphere)
       our_dir[1,*] = !DPI/2 - our_dir[1,*]
    endif
    if (keyword_set(degrees)) then begin
       our_dir[1,*] *= !DTOR
       our_dir[0,*] *= !DTOR
    endif

    Z_U = self->computeReconstructionZU(our_dir)

    A_matrix = real_part(transpose(Z_U) # (*self.A_matrix))

    e = data[(*self.b_high)] / (mean(data[(*self.b_low)]))[0]
    
    Pr = A_matrix # e

    return, Pr

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;
;     mas_odf::prepareBValues
;
; PURPOSE:
;
;     Private method that prepares the b-value array for computation.
;
; ARGUMENTS:
;
;     BVALUES - the array of b-values
;
; KEYWORDS:
;
; MODIFICATION HISTORY:
;
;-

pro mas_odf::prepareBValues, bvalues, bval_thr=bval_thr
    
    b_thr = self.b_thr
    
    high_b_idx = where(bvalues gt b_thr, complement=low_b_idx, count)

    if (count eq 0) then begin
       message, 'There must be a distinct low/high b-value set'
       return
    endif else begin
       ptr_free, [ self.b_high, self.b_low, self.b_values ]
       self.b_high = ptr_new(high_b_idx)
       self.b_low  = ptr_new(low_b_idx)
       self.b_values = ptr_new(bvalues)
    endelse
  
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;
;     mas_odf::prepareGradients
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

pro mas_odf::prepareGradients, gradients, $
                               spherical=spherical, $
                               degrees=degrees

    our_gradients = gradients

    if (keyword_set(spherical)) then begin
       if (self.debug) then print, "mas_odf::prepareGradients: Spherical"
       sz_grad = size(our_gradients, /dimension)
       if (sz_grad[0] eq 2 or sz_grad[0] eq 3) then begin
          if (self.debug) then print, "mas_odf::prepareGradients: Transposing"
          our_gradients = transpose(our_gradients)
       endif

       if (keyword_set(degrees)) then begin
          if (self.debug) then print, "mas_odf::prepareGradients: Degrees"
          our_gradients *= !DTOR
       endif

    endif else begin
       
       our_gradients = cv_coord(from_rect=our_gradients, /to_sphere)
       our_gradients[1,*] = !DPI/2. - our_gradients[1,*]

    endelse

    ptr_free, self.gradients_sph
    self.gradients_sph = ptr_new(our_gradients, /no_copy)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;
;     mas_odf::prepareZQMatrix
;
; PURPOSE:
;
;     Private method to compute the ZQ matrix identified in Ref [1].
;
;-

pro mas_odf::prepareZQMatrix

    N = 1L * (self.l_max+1)*(self.l_max+2)/2

    Z_Q = dcomplexarr(N, n_elements(*self.b_high))
    P   = dblarr(N,N)
    ptr = 0

    lambda_mat = dblarr(N,N)
    lambda = 0.1D

    for el = 0, self.l_max, 2 do begin

       for em = -el, el do begin

          P[ptr,ptr] = legendre(0,el)  
        
          if (self.use_lap_sharp) then begin
             P[ptr,ptr] *= (1.0D + 0.1 * el*(el+1) )
          endif else if (self.use_tikhonov) then begin
             lambda_mat[ptr,ptr] = lambda
          endif else if (self.use_lap_bert) then begin
             lambda_mat[ptr,ptr] = (el^2 * (el+1)^2) * 0.06
          endif
          
          tmp = spher_harm((*self.gradients_sph)[*self.b_high,0], $
                           (*self.gradients_sph)[*self.b_high,1], $
                           el, em)
          
          Z_Q[ptr,*] = tmp

          ptr++

       endfor
       
    endfor

    Z_Q1 = matrix_multiply(Z_Q, conj(transpose(Z_Q))) + lambda_mat
    Z_QI = la_invert(Z_Q1, /double) ;; This is picky
    Z_Q_dagger = matrix_multiply(Z_QI, conj(Z_Q))
    
    ptr_free, [ self.Z_Q_matrix, self.A_matrix, self.Z_Q_dagger_matrix, self.P_lm_matrix ]
    self.Z_Q_matrix        = ptr_new(Z_Q)
    self.A_matrix          = ptr_new(P # Z_Q_dagger)
    self.Z_Q_dagger_matrix = ptr_new(Z_Q_dagger, /no_copy)
    self.P_lm_matrix       = ptr_new(P, /no_copy)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;       mas_odf::setProperty
;
;
; PURPOSE:
;       
;       Change the parameters of the object
;
; KEYWORDS PARAMETERS:
;
;       Each keyword can be set to set its current value in the object

pro mas_odf::setProperty, $
   debug=debug, $
   l_max=l_max, $
   b_values=b_values, $
   gradients=gradients, $
   use_lap_bert=use_lap_bert, $
   use_lap_sharp=use_lap_sharp, $
   use_tikhonov=use_tikhonov

   ;; TODO: implement this method

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;       mas_odf::getProperty
;
;
; PURPOSE:
;       
;       Receive information about the parameters of the object
;
; KEYWORDS PARAMETERS:
;
;       Each keyword can be set to receive its current value

pro mas_odf::getProperty, $
   debug=debug, $
   l_max=l_max, $
   use_lap_bert=use_lap_bert, $
   use_lap_sharp=use_lap_sharp, $
   can_use_weights=can_use_weights, $
   use_tikhonov=use_tikhonov

    debug = self.debug
    l_max = self.l_max
    use_lap_bert = self.use_lap_bert
    use_lap_sharp = self.use_lap_sharp
    use_tikhonov = self.use_tikhonov
    can_use_weights = self.can_use_weights

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;       mas_odf__define
;
; PURPOSE: 
;
;       The MAS_ODF object's structure definition.

pro mas_odf__define

    struct = { MAS_ODF, $
               can_use_weights: 0B, $
               debug: 1B, $
               l_max: 0, $
               b_high: ptr_new(), $
               b_low: ptr_new(), $
               b_thr: float(0), $
               b_values:ptr_new(), $
               Z_Q_matrix: ptr_new(), $
               Z_Q_dagger_matrix: ptr_new(), $
               P_lm_matrix: ptr_new(), $
               A_matrix: ptr_new(), $
               gradients_sph:ptr_new(), $
               use_lap_bert: 0B, $
               use_lap_sharp: 0B, $
               use_tikhonov: 0B }
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mas_odf_bulk_reconstructor__define

    struct = { MAS_ODF_BULK_RECONSTRUCTOR, inherits MAS_ODF, $
               bulk_A_matrix: ptr_new(), $
               reco_vertices: ptr_new(), $
               bootstrap_rhat: ptr_new() }

end

pro mas_odf_bulk_reconstructor::cleanup

    if (self.debug) then print, "mas_odf_bulk_reconstructor::cleanup"
    ptr_free, self.bulk_A_matrix
    ptr_free, self.reco_vertices
    ptr_free, self.bootstrap_rhat
    self->mas_odf::cleanup

end

function mas_odf_bulk_reconstructor::getDisplacementProbability, $
   data=data

    ;;return, self->makeBootstrap(data=data)
    
    if (0 and n_elements(*self.b_low) gt 0) then begin
        divisor = (mean(data[(*self.b_low)]))[0]
    endif else begin
        divisor = 1.0
    endelse

    e = data[(*self.b_high)] / divisor
    fit = *self.bulk_A_matrix # e
    return, fit
        
end

function mas_odf_bulk_reconstructor::makeBootstrap, $
   data=data

    if (0 and n_elements(*self.b_low) gt 0) then begin
        divisor = (mean(data[(*self.b_low)]))[0]
    endif else begin
        divisor = 1.0
    endelse

    e = data[(*self.b_high)];; / divisor
    
    zq = *self.Z_Q_matrix
    zq_t = transpose(zq)
    H = real_part(zq ## la_invert(zq_t ## zq) ## zq_t)
    H_diag = dblarr(n_elements(e))
    for i = 0, n_elements(e)-1 do begin
        H_diag[i] = H[i,i]
    endfor
    
    refit = H ## e
    residuals = real_part(reform(e - refit))
    rhat = residuals/sqrt(1.0 - H_diag)
    
    ptr_free, self.bootstrap_rhat
    self.bootstrap_rhat = ptr_new(rhat)
    
    return, 1
    
end

function mas_odf_bulk_reconstructor::sampleBootstrap, seed, data=data

    e = data[(*self.b_high)]
    n_e = n_elements(e)
    samp = round( randomu(seed, n_e) *n_e)
    e_new = e + (*self.bootstrap_rhat)[samp]
    ;print, (*self.bootstrap_rhat)[samp]
    return, *self.bulk_A_matrix # e_new

end

pro mas_odf_bulk_reconstructor::setReconstructionDirections, $
   directions, $
   spherical=spherical, $
   degrees=degrees

    our_dirs = directions

    if (n_elements(directions) eq 0) then begin
       message, 'Directions must have at least one direction'
    endif
    
    if (not keyword_set(spherical)) then begin
       our_dirs = cv_coord(from_rect=our_dirs, /to_sphere)
       our_dirs[1,*] = !DPI/2 - our_dirs[1,*]
    endif else if (keyword_set(degrees)) then begin
       our_dirs[1,*] *= !DTOR
       our_dirs[0,*] *= !DTOR
    endif

    Z_U = self->computeReconstructionZU(reform(our_dirs[0:1,*]))

    ptr_free, self.bulk_A_matrix

    self.bulk_A_matrix = ptr_new( real_part( transpose(Z_U) # *self.A_matrix ) )


end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro mas_odf_example

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
    bvals = [ replicate(800.0, 46), 0.]

    ;; some diffusion data
    data = [ 94971.8, 96982.2, 68006.5, 82847.6, 60716.6, 74985.8, $
             83197.1, 104232., 96661.7, 84150.5, 84722.3, 94969.2, $
             99566.2, 83487.4, 72094.5, 87053.2, 92669.6, 95042.6, $
             103382., 84024.1, 84522.4, 86291.2, 66306.5, 81414.1, $
             101041., 89636.8, 76524.8, 66982.3, 66651.2, 86889.2, $
             92382.5, 70348.6, 80515.4, 82981.2, 73442.0, 99426.0, $
             79545.3, 71345.7, 74289.4, 85795.4, 101870., 67104.5, $
             59460.6, 76771.4, 77757.9, 85816.7, 150168. ]
    
    ;; odf object
    odf = obj_new('mas_odf', b_values=bvals, $
                  gradients=grad, /spherical, /degrees, $
                  l_max=4, /use_tikhonov)

    ;; get some probabilities
    print, "PR[1,0,0]: ", odf->getDisplacementProbability(data=data, $
                                                          direction=[1.,0,0])
    print, "PR[0,1,0]: ", odf->getDisplacementProbability(data=data, $
                                                          direction=[0,1.,0])
    print, "PR[0,0,1]: ", odf->getDisplacementProbability(data=data, $
                                                          direction=[0,0,1.])

    obj_destroy, odf
    
    ;; special bulk reconstructor object
    odf = obj_new('mas_odf_bulk_reconstructor', b_values=bvals, $
                  gradients=grad, /spherical, /degrees, $
                  l_max=4, /use_tikhonov)

    ;; set the directions
    mesh_obj, 4, vert, poly, fltarr(50,50)+1
    odf->setReconstructionDirections, vert

    ;; get the bulk probabilites
    Pr  = odf->getDisplacementProbability(data=data)

    ;; minmax to enhance shape
    Pr  = (Pr - min(Pr))/(max(Pr)-min(Pr))

    ;; distort sphere to create direction shape
    vert[0,*] *= Pr
    vert[1,*] *= Pr
    vert[2,*] *= Pr

    poly = obj_new('idlgrpolygon', vert, polygon=poly, color=[255,0,0])

    xobjview, poly, background=[0,0,0], /block

    obj_destroy, poly
    obj_destroy, odf

end

