;; $Id$
;;
;*************************************************************************
;+
; NAME:
;       mas_mow__define
;
;
; PURPOSE:
;       
;       Reconstructs HARDI MRI voxel data using the Mixture of Wisharts
;       method outlined in reference [1].
;
; REFERENCES: 
;
;       1. Jian, Bing et al. “A novel tensor distribution model for the 
;          diffusion-weighted MR signal.” NeuroImage 37.1 (2007): 164-176.
;       
;       SEE ALSO:
;
;       2. Jian, B., and B.C. Vemuri. “A Unified Computational Framework 
;          for Deconvolution to Reconstruct Multiple Fibers From Diffusion
;          Weighted MRI.” Medical Imaging, IEEE Transactions on 26.11
;          (2007): 1464-1471.
;       
;       3. C. L. Lawson, R. J. Hanson. "Solving Least Square Problem", 
;          Prentice Hall, Englewood Cliffs NJ, 1974.
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
;        B_MATRIX:     The imaging B-matrix provided as n columns of the row:
;                      [ xx yy zz 2*xy 2*xz 2*yz ], where n is the number of
;                      diffusion-weighted acquisitions. The doubling
;                      of the off-diagonal elements is solely for
;                      compatibility with MAS's b-matrix
;                      storage scheme. The elements are divided by two
;                      before they are used in the computation
;
;        DEBUG:        Set this keyword to have extra debugging messages
;                      sent to the IDL console
;
;        RADIUS:       The desired diffusion radius
;
;        DIFF_TIME:    The diffusion time from the pulse sequence
;                      parameters. 
;
;        COMPUTATION_METHOD: One of the following:
;                      0  Damped least squares. The least squares projection
;                         matrix is regularized using tikhonov regularization
;                         to improve conditioning.
;                      1  Non-negative least squares. Uses method outlined in 
;                         reference [3].
;                      2  SVD pseudo-inverse. Uses SVD to compute a pseudoinverse.
;                         This method is the worst of the three
;                         
;        DECONVOLUTION_LEVEL: Tessellation order to use during weight estimation phase.
;                      0  Regular Icosahedron
;                      1  Level 1 tessellation of #0
;                      2  Level 2 tessellation of #0
;                      ...
;                      
;        EIGENVALUES: Eigenvalues to use for idealized tensor. Default is [1.5, 0.4, 0.4]
;                     from [1]
; OUTPUTS:
;
;
;
; METHODS:
;
;     TODO: Complete this section
;
; PROCEDURE:
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function mas_mow::init, $
   debug=debug, $
   diff_time=diff_time, $
   radius=radius, $
   computation_method=computation_method, $
   bval_thr=bval_thr, $
   b_matrix=b_matrix, $
   eigenvalues=eigenvalues, $
   deconvolution_level=deconvolution_level, $
   nnls_so_location=nnls_so_location

    self.can_use_weights = 1

    if (keyword_set(diff_time)) then begin
       self.diff_time = diff_time
    endif 
    
    if (keyword_set(debug)) then begin
       self.debug =  byte(debug)
    endif
    
    if (keyword_set(bval_thr)) then begin
        self.b_thr = float(bval_thr)
    endif else begin
        self.b_thr = 400. * 1e-3
    endelse
    
    if (keyword_set(radius)) then begin
       self.diff_radius = radius
    endif else begin
       self.diff_radius = 10.0
    endelse
 
    if (keyword_set(b_matrix)) then begin
       if (self.debug) then print, "mas_mow::init: setting b-matrix"
       self->setBMatrix, b_matrix
    endif
    
    if (keyword_set(eigenvalues)) then begin
       self.deco_evals = eigenvalues
    endif else begin
       self.deco_evals = [1.5D, 0.4D, 0.4D]
    endelse

    if (keyword_set(computation_method)) then begin
       self.computation_method = byte(computation_method)
    endif

    if (keyword_set(deconvolution_level)) then begin
       self.deco_level = byte(deconvolution_level)
    endif 

    if (keyword_set(nnls_so_location) && nnls_so_location ne '') then begin
        self.nnls_so_location = nnls_so_location
        self.have_nnls_so = 1
    endif
    
    self.deco_p = 2.0D

    self->prepareDeconvolution

  return, 1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function mas_mow::initializeMCMC, $
    data=data, $
    weights=weights, $
    radius=radius
 
    mcmc = obj_new('mas_mow_MCMC_simulator', $
                   weights=weights, $
                   mow_object=self)
    
    return, mcmc

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function mas_mow::getDisplacementProbability, $
    data=data, $
    weights=weights, $
    radius=radius, $
    direction=direction

    if (not keyword_set(weights)) then begin
       if (n_elements(data) ne 0) then begin
          weights = self->computeWeights(data)
       endif else begin
          message, 'mas_mow::getDisplacementProbability: weights or data needed'
          return, -1
       endelse
    endif 

    if (not keyword_set(radius)) then begin
       radius = self.diff_radius
    endif else begin
       self.diff_radius = radius
    endelse

    w  = weights/sqrt( (4.0*!DPI*self.diff_time)^3 * (*self.deco_D_inverse)[6,*] ) 
    
    n_dirs = n_elements(direction)/3

    rt = dblarr(6,self.deco_nverts)

    sum = dblarr(n_dirs)

    for j = 0, n_dirs-1 do begin

       r = direction[*,j] * radius
       
       rt[0,*] = r[0]*r[0]
       rt[1,*] = r[1]*r[1]
       rt[2,*] = r[2]*r[2]
       rt[3,*] = r[2]*r[1]
       rt[4,*] = r[0]*r[1]
       rt[5,*] = r[0]*r[2]
       
       rt_d_r = total( (*self.deco_D_inverse)[0:5,*] * rt, 1)
       
       exp_arg = -(rt_d_r/(4.0 * self.diff_time))
       
       sum[j] = total(w * exp(exp_arg)) 

    endfor

    return, sum

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function mas_mow::computeWeights, data, residuals=residuals

    if (n_elements(*self.b_low) ne 1) then begin
       low_b_data = (moment(data[*self.b_low]))[0]
    endif else begin
       low_b_data = (moment([data[*self.b_low],data[*self.b_low]]))[0]
    endelse

    high_b_data = double(reform(data[*self.b_high])/low_b_data)
    orig_high_b = high_b_data * low_b_data
    
    case self.computation_method of

       0: begin                 ; DLS
          w = matrix_multiply(*self.DLS_matrix, high_b_data, /atranspose)
       end
       1: begin                 ; NNLS ;; make_dll, 'nnls', 'mas_nnls', input_directory=so, output_Directory=so
          resid = 0d
          if (self.have_nnls_so) then begin
              A = *self.A_matrix
              m = (size(A, /dimensions))[0]; = mda in nnls.c
              n = (size(A, /dimensions))[1]
              w = dblarr(n+1)
              x = w ;; same size
              mode = 0L
              rnorm = 0D
              index = lonarr(n+2)
              zz = dblarr(m+1)
              ok = call_external(self.nnls_so_location, 'nnls_c',$
                                 A, m, m, n, high_b_data, x, rnorm, w, zz, index, mode,$
                                 /cdecl, /auto_glue)
              w = x[0:n-1]
          endif else begin
              mas_nnls, *self.A_matrix, high_b_data, w, resid=resid, resnorm=resnorm
          endelse
          ;w = reform(w)
          ;negs = where(w lt 0, neg_ct)
          ;if (neg_ct gt 0) then begin
          ;   w[negs] = 0.0
          ;endif
       end
       else: begin              ; SVD
          w = transpose(mas_nnls_pinv(*self.A_matrix, junk, /double)) ## high_b_data
          w = transpose(w)
       end

    endcase
    
    if (arg_present(residuals)) then begin
        ;;residuals = high_b_data;; (high_b_data - (*self.A_matrix # w)) * low_b_data
        residuals = orig_high_b - ((*self.A_matrix # w) * low_b_data)
    endif

    return, w

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mas_mow::prepareDeconvolution

    if (self.debug) then begin
       print, "mas_mow::prepareDeconvolution: freeing old pointers"
    endif

    ptr_free, self.deco_vertices
    ptr_free, self.DLS_matrix 
    ptr_free, self.A_matrix 
    ptr_free, self.deco_D_inverse

    if (self.debug) then begin
       print, "mas_mow::prepareDeconvolution: preparing deconvolution"
    endif

    ;; Theoretical values for p and evals (from paper)
    evals = diag_matrix(self.deco_evals * 1.0)
    p     = self.deco_p

    if (self.debug) then begin
       print, "mas_mow::prepareDeconvolution: making vertices, level", $
              self.deco_level
    endif

    mas_tessellator_make_icos, level=self.deco_level, $
                              vertlist=vlist

    nv_deco = (size(vlist, /dimensions))[1]

    self.deco_vertices = ptr_new(vlist, /no_copy)
    self.deco_nverts   = nv_deco

    nb_high = n_elements(*self.b_high)

    D_inverse = dblarr(7, nv_deco)

    A_Mat = dblarr(nb_high, nv_deco)

    pi_frac = 3D
    tx1 = double([ [ cos(!pi/pi_frac), -sin(!pi/pi_frac), 0 ], $
                   [ sin(!pi/pi_frac),  cos(!pi/pi_frac), 0 ], $
                   [                0,                 0, 1 ] ])

    tx2 = double([ [ cos(!pi/pi_frac), 0, sin(!pi/pi_frac) ],$
                   [                0, 1,                0 ], $
                   [-sin(!pi/pi_frac), 0, cos(!pi/pi_frac) ] ])
    tx  = tx1 # tx2

    ;; for each test vertex, compute the entries of A using Eqn 15
    for i = 0, nv_deco-1 do begin
        
        vert = reform((*self.deco_vertices)[*,i])
        vert /= norm(vert)
        tmp = reform(vert) # tx
        tmp1 = crossp(tmp, reform(vert))
        tmp1 = tmp1/norm(tmp1)
        tmp2 = crossp(tmp1, reform(vert))
        tmp2 = tmp2/norm(tmp2)

        ;; orthogonal matrix ``Q''
        Q_mat = transpose([ [reform(vert)] , [tmp1] , [tmp2] ])
        
        ;; this is 'capital sigma' from the paper (= D_i/p)
        Sig_i = (Q_mat # evals # transpose(Q_mat)) / p

        ;; the inverse is used for reconstruction (used in eqn 16)
        D_i_tmp = invert(Sig_i*p)
        D_inverse[0:5,i] = reform([D_i_tmp[0,0],   D_i_tmp[1,1],   D_i_tmp[2,2], $
                                 2*D_i_tmp[2,1], 2*D_i_tmp[1,0], 2*D_i_tmp[2,0] ])

        ;; save the determinant for later should be the same for all
        ;; D's. We should just use product of evals.
        D_inverse[6,i] = determ(Sig_i*p)

        for j = 0, nb_high-1 do begin
            b_ind = (*self.b_high)[j]
            B_jtmp = reform((*self.b_matrix)[*,b_ind])
            B_jtmp[3:5] /= 2.0
            B_j = [ [ B_jtmp[0], B_jtmp[3], B_jtmp[4] ], $
                    [ B_jtmp[3], B_jtmp[1], B_jtmp[5] ], $
                    [ B_jtmp[4], B_jtmp[5], B_jtmp[2] ] ]

            A_mat[j,i] = (1.0 + trace(B_j # Sig_i))^(-p)
        endfor

    endfor

    ;; tikhonov preconditioner for damped least squares
    lambda = 0.6D
    tik_mat = diag_matrix(replicate(lambda,nb_high))
    DLS_mat = invert( (matrix_multiply(a_mat, a_mat, /btranspose) + tik_mat) ) # a_mat
    
    self.DLS_matrix = ptr_new(DLS_mat, /no_copy)
    self.A_matrix   = ptr_new(A_mat, /no_copy)
    self.deco_D_inverse = ptr_new(D_inverse, /no_copy)

    self.deco_prepared = 1B

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mas_mow::setBMatrix, b_matrix

    if (self.debug) then begin
       print, "mas_mow::setBMatrix: b_matrix is", size(b_matrix, /dimensions)
    endif

    ptr_free, self.b_matrix
    ptr_free, self.b_low
    ptr_free, self.b_high

    self.b_matrix = ptr_new(b_matrix)
    
    sz_bmatrix = size(b_matrix, /dimensions)

    bvals = reform(total(b_matrix[0:2, *], 1))

    low_b = where(bvals lt self.b_thr, complement=high_b, count)

    if (count eq 0) then begin
       message, 'There needs to be at least one low b-value acquisition.'
       return
    endif

    if (self.debug) then begin
       print, "mas_mow::setBMatrix: low_idx :", low_b
       print, "mas_mow::setBMatrix: high_idx:", high_b
    endif

    self.b_high   = ptr_new(high_b, /no_copy)
    self.b_low    = ptr_new(low_b, /no_copy)
    self.b_values = ptr_new(bvals, /no_copy)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mas_mow::setProperty, $
    diff_time=diff_time, $
    computation_method=computation_method, $
    deconvolution_level=deconvolution_level, $
    debug_mode=debug_mode, $
    b_matrix=b_matrix

    recompute_decon = 0
  
    if (keyword_set(diff_time)) then $
       self.diff_time = double(diff_time)

    if (keyword_set(debug_mod)) then self.debug = byte(debug_mode)

    if (keyword_set(computation_method)) then $
       self.computation_method = byte(computation_method)
    
    if (keyword_set(deconvolution_level)) then begin
       self.deco_level = deconvolution_level
       recompute_decon = 1
    endif
    
    if (keyword_set(b_matrix)) then begin
       self->setBMatrix, b_matrix
       recompute_decon = 1
    endif

    if (recompute_decon and self.deco_prepared) then begin
       self->prepareDeconvolution
    endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mas_mow::getProperty, $
   diff_time=diff_time, $
   can_use_weights=can_use_weights, $
   computation_method=computation_method, $
   deconvolution_level=deconvolution_level, $
   debug_mode=debug_mode, $
   b_matrix=b_matrix
  
   diff_time = self.diff_time
   computation_method = self.computation_method
   deconvolution_level = self.deco_level
   debug_mode = self.debug
   can_use_weights = self.can_use_weights

   if (ptr_valid(self.b_matrix)) then  b_matrix = *self.b_matrix

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mas_mow::cleanup

    if (self.debug) then print, 'mas_mow::cleanup: cleaning up.'

    ptr_free, self.b_matrix
    ptr_free, self.b_values
    ptr_free, self.b_low
    ptr_free, self.b_high
    ptr_free, self.deco_vertices
    ptr_free, self.deco_D_inverse
    ptr_free, self.DLS_matrix
    ptr_free, self.A_matrix

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mas_mow__define

    struct = { MAS_MOW, $
               have_nnls_so: 0B, $
               nnls_so_location: '', $
               can_use_weights: 0B, $
               debug: 0B, $
               b_thr: float(0), $
               b_matrix: ptr_new(), $
               b_values: ptr_new(), $
               b_low:    ptr_new(), $
               b_high:   ptr_new(), $
               deco_prepared: 0B, $
               deco_level: 2B, $
               deco_vertices: ptr_new(), $
               deco_nverts: 0L, $
               deco_D_inverse: ptr_new(), $
               deco_p: 2D, $
               deco_evals: dblarr(3), $
               diff_time: 0D, $
               diff_radius: 0D, $
               computation_method: 0B, $
               DLS_matrix: ptr_new(), $
               A_matrix: ptr_new() $
             }
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mas_mow_bulk_reconstructor::setReconstructionDirections, $
   directions, $
   spherical=spherical, $
   degrees=degrees
  
    our_dirs = directions

    if (n_elements(directions) eq 0) then begin
       message, 'Directions must have at least one direction'
    endif
    
;    if (not keyword_set(spherical)) then begin
;       our_dirs = cv_coord(from_rect=our_dirs, /to_sphere)
;       our_dirs[1,*] = !DPI/2 - our_dirs[1,*]
;    endif else if (keyword_set(degrees)) then begin
;       our_dirs[1,*] *= !DTOR
;       our_dirs[0,*] *= !DTOR
;    endif

    ptr_free, self.reco_vertices
    self.reco_vertices = ptr_new(our_dirs)

    ;; precompute the arena in which the computation will take place
    ;; to reconstruct the probability values. The "temp" array will be
    ;; multiplied by the deconvolution weights, then summed to get
    ;; the final value. 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    rt = dblarr(6,self.deco_nverts)
    n_dirs = n_elements(our_dirs)/3
    temp = dblarr(self.deco_nverts, n_dirs)

    for j = 0, n_dirs-1 do begin

       r = our_dirs[*,j] * self.diff_radius
       
       rt[0,*] = r[0]*r[0]
       rt[1,*] = r[1]*r[1]
       rt[2,*] = r[2]*r[2]
       rt[3,*] = r[2]*r[1]
       rt[4,*] = r[0]*r[1]
       rt[5,*] = r[0]*r[2]
       
       rt_d_r = total( (*self.deco_D_inverse)[0:5,*] * rt, 1)
       
       exp_arg = -(rt_d_r/(4.0 * self.diff_time))
       temp[*,j] = exp(exp_arg)

    endfor
    ptr_free, self.reco_matrix
    self.reco_matrix = ptr_new(transpose(temp), /no_copy)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function mas_mow_bulk_reconstructor::getDisplacementProbability, $
   data=data, $
   weights=weights, $
   radius=radius

    if (not keyword_set(weights)) then begin
        weights = self->computeWeights(data)
    endif
    
    ;; this computation is undone when the min-max scaling happens, since this is
    ;; just a scalar.
    w  = weights/sqrt( (4.0*!DPI*self.diff_time)^3 * (*self.deco_D_inverse)[6,*] ) 
    ;; we need to "replicate" the weights into an array the same size as the computation
    ;; arena created when the reconstruction directions are set. Congrid allows this to
    ;; happen very quickly (takes an nx1 array and creates an nxm array where each column m
    ;; contains the original nx1 array)
    a = congrid(transpose(w), n_elements(*self.reco_vertices)/3, n_elements(w), interp=0)
    return, total(a * *self.reco_matrix, 2)
    
;    return, self->mas_mow::getDisplacementProbability(data=data, $
;                                                      weights=weights, $
;                                                      radius=radius, $
;                                                      direction=*self.reco_vertices)

end

function mas_mow_bulk_reconstructor::makeBootstrap, $
   data=data

    weights = self->computeWeights(data, residuals=residuals)
    
;;    H_diag = dblarr(n_elements(e))
;;    for i = 0, n_elements(e)-1 do begin
;;        H_diag[i] = H[i,i]
;;    endfor
    
;;    refit = H ## e
;;    residuals = real_part(reform(e - refit))
;;    rhat = residuals/sqrt(1.0 - H_diag)

    ptr_free, self.bootstrap_rhat
    self.bootstrap_rhat = ptr_new(residuals)
    return, 1
    
end

function mas_mow_bulk_reconstructor::sampleBootstrap, seed, data=data

    n_e = n_elements((*self.b_high))
    samp = round( randomu(seed, n_e) * n_e)
    data_new = data
    data_new[*self.b_high] += (*self.bootstrap_rhat)[samp]
    return, self->getDisplacementProbability(data=data_new)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mas_mow_bulk_reconstructor::cleanup

    ptr_free, self.reco_vertices
    ptr_free, self.reco_matrix
    self->mas_mow::cleanup

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mas_mow_bulk_reconstructor__define

    struct = { MAS_MOW_BULK_RECONSTRUCTOR, $
               inherits MAS_MOW, $
               reco_matrix: ptr_new(), $
               reco_vertices: ptr_new(), $
               bootstrap_rhat: ptr_new() }

end

pro mas_mow_example

    ;; Note that this b-matrix was taken from MAS, and so the
    ;; off-diagonal elements are multiplied by two.
    b_matrix = [[ 0.00000, 0.00000, 800.000, 0.00000, 0.00000, 0.00000 ], $
                [ 639.389, 0.00000, 159.847, 0.00000, 639.389, 0.00000 ], $
                [ 60.9408, 579.361, 159.847, -375.802, 197.395, -608.635 ], $
                [ 419.341, 221.341, 159.847, 609.318, -517.805, -376.195 ], $
                [ 419.341, 221.341, 159.847, -609.318, -517.805, 376.195 ], $
                [ 60.9408, 579.361, 159.847, 375.802, 197.395, 608.635 ], $
                [ 104.257, 0.00000, 696.391, 0.00000, 538.901, 0.00000 ], $
                [ 9.85680, 94.1192, 696.391, -60.9168, 165.701, -512.030 ], $
                [ 68.2112, 35.9552, 696.391, 99.0464, -435.898, -316.474 ], $
                [ 68.2112, 35.9552, 696.391, -99.0464, -435.898, 316.474 ], $
                [ 9.85680, 94.1192, 696.391, 60.9168, 165.701, 512.030 ], $
                [ 192.865, 101.959, 505.620, -280.459, 624.552, -454.104 ], $
                [ 28.2752, 266.343, 505.620, 173.562, -239.136, -733.944 ], $
                [ 294.759, 0.00000, 505.620, -0.00000, -772.104, 0.00000 ], $
                [ 28.2752, 266.343, 505.620, -173.562, -239.136, 733.944 ], $
                [ 192.865, 101.959, 505.620, 280.459, 624.552, 454.104 ], $
                [ 362.343, 0.00000, 438.080, 0.00000, 796.832, 0.00000 ], $
                [ 34.6112, 327.680, 438.080, -212.992, 246.272, -757.760 ], $
                [ 236.749, 125.453, 438.080, 344.678, -644.096, -468.864 ], $
                [ 236.749, 125.453, 438.080, -344.678, -644.096, 468.864 ], $
                [ 34.6112, 327.680, 438.080, 212.992, 246.272, 757.760 ], $
                [ 491.725, 94.1192, 213.831, -430.259, 648.525, -283.730 ], $
                [ 259.009, 327.680, 213.831, -582.656, 470.677, -529.408 ], $
                [ 5.64480, 580.723, 213.831, 114.509, -69.4848, -704.774 ], $
                [ 149.991, 436.897, 213.831, 511.979, -358.178, -611.301 ], $
                [ 559.117, 27.0848, 213.831, 246.118, -691.539, -152.205 ], $
                [ 559.117, 27.0848, 213.831, -246.118, -691.539, 152.205 ], $
                [ 149.991, 436.897, 213.831, -511.979, -358.178, 611.301 ], $
                [ 5.64480, 580.723, 213.831, -114.509, -69.4848, 704.774 ], $
                [ 259.009, 327.680, 213.831, 582.656, 470.677, 529.408 ], $
                [ 491.725, 94.1192, 213.831, 430.259, 648.525, 283.730 ], $
                [ 505.620, 266.343, 28.2752, -733.944, 239.136, -173.562 ], $
                [ 73.4472, 697.885, 28.2752, 452.803, -91.1424, -280.947 ], $
                [ 771.459, 0.00000, 28.2752, -0.00000, -295.386, 0.00000 ], $
                [ 73.4472, 697.885, 28.2752, -452.803, -91.1424, 280.947 ], $
                [ 505.620, 266.343, 28.2752, 733.944, 239.136, 173.562 ], $
                [ 744.980, 35.9552, 19.4688, -327.328, 240.864, -52.9152 ], $
                [ 200.000, 580.723, 19.4688, -681.600, 124.800, -212.659 ], $
                [ 7.52720, 773.031, 19.4688, -152.562, 24.2112, -245.357 ], $
                [ 344.269, 436.897, 19.4688, 775.654, -163.738, -184.454 ], $
                [ 655.220, 125.453, 19.4688, 573.408, -225.888, -98.8416 ], $
                [ 655.220, 125.453, 19.4688, -573.408, -225.888, 98.8416 ], $
                [ 344.269, 436.897, 19.4688, -775.654, -163.738, 184.454 ], $
                [ 7.52720, 773.031, 19.4688, 152.562, 24.2112, 245.357 ], $
                [ 200.000, 580.723, 19.4688, 681.600, 124.800, 212.659 ], $
                [ 744.980, 35.9552, 19.4688, 327.328, 240.864, 52.9152 ], $
                [ 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000 ] ]

    b_matrix *= 1e-3

    ;; diffusion time
    dt = 17.2

    ;; some diffusion data
    data = [ 94971.8, 96982.2, 68006.5, 82847.6, 60716.6, 74985.8, $
             83197.1, 104232., 96661.7, 84150.5, 84722.3, 94969.2, $
             99566.2, 83487.4, 72094.5, 87053.2, 92669.6, 95042.6, $
             103382., 84024.1, 84522.4, 86291.2, 66306.5, 81414.1, $
             101041., 89636.8, 76524.8, 66982.3, 66651.2, 86889.2, $
             92382.5, 70348.6, 80515.4, 82981.2, 73442.0, 99426.0, $
             79545.3, 71345.7, 74289.4, 85795.4, 101870., 67104.5, $
             59460.6, 76771.4, 77757.9, 85816.7, 150168. ]

    ;; mow object
    mow = obj_new('mas_mow', b_matrix=b_matrix, debug=1, $
                  diff_time=dt, radius=12.0, computation_method=1, $
                  deconvolution_level=1)

    ;; get some probabilities
    print, "--> Using Data"
    print, "PR[1,0,0]: ", mow->getDisplacementProbability(data=data, $
                                                          direction=[1.,0,0])
    print, "PR[0,1,0]: ", mow->getDisplacementProbability(data=data, $
                                                          direction=[0,1.,0])
    print, "PR[0,0,1]: ", mow->getDisplacementProbability(data=data, $
                                                          direction=[0,0,1.])
    ;; exploit the ease and computational effectiveness of using weights
    print, "--> Using Weights"
    w = mow->computeWeights(data)
    print, "PR[1,0,0]: ", mow->getDisplacementProbability(weights=w, $
                                                          direction=[1.,0,0])
    print, "PR[0,1,0]: ", mow->getDisplacementProbability(weights=w, $
                                                          direction=[0,1.,0])
    print, "PR[0,0,1]: ", mow->getDisplacementProbability(weights=w, $
                                                          direction=[0,0,1.])
    obj_destroy, mow
    
    ;; special bulk reconstructor object
    mow = obj_new('mas_mow_bulk_reconstructor', b_matrix=b_matrix, $
                  diff_time=dt, radius=12.0, computation_method=1, $
                 deconvolution_level=1)

    ;; set the directions
    mesh_obj, 4, vert, poly, fltarr(50,50)+1
    mow->setReconstructionDirections, vert

    ;; get the bulk probabilites
    Pr  = mow->getDisplacementProbability(data=data)

    ;; minmax to enhance shape
    Pr  = (Pr - min(Pr))/(max(Pr)-min(Pr))

    ;; distort sphere to create direction shape
    vert[0,*] *= Pr
    vert[1,*] *= Pr
    vert[2,*] *= Pr

    poly = obj_new('idlgrpolygon', vert, polygon=poly, color=[255,0,0])

    xobjview, poly, background=[0,0,0], /block

    obj_destroy, poly
    obj_destroy, mow
    
end

