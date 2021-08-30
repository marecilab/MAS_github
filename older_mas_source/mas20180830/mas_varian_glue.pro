;+
; :Description:
;    Handles the reading of procpar parameters and setting up
;    the imnd array for a 1D agilent experiment
;
; :Params:
;    opp: mas_varian_procpar object
;
;
;
; :Author: wtriplett
;-
function mas_varian_procpar_to_imnd_1d, opp

    common scan_data, project
    
    next_index = project.ni
    imnd = project.imndarray[next_index]
    
    ;;;;;;;;;;;;;;;;;;;; Compute dimensions of data array ;;;;;;;;;;;;;;;;;;;;;
    names = [ 'np', 'arraydim', $
              'layout', 'sfrq', 'at', 'tn', 'sw', $
              'nt', 'seqfil', 'pslabel', 'oversample']
    error = bytarr(n_elements(names))
    par = 0

    np         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    arraydim   = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    layout     = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    srfq       = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    at         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    tn         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    sw         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_avgs   = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd  
    seqfil     = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    pslabel    = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
   
    imnd.spect_nucleus        = strtrim(tn,2)
    imnd.spect_acq_size       = long(np)/2
    ;;imnd.spect_pulse_length   = float(temp)
    imnd.spect_num_avg        = long(num_avgs[0])
    imnd.spect_bf1            = float(srfq)
   
    if (total(error) ne 0) then begin
        for i=0, n_elements(names)-1 do begin
            if (error[i] eq 1) then begin
                print, "Error: parameter "+names[i]+" not found!"
            endif
        endfor
        print, "procpar file possibly corrupted."
        return, -1
    endif

    imnd.fdim = np/2
    imnd.adim = arraydim
    imnd.n_avg       = num_avgs[0]
    imnd.scan_name   = strupcase(seqfil) ;; 'ONEPULSE' ;; (pslabel eq seqfil) ? pslabel : pslabel+' ('+seqfil+')'
    imnd.orientation = [ "----", "----","----" ]
    scantime         = opp->lookup('time_run')
    ;;1          2         3          4         5         6           7         8         9            10
    ;;2010       07        20     T   14        52        55     X    02        PM        Tue          Jul
    timebits = stregex(scantime, "^([0-9]{4})([0-9]{2})([0-9]{2})T([0-9]{2})([0-9]{2})([0-9]{2})$", $
                       /subexpr, /extract)
             
    imnd.scan_date = timebits[1] + '-' + timebits[2] + '-' + timebits[3]+ ','+timebits[4]+':'+timebits[5] + ':' + timebits[6]
    imnd.dimensions = 1
                             
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; TR CHECK                                                                      ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    
    tr = opp->lookup('tr', is_array=ia, not_found=nf)
    tr *= 1000.0
    imnd.spect_rep_time=float(tr[0])
    if (ia ne 0) then begin
        imnd.rep_time_ptr = ptr_new(tr)
    endif

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; TE CHECK                                                                      ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    te = opp->lookup('te', is_array=ia, not_found=nf)
    te *= 1000.0
    if (ia ne 0) then begin
        imnd.echo_Time_ptr = ptr_new(te)
        imnd.echo_time = te[0]
    endif else begin
        imnd.echo_time = te[0]
    endelse

    if pslabel eq 'cpmgecho' then begin
      imnd.spect_spectral_width = 1000/float(te[0])
      imnd.spect_acq_time       = float(te[0])*1000.0*np/2
    endif else begin
      imnd.spect_spectral_width = float(sw)
      imnd.spect_acq_time       = float(at)*1000.0
    endelse


    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; TI CHECK                                                                      ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ti = opp->lookup('ti', is_array=ia, not_found=nf)
    ti *= 1000.0
    if (ia ne 0) then begin
        imnd.inversion_time_ptr = ptr_new(ti)
    endif
    
    imnd.file_Path    = opp->lookup('procpar_file')
    imnd.display_Name = opp->lookup('parent_dir')
    imnd.image_type   = 99
    
    acq_matrix = ptr_new(diag_matrix(replicate(1.0D,4)))
    imnd.acq_matrix = acq_matrix
    
    return, imnd

end


;+
; :Description:
;    Handles converting / reading an Agilent procpar file and 
;    loading the parameter values into an imnd array.
;
; :Params:
;     opp: mas_varian_procpar object
;
;
; :Author: wtriplett
; Modifications
; - Incorporated VTI processing for EPI acquisitions (Magdoom, 10/18/15)
; - Fixed TE parameter for RARE and fat suppression parameter (Magdoom, 12/24/15)
; - Scan date corrected to time when the scan was actually run rather than the time it was submitted

function mas_varian_procpar_to_imnd, opp

    common scan_data, project
    
    next_index = project.ni
    imnd = project.imndarray[next_index]
    imnd.spect_bf1 = opp->lookup('H1reffrq')  ; MHz
    imnd.B0 = 1e-4*(opp->lookup('B0'))         ; Tesla
    imnd.gcoil =opp->lookup('gcoil')          
    imnd.max_slew_rate = float(opp->lookup('gmax'))/(float(opp->lookup('trise'))*1e6) 
    ;;;;;;;;;;;;;;;;;;;; Compute dimensions of data array ;;;;;;;;;;;;;;;;;;;;;
    names = [ 'seqcon', 'ni', 'np', 'nf', 'arraydim', $
              'ns', 'nD', 'ne', 'nv', 'nv2', 'layout', $
              'nt', 'seqfil', 'pslabel','flip1','fsat','tpe' ]
    error = bytarr(n_elements(names))
    par = 0
    seqcon     = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    ni         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    np         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nf         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    arraydim   = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_slices = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_dims   = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_echos  = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nv         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nv2        = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    layout     = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_avgs   = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd  
    seqfil     = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    pslabel    = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    alpha      = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    fatsup     = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    
    if (total(error) ne 0) then begin
        for i=0, n_elements(names)-1 do begin
            if (error[i] eq 1) then begin
                print, "Error: parameter "+names[i]+" not found!"
            endif
        endfor
        print, "procpar file possibly corrupted."
        return, -1
    endif

    seqcon = string(transpose(byte(seqcon)))    
    data_dims = opp->computeDataDimensions()
    imnd.fdim = data_dims[0]
    imnd.pdim = data_dims[1]
    imnd.sdim = data_dims[2]
    if n_elements(data_dims) gt 4 then imnd.adim = data_dims[3]*data_dims[4] else imnd.adim = data_dims[3]

    oversample = opp->lookup('oversample', not_found=ntfd)
    if (ntfd) then begin
        oversample = 1
    endif
    imnd.fdim *= oversample
    
    ;; For kspace span resolution reduction
    imnd.k_fdim_span_max = (imnd.k_fdim_span = imnd.fdim)
    imnd.k_pdim_span_max = (imnd.k_pdim_span = imnd.pdim)
    imnd.k_sdim_span_max = (imnd.k_sdim_span = imnd.sdim)
    
    imnd.n_echo      = num_echos
    imnd.n_avg       = num_avgs[0]
    imnd.pn_avg      = ptr_new(num_avgs)
    imnd.alpha       = alpha
    imnd.scan_name   = (pslabel eq seqfil) ? pslabel : pslabel+' ('+seqfil+')'
    imnd.orientation = [ opp->lookup('orient'), "----","----" ]
    scantime         = opp->lookup('time_run')
    ;;1          2         3          4         5         6           7         8         9            10
    ;;2010       07        20     T   14        52        55     X    02        PM        Tue          Jul
    timebits = stregex(scantime, "^([0-9]{4})([0-9]{2})([0-9]{2})T([0-9]{2})([0-9]{2})([0-9]{2})", $
                       /subexpr, /extract)
    imnd.scan_date = timebits[1]+'-'+timebits[2]+'-'+timebits[3]+' @ '+$
                     timebits[4]+':'+timebits[5]+':'+timebits[6]

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; EPI Reference Scan number CHECK                                               ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    
    imnd.epi_num_ref_scans = 0L
    if (seqfil eq 'epi' or seqfil eq 'epip') then begin
        ref_type = opp->lookup('epiref_type', not_found=nf)
        if (nf eq 0) then begin
            case ref_type of
                'single': imnd.epi_num_ref_scans = 1L
                'triple': imnd.epi_num_ref_scans = 3L
                'fulltriple' : imnd.epi_num_ref_scans = 3L
                else: begin
                    message, /info, 'Unknown ref scan type: '+ref_type
                end
            endcase
        endif
    endif
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; TR CHECK                                                                      ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    
    tr = opp->lookup('tr', is_array=ia, not_found=nf)
    tr *= 1000.0
    if (ia ne 0) then begin
        imnd.rep_time_ptr = ptr_new(tr)
        imnd.recov_time = tr[0]
    endif else begin
        imnd.recov_time = tr[0]
    endelse

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; TE CHECK                                                                      ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if strmid(seqfil, 0, 3) eq 'fse' then begin
        te = opp->lookup('effte', is_array=ia, not_found=nf) 
        if strcmp(seqfil,'fsemsdw') eq 1 then te = opp->lookup('te1', is_array=ia, not_found=nf)
     endif else if strcmp(layout,'mge3d') eq 1 then te = opp->lookup('TE', is_array=ia, not_found=nf)/1000 else $ 
      te = opp->lookup('te', is_array=ia, not_found=nf)
      
    te *= 1000.0
    if (ia ne 0) then begin
        imnd.echo_Time_ptr = ptr_new(te)
        imnd.echo_time = te[0]
    endif else begin
        imnd.echo_time = te[0]
    endelse

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; TI CHECK                                                                      ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ti = opp->lookup('ti', is_array=ia, not_found=nf)
    ti *= 1000.0
    if seqfil eq 'epi' or seqfil eq 'epip' then begin
      if ref_type eq 'triple' or ref_type eq 'fulltriple' and n_elements(ti) ge 2 then ti = rebin(ti[2:n_elements(ti)-1],(arraydim-2)/2) 
    endif     
    if (ia ne 0) then begin
        imnd.inversion_time_ptr = ptr_new(ti)
    endif
    
    ;; I believe that there is a better way to identify these. Possibly using their
    ;; seqcon, since they both seems to start with 'c'.
    if (strmid(seqfil, 0, 4) eq 'mems' or $
        strmid(seqfil, 0, 5) eq 'mgems') then begin
        te = opp->lookup('TE', is_array=ia)
        if (ia ne 0) then begin
            imnd.echo_Time_ptr = ptr_new(te)
        endif
    endif

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; RARE CHECK                                                                    ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    etl = opp->lookup('etl', not_found=ntfd)
    if ntfd ne 1 then imnd.rare = etl
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; FOV CHECK                                                                     ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; Euler angles
    imnd.theta = ptr_new(opp->lookup('theta'))
    imnd.phi =  ptr_new(opp->lookup('phi'))
    imnd.psi =  ptr_new(opp->lookup('psi'))   
    imnd.f_fov = opp->lookup('lro')
    imnd.p_fov = opp->lookup('lpe')
    imnd.f_fov_offset = opp->lookup('pro')
    imnd.p_fov_offset = opp->lookup('ppe')
    if (num_dims eq 2) then begin
        imnd.s_fov = opp->lookup('lss')
        imnd.s_fov_offset = opp->lookup('pss0')
        if (imnd.s_fov eq 0) then begin
            imnd.s_fov = opp->lookup('thk')*opp->lookup('ns')/10.0
        endif
        if seqfil eq 'epip' then imnd.fdim_shift = round(opp->lookup('pro')/(imnd.f_fov*oversample)*imnd.fdim) 
        imnd.pdim_shift = round(opp->lookup('ppe')*imnd.pdim/imnd.p_fov )
    endif else if (num_dims eq 3) then begin
        imnd.s_fov = opp->lookup('lpe2')
        imnd.s_fov_offset = opp->lookup('ppe2')
        if seqfil eq 'epip' then imnd.fdim_shift = round(opp->lookup('pro')/(imnd.f_fov*oversample)*imnd.fdim)
        imnd.pdim_shift = round(opp->lookup('ppe')/imnd.p_fov * imnd.pdim)
        imnd.sdim_shift = -(round(imnd.sdim/2.0 + (opp->lookup('ppe2')/imnd.s_fov)*imnd.sdim)) mod imnd.sdim
         ;; Note this is n/2 fft correction amount
    endif
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; FLOW CHECK                                                                    ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    imnd.flow = opp->lookup('flow')
    imnd.diff = opp->lookup('diff')
    if strcmp(imnd.flow,'y') eq 1 then begin
      if strcmp(imnd.diff,'y') eq 1 then begin
        imnd.gflow = opp->lookup('gdiff')
        imnd.tflow = opp->lookup('tdelta')
        imnd.bptype = 't'
        imnd.fc = 'n'
        imnd.scheme = 'h'
        imnd.fro = ptr_new(-transpose(opp->lookup('dro')))
        imnd.fpe = ptr_new(-transpose(opp->lookup('dpe')))
        imnd.fsl = ptr_new(-transpose(opp->lookup('dsl')))
      endif else begin
        imnd.gflow = opp->lookup('gflow')
        imnd.tflow = opp->lookup('tflow')
        imnd.bptype = opp->lookup('bptype')
        imnd.fc = opp->lookup('fc')
        imnd.scheme = opp->lookup('scheme')
        if opp->lookup('ncycl') ne -1 then imnd.ncycl = opp->lookup('ncycl')
        if opp->lookup('tDELTA') ne -1 and opp->lookup('tDELTA') ne 0 then imnd.big_delta = ptr_new(opp->lookup('tDELTA') * 1000.0)
        if isa(opp->lookup('pol_r'),/array) then begin
           imnd.fro = ptr_new(transpose(opp->lookup('pol_r')))
           imnd.fpe = ptr_new(transpose(opp->lookup('pol_p')))
           imnd.fsl = ptr_new(transpose(opp->lookup('pol_s')))
        endif
     endelse
   endif 
     
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; EDDY CURRENT MAPPING CHECK                                                    ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if ~isa(opp->lookup('tgraxis'),/number) then begin
      imnd.scheme = opp->lookup('scheme')
      tgrdelay = opp->lookup('tgrdelay')
      tgrdelay += opp->lookup('tpe')     
      if n_elements(tgrdelay) gt 1 then begin
        if imnd.scheme eq 'h' or imnd.scheme eq 'f' then imnd.tgrdelay = ptr_new(1e3*tgrdelay[0:n_elements(tgrdelay)-1:4]) else $
         if imnd.scheme eq 's' then imnd.tgrdelay = ptr_new(1e3*tgrdelay[0:n_elements(tgrdelay)-1:6]) else if n_elements(tgrdelay) mod 2 eq 0 then imnd.tgrdelay = ptr_new(1e3*tgrdelay[0:n_elements(tgrdelay)-1:2]) $
          else imnd.tgrdelay = ptr_new(1e3*tgrdelay[1:n_elements(tgrdelay)-1:2])
      endif
      if isa(opp->lookup('fpdelay')) then begin
        tfpdelay = opp->lookup('fpdelay')
        if n_elements(tfpdelay) gt 1 then begin
          if n_elements(tfpdelay) mod 2 eq 0 then imnd.fpdelay = ptr_new(1e3*tfpdelay[0:n_elements(tfpdelay)-1:2]) $
          else imnd.fpdelay = ptr_new(1e3*tfpdelay[1:n_elements(tfpdelay)-1:2])
          endif
      endif
      
      gradAmp = [[abs((opp->lookup('pol_r'))[0])], [abs((opp->lookup('pol_p'))[0])], [abs((opp->lookup('pol_s'))[0])]]*opp->lookup('tgramp')
      
      imnd.tgraxis  = opp->lookup('tgraxis') 
      imnd.tgrshape = opp->lookup('tgrshape') 
      imnd.tgramp = ptr_new(gradAmp)
      imnd.tgrtime  = 1e3*opp->lookup('tgrtime')      
      case imnd.tgrshape of
        'a' :  begin
               imnd.tasc  = 1e3*opp->lookup('tasc')
               imnd.tdesc = 1e3*opp->lookup('tdesc')
               end
        'f'  : begin
               imnd.BWFresnel = opp->lookup('BWFresnel')
               imnd.betaFresnel = opp->lookup('betaFresnel')
               end
       'i'    : begin
                imnd.BWSinc = opp->lookup('BWsinc')
                imnd.betaSinc = opp->lookup('betasinc')
                if opp->lookup('etasinc') ne -1 then imnd.etaSinc = opp->lookup('etasinc') else imnd.etaSinc = 0.5
                if opp->lookup('scale_sinc') le 0 then imnd.scaleSinc = Si(!DPI) + Si(2*!DPI*imnd.etaSinc*imnd.tgrtime*imnd.BWSinc/float(1000))   $
                  else  imnd.scaleSinc = opp->lookup('scale_sinc')
                end
      else    : break          
      endcase
    endif
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; DIFFUSION CHECK                                                               ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    bval = opp->lookup('bvalue', is_array=ia, not_found=ntfd)
    if (ntfd eq 0 and ia ne 0) then begin
        imnd.bval_array = ptr_new(bval)
        imnd.n_bvals = n_elements(bval)
    endif
    
    if (stregex(opp->lookup('array'), "(d(ro|sl|pe),?){3}", /boolean)) then begin
        ;; note the flipping of pe and sl (y & z in MAS) vectors.
        dro_dpe_dsl_mult = [1, -1, -1]
        
        rot = diag_matrix(dro_dpe_dsl_mult)

        dro = opp->lookup('dro') * dro_dpe_dsl_mult[0]
        dpe = opp->lookup('dpe') * dro_dpe_dsl_mult[1]
        dsl = opp->lookup('dsl') * dro_dpe_dsl_mult[2]
        bvl = opp->lookup('bvalue')
        
        if (opp->paramExists('bvalrr')) then begin
            bmatrix = [ [ opp->lookup('bvalrr') ], $
                        [ opp->lookup('bvalpp') ], $
                        [ opp->lookup('bvalss') ], $
                        [ opp->lookup('bvalrp') ], $
                        [ opp->lookup('bvalrs') ], $
                        [ opp->lookup('bvalsp') ] ]
            bmatrix = transpose(bmatrix)
            num_b = (size(bmatrix, /dimensions))[1]
            ;; This corrects the b-matrix to flip y & z eigenvector component
            ;; so that the diffusion directions will be correct. This correction
            ;; was found to be necessary after some experimentation.
            for i = 0, num_b-1 do begin
                btmp = [ [ bmatrix[0,i], bmatrix[3,i], bmatrix[4,i] ], $
                         [ bmatrix[3,i], bmatrix[1,i], bmatrix[5,i] ], $
                         [ bmatrix[4,i], bmatrix[5,i], bmatrix[2,i] ] ]
                btmp_c = rot ## btmp ## rot
                ; remember: multiplying off-diagonal elements by 2 is a mas convention.
                bmatrix[*,i] = [   btmp_c[0,0],   btmp_c[1,1],   btmp_c[2,2], $
                                 2*btmp_c[1,0], 2*btmp_c[2,0], 2*btmp_c[2,1] ]
            endfor
        endif else begin
            grad = transpose([[dro], [dpe], [dsl]])
            bmatrix = fltarr(6, n_elements(grad)/3)
            for i = 0, n_elements(grad)/3-1 do begin
                bmatrix[0,i] = bvl[i] * grad[0,i] * grad[0,i] ;; xx
                bmatrix[1,i] = bvl[i] * grad[1,i] * grad[1,i] ;; yy
                bmatrix[2,i] = bvl[i] * grad[2,i] * grad[2,i] ;; zz
                bmatrix[3,i] = bvl[i] * grad[0,i] * grad[1,i] * 2.0 ;; xy
                bmatrix[4,i] = bvl[i] * grad[0,i] * grad[2,i] * 2.0 ;; xz
                bmatrix[5,i] = bvl[i] * grad[1,i] * grad[2,i] * 2.0 ;; yz
            endfor
        endelse
        
        ;; handle interleaved reference scans
        image_flags = opp->lookup('image', not_found=nf, is_array=ia)
        if (nf) then begin
            keep = lindgen(n_elements(dro))
        endif else if (ia) then begin
            ;; 'image' = 1 => actual scan; otherwise => some kind of reference scan
            keep = long(where(image_flags eq 1))
        endif else begin
            keep = lindgen(n_elements(dro))
        endelse
        
        imnd.b_matrix = ptr_new(bmatrix[*,keep])
        imnd.bval_array = ptr_new(bvl[keep])
        imnd.n_bvals = n_elements(bvl[keep])
        
        dro = dro[keep]
        dpe = dpe[keep]
        dsl = dsl[keep]
        angle_theta = fltarr(n_elements(dro))
        angle_phi   = fltarr(n_elements(dro))
        
        for v=0, n_elements(dro)-1 do begin
            vec = [ dro[v], dpe[v], dsl[v] ]
            vec_sph = cv_coord(from_rect=vec, /to_sphere, /degrees)
            angle_theta[v] = 90.0 - vec_sph[1]
            angle_phi[v]   = vec_sph[0]
        endfor
        
        imnd.angle_theta = ptr_new(angle_theta)
        imnd.angle_phi   = ptr_new(angle_phi)
        
        imnd.big_delta = ptr_new(opp->lookup('tDELTA') * 1000.0)
        imnd.small_delta = ptr_new(opp->lookup('tdelta') * 1000.0)
        
    endif

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; Dynamics Check                                                                ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    nt = opp->lookup('nt', is_array=ia, not_found=nf) ;nt = number of transients (averages)
    if (ia ne 0) then begin
        na = nt[0]
        tat = (imnd.recov_time)/1000.*imnd.pdim*na ;total acquisition time
        num_pre = imnd.n_pre
        s_nt = size(nt)
        num_scans = s_nt[1]
        scanLengthSeconds = dblarr(num_scans) + 1.0
        scanLengthSeconds *= tat
        time_array = make_array(num_scans, /double, /index)
        time_array *= tat 
        time_array+= tat/2.0
        ;time array formation from the Bruker days (see mas_extract_imnd/extract_Dynamics_time)
        imnd.time.array = ptr_new(time_array)
        imnd.time.length = ptr_new(scanLengthSeconds)
        imnd.time.zero = time_array[num_pre] - scanLengthSeconds[num_pre]/2.0
        imnd.time.min = time_array[0]- scanLengthSeconds[0]/2.0
        imnd.time.max = time_array[num_scans-1]+scanLengthSeconds[num_scans-1]/2
    endif

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; MULTISLICE CHECK                                                              ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    
    imnd.slices  = num_slices
    imnd.thick   = opp->lookup('thk')
    if (num_dims eq 2) then begin
        gap = float(10*opp->lookup('gap'))
        imnd.s_fov = (imnd.thick * imnd.slices + ((imnd.slices-1) * gap))/10.0
        imnd.s_voxsz = imnd.thick
        pss = opp->lookup('pss', not_found=not_found)
        if (not_found eq 0) then begin
            imnd.slice_scheme_ptr = ptr_new(sort(pss))
        endif else begin
            imnd.slice_scheme_ptr = ptr_new(indgen(num_slices))
        endelse
    endif else begin
        gap = 0
        imnd.s_voxsz = imnd.s_fov/float(imnd.sdim)
    endelse

    imnd.f_voxsz = imnd.f_fov/float(imnd.fdim)
    imnd.p_voxsz = imnd.p_fov/float(imnd.pdim)
       
    imnd.dimensions = num_dims
    
    if fatsup eq 'y' then imnd.fat_sup = 1 else imnd.fat_sup = 0
    
    imnd.file_Path    = opp->lookup('procpar_file')
    imnd.display_Name = opp->lookup('parent_dir')
    
    acq_matrix = ptr_new(diag_matrix(replicate(1.0D,4)))
    imnd.acq_matrix = acq_matrix
    
    return, imnd

end

function mas_varian_apply_xfm, pdata, axes

    sz_data = size(*pdata, /dimensions)
    
    i3 = diag_matrix(replicate(1.0, 3))

    ;for i=0, 2 do begin
    ;    if (axes[i] lt 0) then begin
    ;        i3[i,*] *= -1
    ;    endif
    ;endfor
    
    xfm = i3[abs(axes), *]
    
    sz_new = long(abs(sz_data#xfm))
    
    ixfm = invert(xfm)
    
    data_new = fltarr(sz_new)
    
    for x=0, sz_new[0]-1 do begin
        for y=0, sz_new[1]-1 do begin    
            for z=0, sz_new[2]-1 do begin
                xyz = [x,y,z]#ixfm
                data_new[x,y,z] = (*pdata)[xyz[0], xyz[1], xyz[2]]
            endfor
        endfor
    endfor
    
    return, data_new

end

;+
; :Description:
;    Opens up an Agilent FID directory and loads the
;    scan parameters into mas.
;
; :Params:
;    dirname: optional string path to a FID directory.
;
; :Author: wtriplett
; Modifications (Magdoom, 04/06/15)
; - Added capability to read 1D spectroscopy data other than singlepulse

pro mas_varian_open_fid_dir, dirname

    common scan_data, project
    forward_function get_dir_slash
    
    ni = project.ni
    
    sl = (get_dir_slash())[1]
    
    if (n_elements(dirname) eq 0 || dirname eq '') then begin
        dirname = dialog_pickfile(/directory)
        if (dirname eq '') then return
    endif
    
    oprocpar = obj_new('mas_varian_procpar', dirname+sl+'procpar')
    
    fid_file = dirname+"/fid"
    
    if (not file_test(fid_file, /read)) then begin
        void = dialog_message("Could not find FID file in directory.", $
                              /error, /center)
        obj_destroy, oprocpar
        return
    endif
    
;    if (oprocpar->paramExists('nseg', value=nseg)) then begin
;        if (nseg gt 1) then begin
;            void = dialog_message(['Multishot acquisitions are not supported', $ 
;                                   'Please open the Agilent .img folder instead.'], $
;                                   /error, /center)
;            obj_destroy, oprocpar
;            return
;        endif
;    endif
    str_1d = ['im1D*','prescan_ex','imt1puls','imt2puls','imcpmg','adcpuls']
    flag_1d = 0
    for i = 0,n_elements(str_1d)-1 do if strmatch(oprocpar->lookup('apptype'),str_1d[i]) then flag_1d = 1 
    if (flag_1d eq 1) then begin
        imnd = mas_varian_procpar_to_imnd_1d(oprocpar)
    endif else begin
        imnd = mas_varian_procpar_to_imnd(oprocpar)
    endelse
    
    obj_destroy, oprocpar
    
    imnd.state1_load_procedure = 'mas_varian_load_state_1_fid'
    
    project.imndarray[ni] = imnd
    project.current_path = dirname
    project.scan_Open_Flag = 1
    project.ci = project.ni
    project.ni = project.ni + 1
    
    mas_add_current_scan_to_scan_selection
    
    project.procPramArray[ni].signal_type = 0    
    cd, dirname
    
end

;+
; :Description:
;    Opens up an Agilent FDF reconstruction directory and
;    loads the scan parameters into MAS imnd array.
;
; :Params:
;    dirname: optional, a string path to an existing directory
;
;
;
; :Author: wtriplett
;-
pro mas_varian_open_fdf_dir, dirname

    common scan_data, project
    forward_function get_dir_slash
    
    ni = project.ni
    
    sl = (get_dir_slash())[1]
    
    if (n_elements(dirname) eq 0 || dirname eq '') then begin
        dirname = dialog_pickfile(/directory)
        if (dirname eq '') then return
    endif
    
    op = obj_new('mas_varian_procpar', dirname+sl+"procpar")
    
    fdf_files = file_search(dirname+'*.fdf', /fold_case)
    
    if (strlen(fdf_files[0]) eq 0) then begin
        void = dialog_message("Could not find any FDF files in directory.", $
                              /error, /center)
        return
    endif
    
    imnd = mas_varian_procpar_to_imnd(op)
    
    ;; Adjust for oversampling
    if (op->paramExists('oversample', value=oversample)) then begin
        imnd.fdim = long(imnd.fdim/float(oversample))
        imnd.f_voxsz *= oversample
        imnd.k_fdim_span_max = imnd.fdim
    endif
    
    ;; Adjust for reference scans 
    if (op->paramExists('nrefs', value=nrefs)) then begin
        imnd.adim -= nrefs
    endif
    
    obj_destroy, op
    imnd.multi_scan_file_array = ptr_new(fdf_files)
    project.procPramArray[ni].signal_type = 0

    imnd.state1_load_procedure = 'mas_varian_load_state_1_fdf'
    
    project.imndarray[ni] = imnd
    project.scan_Open_Flag = 1
    project.current_path = dirname
    project.ci = project.ni
    project.ni = project.ni + 1
    
    mas_add_current_scan_to_scan_selection
    
end
   
;+
; :Description:
;    state 1 loader for Agilent fid directory
;
;
;
;
;
; :Author: wtriplett
;-
pro mas_varian_load_state_1_fid

    common scan_data, project
    COMMON common_widgets
    forward_function get_dir_slash
    
    ci = project.ci
    
    if (project.procpramarray[ci].state_1 ne 0) then begin
        return
    endif
    
    project.imndarray[ci].image_type = 1
    
    dir_sep = (get_dir_slash())[1]
    
    procpar_path = project.imndarray[ci].file_path
    opp = obj_new('mas_varian_procpar', procpar_path)
    
    seqfil = opp->lookup('seqfil')
    apptype = opp->lookup('apptype')
    
    case project.imndarray[ci].dimensions of 
    
        2: begin

            openr, lun, file_dirname(project.imndarray[ci].file_path) + dir_sep+"fid",$
                   /swap_if_little_endian, /get_lun, error=iserror
            if (iserror) then begin
                void = dialog_message(['Unable to open FID file.', $
                                       !ERROR_STATE.MSG], /error, /center)
                goto, CLEANUP
            endif
            
            if (apptype eq 'imEPI' or apptype eq 'im2Depi') then begin
                if (apptype eq 'imEPI') then   data = mas_varian_read_fid_epi(lun, opp)
                if (apptype eq 'im2Depi') then data = mas_varian_read_fid_epip(lun, opp)
                
                if (opp->paramExists('epiref_type', value=epiref_type)) then begin
                    image_indicator = opp->lookup('image', not_found=ntfd)
                    project.imndarray[project.ci].adim = n_elements(image_indicator)
                    YesNo = dialog_message(['Discard EPI reference scans?', $
                                            '(Correction will still be applied.)'], /question, /center)
                    if (YesNo eq 'Yes') then begin
                        ;; This will remove the reference scans from the final recon array.
                        ;; disabled for now. 
                        actual_images = where(image_indicator eq 1, num_images)
                        if (num_images ne 0) then begin
                            *data = (*data)[*,*,*,actual_images]
                            project.imndarray[project.ci].adim = num_images
                            project.imndarray[project.ci].epi_num_ref_scans = 0
                        endif
                    endif
                endif
            endif else begin
                data = mas_varian_read_fid_data_2d(lun, opp)
            endelse
            close, lun
            project.dataarray[ci].state1 = data

        end
        
        3: begin
        
            openr, lun, file_dirname(project.imndarray[ci].file_path) + dir_sep+"fid",$
                   /swap_if_little_endian, /get_lun, error=iserror
            if (iserror) then begin
                void = dialog_message(['Unable to open FID file.', $
                                       !ERROR_STATE.MSG], /error, /center)
                goto, CLEANUP
            endif
            data = mas_varian_read_fid_data_3d(lun, opp)
            close, lun
            project.dataarray[ci].state1 = data
            
        end
        
        else: begin
            
            ;;void = dialog_message('Varian support only includes 2D or 3D acquisitions.', $
            ;;                      /center, /error)
            openr, lun, file_dirname(project.imndarray[ci].file_path) + dir_sep+"fid",$
                   /swap_if_little_endian, /get_lun, error=iserror
            
            data = mas_varian_read_fid_data_1d(lun, opp)
            close, lun
            project.dataarray[ci].state1 = data
            project.procpramarray[CI].state_1 = 1
            project.imndarray[ci].image_type = 99
            goto, CLEANUP
          
        end
        
    endcase
    
    if project.procPramArray[ci].zpad_flag eq 1 then mas_zeropad, project.dataarray[ci].state1          
        
    mas_kspace_subsample, project.dataarray[ci].state1
    mas_kspace_shift, project.dataArray[CI].state1

    IF (project.procPramArray[ci].filterflag eq 'Frequency') THEN frequency_filter
    
    signal_type = project.procPramArray[CI].signal_type
        
    IF signal_type LT 5 OR signal_type GT 7 THEN BEGIN
        ;; data that needs FFT, if motion correct, apply
        ;; now
        mas_fft,   project.dataArray[CI].state1
        mas_shift, project.dataArray[CI].state1
    ENDIF
        
    mas_signal_selection , project.dataArray[CI].state1
    
    ;; slice-direction part for 2D multislice "kspace" span reduction.
    mas_kspace_subsample_multislice, project.dataarray[CI].state1

    ;; get rid of oversampled part
    oversample = opp->lookup('oversample', not_found=ntfd)
    if (1 and ntfd eq 0) then begin
        data = project.dataArray[CI].state1
        old_fdim = (size(*data, /dimensions))[0]
        new_fdim = old_fdim / oversample
        *data = (*data)[(old_fdim-new_fdim)/2:(old_fdim+new_fdim)/2-1,*,*,*]
        if project.procPramArray[ci].zpad_flag eq 1 then project.imndarray[ci].fdim = new_fdim/2 else project.imndarray[ci].fdim = new_fdim
        project.imndarray[ci].f_voxsz *= oversample
    endif

    IF project.procPramArray[CI].mc_enable NE 0 THEN BEGIN
        mas_motion_correct
    ENDIF
    
    mas_smoothing , project.dataArray[CI].state1
        
    IF (project.procPramArray[ci].image_domain_filterflag eq 'Image Domain') THEN BEGIN
        image_domain_filter
    ENDIF
    
    if (project.procpramarray[ci].freq_interp ne 1 or $
        project.procpramarray[ci].phase_interp ne 1 or $
        project.procpramarray[ci].slice_interp ne 1) then begin
        
        temp = mas_interpolate(project.dataarray[CI].state1)
        ptr_free, project.dataarray[CI].state1
        project.dataarray[CI].state1 = temp
    endif
    
    project.procpramarray[CI].state_1 = 1

    CLEANUP:
    
    obj_destroy, opp
    
end

;+
; :Description:
;    state 1 loader for Agilent FDF reconstruction directories.
;
;
; :Author: wtriplett
;-
pro mas_varian_load_state_1_fdf

    common scan_data, project
    
    ci = project.ci
    
    fdf_files = *project.imndarray[ci].multi_scan_file_array
    scan_name = project.imndarray[ci].scan_name
    scan_indicator = strmid(scan_name, 0, 4)
    
    fdim = project.imndarray[ci].fdim
    pdim = project.imndarray[ci].pdim
    sdim = project.imndarray[ci].sdim
    adim = project.imndarray[ci].adim
    
    num_echoes = project.imndarray[ci].n_echo
    
    data = fltarr(fdim, pdim, sdim, adim)
    nvols = 0
    
    progressbar = obj_new('progressbar', title='Loading Files...', text='Loading .fdf files...', /nocancel)
    progressbar->Start
    
    ;; This is required so that imagery read from
    ;; .fdf data will have the proper read/freq, phase axes
    ;; in mas (read -> image-x, phase -> image-y).
    rotate_amount = 3
      
    for f = 0, n_elements(fdf_files)-1 do begin
    
        tmpfdf = mas_varian_read_fdf_file(fdf_files[f])
        if (f eq 0) then begin
            ;; At this time we will check for any zero-filling that has
            ;; happened in recon and adjust the matrix sizes and voxel
            ;; dims accordingly.
            img_size = size(*tmpfdf.img_data, /dimension)
            if (fdim ne img_size[0] or pdim ne img_size[1]) then begin
                data = fltarr(img_size[0], img_size[1], sdim, adim)
                project.imndarray[ci].fdim = img_size[0]
                project.imndarray[ci].pdim = img_size[1]
                project.imndarray[ci].f_voxsz = project.imndarray[ci].f_fov/img_size[0]
                project.imndarray[ci].p_voxsz = project.imndarray[ci].p_fov/img_size[1]
            endif
        endif
        ;; must be a slice or slab
        if (strpos(strupcase(fdf_files[f]), 'SLICE') ne -1 and $
            strpos(strupcase(fdf_files[f]), 'SLAB') ne -1) then continue
            
        if (strpos(strupcase(fdf_files[f]), 'SLICE') ne -1) then begin
            
            ;; if its a slice then it will have these paramets
            slice = tmpfdf.slice_no-1
            ;;if (scan_indicator eq 'mems' or scan_indicator eq 'tagc') then begin
            if (num_echoes gt 1) then begin
                aidx  = tmpfdf.echo_no-1
            endif else begin
                aidx  = tmpfdf.array_index-1
            endelse
            data[*,*,slice,aidx] = rotate(*tmpfdf.img_data, rotate_amount)
                
        endif else begin
        
            data[*,*,*,nvols++] = rotate(*tmpfdf.img_data, rotate_amount)
            if (nvols eq adim) then break
            
        endelse
        
        ptr_free, tmpfdf.img_data
        
        progressbar->Update, float(f)/n_elements(fdf_files) * 100.0
        
    endfor
    progressbar->Destroy
    
    project.dataarray[ci].state1 = ptr_new(data, /no_copy)
    
    if (project.procpramarray[ci].freq_interp ne 1 or $
        project.procpramarray[ci].phase_interp ne 1 or $
        project.procpramarray[ci].slice_interp ne 1) then begin
        
        temp = mas_interpolate(project.dataarray[CI].state1)
        ptr_free, project.dataarray[CI].state1
        project.dataarray[CI].state1 = temp
    endif

    project.procpramarray[ci].state_1 = 1
    
end

;+
; :Description:
;    Event handles for Agilent "open experiment directory" GUI.
;
; :Params:
;    event: IDL generated event
;
;
;
; :Author: wtriplett
;-
pro mas_varian_summarize_experiment_dir_event, event

    widget_control, event.top, get_uvalue=state
    ui = (*state).ui
    
    case event.id of
    
        ui.table: begin
            if ((event.type eq 9 or event.type eq 4) && event.sel_top ge 0) then begin
        
                ;; extract just the rows that are selected, regardless
                ;; of whatever columns are selected. selection will be
                ;; expanded to cover all columns in each row.
                select = widget_info(event.id, /table_select)
                
                n_cols = n_elements(widget_info(event.id, /column_widths))
                select = reform(select[1,*])
                select = select[uniq(select[sort(select)])]
                ;; the trick is to adjust the selection and then to keep the view the same.
                table_view = widget_info(event.id, /table_view)
                
                if (event.type eq 9) then begin
                    exclude = lindgen((event.sel_bottom - event.sel_top + 1)) + event.sel_bottom
                    ;;print, exclude 
                endif
                
                for c = 0, n_elements(select)-1 do begin
                
                    if (event.type eq 9 && select[c] eq exclude) then continue
                    
                    ;; fill the selection out so that entire rows appear to be selected.
                    if ( n_elements(new_select) eq 0 ) then begin
                        new_select = transpose([[lindgen(n_cols)], [replicate(select[c],n_cols)]])
                    endif else begin
                        new_select = [ [new_select], [transpose([[lindgen(n_cols)], [replicate(select[c],n_cols)]])] ]
                    endelse
                endfor
        
                widget_control, event.id, set_table_select=new_select
                widget_control, event.id, set_table_view=table_view
            
            endif
         end
         
         ui.btn_open_selected: begin
             select = widget_info(ui.table, /table_select)
             select = reform(select[1,*])
             select = select[uniq(reform(select))]
             for s = 0, n_elements(select)-1 do begin
                mas_varian_open_fid_dir, (*state).data[select[s]].(0)
             end
         end
         
;         ui.btn_export_csv: begin
;         
;         end
         
         ui.btn_close: begin
            widget_control, event.top, /destroy
         end
         
     endcase
     
    
end

;+
; :Description:
;    Reads all of the FID directories in a study folder and 
;    presents the user with a tabular GUI for selecting which
;    ones to open.
;
; :Params:
;    dir: optional string path of an experiment/study folder
;
; :Author: wtriplett
;-
pro mas_varian_summarize_experiment_dir, dir

    common scan_data, project
    common common_widgets
    
    cd, current=cwd
    
    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        cd, cwd
        void = dialog_message(['An error occurred while browsing varian directory.', $
                               !ERROR_STATE.MSG], /error, /center)
        if (obj_valid(pbar)) then pbar->Destroy
        if (n_elements(table) ne 0 && widget_info(table, /valid_id)) then begin
            widget_control, table, /destroy
        endif
        if (ptr_valid(state)) then ptr_free, state
        return
    endif
    catch, /cancel
    
    if (n_elements(dir) eq 0) then begin
        dir = dialog_pickfile(path=project.current_path, title='Please select a Varian experiment directory.', $
                              /dir, /must_exist)
        if (dir eq '') then return
    endif
    
    cd, dir
    
    pbar = obj_new('progressbar', title='Please Wait', text='Reading directory...')
    pbar->Start
    
    files = file_search('.', 'procpar', /fully_qualify_path)
    
    if (n_elements(files) eq 0) then begin
        void = dialog_message('No valid varian .fid experiments found.', /center, /error)
        return
    endif
    
    imnd_names = strupcase(['file_path','display_name','scan_name','scan_date','fdim','pdim','sdim','f_fov','p_fov','s_fov', $
                            'f_voxsz','p_voxsz','s_voxsz','recov_time','echo_time','n_avg','adim','epi_num_ref_scans'])
                    
    head_names = ['Full Path','Scan Dir','Scan Name','Scan Date', 'x dim', 'y dim', 'z dim', 'x fov (cm)', 'y fov (cm)', 'z fov (cm)', $
                  'x vox dim (cm)', 'y vox dim (cm)', 'z vox dim (cm)', 'TR (ms)', 'TE (ms)', 'Num Avg', 'Num Reps', $
                  'Num Ref Scans']
    col_widths = [70,85,90,70,45,45,45,72,72,72,80,80,80,60,60,60,60,75]              
    header_str = '"'+strjoin(head_names, '","')+'"'
    
    tags = tag_names(project.imndarray[0])
    
;    openw, lun, dialog_pickfile(default_extension='csv'), /get_lun
    
;    printf, lun, header_str
    
    if (!VERSION.OS_FAMILY eq 'Windows') then begin
        ;; double path sep to escape backslash (\\)
        ps = path_sep()+path_sep()
    endif else begin
        ps = path_sep()
    endelse
    
    for f = 0, n_elements(files)-1 do begin
    
        if (not strmatch(files[f], '*.fid'+ps+'procpar')) then continue
        
        opp = obj_new('mas_varian_procpar', files[f])

        imnd = mas_varian_procpar_to_imnd(opp)
        
        if (size(imnd, /type) ne 8) then begin
            message, /info, files[f]+' possibly corrupted. SKipping.'
            continue
        endif
        
        str = ''
        
        for i = 0, n_elements(imnd_names)-1 do begin
            
            val = imnd.(where(strmatch(tags, imnd_names[i])))
            
            if (n_elements(tmpstruct) eq 0) then begin
                tmpstruct = create_struct(imnd_names[i], val)
            endif else begin
                tmpstruct = create_struct(tmpstruct, imnd_names[i], val)
            endelse
            
            if (1 or i lt 2) then begin
                str += '"'+strtrim(val,2)+'"'
            endif else begin
                str += strtrim(val,2)
            endelse
            str += ','
            
        endfor
        tmpstruct.(0) = file_dirname(tmpstruct.(0))
        tmpstruct.(1) = file_basename(tmpstruct.(1))
        
        if (n_elements(data) eq 0) then begin
            data = temporary(tmpstruct)
        endif else begin
            data = [ data, temporary(tmpstruct) ]
        endelse
        
        ;printf, lun, str
        
        obj_destroy, opp

        pbar->Update, float(f)/(n_elements(files)-1) * 100.0
        
    endfor

    pbar->Destroy
    
    if (n_elements(data) eq 0) then begin
        void = dialog_message('No eligible experiment directories found.', /error, /center)
        return
    endif
    
    base = widget_base(title='Varian Experiment Directory', /column, group_leader=WID_BASE_MAIN)
    table = widget_table(base, value=data, column_widths=col_widths, column_labels=head_names, $
                         /all_events, /resizeable_columns, /disjoint_selection, /scroll, $
                         scr_xsize=900, scr_ysize=500, /no_row_headers)
            
    b = widget_base(base, /row, /base_align_center, /align_center)
    btn_open_selected = widget_button(b, value='Open Selected Scans')
    ;;btn_export_csv  = widget_button(b, value='Export Table as .csv')
    btn_close = widget_button(b, value='Close Window')
    
    state = ptr_new({ ui: { tlb: base, $
                            table: table, $
                            btn_open_selected: btn_open_selected, $
                            $;btn_export_csv: btn_export_csv, $
                            btn_close: btn_close }, $
                      data: data $
                    })
    
    widget_control, base, /realize, set_uvalue=state
    
    xmanager, 'mas_varian_summarize_experiment_dir', base
    
    ptr_free, state
    
   ; close, lun & free_lun, lun
    
end
