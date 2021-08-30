;; $Id$
;;
;; Subroutine name: mc_compute_mi
;;
;; Created by: BT 2008-03-10
;;
;; Calling Information:
;;   im_a, im_b - images to compare
;;           bs - bin size to use in histogramming
;;
;; Bugs or Important Comments to Developers:
;;
;;
;; Purpose of subroutine:
;;  computes mutual information of two images
;;
;; Editing Information:
;;

function mc_compute_mi, im_a, im_b, bs, normalized=normalized

    ;max_A = max(im_a)
    ;max_B = max(im_b)
    
    ;A = im_a/max([max_A,max_B])
    ;B = im_b/max([max_A,max_B])

    A = bytscl(im_A)
    B = bytscl(im_B)

    sz_a = size(A)
    sz_b = size(B)

    pab = float(hist_2d(A, B, bin1=bs, bin2=bs));, max1=250, max2=250, min1=5, min2=5))

    hist_a = float(histogram(A, binsize=bs));, max=250, min=5))
    hist_a = hist_a/total(hist_a)

    hist_b = float(histogram(B, binsize=bs));, max=250, min=5))
    hist_b = hist_b/total(hist_b)        
    
    Ha = 0.0
    Hb = 1.0
    if (0 and keyword_set(normalized)) then begin
        normalized = 1
        if (normalized eq 1) then begin
            n_save = where(hist_a ne 0)
            Ha = total((-hist_a(n_save)*alog(hist_a[n_save])))/alog(2)
            n_save = where(hist_b ne 0)
            Hb = total((-hist_b(n_save)*alog(hist_b[n_save])))/alog(2)
        endif
    endif

    pab = pab / total(pab)
    pa_pb = matrix_multiply(hist_a, hist_b, /btranspose);;  transpose(hist_b)
    save = where(pa_pb gt 1e-12)
    u = pab(save)
    v = pa_pb(save)
    t = where(u gt 1e-12)
    y = u[t] * (alog(u[t]/v[t])/alog(2))

    return, total(y)/(Ha+Hb)

end

;; Subroutine name: mc_objective_function
;;
;; Created by: BT 2008-03-14
;;
;; Calling Information:
;;  steps - a 2x1 vector of shift coordinates
;;
;; Bugs or Important Comments to Developers:
;;  The usefulness of this hasn't been established
;;
;; Purpose of subroutine:
;;   returns the _negative_ of the mutual information, for use with
;;   IDL's built-in optimization (minimization) procedures (like
;;   POWELL)
;;
;; Editing Information:
;;

function mc_objective_function, steps

    common _mc_working_set, working_set, _mc_working_set

  ;;  print, steps

    txd = mc_transform_image(working_set.tar_image, steps)

    m = mc_compute_mi(float(working_set.ref_image), $
                      float(txd),$
                      working_set.bs,$
                      normalized=working_set.normalized)

    return, -m
end

;; Subroutine name: mas_motion_correct_mi
;;
;; Created by: BT 2008-03-10
;;
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine: motion correct a series of images
;;
;; Editing Information:
;;  Created 2008-03-05
;;
pro mas_motion_correct_mi, ref_idx, range_low, range_hi,  $
                           start_mesh, n_passes, binsize, $
                           cor_slice, TOL=tol
    ;; common block for scan data
    common scan_data

    ;; anonymous structure used to pass image data of the current
    ;; working set back and forth between this procedure and the
    ;; objective function. Used by POWELL's method of optimization
    common _mc_working_set, working_set, _mc_working_set

    if not keyword_set(tol) then begin
        tol = 0.0001
    endif

    print, "Starting motion correction with tolerance = "+string(tol, FORMAT='(D0.10)')

    CI = project.ci

    fdim = project.imndArray[CI].fdim
    pdim = project.imndArray[CI].pdim
    sdim = project.imndArray[CI].sdim
    adim = project.imndArray[CI].adim

    ;; Set of the images in intensity mode
    image_set = ptr_new(dblarr(fdim,pdim,adim))

    data_file = project.imndarray[ci].file_path+'.motion_correction.dat'
    save_ok = 0

    if (project.imndArray[ci].image_type eq 11) then begin
        ;; PAR/REC scans don't have parent dirs, per se. check to make
        ;; sure that the dir they are in is writeable
        fs_test = file_dirname(project.imndarray[ci].file_path)
    endif else begin
        ;; scans that have their own special scan dir structure
        fs_test = project.imndarray[ci].file_path
    endelse

    if (file_test(fs_test, /WRITE)) then begin 
    
        ;; Save shift values here
        openw, lun, data_file, /GET_LUN, BUFSIZE=0, ERROR=err

        if (err eq 0) then begin
            
            settings     = fltarr(16)
            settings[0]  = range_low
            settings[1]  = range_hi
            settings[2]  = cor_slice
            settings[3]  = ref_idx
            settings[4]  = binsize
            settings[5]  = start_mesh
            settings[6]  = n_passes
            settings[7]  = project.procPramArray[ci].mc_mi_min_method
            settings[8]  = tol
            settings[9]  = project.procPramArray[ci].mc_mi_use_translation
            settings[10] = project.procPramArray[ci].mc_mi_use_rotation
            settings[11] = project.procPramArray[ci].mc_mi_use_shear
            settings[12] = project.procPramArray[ci].mc_mi_use_dilation
            settings[13] = project.procPramArray[ci].mc_mi_use_dilation_h
            settings[14] = project.procPramArray[ci].mc_mi_use_dilation_v
            settings[15] = project.procPramArray[ci].mc_mi_normalized
            printf, lun, '0.4'
            printf, lun, settings
            
            save_ok = 1

        endif
    endif
    
    if save_ok ne 1 then begin
        
        junk = dialog_message('Unable to write to scan directory. '+$
                              'MC data will not be saved. Continue?', /QUESTION, /center)
        if junk eq 'No' then return
        
    endif

    ;; For the progress bar
    if (cor_slice eq -1) then begin
        max_iter = 0 ;(sdim+fdim+pdim)*(adim)
        if (project.procpramarray[ci].mc_mi_correct_dir_p eq 1) then max_iter += pdim
        if (project.procpramarray[ci].mc_mi_correct_dir_r eq 1) then max_iter += fdim
        if (project.procpramarray[ci].mc_mi_correct_dir_s eq 1) then max_iter += sdim
        max_iter *= adim
    endif else begin
        max_iter = adim
    endelse
    cur_iter = 0

    progressbar = OBJ_NEW('progressbar', Color='red', $
                          Title='Performing MI Motion Correction . . .')
    progressbar -> Start

    ;; grab the current time to measure elapsed
    st_time = systime(1)
    

    for dimension = 0, 2 do begin

        if (dimension eq 2 and project.procpramarray[ci].mc_mi_correct_dir_p eq 0) then begin
            print, "Skipping phase direction"
            ;;cur_iter += pdim*adim
            continue
        endif
        if (dimension eq 1 and project.procpramarray[ci].mc_mi_correct_dir_r eq 0) then begin
            print, "Skipping freq. direction"
            ;;cur_iter += fdim*adim
            continue
        endif
        if (dimension eq 0 and project.procpramarray[ci].mc_mi_correct_dir_s eq 0) then begin
            ;;cur_iter += sdim*adim
            print, "Skipping slice direction"
            continue
	endif

        case dimension of
            
            0: begin
                sl_dir = 's'
                dim = sdim
                dim_start = 0
                dim_end = sdim-1
                image_set = ptr_new(dblarr(fdim,pdim,adim))
            end
            1: begin
                sl_dir = 'r'
                dim = fdim
                dim_start = 0
                dim_end = fdim-1
                image_set = ptr_new(dblarr(pdim,sdim,adim))
            end
            2: begin
                sl_dir = 'p'
                dim = pdim
                dim_start =  0
                dim_end = pdim-1
                image_set = ptr_new(dblarr(fdim,sdim,adim))
            end

        endcase

        ;; eliminate no-content slices
        mn = fltarr(dim_end - dim_start + 1)
        ;;for i=dim_start, dim_end do begin
        ;;    case dimension of
        ;;        0: slc = (*project.dataarray[ci].state1)[*,*,i,ref_idx]
        ;;        1: slc = (*project.dataarray[ci].state1)[i,*,*,ref_idx]
        ;;        2: slc = (*project.dataarray[ci].state1)[*,i,*,ref_idx]
        ;;    endcase
        ;;    mn[i] = mean(slc)
        ;;endfor
        
        ;; determine which slices are to be skipped
        ;std = stdev(mn)
        ;skip_slices = where(mn - .25*std le 0)
        
        ;; Now consider each slice 
        for slice = dim_start, dim_end do begin 
            
            progressbar->setProperty, TITLE='MI Motion Corr. ('+sl_dir+')'
            
            ;; True if user selects only one slice to correct
            if slice ne cor_slice and cor_slice ne -1 then begin
                ;; leave this slice alone
                continue
            endif
            
            ;;if (where(skip_slices eq slice) ne -1) then begin
            ;;    print, "Skipping "+sl_dir+string(slice)
            ;;    cur_iter += adim
 ;          ;;     continue
            ;;endif

            ;; *image_set contains the data to be operated on.
            ;; it contains the _intensity_ data, not kspace (FFT)!
            
            case dimension of 

                0: (*image_set)[*,*,*] = REFORM( (*project.dataArray[ci].state1)[*,*,slice,*] )
                1: (*image_set)[*,*,*] = REFORM( (*project.dataArray[ci].state1)[slice,*,*,*] )
                2: (*image_set)[*,*,*] = REFORM( (*project.dataArray[ci].state1)[*,slice,*,*] )
                
            endcase

            ref_image = REFORM((*image_set)[*,*,ref_idx])
            ref_image = mc_pad_image(ref_image, out_x=pd_x, out_y=pd_y)
            
            ;; Ok, now consider each slice's array data
            for aidx = 0, adim-1 do begin
                
                cur_iter++
                progressBar -> SetProperty, $
                  Text=strcompress('Correcting slice: '+ $
                                   string(slice)+'/' + $
                                   string(dim_end-dim_start)+', volume: ' + $
                                   string(aidx+1)+'/'+string(adim))
                progressBar -> Update, (float(cur_iter)/float(max_iter))*100.0
                
                ;; Don't correct the reference image
                if aidx eq ref_idx then  continue
                
                ;; The working magnitude mode image
                tar_image = mc_pad_image(REFORM((*image_set)[*,*,aidx]))
                
                ;; this is a common block to share current working
                ;; image data with the objective function. Pointers
                ;; should reduce memory usage and keep the footprint
                ;; of the common block to a mimimum
                
                working_set = { ref_image: ref_image, $
                                tar_image: tar_image, $
                                bs:        binsize,   $
                                normalized: project.procPramArray[ci].mc_mi_normalized }
                max_mi     = 0
                best_mi    = 0
                scls       = [0,-2,2,4,-4,-6,6]
                start_pt   = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                best_ever  = start_pt
                best_shift = [ start_pt[0], start_pt[1] ]
                
                use_trl = float(project.procpramarray[ci].mc_mi_use_translation)
                use_rot = float(project.procpramarray[ci].mc_mi_use_rotation)
                use_shr = float(project.procpramarray[ci].mc_mi_use_shear)
                use_scl = 0 ;float(project.procpramarray[0].mc_mi_use_dilation)
                
                optvec = [use_trl, use_trl, use_shr, use_shr, use_rot, use_scl, use_scl]

                if (1) then begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                   print, "******** Using Cascade Method ********"
                   dir_vec = fltarr(n_elements(optvec))
                   iterations = total([use_trl, use_shr, use_rot])
                   n_iter = 0
                   
                   for tx = 0, 3-1 do begin ;;iterations-1 do begin
                   
                      if (([use_trl, use_shr, use_rot])[tx] eq 0) then continue
                      
                      dir_vec_wk = dir_vec
                      case tx of
                         0: dir_vec_wk[0:1] = optvec[0:1]
                         1: dir_vec_wk[2:3] = optvec[2:3]
                         2: dir_vec_wk[4]   = optvec[4]
                      endcase
                      
                      dir_vec_wk = diag_matrix(dir_vec_wk)
                      
                      powell, start_pt, dir_vec_wk, tol, max_mi,$
                              'mc_objective_function', ITER=n_iter
                      
                      best_mi = -max_mi
                      best_ever = start_pt
                      
                      print, start_pt
                      
                   endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                endif else begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                   print, "******** Using Iteration Method ********"
                   directions = diag_matrix(float(optvec))
                   
                   scale_first = 0
                   powell, start_pt, directions, tol, max_mi,$
                           'mc_objective_function', ITER=n_iter
                   
                   best_mi = -max_mi
                   best_ever = start_pt                   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                endelse

                scx = 0
                scy = 0
                
                nscls_x = n_elements(scls)
                nscls_y = n_elements(scls)
                
                if (project.procpramarray[ci].mc_mi_use_dilation_h eq 0) then begin
                    nscls_x = 1
                endif
                if (project.procpramarray[ci].mc_mi_use_dilation_v eq 0) then begin
                    nscls_y = 1
                endif
                
                if (project.procpramarray[ci].mc_mi_use_dilation) then begin
                    for scx = 0, nscls_x-1 do begin                
                        for scy = 0, nscls_y-1 do begin                
                            
                            curr_tx = [ best_ever[0:4], scls[scx], scls[scy] ]
                            tmp = mc_transform_image(tar_image, curr_tx)
                            
                            cur_mi = mc_compute_mi(float(working_set.ref_image), $
                                                   float(tmp), $
                                                   working_set.bs, normalized=working_set.normalized)
                            
                            if (cur_mi gt best_mi) then begin
                                best_mi = cur_mi
                                best_ever = curr_tx
                            endif
                            
                        endfor
                    endfor
                endif
                
                best_tx = best_ever
                
                print, strcompress('Found ('+sl_dir+':'+string(slice)+','+string(aidx)+') '+string(max_mi)+$
                                   ' trl'+array2string(best_tx[0:1])+ $
                                   ' shr'+array2string(best_tx[2:3])+ $
                                   ' rot['+string(best_tx[4])+$
                                   '] scl'+array2string(best_tx[5:6])+$
                                   ' in '+string(n_iter)+' iterations.')
                
                ;; User Cancelled correction, clean up and then
                ;; reload state1
                if (progressBar->CheckCancel()) then begin
                    print, "Motion Correction cancelled."
                    ptr_free, image_set
                    working_set = 0
                    if (save_ok) then begin
                        close, lun
                        free_lun, lun
;                        file_delete, data_file
                    endif
                    progressBar->destroy
                    project.procpramarray[project.ci].mc_enable = 0
                    mas_redraw
                    mas_load_state_1
                    return
                endif
                
                if (save_ok) then begin
                    ;; write parameters to save file
                    printf, lun, strcompress(sl_dir+' '+string(slice)+' '+string(aidx)+' '+$
                                             string(best_tx[0])+' '+string(best_tx[1])+' '+$
                                             string(best_tx[2])+' '+string(best_tx[3])+' '+$
                                             string(best_tx[4])+' '+string(best_tx[5])+' '+$
                                             string(best_tx[6]))
                endif
                
                tmp = mc_transform_image(tar_image, best_tx, scale_first=scale_first)
                
                case dimension of 
                    
                    0: (*project.dataArray[ci].state1)[*,*,slice, aidx] = mc_unpad_image(tmp, pd_x, pd_y)
                    1: (*project.dataArray[ci].state1)[slice,*,*, aidx] = mc_unpad_image(tmp, pd_x, pd_y)
                    2: (*project.dataArray[ci].state1)[*,slice,*, aidx] = mc_unpad_image(tmp, pd_x, pd_y)
                    
                endcase

            endfor ;; adim
            
        endfor ;; slice
        
    endfor ;; dimension (Freq/Phase, Phase/Slice)
    
    if (save_ok) then begin
        close, lun
        free_lun, lun
    endif
    
    print, strcompress('Total Elapsed Time:'+string(systime(1)-st_time)+' seconds.')
    progressbar -> Destroy
    ptr_free, image_set
    
    HEAP_GC
    
end

;; Subroutine name: mc_apply_shifts
;;
;; Created by: BT 2008-03-10
;;
;; Calling Information:
;;  shift_values - an array indexed by slice/arrayindex containing
;;                 shift values
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;  Applied bulk shifts
;;
;; Editing Information:
;;  Created 2008-03-05
;;
pro mc_apply_shifts, shift_values, TARGET=target, KSPACE=kspace

    common scan_data

    update_status_bar, 'Applying motion correction...'

    shift_ar_sz = size(shift_values)
    
    sdim = shift_ar_sz[1]
    adim = shift_ar_sz[2]

    ;; the data may come from elsewhere
    if (keyword_set(TARGET)) then begin
        tar = target
    endif else begin
        tar = project.dataArray[project.ci].state1
    endelse

    nofft = (keyword_set(kspace)) ? 1 : 0

    progressbar = OBJ_NEW('progressbar', /fast_loop, $
                          Title='Applying Saved Motion Correction...')
    progressBar -> SetProperty, $
      Text=strcompress('Applying saved motion correction...')
    progressbar -> Start
    
    cur_iter    = 0
    max_iter    = sdim*adim
    scale_first = 0

    for slice = 0, sdim-1 do begin 

        for aidx = 0, adim-1 do begin

            cur_iter = cur_iter + 1
            progressBar -> Update, (float(cur_iter)/float(max_iter))*100.0

            mc_data = shift_values[slice, aidx, *]
            
            ;; no need to go through the rigamorole for zero shift
            if (max(abs(mc_data)) eq 0) then begin ;;(bs_x eq 0.0 and bs_y eq 0.0) then begin
                continue
            endif

            tmp = mc_pad_image( (*tar)[*,*,slice,aidx], out_x=pd_x, out_y=pd_y)
            tmp = mc_transform_image(tmp, mc_data)
            (*tar)[*,*,slice,aidx] = mc_unpad_image( tmp, pd_x, pd_y )

            if (progressBar->checkcancel()) then begin

                print, "Motion Correction cancelled."
                ptr_free, tar

                progressBar->destroy
                project.procpramarray[project.ci].mc_enable = 0
                mas_redraw
                mas_load_state_1
                return

            endif

        endfor

    endfor

    progressBar->Destroy
    update_status_bar, ''

end

;; Subroutine name: mc_load_saved_values
;;
;; Created by: BT 2008-03-10
;;
;; Calling Information:
;;  data_file - file containing saved shift values and settings
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine: 
;;  Load in a data file of shift values
;;
;; Editing Information:
;;  Created 2008-03-05
;;
pro mc_load_saved_values, data_file

    common scan_data

    ci = project.ci
    sdim = project.imndarray[ci].sdim
    adim = project.imndarray[ci].adim
    fdim = project.imndarray[ci].fdim
    pdim = project.imndarray[ci].pdim

    ;; Holds the loaded shift values
    mc_data = fltarr(sdim, adim, 7)
    
    ;; it exists or we wouldn't be here.
    n_lines = file_lines(data_file) ;; number of iterations
    openr, fileid, data_file, /GET_LUN

    ;; read file format version
    readf, fileid, version

    if (version eq 0.1) then begin
        ;; for compatibility with some testing versions of the motion
        ;; correction we won't be able to restore the settings, but
        ;; we can salvage the parameters
        settings = fltarr(7)

        readf, fileid, settings

        project.procPramArray[ci].mc_mi_rng_low    = settings[0] 
        project.procPramArray[ci].mc_mi_rng_hi     = settings[1]
        project.procPramArray[ci].mc_mi_slice      = settings[2]
        project.procPramArray[ci].mc_mi_refidx     = settings[3]
        project.procPramArray[ci].mc_mi_binsize    = alog(settings[4])/alog(2)
        project.procPramArray[ci].mc_mi_start_mesh = settings[5]
        project.procPramArray[ci].mc_mi_n_passes   = settings[6]
        project.procPramArray[ci].mc_mi_use_translation = 1
        project.procPramArray[ci].mc_mi_use_rotation = 0
        project.procPramArray[ci].mc_mi_use_shear = 0
        project.procPramArray[ci].mc_mi_use_dilation = 0
        project.procPramArray[ci].mc_mi_normalized = 0

        n_lines -= 2

    endif else if (version eq 0.2) then begin

        settings = fltarr(9)

        readf, fileid, settings

        project.procPramArray[ci].mc_mi_rng_low    = settings[0] 
        project.procPramArray[ci].mc_mi_rng_hi     = settings[1]
        project.procPramArray[ci].mc_mi_slice      = settings[2]
        project.procPramArray[ci].mc_mi_refidx     = settings[3]
        project.procPramArray[ci].mc_mi_binsize    = alog(settings[4])/alog(2)
        project.procPramArray[ci].mc_mi_start_mesh = settings[5]
        project.procPramArray[ci].mc_mi_n_passes   = settings[6]
        project.procPramArray[ci].mc_mi_min_method = settings[7]
        project.procPramArray[ci].mc_mi_search_tolerance = round(alog(settings[8])/alog(10))
        project.procPramArray[ci].mc_mi_use_translation = 1
        project.procPramArray[ci].mc_mi_use_rotation = 0
        project.procPramArray[ci].mc_mi_use_shear = 0
        project.procPramArray[ci].mc_mi_use_dilation = 0
        project.procPramArray[ci].mc_mi_normalized = 0

        n_lines -= 2

    endif else if (version ge 0.3) then begin

        settings = fltarr(15)

        readf, fileid, settings

        project.procPramArray[ci].mc_mi_rng_low    = settings[0] 
        project.procPramArray[ci].mc_mi_rng_hi     = settings[1]
        project.procPramArray[ci].mc_mi_slice      = settings[2]
        project.procPramArray[ci].mc_mi_refidx     = settings[3]
        project.procPramArray[ci].mc_mi_binsize    = alog(settings[4])/alog(2)
        project.procPramArray[ci].mc_mi_start_mesh = settings[5]
        project.procPramArray[ci].mc_mi_n_passes   = settings[6]
        project.procPramArray[ci].mc_mi_min_method = settings[7]
        project.procPramArray[ci].mc_mi_search_tolerance = round(alog(settings[8])/alog(10))
        project.procPramArray[ci].mc_mi_use_translation = settings[9]
        project.procPramArray[ci].mc_mi_use_rotation = settings[10]
        project.procPramArray[ci].mc_mi_use_shear = settings[11]
        project.procPramArray[ci].mc_mi_use_dilation = settings[12]
        project.procPramArray[ci].mc_mi_normalized = settings[13]

        n_lines -= 2

    endif else begin
        ;; prerelease versions contained no header, just settings
        ;; rewind to get the setting that was popped off when the
        ;; version string was read
        point_lun, fileid, 0

    endelse

    if (version le 0.2) then begin
        
        fmt = '(2I0,2F0)' 
        
        ;; read in the shifts
        while ~ eof(fileid) do begin
            
            readf, fileid, FORMAT=fmt, $
              sidx, aidx, shft_x, shft_y
            
            mc_data[sidx, aidx, 0] = shft_x
            mc_data[sidx, aidx, 1] = shft_y
            
        endwhile
        
    endif else if (version eq 0.3) then begin

        print, "Version 3"

        tmp = fltarr(7)
        fmt = '(2I0,7F0)' 
        
        ;; read in the shifts
        while ~ eof(fileid) do begin
            
            readf, fileid, FORMAT=fmt, sidx, aidx, tmp
            
            mc_data[sidx, aidx, *] = tmp
            
        endwhile

    endif else begin

        print, "Version 4"

        tform = fltarr(7)
        dim = ''

        fmt = '(A1,2I0,7F0)' 
        tar = project.dataArray[project.ci].state1

        progressBar = obj_new('progressbar', color='red', /fast_loop, $
                              title='Progress Bar')
        progressBar -> SetProperty, text='Applying Motion Correction'
        
        progressbar->start
        total_iter = float(n_lines) ;;(fdim+pdim+sdim)*adim*1.0
        curr_iter = float(0.0)
        ;; read in the shifts
        while ~ eof(fileid) do begin
            
            readf, fileid, FORMAT=fmt, dim, sidx, aidx, tform
            
            if (max(abs(tform)) eq 0) then continue

            case dim of 
                's': begin
                    img = mc_pad_image( reform( (*tar)[*,*,sidx,aidx] ), out_x=pd_x, out_y=pd_y)
                    img = mc_transform_image(img, tform)
                    (*tar)[*,*,sidx,aidx] = mc_unpad_image( img, pd_x, pd_y )

                end
                
                'r': begin
                    img = mc_pad_image( reform( (*tar)[sidx,*,*,aidx] ), out_x=pd_x, out_y=pd_y) 
                    img = mc_transform_image(img, tform)
                    (*tar)[sidx,*,*,aidx] = mc_unpad_image( img, pd_x, pd_y )
                end
                
                'p': begin
                    img = mc_pad_image( reform( (*tar)[*,sidx,*,aidx] ), out_x=pd_x, out_y=pd_y) 
                    img = mc_transform_image(img, tform)
                    (*tar)[*,sidx,*,aidx] = mc_unpad_image( img, pd_x, pd_y )
                end
                
                else:
             endcase

            if (curr_iter++ mod 30. eq 0) then $
              progressbar->update, curr_iter/total_iter * 100.0
            
            if (progressBar->checkcancel()) then begin
                
                print, "Motion Correction cancelled."
                ptr_free, tar
                close, fileid
                free_lun, fileid
                progressBar->destroy
                project.procpramarray[project.ci].mc_enable = 0
                mas_redraw
                mas_load_state_1
                return
                
            endif

        endwhile
        
        close, fileid
        free_lun, fileid
        progressbar->update, 100.0
        progressbar->destroy

        return

    endelse

    close, fileid
    free_lun, fileid

    ;; a separate function for application, since it might be useful
    ;; to have it separated from the actual file i/o
    
    if project.imndArray[ci].image_type eq 11 then begin
        ;; PAR/REC is already in magnitude mode
        mc_apply_shifts, mc_data
    endif else begin
        mc_apply_shifts, mc_data, /KSPACE
    endelse

end

;; Subroutine name: mas_mc_do_pca
;;
;; Created by: BT 2008-11-06
;;
;; Calling Information:
;;  scan_id      - IN  - The ID of the scan analyze indexed from zero
;;                       as appearing in the MAS main window list
;;  slicenum     - IN  - The slice to analyze (in slice direction)
;;                       The default is the sdim_start
;;  show         - IN  - If set, the "result" image will be displayed
;;
;;  result       - OUT - Contains the PCA imagery
;;  eigenvalues  - OUT - Contains the eigenvalues cooresponding to nth PC
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine: 
;;  Load in a data file of shift values
;;
;; Editing Information:
;;  Created 2008-11-06
;;
pro mas_mc_do_pca, scan_id, slicenum=slicenum, result=result, $
                   eigenvalues=eigenvalues, show=show

    common scan_data

    ;; make sure that the scan exists
    if (n_elements(scan_id) ne 0) then begin
       if (scan_id ge project.ni or scan_id lt 0) then begin
          junk = dialog_message('The requested scan does not exist', /error, /center)
          return
       endif else begin
          ci = scan_id
       endelse
    endif else begin
       ci = project.ci
    endelse
    
    data = project.dataarray[ci].state1
    if not ptr_valid(data) then return
    
    slicenum = keyword_set(slicenum) ? slicenum : project.procpramarray[ci].sdim_start
    slicedata = reform((*data)[*,*,slicenum,*])
    dims = size(slicedata, /dimensions)
    pca_data = fltarr(dims[0]*dims[1], dims[2])

    ;; set up the data to be analyize. 2D images -> 1D arrays
    for ar = 0, dims[2]-1 do begin

       pca_data[*,ar] = (reform(slicedata[*,*,ar]))[*]

    endfor

    result = pcomp(transpose(temporary(pca_data)), eigenvalues=eigenvalues, $
                   /standardize, /covariance)
    result = temporary(transpose(result))

    ;; now convert the 1D result back to 2D images
    for ar = 0, dims[2]-1 do begin

       slicedata[*,*,ar] = (reform(result[*,ar]))[*]

    endfor

    result = temporary(slicedata)

    ;; if the user wants to see the image, then display it
    if (keyword_set(show)) then begin
       
       bytimage = ptr_new(bytarr(dims), /no_copy)

       ;; bytscl'ing
       for ar = 0, dims[2]-1 do begin

          (*bytimage)[*,*,ar] = bytscl(result[*,*,ar])

       endfor
       
       ;; make them into a "sheet"
       img_struct = mas_make_sheet_display(bytimage, 1, 1)
       image = *img_struct.image

       ptr_free, img_struct.image
       ptr_free, bytimage

       mas_display_multi, image, $
                          tab_title=strcompress("PCA Imagery: Scan "+string(ci+1)+$
                                                ", Slice "+string(slicenum+1))

    endif

end

pro mas_mc_pca_analysis, first_scan=first_scan, second_scan=second_scan, $
                         slicenum=slicenum, show=show, $
                         start_ev=start_ev, end_ev=end_ev
   
    if (not keyword_set(slicenum)) then begin
       slicenum = 0
    endif

    if (arg_present(second_scan)) then begin

       mas_mc_do_pca, first_scan,  eigenvalues=ev_1, slicenum=slicenum, show=show
       mas_mc_do_pca, second_scan, eigenvalues=ev_2, slicenum=slicenum, show=show
       plot_title = strcompress('PCA Comparison, Scan '+string(first_scan+1)+' vs. '+string(second_scan+1))
       
    endif else begin
       
       plot_title = strcompress('PCA Analysis, Scan '+string(first_scan+1))
       mas_mc_do_pca, first_scan,  eigenvalues=ev_1, slicenum=slicenum, show=show
       
    endelse
    
    if (n_elements(ev_1) eq 0) then begin
       junk = dialog_message('PCA could not be computed.', /error, /center)
       return
    endif
    
    start_ev = keyword_set(start_ev) ? start_ev : 0
    end_ev   = keyword_set(end_ev)   ? end_ev   : n_elements(ev_1)

    end_ev   = (end_ev > 0)   < n_elements(ev_1)
    start_ev = (start_ev > 0) < end_ev
    
    pca_indices = indgen((end_ev - start_ev+1))+start_ev
    ;;pca_indices = indgen(n_elements(ev_1)-1)+1
    
    iplot, pca_indices, ev_1[pca_indices], color=[255,0,0], sym_index=2, $
           linestyle=2, name=strcompress('Scan #'+string(first_scan+1)), $
           identifier=foo, /insert_legend, title=plot_title, view_title=plot_title, $
           xtitle='nth Principal Component', ytitle='Eigenvalue'
    
    if (n_elements(ev_2) ne 0) then begin
       iplot, pca_indices, ev_2[pca_indices], color=[0,0,255], sym_index=5, $
              linestyle=5, name=strcompress('Scan #'+string(second_scan+1)),$
              overplot=foo, /insert_legend
    endif

end

pro mas_mc_pca_gui_event, event

    uname = widget_info(event.id, /uname)

    if (uname eq 'btn_cancel') then begin
       widget_control, event.top, /destroy
       return
    endif else if (uname ne 'btn_run') then begin
       ;; we ignore the other
       return
    endif

    widget_control, widget_info(event.top, find_by_uname='first_txt'), get_value=first
    first = fix(first)-1

    widget_control, widget_info(event.top, find_by_uname='second_txt'), get_value=second
    second = fix(second)-1

    widget_control, widget_info(event.top, find_by_uname='slice_sl'), get_value=slice
    slice = fix(slice)-1
    
    show = widget_info( widget_info(event.top, find_by_uname='show_btn'), /button_set )

    widget_control, event.top, /destroy

    if (second gt 0) then begin
        mas_mc_pca_analysis, first_scan=first, second_scan=second, slicenum=slice, show=show, start_ev=1, end_ev=20
     endif else begin
        mas_mc_pca_analysis, first_scan=first, slicenum=slice, show=show, start_ev=1, end_ev=20
     endelse

end


pro mas_mc_pca_gui

    common common_widgets
    common scan_data

    if (ptr_valid(project.dataarray[project.ci].state1) eq 0) then begin
       junk = dialog_message(['No scans found to run analysis.', $
                              'Please make sure that a scan is loaded', $
                              'in the MAS main window.'] , /error, /center)
       return
    endif

    base = widget_base(group_leader=wid_base_main, title='PCA Comparison', /column, /modal, /align_center)

    first_base = widget_base(base, /row)
    first_lbl = widget_label(first_base, value='First Scan #', /align_right)
    first_txt = widget_text(first_base, xsize=3, uname='first_txt', /edit)

    second_base = widget_base(base, /row)
    second_lbl = widget_label(second_base, value='Second Scan #', /align_right)
    second_txt = widget_text(second_base, xsize=3, uname='second_txt', /edit)

    show_base = widget_base(base, /row, /nonexclusive)
    show_btn = widget_button(show_base, value="Show PCA Imagery", uname='show_btn')

    slice_base = widget_base(base, /column)
    slice_sl = widget_slider(slice_base, title='Slice', value=project.procpramarray[project.ci].sdim_start+1, $
                             uname='slice_sl', $
                             max=(size(*project.dataarray[project.ci].state1, /dimensions))[2], $
                             min=1)
    
    btn = widget_button(base, value="Run Analysis", uname='btn_run')
    btn = widget_button(base, value="Cancel", uname='btn_cancel')

    widget_control, base, /realize
    xmanager, 'mas_mc_pca_gui', base

end
