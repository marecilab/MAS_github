;; $Id$
;;
 function philips_transform_Tang, apfhrl, inv=inv

     gradient = apfhrl
     ap = gradient[0]
     fh = gradient[1]
     rl = gradient[2]
    
     Tang_rl = float([ [1, 0      , 0       ], $
                       [0, cos(rl), -sin(rl)], $
                       [0, sin(rl),  cos(rl)] ])
    
     Tang_ap = float([ [cos(ap), 0, sin(ap)], $
                       [      0, 1, 0      ], $
                       [-sin(ap),0, cos(ap)] ])
    
     Tang_fh = float([ [cos(fh), -sin(fh), 0], $
                       [sin(fh),  cos(fh), 0], $
                       [      0,        0, 1] ])
    
;     Tang_outtmp = matrix_multiply(Tang_rl,Tang_ap)
;     Tang_out    = matrix_multiply(Tang_outtmp, Tang_fh)

     Tang_outtmp = Tang_rl     ## Tang_ap
     Tang_out    = Tang_outtmp ## Tang_fh

     tiny = where(abs(Tang_out) lt 1e-6, ct)
     if (ct gt 0) then begin
        Tang_out[tiny] = 0.0
     endif

     return, keyword_set(inv) ? invert(Tang_out) : Tang_out

 end

 function philips_get_tx_matrix, apfhrl, $
                                 orientation, $
                                 fat_shift, $
                                 fold_over, $
                                 user_defined=user_defined, overplus=overplus
                                 

     ident_3 = diag_matrix(replicate(1.0, 3))

 ;;;;;;;;;;;;;; TPO ;;;;;;;;;;;;;;;;;;;;;;

     Tpo_supine = float([ [ 1.0,  0.0,  0.0], $
                          [ 0.0,  1.0,  0.0], $
                          [ 0.0,  0.0,  1.0] ])

     Tpo_prone  = float([ [-1.0,  0.0,  0.0], $
                          [ 0.0, -1.0,  0.0], $
                          [ 0.0,  0.0,  1.0] ])

     Tpo_right_dec = float([ [ 0.0, -1.0,  0.0], $
                             [ 1.0,  0.0,  0.0], $
                             [ 0.0,  0.0,  1.0] ])
    
     Tpo_left_dec = float([ [ 0.0,  1.0,  0.0], $
                            [-1.0,  0.0,  0.0], $
                            [ 0.0,  0.0,  1.0] ])

 ;;;;;;;;;;;;;; TPP ;;;;;;;;;;;;;;;;;;;;;;

     Tpp_feetfirst = float([ [ 0.0, -1.0,  0.0], $
                             [-1.0,  0.0,  0.0], $
                             [ 0.0,  0.0,  1.0] ])
;;     Tpp_feetfirst = transpose(Tpp_feetfirst)

     Tpp_headfirst = float([ [ 0.0,  1.0,  0.0], $
                             [-1.0,  0.0,  0.0], $
                             [ 0.0,  0.0, -1.0] ])
;;     Tpp_headfirst = transpose(Tpp_headfirst)

     Tpom = Tpo_supine ## Tpp_headfirst

 ;;;;;;;;;;;;;; TSOM ;;;;;;;;;;;;;;;;;;;;;;

     Tsom_sag = float([ [ 0.0,  0.0, -1.0], $
                        [ 0.0, -1.0,  0.0], $
                        [ 1.0,  0.0,  0.0] ])

     Tsom_cor = float([ [ 0.0, -1.0,  0.0], $
                        [ 0.0,  0.0,  1.0], $
                        [ 1.0,  0.0,  0.0] ])

     Tsom_tra = float([ [ 0.0, -1.0,  0.0], $
                        [-1.0,  0.0,  0.0], $
                        [ 0.0,  0.0,  1.0] ])
 
 ;;;;;;;;;;;;;; TPREP ;;;;;;;;;;;;;;;;;;;;;;
   
     Tprep_par = ident_3

     Tprep_per = float([ [ 0.0, -1.0,  0.0], $
                         [ 1.0,  0.0,  0.0], $
                         [ 0.0,  0.0,  1.0] ])

 ;;;;;;;;;;;;;; TFSD ;;;;;;;;;;;;;;;;;;;;;;

     Tfsd_m = float([ [-1.0,  0.0,  0.0], $
                      [ 0.0,  1.0,  0.0], $
                      [ 0.0,  0.0,  1.0] ])

     Tfsd_p = float([ [ 1.0,  0.0,  0.0], $
                      [ 0.0, -1.0,  0.0], $
                      [ 0.0,  0.0,  1.0] ])

     Tfsd_s = float([ [-1.0,  0.0,  0.0], $
                      [ 0.0,  1.0,  0.0], $
                      [ 0.0,  0.0, -1.0] ])

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

     case orientation of

         1: begin                      ;; TRA
             Tsom = Tsom_tra
             print, "ORIENTATION: TRANSVERSE"
             case fold_over of
                 1: begin 
                    print, "FOLD_OVER: AP"
                     Tprep = Tprep_par  ;; AP
                     if (fat_shift eq 'L' or fat_shift eq 'P') then Tfsd = Tfsd_p
                     if (fat_shift eq 'R' or fat_shift eq 'F' or fat_shift eq 'A') then Tfsd = Tfsd_m
                     if (fat_shift eq 'H') then Tfsd = Tfsd_s
                 end
                 2: begin
                    print, "FOLD_OVER: RL"
                     Tprep = Tprep_per  ;; RL
                     if (fat_shift eq 'R' or fat_shift eq 'P') then Tfsd = Tfsd_p
                     if (fat_shift eq 'L' or fat_shift eq 'F' or fat_shift eq 'A') then Tfsd = Tfsd_m
                     if (fat_shift eq 'H') then Tfsd = Tfsd_s
                 end
             endcase
         end
         2: begin                      ;; SAG
             Tsom = Tsom_sag
             case fold_over of
                 1: begin
                     Tprep = Tprep_par  ;; AP
                     if (fat_shift eq 'F' or fat_shift eq 'A') then Tfsd = Tfsd_p
                     if (fat_shift eq 'L' or fat_shift eq 'H' or fat_shift eq 'P') then Tfsd = Tfsd_m
                     if (fat_shift eq 'R') then Tfsd = Tfsd_s
                 end
                 3: begin
                     Tprep = Tprep_per  ;; FH
                     if (fat_shift eq 'F' or fat_shift eq 'P') then Tfsd = Tfsd_p
                     if (fat_shift eq 'L' or fat_shift eq 'H' or fat_shift eq 'A') then Tfsd = Tfsd_m
                     if (fat_shift eq 'R') then Tfsd = Tfsd_s
                 end
             endcase
         end
         3: begin                      ;; COR
             Tsom = Tsom_cor
             case fold_over of
                 2: begin
                     Tprep = Tprep_par  ;; RL
                     if (fat_shift eq 'F' or fat_shift eq 'R') then Tfsd = Tfsd_p
                     if (fat_shift eq 'L' or fat_shift eq 'H' or fat_shift eq 'A') then Tfsd = Tfsd_m
                     if (fat_shift eq 'P') then Tfsd = Tfsd_s
                 end
                 3: begin
                     Tprep = Tprep_per  ;; FH
                     if (fat_shift eq 'F' or fat_shift eq 'L') then Tfsd = Tfsd_p
                     if (fat_shift eq 'R' or fat_shift eq 'H' or fat_shift eq 'A') then Tfsd = Tfsd_m
                     if (fat_shift eq 'P') then Tfsd = Tfsd_s
                 end
             endcase
         end
     endcase

     if (keyword_set(user_defined)) then begin
         frame='LPH'
     endif else if keyword_set(overplus) then begin
         frame='LPH'
     endif else begin
         frame='MPS'
     endelse

     print, "FRAME: "+frame
     case frame of

         'MPS': begin
             Txyz = Tprep ## Tfsd ;;matrix_multiply(Tprep, Tfsd)
             return, Txyz
         end

         'LPH': begin           
            Tsom = invert(Tsom)
            Ta = philips_transform_Tang(apfhrl, /inv)
            return, Ta ## Tsom   ;;matrix_multiply(Ta, Tsom)
          end
      endcase
     
  end


; Subroutine name: process_Philips_image
; Created by: CD 9/27/07
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; this routine will rescale all Philips data so values proportional to MR signal
; and will reorient the image by flipping the image in the Y direction
; the images must be flipped in the Y direction as we define the 0,0 coordinate
; as lower left hand corner while Philips defines this as top left hand corner
; since Philips zero pads images, the final fdim and pdim are not known until
; after the zero-padding is removed by this function, these values are returned
; so that MAS can be updated appropriately
;
; Editing Information:

FUNCTION process_Philips_image, imgdata_ptr, rs_slope_ptr, rs_intercept_ptr, sc_slope_ptr $
                                ,FINAL_RDIM=final_fdim, FINAL_PDIM=final_pdim

    DEBUG = 0
    
    progressbar = OBJ_NEW('progressbar', Color='red', Text='Processing data from Philips file.',/NOCANCEL)
    progressbar -> Start
    
    img_sz = SIZE(*imgdata_ptr)
    
    fdim = img_sz[1]
    pdim = img_sz[2]
    sdim = 1
    adim = 1
    
    CASE img_sz[0] OF
        4: BEGIN
            sdim = img_sz[3]
            adim = img_sz[4]
        END
        
        3: sdim = img_sz[3]
        ELSE:
    ENDCASE
    
    IF adim EQ 1 THEN BEGIN
        IF sdim EQ 1 THEN img_sz[0] = 2 ELSE img_sz[0] = 3
    ENDIF
    
    IF N_ELEMENTS(*rs_slope_ptr) EQ 0 THEN BEGIN
        rs_slope_ptr = PTR_NEW(FLTARR(sdim*adim))
        (*rs_slope_ptr)[*] = 1
    ENDIF
    
    IF N_ELEMENTS(*rs_intercept_ptr) EQ 0 THEN BEGIN
        rs_intercept_ptr = PTR_NEW(FLTARR(sdim*adim))
        (*rs_intercept_ptr)[*] = 0
    ENDIF
    
    IF N_ELEMENTS(*sc_slope_ptr) EQ 0 THEN BEGIN
        sc_slope_ptr = PTR_NEW(FLTARR(sdim*adim))
        (*sc_slope_ptr)[*] = 1
    ENDIF
    
    
    imgdata = FLTARR(fdim,pdim,sdim,adim)
    
    ;; Philips scaling is very mysterious
    FOR jj = 0, adim-1 DO BEGIN
        FOR ii = 0, sdim-1 DO BEGIN
            rs_iter = (*rs_slope_ptr)[jj*sdim+ii]
            ss_iter = (*sc_slope_ptr)[jj*sdim+ii]
            ri_iter = (*rs_intercept_ptr)[jj*sdim+ii]
            
            imgdata[*,*,ii,jj] = ((REVERSE(ROTATE((*imgdata_ptr)[*,*,ii,jj],2)) * rs_iter) + ri_iter) / (rs_iter * ss_iter)

;            (img * rs) + ri    img   +  ri
;            --------------- =  ---     ------
;               (rs * ss)       ss      rs * ss

        ENDFOR
        
        progressBar -> Update, (float(jj)/float(adim-1))*100.0
    ENDFOR
    
    PTR_FREE, imgdata_ptr
    
    ;; this code will get rid of the zero-padding applied by Phillips to ensure
    ;; that the image is square, this padding is added in excess to the image so
    ;; its removal does not alter the proportions of the image, due to the way the
    ;; data is reported by Philips, we do not know the fdim or pdim of the actual
    ;; image prior to the removal of the padding
    
    ;;first sum the adim and sdim
    CASE img_sz[0] OF
        4: test_sum = TOTAL(TOTAL(imgdata,4),3)
        3: test_sum = TOTAL(imgdata,3)
        2: test_sum = imgdata
    ENDCASE
    
                                ;sum along the y direction
    test_x   = TOTAL(test_sum,2)
    data_x   = WHERE(test_x NE 0)
    
                                ;sum along the x direction
    test_y   = TOTAL(test_sum,1)
    data_y   = WHERE(test_y NE 0)
    
    progressbar -> Destroy
    
    IF DEBUG EQ 1 THEN BEGIN
        final_fdim = fdim
        final_pdim = pdim
        
        print, 'DEBUGGING PHILIPS DATA'
        
        IF (WHERE(test_x EQ 0))[0] NE -1 THEN BEGIN
            imgdata[WHERE(test_x EQ 0),*,*,*] = 99999
        ENDIF
        
        IF (WHERE(test_y EQ 0))[0] NE -1 THEN BEGIN
            imgdata[*,WHERE(test_y EQ 0),*,*] = 99999
        ENDIF
        
        RETURN, PTR_NEW(imgdata)
    ENDIF ELSE BEGIN
        final_fdim = fdim;(SIZE(data_x))[1]
        final_pdim = pdim;(SIZE(data_y))[1]
        
        RETURN, PTR_NEW(imgdata);PTR_NEW(imgdata[data_x,data_y,*,*])
    ENDELSE
END


; Subroutine name: add_Philips_to_MAS
; Created by: CD 2/15/07
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:
; CD : moved to readRECdata.pro to avoid compilation errors

PRO add_Philips_to_MAS, data

    COMMON scan_data
    
    ni = project.ni
    
    project.imndArray[ni].fdim = data.fdim
    project.imndArray[ni].pdim = data.pdim
    project.imndArray[ni].sdim = data.sdim
    project.imndArray[ni].adim = data.adim
    
    project.imndArray[ni].scan_date = data.scan_date
    project.imndArray[ni].scan_name = data.scan_name
    project.imndArray[ni].recov_time = data.recov_time
    project.imndArray[ni].echo_time = data.echo_time
    project.imndArray[ni].n_echo = data.n_echo
    project.imndArray[ni].n_avg = data.n_avg
    
    project.imndArray[ni].pn_avg = PTR_NEW(data.pn_avg)
    
    project.imndArray[ni].f_fov = data.f_fov
    project.imndArray[ni].p_fov = data.p_fov
    project.imndArray[ni].s_fov = data.s_fov
    
    project.procPramArray[ni].orientation_Imin = data.Imin
    project.procPramArray[ni].orientation_Jmin = data.Jmin
    project.procPramArray[ni].orientation_Kmin = data.Kmin
    project.procPramArray[ni].have_orientation = 1

    project.imndArray[ni].slices = data.slices
    project.imndArray[ni].thick = data.thick

    project.imndarray[ni].f_voxsz = data.f_voxsz
    project.imndarray[ni].p_voxsz = data.p_voxsz
    project.imndarray[ni].s_voxsz = data.s_voxsz

    project.imndArray[ni].dimensions = data.dimensions
    
    project.imndArray[ni].fat_sup = data.fat_sup
    
    project.imndArray[ni].file_Path = data.file_Path
    
    IF data.DICOMtype NE 4 THEN BEGIN
        data.display_Name = STRMID(data.display_Name,0,STRPOS(data.display_Name,'\',/REVERSE_SEARCH)+1)
    ENDIF
    
    project.imndArray[ni].display_Name = data.display_Name
    
    project.imndArray[ni].orientation = data.orientation
    
    IF data.diffimage EQ 1 THEN BEGIN
        project.imndArray[ni].bval_Array  = data.diffdata.bvals
        project.imndArray[ni].b_matrix    = data.diffdata.bmatrix
        project.imndArray[ni].angle_theta = data.diffdata.theta
        project.imndArray[ni].angle_phi   = data.diffdata.phi
        project.imndArray[ni].n_bvals	  = n_elements(*data.diffdata.bvals)
        project.imndArray[ni].small_delta = data.diffdata.sdeltas
        project.imndArray[ni].big_delta   = data.diffdata.bdeltas
    ENDIF

    project.imndArray[ni].acq_matrix  = data.acq_matrix
    
    project.imndArray[ni].PARRECtype = data.PARRECtype
    
    IF data.PARRECtype NE 0 THEN BEGIN
        project.imndArray[ni].REC_rs_int = data.REC_rs_int
        project.imndArray[ni].REC_rs_slp = data.REC_rs_slp
        project.imndArray[ni].REC_sc_slp = data.REC_sc_slp
    ENDIF
    
    IF data.Philips_Data EQ 1 THEN BEGIN
        project.imndArray[ni].Philips_init_fdim = data.init_fdim
        project.imndArray[ni].Philips_init_pdim = data.init_pdim
    ENDIF
    
    project.imndArray[ni].image_type = 11
    
    project.imndArray[ni].DICOMtype  = data.DICOMtype
    project.imndArray[ni].DICOMfiles = data.DICOMfiles
    
    project.imndarray[ni].k_fdim_span = data.init_fdim
    project.imndarray[ni].k_pdim_span = data.init_pdim
    project.imndarray[ni].k_sdim_span = data.sdim
    project.imndarray[ni].k_fdim_span_max = data.init_fdim
    project.imndarray[ni].k_pdim_span_max = data.init_pdim
    project.imndarray[ni].k_sdim_span_max = data.sdim
    
    project.imndArray[ni].state1_load_procedure = 'mas_read_philips_parrec_load_state1'
    
    project.current_path = file_dirname(data.file_path)

    project.scan_Open_Flag = 1
    
    project.ci = project.ni
    
    project.ni = project.ni + 1
    
    mas_add_current_scan_to_scan_selection
    
    mas_redraw_GUI
    
    update_status_bar, ''
    
END


; Subroutine name: readRECdata
; Created by: CD 9/24/07
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;	reads in all data from Philips REC file
;
; Editing Information:
;   BT, 2008-04-04 - make reading version 4.2 PAR/REC more
;                    like version 4.1. There are only minor
;                    differences. 

FUNCTION readRECdata

    COMMON scan_data
    
    CI = project.ci
    
    rec_file = project.imndArray[CI].file_Path + '.REC'
    
    fdim = project.imndArray[CI].Philips_init_fdim
    pdim = project.imndArray[CI].Philips_init_pdim
    sdim = project.imndArray[CI].sdim
    adim = project.imndArray[CI].adim


    CASE project.imndArray[CI].PARRECtype OF
        1: BEGIN
           ;; FIXME - endian handling should be delegated to the
           ;;         user. This is a quick-fix for ppc Mac, and
           ;;         assumes that _all_ data needs swapping on ppc,
           ;;         which may not be the case.
            rawd = READ_BINARY(rec_file,DATA_TYPE=12, ENDIAN=(project.big_endian eq 1) ? 'little' : 'native')

            ;; data organized, first all sdim then next adim
            imgdata_ptr = PTR_NEW(REFORM(TEMPORARY(rawd),fdim,pdim,sdim,adim))
        END
        
        2: BEGIN
            ;; data organized, first all adim then next sdim
            ;; can utilize a fdim x pdim x (X) where X = adim*sdim matrix so
            ;; that when give to process_Philips_image, the scaling parameters
            ;; will index with the correct image.
            ;; (Not sure that the above still applies (BT 20080404)

            rawd = READ_BINARY(rec_file,DATA_TYPE=12, data_dims=[fdim,pdim,sdim,adim], ENDIAN=(project.big_endian eq 1) ? 'little' : 'native')
            imgdata_ptr = PTR_NEW(temporary(rawd))

        END
    ENDCASE
    
    imgdata_ptr = process_Philips_image(imgdata_ptr,                      $
                                        project.imndArray[CI].REC_rs_slp, $
                                        project.imndArray[CI].REC_rs_int, $
                                        project.imndArray[CI].REC_sc_slp, $
                                        FINAL_RDIM=final_fdim,            $
                                        FINAL_PDIM=final_pdim)

    ;; commented, BT 2008-04-04
    ;; reorganize images into fdim x pdim x sdim x adim
    ;;IF project.imndArray[CI].PARRECtype EQ 2 THEN BEGIN
    ;;
    ;;    imgdata = FLTARR(final_fdim,final_pdim,sdim,adim)
    ;;    
    ;;    FOR ii = 0, sdim-1 DO BEGIN
    ;;        FOR jj = 0, adim-1 DO BEGIN
    ;;            imgdata[*,*,ii,jj] = REVERSE(ROTATE((*imgdata_ptr)[*,*,ii,jj],2))
    ;;        ENDFOR
    ;;    ENDFOR
    ;;
    ;;    PTR_FREE, imgdata_ptr
    ;;    imgdata_ptr = PTR_NEW(imgdata)
    ;;ENDIF

    project.imndArray[CI].fdim = final_fdim
    project.imndArray[CI].pdim = final_pdim
    mas_redraw_GUI

    RETURN, imgdata_ptr
END


; Subroutine name: readPARdata
; Created by: CD 9/24/07
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;	reads in all data from Philips PAR file
;
; Editing Information:
;  Modified BT 3/12/08 to handle missing PAR/REC files gracefully

PRO readPARdata

    COMMON scan_data

    path_pattern = DIALOG_PICKFILE( Title = 'Select either of the PAR/REC files to open.' ,$
                                    path=project.current_path)
    
    IF (path_pattern EQ '') THEN RETURN
    
    pos = STRPOS(STRMID(path_pattern,0,STRLEN(path_pattern)),'.',/REVERSE_SEARCH)
    path_pattern = STRMID(path_pattern,0,pos)
    
    par_file = path_pattern + '.PAR'
    rec_file = path_pattern + '.REC'
    
    if (file_test(par_file, /READ) eq 0) then begin
        junk = dialog_message('Unable to read PAR data from location.', /error, /center)
        return
    endif else if (file_test(rec_file, /READ) eq 0) then begin
        junk = dialog_message('Unable to read REC data from location.', /error, /center)
        return
    endif
    
    textline = ''
    linenum  = 0
    exitloop = 0
    headerdone = 0

    deg2rad = !PI / 180
    
    OPENR, fid, par_file , /GET_LUN
    
    WHILE (EOF(fid) NE 1) DO BEGIN
        
        READF, fid, textline
        linenum = linenum + 1

        CASE linenum OF
            8:  BEGIN
                pos = STRPOS(textline,'V4.1',/REVERSE_SEARCH)
                IF pos NE -1 THEN BEGIN
                    PARRECtype = 1
                    dtstart = 98
                ENDIF ELSE BEGIN
                    pos = STRPOS(textline,'V4.2',/REVERSE_SEARCH)
                    IF pos NE -1 THEN BEGIN
                        PARRECtype = 2
                        dtstart = 100
                    ENDIF ELSE BEGIN
                        void = DIALOG_MESSAGE('PARREC version not supported.',/ERROR, /center)
                        RETURN
                    ENDELSE
                ENDELSE
            END
            
            13: BEGIN
                scan_name = STRTRIM(STRMID(textline,STRPOS(textline,':',/REVERSE_SEARCH)+2),2)
            END
            
            14: BEGIN
                sequence_name = STRTRIM(STRMID(textline,STRPOS(textline,':',/REVERSE_SEARCH)+2),2)
            END
            
            15:	BEGIN
               months = ['Nul','Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', $
                         'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
                scan_date = STRTRIM(STRMID(textline,STRPOS(textline,':')+2),2)
                temp = stregex(scan_date, '^([0-9]{4})\.([0-9]{2})\.([0-9]{2}) / ([0-9]{2}):([0-9]{2}):([0-9]{2})', $
                               /extract, /subexp)
                scan_date = months[fix(temp[2])]+' '+temp[3]+ ' '+temp[1]+' '+temp[4]+':'+temp[5]+':'+temp[6]
            END
            
            21: BEGIN
                n_echo = FIX(STRMID(textline,STRPOS(textline,':',/REVERSE_SEARCH)+2))
            END
            
            22: BEGIN
                sdim = FIX(STRMID(textline,STRPOS(textline,':',/REVERSE_SEARCH)+2))
                slices = sdim
            END
            
            25: BEGIN
                orientation = STRTRIM(STRMID(textline,STRPOS(textline,':',/REVERSE_SEARCH)+2),2)
            END
            
            26: BEGIN
                prep_str = STRTRIM(STRMID(textline,STRPOS(textline,':',/REVERSE_SEARCH)+2),2)
                
                CASE prep_str OF
                    'Anterior-Posterior': prep_dir = 1
                    'Right-Left'        : prep_dir = 2
                    'Foot-Head'         : prep_dir = 3
                    'Feet-Head'         : prep_dir = 3
                ENDCASE
            END
            
            30: BEGIN
                recov_time = FLOAT(STRMID(textline,STRPOS(textline,':',/REVERSE_SEARCH)+2))
            END
            
            31: BEGIN
                str1 = STRTRIM(STRMID(textline,STRPOS(textline,':',/REVERSE_SEARCH)+2),2)
                brk1 = STRPOS(str1,' ')
                brk2 = STRPOS(str1,' ',brk1+2)
                
                ap_fov = FLOAT(STRMID(str1,0,brk1)) / 10
                fh_fov = FLOAT(STRMID(str1,brk1+1,brk2)) / 10
                rl_fov = FLOAT(STRMID(str1,brk2+1)) / 10
            END
            
            33: BEGIN
                str1 = STRTRIM(STRMID(textline,STRPOS(textline,':',/REVERSE_SEARCH)+2),2)
                brk1 = STRPOS(str1,' ')
                brk2 = STRPOS(str1,' ',brk1+2)
                
                ap_ang = FLOAT(STRMID(str1,0,brk1)) * deg2rad
                fh_ang = FLOAT(STRMID(str1,brk1+1,brk2)) * deg2rad
                rl_ang = FLOAT(STRMID(str1,brk2+1)) * deg2rad
                apfhrl_ang =  [ap_ang, fh_ang, rl_ang]
            END
            
            39: BEGIN
                ;;check for fat suppression, we are actually checking specifically for SPIR
                ;;so if the SPIR flag is set to 1, then set fat_sup to 2, so that in MAS the
                ;;label is set to 'SPIR' instead of just 'Yes'
                fat_sup = FIX(STRMID(textline,STRPOS(textline,':',/REVERSE_SEARCH)+2))
                IF fat_sup EQ 1 THEN fat_sup = 2
            END
            
            42: BEGIN
                ;; Check to see if this is a diffusion scan. If it is,
                ;; then we need to account for the extra philips
                ;; volume
                diffusion = FIX(STRMID(textline,STRPOS(textline,':',/REVERSE_SEARCH)+2))
                diffusion = (diffusion ne 0) ? 1 : 0
            END

            45: BEGIN
                ;; there is sometimes an extra adim when using Philips defined gradient directions that
                ;; is not present when using user defined gradient directions
                ;;adim = FIX(STRMID(textline,STRPOS(textline,':',/REVERSE_SEARCH)+2)) + 1
                n_diff_volumes = FIX(STRMID(textline,STRPOS(textline,':',/REVERSE_SEARCH)+2))
                adim = n_diff_volumes
                
                ;; later in this process, we'll compare the image
                ;; count in the .par file with this value and adjust
                ;; the DTI parameters so that the extra volume is not
                ;; factored into the ADT processing

            END
            
            46: BEGIN
                exitloop = 1
                headerdone = 1
            END
            
            ELSE :
        ENDCASE
        
        IF exitloop EQ 1 THEN BREAK
        
    ENDWHILE
    
    CLOSE, fid
    
;; >>>>>>>> BT    
;;    CASE PARRECtype OF
;;        1: BEGIN
;;            dstart = 98
;;            adim_idx = INTARR(adim)
;;            adim_idx[0] = 0
;;            FOR ii = 1, adim-1 DO adim_idx[ii] = adim_idx[ii-1]+sdim
;;        END
;;        
;;        2: BEGIN
;;            adim_idx = INDGEN(adim)
;;            dstart = 100
;;        END
;;    ENDCASE
;; <<<<<<<<< BT

    CASE PARRECtype OF
        1: BEGIN
            dstart = 98
        END
        
        2: BEGIN
            dstart = 100
        END

    ENDCASE

    n_diff_volumes_idx = INTARR(n_diff_volumes)
    n_diff_volumes_idx[0] = 0
    FOR ii = 1, n_diff_volumes-1 DO n_diff_volumes_idx[ii] = n_diff_volumes_idx[ii-1]+sdim

    img_params = READ_ASCII(par_file, count=num_img, DATA_START=dstart)
    
    num_img = num_img - 1
    
; do to zero padding of the image by philips to always have a square image
; these are the dimensions of the zero-padded square, the real dimensions
; are not known prior to reading in the image as part of the mas_load_state1
; routine


    ;; total count of images (zero-indexed) divided by number of slices.
    adim       = (max(img_params.field01[6,*])+1)/sdim

    init_fdim  = img_params.field01[9,0]
    init_pdim  = img_params.field01[10,0]
    
    thick      = img_params.field01[22,0]
    gap        = img_params.field01[23,0]
    spacing_x  = img_params.field01[28,0]
    spacing_y  = img_params.field01[29,0]

    voxdim = [spacing_x, spacing_y, thick+gap]/10. ; mm => cm
    ;;voxdim /= min(voxdim)

    n_avg      = img_params.field01[34,0]
    dimensions = 2
    echo_time  = img_params.field01[30,0]
    
    img_orient = img_params.field01[25,0]
    
    slice_idx  = img_params.field01[0,0:num_img-1]
    
; save these parameters, they will be used when reading in the REC (image) file
    rescale_intercept = img_params.field01[11,*]
    rescale_slope     = img_params.field01[12,*]
    scale_slope       = img_params.field01[13,*]
    
    pn_avg = REFORM(TRANSPOSE(img_params.field01[34,n_diff_volumes_idx]))
    
    bvals = img_params.field01[33,n_diff_volumes_idx]

    APcos = img_params.field01[45,n_diff_volumes_idx]
    FHcos = img_params.field01[46,n_diff_volumes_idx]
    RLcos = img_params.field01[47,n_diff_volumes_idx]
    
    bmatrix = FLTARR(6,n_diff_volumes)
    theta   = FLTARR(n_diff_volumes)
    phi     = FLTARR(n_diff_volumes)

    if (diffusion gt 0) then begin

        case img_orient of
            1: begin ;; Transverse
                acq_matrix = double([ [ 1, 0, 0, 0], $ 
                                      [ 0,-1, 0, 0], $
                                      [ 0, 0, 1, 0], $
                                      [ 0, 0, 0, 1] ])
            end
            2: begin ;; Sagittal
                acq_matrix = double([ [ 0, 1, 0, 0], $
                                      [ 0, 0, 1, 0], $
                                      [ 1, 0, 0, 0], $
                                      [ 0, 0, 0, 1] ])
            end
            3: begin ;; Coronal
                acq_matrix = double([ [ 1, 0, 0, 0], $
                                      [ 0, 0, 1, 0], $
                                      [ 0, 1, 0, 0], $
                                      [ 0, 0, 0, 1] ])
            end
        endcase

;; >>>>>>>>>>>>>>
;; DO NOT USE
;;         pat_ori       = 0
;;         pat_prep      = prep_dir
;;         slc_ori       = img_orient
;;         fat_shift_dir = 'P'
;;         fold_over     = prep_dir
;;         user_defined  = 1
;;         overplus      = 0
        
;;         phil_tx = fltarr(4,4)
;;         phil_tx[0:2,0:2]  =  philips_get_tx_matrix(apfhrl_ang, $
;;                                                    slc_ori, $
;;                                                    fat_shift_dir, $
;;                                                    fold_over, $
;;                                                    user_defined=user_defined, $
;;                                                    overplus=overplus)
;;         phil_tx[3,3] = 1.0

;;         acq_matrix = phil_tx

;;         print, 'ACQ_MATRIX:' & print, acq_matrix

;;         grad_dirs_rect = [ APcos, FHcos, RLcos ]
;;         grad_dirs_orig = grad_dirs_rect

;;         grad_dirs = vert_t3d(grad_dirs_rect, matrix=acq_matrix, /no_copy)

;;         AP_cos  = grad_dirs[0,*]
;;         FH_cos  = grad_dirs[1,*]
;;         RL_cos  = grad_dirs[2,*]


;;         APcos_b = grad_dirs[0,*]
;;         FHcos_b = grad_dirs[1,*]
;;         RLcos_b = grad_dirs[2,*]

;; <<<<<<<<<<<<<<<<

;orig
      grad_dirs_rect = [ RLcos, APcos, FHcos ]
      grad_dirs_orig = [ RLcos, APcos, FHcos ]
       
      grad_dirs = vert_t3d(grad_dirs_rect, matrix=acq_matrix, /no_copy)
 
      RLcos_b = grad_dirs[0,*]
      APcos_b = grad_dirs[1,*]
      FHcos_b = grad_dirs[2,*]
      
      ;;RLcos = grad_dirs[0,*]
      ;;APcos = grad_dirs[1,*]
      ;;FHcos = grad_dirs[2,*]
     
    endif else begin

        RLcos_b = RLcos
        APcos_b = APcos
        FHcos_b = FHcos

    endelse


    ;; reconstruct B-matrix from gradient cosines in PAR file
    ;; Idea:
    ;;  Sum_ij(b_{ij} D_{ij}) = b g^t D g.
    ;;  Expanding the right-hand side:
    ;;    g^t D g = Sum_ij(g_i g_j D_{ij}).
    ;;  Hence:
    ;;     b_{ij} = b g_i g_j, where
    ;;   g_x = RJcos )
    ;;   g_y = APcos )  <-- These come from 
    ;;   g_z = FHcos )      the .PAR file
    ;;   b   = bvals )

    FOR jj = 0, n_diff_volumes-1 DO BEGIN

        bval_iter  = bvals[jj]
        
;; ORIG
        xc = RLcos_b[jj] 
        yc = APcos_b[jj]
        zc = FHcos_b[jj]

;; This is the order reported in the par file.
;;        xc = APcos_b[jj]
;;        yc = FHcos_b[jj]
;;        zc = RLcos_b[jj]

        xx = bval_iter * xc * xc    ;;RLcos_iter^2
        yy = bval_iter * yc * yc    ;;APcos_iter^2
        zz = bval_iter * zc * zc    ;;FHcos_iter^2
        
        xy = bval_iter * xc * yc * 2 ;;RLcos_iter * APcos_iter * 2
        xz = bval_iter * xc * zc * 2 ;;RLcos_iter * FHcos_iter * 2
        yz = bval_iter * yc * zc * 2 ;;APcos_iter * FHcos_iter * 2
        
        bmatrix[*,jj] = [xx, yy, zz, xy, xz, yz]
        
                                ; convert to / save spherical
                                ; coordinates
        ;; Use the original values from par file to create the
        ;; gradient angles. This is to be consistent with the was
        ;; the imnd file reports the angles/b-matrix for bruker data
;;orig        
         xc = RLcos[jj]  ;; RLcos_iter
         yc = APcos[jj]  ;; APcos_iter
         zc = FHcos[jj]  ;; FHcos_iter
        
        pval = SQRT(xc^2 + yc^2 + zc^2)
        
                                ; theta is angle with respect to z-axis
        IF pval EQ 0 THEN theta[jj] = 0 $
        ELSE theta[jj] = ACOS(zc/pval)
        
        phi[jj] = ATAN(yc,xc)
        
                                ; convert to degrees from radians
        phi[jj] = phi[jj] * 180 / !PI
        theta[jj] = theta[jj] * 180 / !PI
        
                                ; rotate negative angles by 360 degrees
        idx = WHERE((phi LT 0),zcount)
        IF zcount NE 0 THEN phi[idx] = phi[idx] + 360
        
    ENDFOR
    
    print, 'WARNING : USING HARDCODED PHILIPS SMALL/BIG DELTA VALUES'
    bdeltas = FLTARR((SIZE(bvals))[1])
    sdeltas = FLTARR((SIZE(bvals))[1])
    
    bdeltas[*] = 29.73          ;ms
    sdeltas[*] = 16.91          ;ms
    
    diffdata = { bmatrix:PTR_NEW(bmatrix), $
                 acq_matrix:ptr_new(acq_matrix), $
                 bvals:PTR_NEW(bvals),     $
                 theta:PTR_NEW(theta),     $
                 phi:PTR_NEW(phi),         $
                 bdeltas:PTR_NEW(bdeltas), $
                 sdeltas:PTR_NEW(sdeltas) }
    
; adjust parameters for slice orientation
; img_orient, 1=TRA 2=SAG 3=COR
; prep_dir (phase dir), 1=AP 2=RL 3=FH
; images appear to be oriented regardless of phase direction, not sure if this affects fov
; prior to rotate the images during read in, the image orientation uses philips ap, rl
; fh is defined in reverse (i.e. hf)

;; Below, (i,j,k) refer to state1 data array indices: (*state1)[i,j,k,...]
    CASE img_orient OF
        1: BEGIN ;; Transverse
            ;; in-file i(0,...,end) => scanner R->L
            ;; in-file j(0,...,end) => scanner P->A
            ;; in-file k(0,...,end) == scanner F->H
            f_fov = rl_fov
            p_fov = ap_fov
            s_fov = fh_fov
            Imin = 4 ;; R
            Jmin = 3 ;; P
            Kmin = 2 ;; I
            CASE prep_dir OF
                1: orientation = ['Transverse','RL','(PA)']
                2: orientation = ['Transverse','(RL)','PA']
            ENDCASE
        END
        
        2: BEGIN ;; Sagittal
            ;; in-file i(0,...,end) == scanner A->P
            ;; in-file j(0,...,end) == scanner F->H
            ;; in-file k(0,...,end) == scanner L->R
            f_fov = ap_fov
            p_fov = fh_fov
            s_fov = rl_fov
            Imin = 0 ;; 'A'
            Jmin = 2 ;; 'I'
            Kmin = 1 ;; 'L'
            CASE prep_dir OF
                1: orientation = ['Sagittal','(AP)','FH']
                3: orientation = ['Sagittal','AP','(FH)']
            ENDCASE
        END
        
        3: BEGIN ;; Coronal
            ;; in-file i(0,...,end) == scanner R->L
            ;; in-file j(0,...,end) == scanner F->H
            ;; in-file k(0,...,end) == scanner A->P
            f_fov = rl_fov
            p_fov = fh_fov
            s_fov = ap_fov
            Imin = 4 ;; 'R'
            Jmin = 2 ;; 'I'
            Kmin = 0 ;; 'A'
            CASE prep_dir OF
                2: orientation = ['Coronal','(RL)','FH']
                3: orientation = ['Coronal','RL','(FH)']
            ENDCASE
        END
    ENDCASE
    
    data = { image_Number:INDGEN(num_img) $
             , fdim:0 $
             , pdim:0 $
             , init_fdim:init_fdim $
             , init_pdim:init_pdim $
             , sdim:sdim $
             , adim:adim $
             , sequence_name:sequence_name $
             , scan_date:scan_date $
             , scan_name:sequence_name $
             , recov_time:recov_time $
             , echo_time:echo_time $
             , n_echo:n_echo $
             , n_avg:n_avg $
             , pn_avg:pn_avg $
             , f_fov:f_fov $
             , p_fov:p_fov $
             , s_fov:s_fov $
             , Imin:Imin $
             , Jmin:Jmin $
             , Kmin:Kmin $
             , slices:slices $
             , thick:thick $
             , f_voxsz: voxdim[0] $
             , p_voxsz: voxdim[1] $
             , s_voxsz: voxdim[2] $
             , dimensions:dimensions $
             , file_Path:path_pattern $
             , display_Name:path_pattern $
             , orientation:orientation $
             , diffimage:diffusion $
             , diffdata:diffdata $
             , DICOMtype:4 $
             , DICOMfiles:PTR_NEW(par_file) $
             , chkval:1 $
             , PARRECtype:PARRECtype $
             , REC_rs_int:PTR_NEW(rescale_intercept) $
             , REC_rs_slp:PTR_NEW(rescale_slope) $
             , REC_sc_slp:PTR_NEW(scale_slope) $
             , fat_sup:fat_sup $
             , Philips_Data:1 $
             , acq_matrix: ptr_new(acq_matrix) $
           }
    
    add_Philips_to_MAS, data
    
END

pro mas_read_philips_parrec_load_state1

    common scan_data, project
    forward_function mas_extract_data
    forward_function DICOM_load_image
    
    ci = project.ci
    
    ;; if .dcm file extension then image extraction still neccessary
    ;; otherwise multi-frame dicom is already extracted
    dcmfound = STRPOS(project.imndArray[ci].file_Path,'.dcm')
    
    IF (dcmfound NE -1) THEN BEGIN
    
        project.dataArray[CI].state1 = mas_extract_data()
        
    ENDIF ELSE BEGIN
    
        ;; reading of the REC file happens here.
        project.dataArray[CI].state1 =  $
            DICOM_load_image( project.imndArray[CI].DICOMtype, $
            project.imndArray[CI].DICOMfiles)
            
    ENDELSE
   
    ;; only MI based MC works for these images
    ;; ideally we need to display a message
    IF project.procPramArray[CI].mc_enable eq 2 then begin
        mas_motion_correct
    endif
    
    project.dataArray[CI].state1 = temporary(mas_interpolate( project.dataArray[CI].state1 ))
    
    ;;mas_shift,            project.dataArray[CI].state1, ISDICOM=1
    ;;mas_signal_selection, project.dataArray[CI].state1
    mas_smoothing,        project.dataArray[CI].state1
    
    IF (project.procPramArray[ci].image_domain_filterflag eq 'Image Domain') THEN BEGIN
        image_domain_filter
    ENDIF
    
    project.procpramarray[ci].state_1 = 1
    
end

pro mas_read_philips_parrec


end
