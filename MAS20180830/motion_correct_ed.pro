;; $Id$
;;

FUNCTION mas_2D_fft_shift, imgdata, rshift, pshift
    COMMON scan_data

    CI = project.ci
    
    fdim = project.imndArray[CI].fdim
    pdim = project.imndArray[CI].pdim
    sdim = project.imndArray[CI].sdim
    adim = project.imndArray[CI].adim

    r_dk = (1/project.imndArray[CI].f_fov)
    p_dk = (1/project.imndArray[CI].p_fov)
    
    svec = COMPLEXARR(fdim,pdim)
    
    r_val = 2*!PI*r_dk*rshift
    p_val = 2*!PI*p_dk*pshift
    
    FOR idx=0,fdim-1 DO BEGIN
        svec[idx,*] = REPLICATE(EXP(-COMPLEX(0.,idx*r_val)),pdim)
    ENDFOR
    
    FOR idx=0,pdim-1 DO BEGIN
        svec[*,idx] = svec[*,idx] * REPLICATE(EXP(-COMPLEX(0.,idx*p_val)),fdim)
    ENDFOR


    ;;imgdata = TEMPORARY(FFT(imgdata * svec))

    imgdata = imgdata * svec
    
    RETURN, imgdata
END

PRO mas_motion_correct_der, REF_IDX, NUM_PIX_CHK, SHFTS_PER_VXL, pathtosave

    mas_motion_correct_mi, REF_IDX, -1, 1, .125, 2 ;;NUM_PIX_CHK

END

PRO mas_motion_correct_ed, REF_IDX, NUM_PIX_CHK, SHFTS_PER_VXL, pathtosave

    COMMON scan_data

    CI = project.ci
    
    fdim = project.imndArray[CI].fdim
    pdim = project.imndArray[CI].pdim
    sdim = project.imndArray[CI].sdim
    adim = project.imndArray[CI].adim
    
    ffov = project.imndArray[CI].f_fov
    pfov = project.imndArray[CI].p_fov
    
    ;; mc_res = FLTARR(fdim,pdim,sdim,adim)
    mc_res = COMPLEXARR(fdim,pdim,sdim,adim)
    mc_kspace = COMPLEXARR(fdim,pdim,sdim,adim)
    
    rstp = ffov / (fdim * SHFTS_PER_VXL)
    pstp = pfov / (pdim * SHFTS_PER_VXL)
    
    num_stps = (SHFTS_PER_VXL*NUM_PIX_CHK)*2+1
    
    progressbar = OBJ_NEW('progressbar', Color='red', $
                          Text='Performing Motion Correction . . .', $
                          /NOCANCEL)
    progressbar -> Start

    ;; Do for every slice in the array
    FOR k = 0, sdim-1 DO BEGIN
        ;; FOR k = 5,5 DO BEGIN
        
        ref_img = mas_2D_fft_shift(REFORM((*project.dataArray[CI].state1)[*,*,k,REF_IDX]),0,0)
        ref_img = TEMPORARY(FFT(ref_img))
        
        ref_img = SOBEL(ABS(SHIFT(ref_img,fdim/2,pdim/2)))
        
        print, 'MOTION CORRECTION : SLICE ' + STRTRIM(STRING(k),2) + $
          ' output: adim, xshift (vxls), yshift (vxls)'
        
        ;; align each adim
        FOR aidx=0,adim-1 DO BEGIN
            
            IF aidx NE REF_IDX THEN BEGIN
                test_array = FLTARR(num_stps,num_stps,3)
                
                rshift = -rstp * SHFTS_PER_VXL * NUM_PIX_CHK
                FOR i = 0, num_stps-1 DO BEGIN
                    
                    pshift = -pstp * SHFTS_PER_VXL * NUM_PIX_CHK
                    FOR j = 0, num_stps-1 DO BEGIN
                        
                        test_img = mas_2D_fft_shift(REFORM((*project.dataArray[CI].state1)[*,*,k,aidx]),rshift,pshift)
                        test_img = TEMPORARY(FFT(test_img))
                        
                        test_img = SOBEL(ABS(SHIFT(test_img,fdim/2,pdim/2)))
                        
                        test_array[i,j,*] = [TOTAL((test_img - ref_img)^2),rshift,pshift]

                        pshift = pshift + pstp
                    ENDFOR
                    rshift = rshift + rstp
                ENDFOR
                
                min_idx = 0
                min_val = MIN(REFORM(test_array[*,*,0]),min_idx)
                
                f_shft = REFORM((test_array[*,*,1])[*])
                p_shft = REFORM((test_array[*,*,2])[*])
                
                f_shft[WHERE(ABS(f_shft) LT rstp)] = 0
                p_shft[WHERE(ABS(p_shft) LT pstp)] = 0
                
                print, aidx, f_shft[min_idx] / (rstp * SHFTS_PER_VXL), $
                  p_shft[min_idx] / (pstp * SHFTS_PER_VXL)
                
                final_img = mas_2D_fft_shift(REFORM((*project.dataArray[CI].state1)[*,*,k,aidx]),f_shft[min_idx],p_shft[min_idx])
                mc_kspace[*,*,k,aidx] = final_img
                
                final_img = TEMPORARY(FFT(final_img))
                
                ;; for debugging
                ;;mc_res[*,*,k,aidx] = ABS(SHIFT(final_img,fdim/2,pdim/2))
                mc_res[*,*,k,aidx] = final_img
            ENDIF ELSE BEGIN
                ;; fft reference image without shifting
                final_img = $
                  mas_2D_fft_shift(REFORM((*project.dataArray[CI].state1)[*,*,k,aidx]),0,0)
                
                mc_kspace[*,*,k,aidx] = final_img
                
                final_img = TEMPORARY(FFT(final_img))
                
                ;; for debugging
                ;;mc_res[*,*,k,aidx] = ABS(SHIFT(final_img,fdim/2,pdim/2))
                
                mc_res[*,*,k,aidx] = final_img
            ENDELSE
        ENDFOR
        progressBar -> Update, (float(k)/float(sdim-1))*100.0
    ENDFOR
    
    progressbar -> Destroy
    
    project.dataArray[CI].state1 = PTR_NEW(mc_kspace);;res)
    
    ;; writing the raw corrected images to file
    
    sz_data = size(mc_kspace)
    save_path = pathtosave
    
    ;; Now we write the header and the data to the desired path.
    OpenW, lun0, save_path+'.raw', /SWAP_IF_LITTLE_ENDIAN, /get_lun
    WriteU, lun0, sz_data, mc_kspace
    Free_LUN, lun0
    
    ;; Write the GEO_INFO filae which contains the minimum amount of information to
    ;; read in and display the data.
    write_geo_info_file, save_path, sz_data
END


PRO mas_motion_correct_ed_GUI

    HEAP_GC

    COMMON MC_GUI_DATA, MC_WINDOW_BASE, MC_button, MC_refidx, MC_shftvxl, $
      MC_pixchk, MC_filepath, MC_changefolder
    COMMON scan_data

    title = STRING('MAS Motion Correction Parameters (ED)')

    MC_WINDOW_BASE = widget_base(TITLE=title, ROW=2)
    
    COL_ONE = WIDGET_BASE(MC_WINDOW_BASE, COLUMN=3, /ALIGN_RIGHT)    
    ROW_ONE = WIDGET_BASE(COL_ONE, ROW=2, /ALIGN_RIGHT)
    
    MC_refidx = WIDGET_LABEL(ROW_ONE, VALUE='Reference Image Index :')
    MC_refidx = WIDGET_TEXT(ROW_ONE, UNAME='MC_refidx', VALUE='0', /EDITABLE, /ALL_EVENTS )
    
    ROW_TWO = WIDGET_BASE(COL_ONE, ROW=2, /ALIGN_RIGHT)
    
    MC_shftvxl = WIDGET_LABEL(ROW_TWO, VALUE='Shifts Per Voxel :')
    MC_shftvxl = WIDGET_TEXT(ROW_TWO, UNAME='MC_shftvxl', VALUE='5', /EDITABLE, /ALL_EVENTS )
    
    ROW_THR = WIDGET_BASE(COL_ONE, ROW=2, /ALIGN_RIGHT)
    
    MC_pixchk = WIDGET_LABEL(ROW_THR, VALUE='Voxels to Check :')
    MC_pixchk = WIDGET_TEXT(ROW_THR, UNAME='MC_pixchk', VALUE='2', /EDITABLE, /ALL_EVENTS )
    
    COL_TWO = WIDGET_BASE(MC_WINDOW_BASE, COLUMN=1, /ALIGN_RIGHT)
    
    MC_changefolder =WIDGET_BUTTON(COL_TWO, VALUE ='Save Folder', /ALIGN_LEFT)
    
    MC_button = WIDGET_BUTTON(COL_TWO, VALUE ='Start Motion Correction', /ALIGN_LEFT)
    
    COL_THREE = WIDGET_BASE(COL_TWO, COLUMN=1, /ALIGN_RIGHT)

    MC_filepath = WIDGET_TEXT(COL_THREE, UNAME='MC_filepath',$
                              VALUE = project.current_path , $
                              /EDITABLE, /ALL_EVENTS )


    widget_control, MC_WINDOW_BASE, /realize
    xmanager, 'mas_motion_correct_ed_GUI', MC_WINDOW_BASE, Event_Handler='mas_motion_correct_ed_GUI_event'
END

PRO mas_motion_correct_ed_GUI_event, Event

    COMMON MC_GUI_DATA
    common scan_data

    ci = project.ci

    CASE Event.id OF
        
        MC_button: BEGIN
            WIDGET_CONTROL, MC_refidx, GET_VALUE=TEMP
            REF_IDX = (FIX(TEMP))[0]
            project.procPramArray[ci].mc_ed_refidx = REF_IDX

            WIDGET_CONTROL, MC_shftvxl, GET_VALUE=TEMP
            SHFTS_PER_VXL = (FIX(TEMP))[0]
            project.procPramArray[ci].mc_ed_spv = SHFTS_PER_VXL

            WIDGET_CONTROL, MC_pixchk, GET_VALUE=TEMP
            NUM_PIX_CHK = (FIX(TEMP))[0]
            project.procPramArray[ci].mc_ed_vxl_chk = NUM_PIX_CHK

            WIDGET_CONTROL, MC_filepath, GET_VALUE=TEMP
            pathtosave = TEMP
            project.procPramArray[ci].mc_ed_save_fldr = pathtosave

            mas_motion_correct_ed, REF_IDX, NUM_PIX_CHK, SHFTS_PER_VXL, pathtosave
            
            WIDGET_CONTROL, MC_WINDOW_BASE, /DESTROY
        END
        
        MC_changefolder: BEGIN
            file = DIALOG_PICKFILE(/WRITE)
            WIDGET_CONTROL, MC_filepath, SET_VALUE=file   
         
        END
        
        ELSE:
        
;;        mas_fft, project.dataArray[CI].state1
        
    ENDCASE
END
