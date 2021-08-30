;; $Id$
;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; DIFFUSION WEIGHTING EPI NOTES
;
;Diffusion weighting with the simple, two gradient-pulse spin echo
;(appropriate for EPI)

;For the simple two gradient-pulse spin echo diffusion weighting,
;the detailed b-matrix formulation (with matrix elements, bij) of the
;diffusion weighting can be equated a signal b-value for the diffusion
;weighting with the direction specified by the gradient unit vector,
;g = (g1 g2, g3), where gi is the direction cosine of the unit vector
;in the ith direction.
;
;Sum i,j (bij Dij) = b g^tDg, for i,j = 1,2,3                      [1]
;
;The vector-matrix product on the right-hand side of equation [1], can
;be expanded to the following (where g^t is the transpose of g, and
;D is the rank-2 diffusion tensor):
;
;g^tDg = g1 D11 G1 + g1 D12 g2 + g1 D13 g3
;      + g2 D21 G1 + g2 D22 g2 + g2 D23 g3
;      + g3 D31 G1 + g3 D32 g2 + g3 D33 g3
;
;     = Sum i,j (gi gj Dij), for i,j = 1,2,3                      [2]
;
;Therefore the b-matrix elements can be equated to the product of the
;direction cosines of the gradient unit vector.
;
;bij = b gi gj                                                    [3]
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; PROGRAMING EXCEPTIONS / NOTES

; PHILIPS PRIVATE CLASSES USED
;
; NAME				DICOM CODE
;
; Slice Index		2001	100A
; # Slices			2001	1018
; bval				2001	1003
; RLcos				2005	10B0
; APcos				2005	10B1
; FHcos				2005	10B2
; scale slope		2005	100E

; OTHER NOTES
; Dimensions parameter set = 2 if not defined in DICOM (0018,0023)
;
; Processing Note : if a DICOM file has a valid frame (0028,0008)
; then file is assumed to be multi-frame(image), if not defined and
; the number of slices (2001,1018) is defined, then assumed to be
; a multi-image multi-file dataset (one per file), else the DICOM file
; must be a single-image single-file dataset

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; DICOM types					;
; 1 = single file, single image ;
; 2 = multi-file, multi-frame	;
; 3 = single-file, multi-frame	;
; 4 = PAR/RED					;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; Subroutine name: DICOM_load_image
; Created by: CD 3/22/07
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

FUNCTION DICOM_load_image, DICOMtype, FLIST

   COMMON scan_data

   CI = project.ci
   
   FLIST = *FLIST
   num_files = N_ELEMENTS(FLIST)

   
   CASE DICOMtype OF
      
      1:	BEGIN           ;SINGLE FILE - SINGLE IMAGE
         oDICOM = OBJ_NEW( 'IDLffDICOM' )
         IF oDICOM->READ(FLIST) NE 1 THEN BEGIN
            void = DIALOG_MESSAGE('DICOM : File Reloading Error.',/ERROR)
            RETURN, -1
         ENDIF
         
         rs_intercept_ptr = ptr_new(0.0)                ; PTR_NEW(FLOAT(*(oDICOM->GetValue('0028'x,'1052'x,/NO_COPY))[0]))
         rs_slope_ptr     = ptr_new(1.0)                ; PTR_NEW(FLOAT(*(oDICOM->GetValue('0028'x,'1053'x,/NO_COPY))[0]))
         sc_slope_ptr     = ptr_new(1.0)                ;PTR_NEW(FLOAT(*(oDICOM->GetValue('2005'x,'100E'x,/NO_COPY))[0]))
         imgdata_ptr      = PTR_NEW(oDICOM->GetValue('7fe0'x,'0010'x,/NO_COPY))
         
         imgdata_ptr = process_Philips_image(imgdata_ptr, rs_slope_ptr,rs_intercept_ptr, sc_slope_ptr $
                                             ,FINAL_RDIM=final_fdim, FINAL_PDIM=final_pdim)
         
         project.imndArray[CI].fdim = final_fdim
         project.imndArray[CI].pdim = final_pdim
         mas_redraw_GUI
      END
      
      2: 	BEGIN           ;MULTI_FILE - MULTI FRAME
         oDICOM = OBJ_NEW( 'IDLffDICOM' )
         progressbar = OBJ_NEW('progressbar', Color='red', Text='Extracting DICOM data from multiple files',/NOCANCEL)
         progressbar -> Start
         
         slices = project.imndArray[CI].slices
         
         fdim = project.imndArray[CI].Philips_init_fdim
         pdim = project.imndArray[CI].Philips_init_pdim
         sdim = project.imndArray[CI].sdim
         adim = project.imndArray[CI].adim
         
         raw_img      = FLTARR(fdim,pdim,num_files)
         slice_idx    = FLTARR(num_files)
         inst_idx     = FLTARR(num_files)
         rs_intercept = FLTARR(num_files)
         rs_slope     = FLTARR(num_files)
         sc_slope     = FLTARR(num_files)
         
         FOR ii = 0, num_files-1 DO BEGIN
                                ; check is DICOM readable, else discard
            IF oDICOM->READ(FLIST[ii]) EQ 1 THEN BEGIN
                                ; slice index, scale_slope USES PHILIPS PRIVATE CLASS
               slice_idx[ii]    = FLOAT(*(oDICOM->GetValue('2001'x,'100A'x,/NO_COPY))[0])
               inst_idx[ii]     = FLOAT(*(oDICOM->GetValue('0020'x,'0013'x,/NO_COPY))[1])
               rs_intercept[ii] = 0.0                          ; FLOAT(*(oDICOM->GetValue('0028'x,'1052'x,/NO_COPY))[0])
               rs_slope[ii]     = 1.0                          ; FLOAT(*(oDICOM->GetValue('0028'x,'1053'x,/NO_COPY))[0])
               sc_slope[ii]     = 1.0                          ; FLOAT(*(oDICOM->GetValue('2005'x,'100E'x,/NO_COPY))[0])
               raw_img[*,*,ii]  = *(oDICOM->GetValue('7fe0'x,'0010'x,/NO_COPY))[0]
            ENDIF ELSE BEGIN
               void = DIALOG_MESSAGE('DICOM : File Reloading Error.',/ERROR)
               RETURN, -1
            ENDELSE
            
            progressBar -> Update, (float(ii)/float(num_files-1))*100.0
            
         ENDFOR
         
         progressbar -> Destroy
         
         rawdata_ptr = process_Philips_image(PTR_NEW(raw_img),PTR_NEW(rs_slope),PTR_NEW(rs_intercept), $
                                             PTR_NEW(sc_slope),FINAL_RDIM=final_fdim, FINAL_PDIM=final_pdim)
         
         unq_slice = slice_idx(UNIQ(slice_idx))
         
         IF ((SIZE(unq_slice))[1]) NE slices THEN BEGIN
            void = DIALOG_MESSAGE('DICOM : File Reloading Error.'/ERROR)
            RETURN, -1
         ENDIF
         
                                ; ALIGN DATA ALONG (fdim,pdim,sdim,adim)
         imgdata = FLTARR(final_fdim,final_pdim,sdim,adim)
                                ; process all slices
         FOR ii=0,slices-1 DO BEGIN
            curr_idx = WHERE((slice_idx EQ unq_slice[ii]))
            
                                ;get instances
            sort_idx = SORT(inst_idx(curr_idx))
            curr_idx = curr_idx(sort_idx)
            
            imgdata[*,*,ii,*] = (*rawdata_ptr)[*,*,curr_idx]
         ENDFOR
         
         project.imndArray[CI].fdim = final_fdim
         project.imndArray[CI].pdim = final_pdim
         mas_redraw_GUI
         
         imgdata_ptr = PTR_NEW(imgdata)
      END
      
      3:	BEGIN           ;SINGLE FILE - MULTI-FRAME
         oDICOM = OBJ_NEW( 'IDLffDICOM' )
         IF oDICOM->READ(FLIST) NE 1 THEN BEGIN
            void = DIALOG_MESSAGE('DICOM : File Reloading Error.',/ERROR)
            RETURN, -1
         ENDIF
         
         fdim = project.imndArray[CI].Philips_init_fdim
         pdim = project.imndArray[CI].Philips_init_pdim
         sdim = project.imndArray[CI].sdim
         adim = project.imndArray[CI].adim
         
         raw_img_ptr = oDICOM->GetValue('7fe0'x,'0010'x,/NO_COPY)
         
         rs_intercept_ptr = oDICOM->GetValue('0028'x,'1052'x,/NO_COPY)
         rs_slope_ptr     = oDICOM->GetValue('0028'x,'1053'x,/NO_COPY)
         sc_slope_ptr     = oDICOM->GetValue('2005'x,'100E'x,/NO_COPY)
         
         num_slopes = sdim*adim
         
         rs_intercept = FLTARR(num_slopes)
         rs_slope     = FLTARR(num_slopes)
         sc_slope     = FLTARR(num_slopes)
         raw_img      = FLTARR(fdim,pdim,num_slopes)
         
         FOR ii = 0, num_slopes -1 DO BEGIN
            rs_intercept[ii] = *(rs_intercept_ptr[ii])
            rs_slope[ii]     = *(rs_slope_ptr[ii])
            sc_slope[ii]     = *(sc_slope_ptr[ii])
            raw_img[*,*,ii]  = *(raw_img_ptr[ii])
         ENDFOR
         
         rawdata_ptr = process_Philips_image(PTR_NEW(raw_img),PTR_NEW(rs_slope),PTR_NEW(rs_intercept), $
                                             PTR_NEW(sc_slope),FINAL_RDIM=final_fdim, FINAL_PDIM=final_pdim)
         
         imgdata = FLTARR(final_fdim,final_pdim,sdim,adim)
         
         FOR ii = 0, sdim-1 DO BEGIN
            FOR jj = 0, adim-1 DO BEGIN
               imgdata[*,*,ii,jj] = REVERSE(ROTATE((*rawdata_ptr)[*,*,ii*adim+jj],2))
            ENDFOR
         ENDFOR
         
         project.imndArray[CI].fdim = final_fdim
         project.imndArray[CI].pdim = final_pdim
         mas_redraw_GUI
         
         imgdata_ptr = PTR_NEW(imgdata)
      END
      
      4: BEGIN                  ;READ PAR_REC IMAGE
         imgdata_ptr = readRECdata()
      END
      
   ENDCASE
   
   if (obj_valid(oDICOM)) then OBJ_DESTROY, oDICOM
   HEAP_GC
   
   RETURN, imgdata_ptr
END


; Subroutine name: DICOM_calc_DTI_tensor
; Created by: CD 3/07/07
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

FUNCTION DICOM_calc_DTI_tensor, adim, ibvals, iRLcos, iAPcos, iFHcos, img_orient=img_orient
  
    bmatrix = FLTARR(6,adim)
    theta   = FLTARR(adim)
    phi     = FLTARR(adim)
    
    acq_matrix = double([ [ 1, 0, 0, 0], $
                          [ 0, 1, 0, 0], $
                          [ 0, 0, 1, 0], $
                          [ 0, 0, 0, 1] ])
 
    if (keyword_set(img_orient)) then begin
       case img_orient of
          'Transverse': begin ;; Transverse
             acq_matrix = double([ [ 1, 0, 0, 0], $
                                   [ 0,-1, 0, 0], $
                                   [ 0, 0, 1, 0], $
                                   [ 0, 0, 0, 1] ])
          end
          'Sagittal': begin ;; Saggital
             acq_matrix = double([ [ 0, 1, 0, 0], $
                                   [ 0, 0, 1, 0], $
                                   [ 1, 0, 0, 0], $
                                   [ 0, 0, 0, 1] ])
          end
          'Coronal': begin ;; Coronal
             acq_matrix = double([ [ 1, 0, 0, 0], $
                                   [ 0, 0, 1, 0], $
                                   [ 0, 1, 0, 0], $
                                   [ 0, 0, 0, 1] ])
          end
          else:
       endcase
    endif
    
    grad_dirs_rect = transpose([ [iRLcos], [iAPcos], [iFHcos] ])
    grad_dirs_orig = transpose([ [iRLcos], [iAPcos], [iFHcos] ])
    
    grad_dirs = vert_t3d(grad_dirs_rect, matrix=acq_matrix)
    
    RLcos_b = grad_dirs[0,*]
    APcos_b = grad_dirs[1,*]
    FHcos_b = grad_dirs[2,*]

                                ; set values for b0 image
    bmatrix(*,0) = [0,0,0,0,0,0]
    theta(0) = 0
    phi(0) = 0
    
                                ; skip b0 image
    FOR jj = 1, adim-1 DO BEGIN
                                ; extract data
       
       bval  = ibvals[jj]
       RLcos = RLcos_b[jj]          ; X
       APcos = APcos_b[jj]          ; Y
       FHcos = FHcos_b[jj]          ; Z
       
                                ; RLcos = iAPcos[jj]
                                ; APcos = iRLcos[jj]
       
                                ; bmatrix
       
       xx = bval * RLcos^2
       yy = bval * APcos^2
       zz = bval * FHcos^2
       
;		 double the off-diagonal terms
       xy = bval * RLcos * APcos * 2.
       xz = bval * RLcos * FHcos * 2.
       yz = bval * APcos * FHcos * 2.
       
                                ; convert to / save spherical coordinates
       xc = iRLcos[jj]
       yc = iAPcos[jj]
       zc = iFHcos[jj]
       
       pval = SQRT(xc^2 + yc^2 + zc^2)
       
                                ; theta is angle with respect to z-axis
       IF pval EQ 0 THEN theta[jj] = 0 $
       ELSE theta[jj] = ACOS(zc/pval)
       
       phi[jj] = ATAN(yc,xc)
       
       bmatrix[*,jj] = [xx, yy, zz, xy, xz, yz]
       
                                ; convert to degrees from radians
       phi[jj] = phi[jj] * 180 / !PI
       theta[jj] = theta[jj] * 180 / !PI
       
                                ; rotate negative angles by 360 degrees
       idx = WHERE((phi LT 0),zcount)
       IF zcount NE 0 THEN phi[idx] = phi[idx] + 360
       
    ENDFOR
    
    print, 'WARNING : USING HARDCODED PHILIPS SMALL/BIG DELTA VALUES'
    bdeltas = FLTARR((SIZE(ibvals))[1])
    sdeltas = FLTARR((SIZE(ibvals))[1])
    
    bdeltas[*] = 29.73          ;ms
    sdeltas[*] = 16.91          ;ms
    
    data = {bmatrix:PTR_NEW(bmatrix), $
            acq_matrix:ptr_new(acq_matrix), $
            bvals:PTR_NEW(ibvals), $
            theta:PTR_NEW(theta), $
            phi:PTR_NEW(phi), $
            bdeltas:PTR_NEW(bdeltas), $
            sdeltas:PTR_NEW(sdeltas) }
    
    RETURN, data
END


; Subroutine name: get_basic_dicom_params
; Created by: CD 2/15/07
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

FUNCTION get_basic_dicom_params, oDICOM
	; READ IN BASIC PARAMETERS

    ;get the freq. dimension
    temp = oDICOM->GetValue('0028'x,'0010'x,/NO_COPY)
    fdim = *temp[0]
   
                                ;get the phase dimension
    temp = oDICOM->GetValue('0028'x,'0011'x,/NO_COPY)
    pdim = *temp[0]
    
                                ;get the slice thickness
    temp = oDICOM->GetValue('0018'x,'0050'x,/NO_COPY)
    thick = *temp[0]
             
    temp = oDICOM->GetValue('0028'x,'0030'x,/NO_COPY)
    rp_pixdim = float(strsplit(*temp[0], '\', /extract))

    voxel_dims = [ rp_pixdim[0], rp_pixdim[1], float(thick) ]/10. ;;mm => cm

                       ;store pointer to # of frames (may be invalid)
    frame_ptr = oDICOM->GetValue('0028'x,'0008'x,/NO_COPY)
    
                                ;store pointer to # of slices (may be invalid) / USES PHILIPS PRIVATE CLASS
    slice_ptr = oDICOM->GetValue('2001'x,'1018'x,/NO_COPY)
    
                                ;check for SPIR fat suppresion
    SPIR = oDICOM->GetValue('2001'x,'1021'x,/NO_COPY)
    
                                ; set fat_sup to 2, so that the label 'SPIR' is displayed and not just a 'Yes'
    IF STRCMP(*SPIR[0],'Y') EQ 1 THEN fat_sup =2 ELSE fat_sup = 0
    
    data = {fdim:fdim, pdim:pdim, thick:thick, frame_ptr:frame_ptr, $
            slice_ptr:slice_ptr, fat_sup:fat_sup, voxel_dims:voxel_dims }
    
    RETURN, data
END


; Subroutine name: get_other_dicom_params
; Created by: CD 2/15/07
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

FUNCTION get_other_dicom_params, oDICOM, bparams
	; READ IN OTHER PARAMETERS

                                ;get the image number
    temp = oDICOM->GetValue('0020'x,'0013'x,/NO_COPY)
    image_Number = (*temp[0]) - 1
    
                                ;get the pixel spacing in mm (format: y/x)
                                ;calculate field of view display in cm units
    temp = oDICOM->GetValue('0028'x,'0030'x,/NO_COPY)
    slashpos = STRPOS(*temp[0],(get_dir_slash())[0]) ;
    f_fov = bparams.fdim * FLOAT(STRMID(*temp[0],0,slashpos-1)) / 10
    p_fov = bparams.pdim * FLOAT(STRMID(*temp[0],slashpos+1,STRLEN(*temp[0]))) / 10
    
                                ;get the dimensions
    temp = oDICOM->GetValue('0018'x,'0023'x,/NO_COPY)
    
                                ;if not defined, assume is 2D
    IF STRCMP(*temp[0],'1D') EQ 1 THEN dimensions = 1 ELSE $
       IF STRCMP(*temp[0],'2D') EQ 1 THEN dimensions = 2 ELSE $
          IF STRCMP(*temp[0],'3D') EQ 1 THEN dimensions = 3 ELSE dimensions = 2
    
                                ;the orientation is usually HFS

    temp = oDICOM->GetValue('2001'x,'100B'x, /no_copy)
    temp = strcompress(*temp[0], /remove_all)
    case temp of 
       'TRANSVERSAL': img_orient = 'Transverse'
       'SAGITTAL': img_orient = 'Sagittal'
       'CORONAL': img_orient = 'Coronal'
       else: img_orient = 'Unknown'
    endcase

    temp = oDICOM->GetValue('0018'x,'5100'x,/NO_COPY)
    orientation = [img_orient, *temp[0]]
    
                                ;get the sequence name
    temp = oDICOM->GetValue('0018'x,'9005'x,/NO_COPY)
                                ;if not available, display NA
    IF PTR_VALID(temp) EQ 1 THEN  sequence_name =  *temp[0] $
    ELSE sequence_name = 'NA'
    
                                ;get the scan date
    temp = oDICOM->GetValue('0008'x,'0023'x,/NO_COPY)
    scan_date = *temp[0]
    
                                ;get the scan name
    temp = oDICOM->GetValue('0008'x,'103E'x,/NO_COPY)
    scan_name = *temp[0]
    
                                ;get the reocovery time
    temp = oDICOM->GetValue('0018'x,'0080'x,/NO_COPY)
    recov_time = *temp[0]
    
                                ;get the echo time
    temp = oDICOM->GetValue('0018'x,'0081'x,/NO_COPY)
    echo_time = *temp[0]
    
                                ;get the number of echos
    temp = oDICOM->GetValue('0018'x,'0086'x,/NO_COPY)
    n_echo = *temp[0]
    
                                ;get the number of averages
    temp = oDICOM->GetValue('0018'x,'0083'x,/NO_COPY)
    n_avg = *temp[0]
    
    data = {image_Number:image_Number, f_fov:f_fov, p_fov:p_fov, $
            dimensions:dimensions, orientation:orientation, sequence_name:sequence_name, $
            scan_date:scan_date, scan_name:scan_name, recov_time:recov_time, $
            echo_time:echo_time, n_echo:n_echo, n_avg:n_avg}
    
    RETURN, data
    
END


; Subroutine name: exam_card_GUI_event
; Created by: CD 2/14/07
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

PRO exam_card_GUI_event, event
	IF event.clicks EQ 2 THEN BEGIN
		WIDGET_CONTROL, event.top, GET_UVALUE=info, /No_Copy

		acq_id = info.unq_acq(event.index)
		flist = info.FLIST(WHERE((info.acq_idx EQ acq_id)))

		; load selected protocol
		load_DICOM_multi_file, flist, info.oDICOM, info.location

		WIDGET_CONTROL, event.top, SET_UVALUE=info, /No_Copy
	END
END


; Subroutine name: load_exam_card_GUI
; Created by: CD 2/15/07
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

PRO exam_card_GUI_create, location, oDICOM

    progressbar = OBJ_NEW('progressbar', Color='red', Text='Extracting available protocols',/NOCANCEL)
    progressbar -> Start

    ;;retrieve the directory from filename
    lastslash = STRPOS(location,(get_dir_slash())[0],REVERSE_SEARCH=1)
    location = STRMID(location,0,lastslash+1)

    FLIST = FILE_SEARCH(location,'IM_*',count=num_files)
    
    acq_idx  = STRARR(num_files)
    seq_idx  = STRARR(num_files)
    date_idx = STRARR(num_files)
    name_idx = STRARR(num_files)
    
    ;; dicoms come in sorted strangely, so we will sort
    
    image_numbers = lonarr(num_files)
    flist_sorted = strarr(num_files)
    
    FOR ii = 0, num_files-1 DO BEGIN
       IF oDICOM->Read(FLIST[ii]) EQ 1 THEN BEGIN
          
          temp = oDICOM->GetValue('0020'x,'000E'x,/NO_COPY)
          acq_idx[ii] = *temp[0]
          
          ;;get the sequence name
          temp = oDICOM->GetValue('0018'x,'9005'x,/NO_COPY)
          ;;if not available, display NA
          IF PTR_VALID(temp) EQ 1 THEN seq_idx[ii] = *temp[0] $
          ELSE seq_idx[ii] = 'NA'
          
          ;;get the scan date
          temp = oDICOM->GetValue('0008'x,'0023'x,/NO_COPY)
          date_idx[ii] = *temp[0]
          
          ;;get the scan name
          temp = oDICOM->GetValue('0008'x,'103E'x,/NO_COPY)
          name_idx[ii] = *temp[0]
          
          ;;get the image_number for sorting
          temp = oDICOM->GetValue('0020'x,'0013'x,/no_copy)
          image_numbers[ii] = *temp[1]
          
       ENDIF
       
       progressBar -> Update, (float(ii)/float(num_files-1))*100.0
       
    ENDFOR
    
    progressbar -> Destroy
    
    sorted_images = sort(image_numbers)
    flist = flist[sorted_images]
    
    ;; create selection screen
    ;; get unique protocals
    
    unq_idx = UNIQ(acq_idx)
    
    unq_seq  = seq_idx(unq_idx)
    unq_acq  = acq_idx(unq_idx)
    unq_name = name_idx(unq_idx)
    unq_date = date_idx(unq_idx)
    
    display_title = 'Please double click to select protocol.'
    baseID = WIDGET_BASE(Row=1, Title=display_title, TLB_Size_Events=1)
    
    listID = WIDGET_LIST(baseID , FRAME=1 , SENSITIVE=1, $
                         SCR_XSIZE=100 ,SCR_YSIZE=250 ,XSIZE=11 ,YSIZE=0)
    
    seqID = WIDGET_LIST(baseID , FRAME=1 , SENSITIVE=1, $
                        SCR_XSIZE=100 ,SCR_YSIZE=250 ,XSIZE=11 ,YSIZE=0)
    
    dateID = WIDGET_LIST(baseID , FRAME=1 , SENSITIVE=1, $
                         SCR_XSIZE=100 ,SCR_YSIZE=250 ,XSIZE=11 ,YSIZE=0)
    
    WIDGET_CONTROL, listID , SET_VALUE= unq_name ,SENSITIVE=1
    WIDGET_CONTROL, seqID , SET_VALUE= unq_seq ,SENSITIVE=0
    WIDGET_CONTROL, dateID , SET_VALUE= unq_date ,SENSITIVE=0
    
    WIDGET_CONTROL, baseID, /Realize
    
    info = {FLIST:FLIST, unq_acq:unq_acq, acq_idx:acq_idx, oDICOM:oDICOM, $
            num_files:num_files, location:location}
    
    WIDGET_CONTROL, baseID, Set_UValue=info, /No_Copy
    
    XMANAGER, 'get_dicom_header_and_image', baseID, Event_Handler='exam_card_GUI_event', /NO_BLOCK

END


; Subroutine name: load_DICOM_multi_file
; Created by: CD 2/15/07
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

PRO load_DICOM_multi_file, FLIST, oDICOM, location

    num_files = ((SIZE(FLIST))[1])

    ;; read in first file in list then retrieve bparams
    tsh = oDICOM->READ(FLIST[0])
    
    bparams = get_basic_dicom_params(oDICOM)
    oparams = get_other_dicom_params(oDICOM,bparams)
        
    slice_ptr = oDICOM->GetValue('2001'x,'1018'x,/NO_COPY)
    
    sdim   = *slice_ptr[0]
    slices = *slice_ptr[0]
    
    ;;calculate s_fov
    ;; retrieve space between slices, convert to cm from mm
    temp = oDICOM->GetValue('0018'x,'0088'x,/NO_COPY)
    spacing = FLOAT(*temp[0]) / 10
    thick   = FLOAT(bparams.thick) / 10

    IF spacing EQ thick THEN s_fov = thick * slices $
    ELSE IF spacing LT thick THEN s_fov = thick + spacing * (slices-1) $
    ELSE s_fov = thick * slices + (spacing-thick) * (slices-1)
    
    adim = num_files / slices
    
    ;; if adim > 1, need to check for valid bvals, if found, then is diffusion image
    IF adim GT 1 THEN BEGIN
       
       progressbar = OBJ_NEW('progressbar', Color='red', $
                             Text='Checking for DICOM diffusion data from multiple files',/NOCANCEL)
       progressbar -> Start
       
       slice_idx = FLTARR(num_files)
       inst_idx  = FLTARR(num_files)
       
       bval  = FLTARR(num_files)
       RLcos = FLTARR(num_files)
       APcos = FLTARR(num_files)
       FHcos = FLTARR(num_files)
       
       FOR ii = 0, num_files-1 DO BEGIN
          
          ;; check is DICOM readable, else discard
          IF oDICOM->READ(FLIST[ii]) EQ 1 THEN BEGIN
             
             ;; slice index, USES PHILIPS PRIVATE CLASS
             temp = oDICOM->GetValue('2001'x,'100A'x,/NO_COPY)
             slice_idx[ii] = FLOAT(*temp[0])
             
             ;; instance number
             temp = oDICOM->GetValue('0020'x,'0013'x,/NO_COPY)
             inst_idx[ii] = FLOAT(*temp[1])
             
             ;;read in params for processing tensor
             temp = oDICOM->GetValue('2001'x,'1003'x,/NO_COPY)
             bval[ii] = *temp[0]
             
             temp = oDICOM->GetValue('2005'x,'10B0'x,/NO_COPY)
             RLcos[ii] = *temp[0]
             
             temp = oDICOM->GetValue('2005'x,'10B1'x,/NO_COPY)
             APcos[ii] = *temp[0]
             
             temp = oDICOM->GetValue('2005'x,'10B2'x,/NO_COPY)
             FHcos[ii] = *temp[0]
          ENDIF ELSE BEGIN
             void = DIALOG_MESSAGE('Error loading multi-file DICOM protocol.',/ERROR)
             RETURN
          ENDELSE
          
          progressBar -> Update, (float(ii)/float(num_files-1))*100.0
       ENDFOR
       
       progressbar -> Destroy
       
       unq_slice = slice_idx(UNIQ(slice_idx))
       
       IF ((SIZE(unq_slice))[1]) NE slices THEN BEGIN
          void = DIALOG_MESSAGE('Error loading multi-file DICOM protocol : ' $
                                + 'slice count mismatch.',/ERROR)
          RETURN
       ENDIF
       
       ;; process first slice only as this data will be sent to calculate tensor
       ;; if image is DTI, other image data loaded in DICOM_load_image
       
       curr_idx = WHERE((slice_idx EQ unq_slice[0]))
       
       ;;get instances
       sort_idx = SORT(inst_idx(curr_idx))
       
       curr_idx = curr_idx(sort_idx)
       
       bval = bval(curr_idx)
       
       valid_bval = 0
       FOR ii = 0, N_ELEMENTS(bval)-1 DO BEGIN
          IF bval[ii] NE 0 THEN BEGIN
             valid_bval = 1
             BREAK
          ENDIF
       ENDFOR
       
       IF valid_bval EQ 1 THEN diffimage = 1 $
       ELSE diffimage = 0
       
       IF diffimage EQ 1 THEN BEGIN
        
          diffdata = DICOM_calc_DTI_tensor(adim,bval, $
                                         RLcos(curr_idx), $
                                         APcos(curr_idx), $
                                         FHcos(curr_idx), $
                                         img_orient=oparams.orientation[0])
          acq_matrix = *diffdata.acq_matrix
          
       ENDIF
    ENDIF ELSE BEGIN
       diffimage = 0
       diffdata = {nothing:1}
       acq_matrix = diag_matrix([1.0,1,1,1])
    ENDELSE
    
    ;;the file path was passed in
    display_Name = oparams.scan_name + ' : ' + oparams.sequence_name + $
                   ' : ' + oparams.scan_date
    
    ;;the images location
    file_Path = STRMID(location,0,STRPOS(location,(get_dir_slash())[0],/REVERSE_SEARCH)+1)
    
    ;; on philips data, number of averages is same for all images in a scan
    pn_avg = REPLICATE(oparams.n_avg,sdim)
    
    dicom = { image_Number:oparams.image_Number $
              , fdim:0 $
              , pdim:0 $
              , init_fdim:bparams.fdim $
              , init_pdim:bparams.pdim $
              , sdim:sdim $
              , adim:adim $
              , sequence_name:oparams.sequence_name $
              , acq_matrix: ptr_new(acq_matrix) $
              , scan_date:oparams.scan_date $
              , scan_name:oparams.scan_name $
              , recov_time:oparams.recov_time $
              , echo_time:oparams.echo_time $
              , n_echo:oparams.n_echo $
              , n_avg:oparams.n_avg $
              , pn_avg:pn_avg $
              , f_fov:oparams.f_fov $
              , p_fov:oparams.p_fov $
              , s_fov:s_fov $
              , slices:slices $
              , thick:bparams.thick $
              , f_voxsz: bparams.voxel_dims[0]  $
              , p_voxsz: bparams.voxel_dims[1]  $
              , s_voxsz: bparams.voxel_dims[2]  $
              , Imin: 4  $
              , Jmin: 2 $
              , Kmin: 0 $
              , dimensions:oparams.dimensions $
              , file_Path:file_Path $
              , display_Name:display_Name $
              , orientation:oparams.orientation $
              , diffimage:diffimage $
              , diffdata:diffdata $
              , DICOMtype:2 $
              , DICOMfiles:PTR_NEW(FLIST) $
              , chkval:1 $
              , fat_sup:bparams.fat_sup $
              , PARRECtype:0 $
              , Philips_Data:1 $
            }
    
    add_Philips_to_MAS, dicom
    
END


; Subroutine name: load_DICOM_multiframe
; Created by: CD 2/15/07
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:
; when we get a diffusion multi-frame, check that read in order is adim then sdim

FUNCTION load_DICOM_multiframe, bparams, oDICOM

	frames = *bparams.frame_ptr[0]

	; if slice_ptr is defined, then image has an adim
	IF PTR_VALID(bparams.slice_ptr) EQ 1 THEN BEGIN
		sdim = *bparams.slice_ptr[0]
		slices = *bparams.slice_ptr[0]
		adim = frames / slices
	ENDIF ELSE BEGIN
		sdim = frames
		slices = frames
		adim = 1
	ENDELSE

	; if adim > 1, check to see if valid bvals, if so then diffusion image
	IF adim GT 1 THEN BEGIN

		;get the image type, check if bvals exist
		bval_ptr = oDICOM->GetValue('2001'x,'1003'x,/NO_COPY)

		IF PTR_VALID(bval_ptr[0]) EQ 1 THEN BEGIN

			valid_bval = 0
			FOR ii = 0, N_ELEMENTS(bval_ptr)-1 DO BEGIN
				IF *bval_ptr[ii] NE 0 THEN BEGIN
					valid_bval = 1
					BREAK
				ENDIF
			ENDFOR

			IF valid_bval EQ 1 THEN diffimage = 1 $
			ELSE diffimage = 0
		ENDIF ELSE diffimage = 0
	ENDIF ELSE diffimage = 0

	;calculate s_fov
	; retrieve space between slices, convert to cm from mm
        temp = oDICOM->GetValue('0018'x,'0088'x,/NO_COPY)
        spacing = FLOAT(*temp[0]) / 10
        thick   = bparams.thick / 10

	IF spacing EQ thick THEN s_fov = thick * slices $
	ELSE IF spacing LT thick THEN s_fov = thick + spacing * (slices-1) $
	ELSE s_fov = thick * slices + (spacing-thick) * (slices-1)

	;if diffusion image, calculate tensor
	IF diffimage EQ 1 THEN BEGIN

		bval  = FLTARR(frames)
		RLcos = FLTARR(frames)
		APcos = FLTARR(frames)
		FHcos = FLTARR(frames)

		temp_bv = oDICOM->GetValue('2001'x,'1003'x,/NO_COPY)
		temp_rl = oDICOM->GetValue('2005'x,'10B0'x,/NO_COPY)
		temp_ap = oDICOM->GetValue('2005'x,'10B1'x,/NO_COPY)
		temp_fh = oDICOM->GetValue('2005'x,'10B2'x,/NO_COPY)

		FOR jj = 0, adim-1 DO BEGIN
			; extract data

			bval[jj]  = *temp_bv[jj]
			RLcos[jj] = *temp_rl[jj] ; X
			APcos[jj] = *temp_ap[jj] ; Y
			FHcos[jj] = *temp_fh[jj] ; Z
		ENDFOR

		diffdata = DICOM_calc_DTI_tensor(adim,bval,RLcos,APcos,FHcos)

	ENDIF ELSE diffdata = {nothing:1}

	data = {sdim:sdim, slices:slices, adim:adim, s_fov:s_fov, $
			diffdata:diffdata, diffimage:diffimage}

	RETURN, data
END


; Subroutine name: load_DICOM_single_image
; Created by: CD 2/15/07
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

FUNCTION load_DICOM_single_image, oDICOM, bparams

	;assume image is a single image DICOM file
	adim = 1
	sdim = 1
	slices = 1
	s_fov = FLOAT(bparams.thick) / 10

	data = {sdim:sdim, slices:slices, adim:adim, s_fov:s_fov, $
			diffdata:{nothing:1}, diffimage:0}

	RETURN, data
END


; Subroutine name: get_dicom_header_and_image
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; This fuction will extract as much imnd like info from the
; DICOM header as possible and return a structure, described below,
; that represents the scan pointed to by the location variable
;
; Editing Information:
    ;Edited by HS 2006/10/05.
    ;Fix spelling mistakes and commenting
    ;Edited by CD 2007/01/10
    ;Add additional DICOM parameters

FUNCTION get_dicom_header_and_image, location

	; dicom.chkval = 1 indicates successful file read, -1 indicates error

	oDICOM = OBJ_NEW( 'IDLffDICOM' )

	; Open the file
    IF oDICOM->READ(location) NE 1 THEN BEGIN
		void = DIALOG_MESSAGE('This file is not in a supported DICOM format.',/ERROR)
        RETURN, {chkval:-1}
    ENDIF

	bparams = get_basic_dicom_params(oDICOM)

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; DETERMINE NUMBER OF IMAGES IN FILE, PROCESS ACCORDINGLY ;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;check the number of frames, if defined, process as multiframe DICOM file

	IF PTR_VALID(bparams.frame_ptr) EQ 1 THEN BEGIN

		fparams = load_DICOM_multiframe(bparams,oDICOM)
		DICOMtype = 3

	ENDIF ELSE BEGIN

		;if defined image is part of a multi file image set
		;assumes all files are from same imaging exam card
		IF PTR_VALID(bparams.slice_ptr) EQ 1 THEN BEGIN

			; display exam card selection
			exam_card_GUI_create, location, oDICOM

			; exit routine after loading GUI
			; GUI events will handle protocol loading
			RETURN, {chkval:-1}

		ENDIF ELSE BEGIN

			fparams = load_DICOM_single_image(oDICOM,bparams)
			DICOMtype = 1

		ENDELSE
	END

	oparams = get_other_dicom_params(oDICOM,bparams)

    ;the file path was passed in
    file_Path = STRMID(location,0,STRPOS(location,'\',/REVERSE_SEARCH)+1)


    ;the display name is usually the images location
    display_Name = location

	pn_avg = INTARR(fparams.sdim)
	pn_avg[*] = oparams.n_avg

    dicom = { image_Number:oparams.image_Number $
              , fdim:0 $
              , pdim:0 $
              , init_fdim:bparams.fdim $
              , init_pdim:bparams.pdim $
              , sdim:fparams.sdim $
              , adim:fparams.adim $
              , sequence_name:oparams.sequence_name $
              , scan_date:oparams.scan_date $
              , scan_name:oparams.scan_name $
              , recov_time:oparams.recov_time $
              , echo_time:oparams.echo_time $
              , n_echo:oparams.n_echo $
              , n_avg:oparams.n_avg $
              , pn_avg:pn_avg $
              , f_fov:oparams.f_fov $
              , p_fov:oparams.p_fov $
              , s_fov:fparams.s_fov $
              , slices:fparams.slices $
              , thick:bparams.thick $
              , f_voxsz: bparams.voxel_dims[0]  $
              , p_voxsz: bparams.voxel_dims[1]  $
              , s_voxsz: bparams.voxel_dims[2]  $
              , dimensions:oparams.dimensions $
              , file_Path:file_Path $
              , display_Name:display_Name $
              , orientation:oparams.orientation $
              , diffimage:fparams.diffimage $
              , diffdata:fparams.diffdata $
              , DICOMfiles:PTR_NEW(location) $
              , DICOMtype:DICOMtype $
              , chkval:1 $
              , fat_sup:bparams.fat_sup $
              , PARRECtype:0 $
              , Philips_Data:1 $
            }

    OBJ_DESTROY, oDICOM
    HEAP_GC

    RETURN, dicom
END


; Subroutine name: mas_open_multi_dicom
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
; CD 2007/02/15. OLD VERSION OF CODE in mas_open
;
; Purpose of subroutine:
; HS 2006/09/29. Opens DICOM files.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting
    ;Edited by CD 2007/01/10
    ;Add multi-frame generic dicom reader

PRO mas_open_multi_dicom

    COMMON scan_data

	ni = project.ni

    junk = dialog_message(["This function is no longer in use.", $
                           "Use the 'Generic DICOM Importer' instead."], $ 
                           /error, /center)
    return
     
;    IF ni EQ 50 THEN BEGIN
;       update_status_bar, ' Reached max limit of scans able to open. 50 '
;       RETURN
;    END
;
;    dirPath = DIALOG_PICKFILE( Title = 'Select scan file to open' $
;    	,path=project.current_Path)
;
;    IF dirPath EQ '' THEN RETURN
;
;	;check for .dcm file extension, if found, process according to old dicom
;	;processing, if not found, assume is a generic multi-frame dicom file
;	;and process accordingly
;
;	dcmfound = STRPOS(dirpath,'.dcm')
;
;   	IF (dcmfound NE -1) THEN BEGIN
;   		return
;	ENDIF ELSE BEGIN
;	    dicom = get_dicom_header_and_image(dirPath)
;
;		IF dicom.chkval EQ -1 THEN RETURN
;
;		add_Philips_to_MAS, dicom
;
;	ENDELSE

END
