;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

; Subroutine name: erosion_dilation
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This function will return the mask eroded and then dilated with the radious specified.

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

FUNCTION erosion_dilation, p_mask, radius

    dims = size((*p_mask),/dimension)

    ;IF the image is made up of all one pixel value then return a white image.
;    histo = histogram(mask)
;    IF histo[1] eq (dims[0]*dims[1]) or histo[0] eq (dims[0]*dims[1]) THEN $
;       RETURN, REPLICATE(1B, dims[0], dims[1] )

    strucElem = shift(DIST(2*radius+1), radius, radius) LE radius
    
    eros_dila = replicate(0B, dims[0]+radius, dims[1]+radius)
    eros_dila [radius/2,radius/2]  = (*p_mask)
    eros_dila = ERODE(eros_dila, strucElem )
    eros_dila = DILATE(eros_dila, strucElem )
    
    RETURN, ptr_new(eros_dila[radius/2:dims[0],radius/2:dims[1]])
END


; Subroutine name: auto_select_region
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This function will select the proper region based on the maks and image provided
; returning the region within the largest varience from image statistics.

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting


FUNCTION auto_select_region, p_image, p_mask

;    CATCH, Error_status
;    IF Error_status NE 0 THEN BEGIN
;
;
;       IF Error_status eq -161 THEN BEGIN
;         CATCH, /CANCEL
;         PRINT, 'Region selection error in AUTO_SELECT_REGION.'
;         PRINT, 'Setting the whole mask to 1 and returning.'
;         RETURN, PTR_NEW(REPLICATE(1B, dims[0], dims[1] ))
;
;      END
;
;      PRINT, 'inside error message of AUTO_SELECT_REGION'
;       PRINT, 'Error index: ', Error_status
;       PRINT, !ERROR_STATE.msg
;      CATCH, /CANCEL
;    END


    dims = size((*p_image),/dimension)

    regions = LABEL_REGION((*p_mask))

;   regions[0,0], '     *0', $
;     regions[dims[0]-1,0], '     0*', $
;     regions[0,dims[1]-1], '     **', $
;     regions[dims[0]-1,dims[1]-1]



    ;how about this idead the regions that are in the corners dont use them.

    regions[where(regions EQ 0 )] = 0b

;    hist = histogram(regions)
;
;
;    sort_hist = hist[reverse(sort(hist))]
;
;    ;now find the maxium pixel population
;
;    max_hist = max(hist, largest_region)
;
;
;
;    ;calculate image statistics for the first to regions
;    ;then choose which ever has the larger sigma variance.
;
;    mask2 = BYTARR(dims[0], dims[1], /NOZERO)
;    mask2[*,*] =1b
;    mask2[where(regions ne (reverse(sort(hist)))[0])] = 0b
;
;    IMAGE_STATISTICS, (*p_image), mask=mask2, VARIANCE=region0Var
;
;
;
;    mask2 = BYTARR(dims[0], dims[1], /NOZERO)
;    mask2[*,*] =1b
;    mask2[where(regions ne (reverse(sort(hist)))[1])] = 0b
;
;    IMAGE_STATISTICS, (*p_image), mask=mask2, VARIANCE=region1Var
;
;    ;use either the the first or the second largest region with the highest variance.
;    IF region0Var gt region1Var then useRegion = 0 else useRegion = 1
;
;    regions[where(regions ne (reverse(sort(hist)))[useRegion])] = 0b
    p_region = ptr_new(byte(regions))

    RETURN, p_region
END


; Subroutine name: dilation_erosion
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This function will dilate a mask then erode.

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting


FUNCTION dilation_erosion, p_mask, radius

    dims = size((*p_mask),/dimension)
;    histo = histogram(mask)
;    ;IF the image is made up of all one pixel value then return a white image.
;    IF histo[1] eq (dims[0]*dims[1]) or histo[0] eq (dims[0]*dims[1]) THEN $
;       RETURN, REPLICATE(1B, dims[0], dims[1] )

    strucElem = shIFt(DIST(2*radius+1), radius, radius) LE radius

    eros_dila = replicate(0B, dims[0]+radius, dims[1]+radius)
    eros_dila [radius/2,radius/2]  = (*p_mask)
    eros_dila = DILATE(eros_dila, strucElem )
    eros_dila = ERODE(eros_dila, strucElem )

    RETURN, ptr_new(eros_dila[radius/2:dims[0],radius/2:dims[1]])
END


; Subroutine name: background_threshold_mask
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This function will generate a mask based on the image passed in and the parameters for
; threshold. IF there is a problem the don't exclude any data.

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting


function background_threshold_mask, p_image, threshold

    lo_mask = (*p_image) GE threshold[0]
    hi_mask = (*p_image) LE threshold[1]

    return_mask = lo_mask AND hi_mask

    dims = size((*p_image),/dimension)
    ;IF the image is made up of all one pixel value then return a white image.
    IF (histogram(return_mask))[1] eq (dims[0]*dims[1]) THEN $
       RETURN, PTR_NEW(REPLICATE(1B, dims[0], dims[1] ))

    return, ptr_new(return_mask, /NO_COPY)
END


; Subroutine name: update_threshold_display
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This procedure will update the threshold display image data.

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

pro update_threshold_display
    COMMON scan_data
    CI = project.ci
    
    key_adim = project.procPramArray[ci].key_adim
    key_slice = project.procPramArray[ci].key_slice
    threshold = project.procPramArray[ci].threshold
    
    p_image = ptr_new((*project.dataArray[ci].state1)[*,*,key_slice,key_adim])
    
    project.procPramArray[ci].remove_background_thres_image = $
      background_threshold_mask( p_image , threshold )

END

; Subroutine name: generate_background_mask
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This procedure will generate the background filter removial and store it at the
; end of a pointer.

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

FUNCTION generate_background_mask, p_image
    COMPILE_OPT IDL2

    ;bring in the appropriate global variables
    COMMON scan_data
    CI = project.ci


    dims = size(*p_image, /dimension)
    ;PRINT, 'Dimension of data', dims

    ;make a mask array the same size as the state1 array but only the first 3 dimesions.
    masks = bytarr(dims[0], dims[1])

    ;create the mask desired.


    p_mask = background_threshold_mask( p_image , project.procPramArray[ci].threshold )


    p_mask = erosion_dilation( p_mask, 2)

    p_mask = AUTO_SELECT_REGION( p_image , p_mask)

    p_mask = dilation_erosion( p_mask, 2)

    return,p_mask

END

function mas_background_filter_getmask, pdata

  common scan_data
  ci = project.ci

  dims = size(*pdata, /dimensions)
  data_sz = size(*pdata)
  key_adim  =  project.procPramArray[ci].key_adim
  mask = bytarr(dims)
  work = ptr_new(fltarr(dims[0], dims[1]))

  if (data_sz[0] eq 4) then begin
     
     for sl = 0, dims[2]-1 do begin
        
        for ar = 0, dims[3]-1 do begin
           *work = reform((*pdata)[*,*,sl,key_adim])
           pmask = generate_background_mask(work)
           mask[*,*,sl,ar] = *pmask
           ptr_free, pmask
        endfor

     endfor

     ptr_free, work

     return, ptr_new(mask, /no_copy)
     
  endif else if (data_sz[0] eq 3) then begin
     
     for sl = 0, dims[2]-1 do begin
        
        *work = reform((*pdata)[*,*,sl])
        pmask = generate_background_mask(work)
        mask[*,*,sl] = *pmask
        ptr_free, pmask
        
     endfor
     
     ptr_free, work
     
     return, ptr_new(mask, /no_copy)

  endif else begin

     return, ptr_new()

  endelse

end

; Subroutine name: mas_background_filter_load
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:

pro mas_background_filter_load, junk

    common scan_data
    forward_function mas_read_nifti
    forward_function mas_interpolate
    
    mas_load_state_1
    
    adim = project.imndarray[project.ci].adim
    
    ;file = dialog_pickfile(path=project.current_path, /read, /file, filter=['*.nii'], $
    ;                         title="Choose a NIFTI mask file:")
    ;if (file eq '' or not file_test(file, /read)) then return
    
    nif = mas_read_nifti(nifti_filename=file, read_status=rs)

    if (rs eq 1) then begin

        nif.voxel_data = mas_interpolate(nif.voxel_data)
        vox_ge = where(*nif.voxel_data ge 0.5, n_ge, complement=vox_lt, ncomplement=n_lt)
        if (n_ge gt 0) then begin
            (*nif.voxel_data)[vox_ge] = 1
        endif
        if (n_lt gt 0) then begin
            (*nif.voxel_data)[vox_lt] = 0
        endif
        mask_dims = size(*nif.voxel_data, /dimensions)
        n_mask_dims = n_elements(mask_dims)
        data_dims = size(*project.dataarray[project.ci].state1, /dimensions)
        n_data_dims = n_elements(data_dims)
        
        if (n_data_dims eq 3) then begin
            ;; mask_dims must equal data_dims
            if (n_mask_dims ne 3) then begin
                void = dialog_message('Data Dimensions does not equal Mask Dimensions.', /error, /center)
                ptr_free, nif.voxel_data
                return
            endif else begin
                dim_diff = total(abs(data_dims - mask_dims))
                if (dim_diff ne 0) then begin
                    void = dialog_message('Data Dimensions does not equal Mask Dimensions.', /error, /center)
                    ptr_free, nif.voxel_data
                    return
                endif
            endelse
         endif else if (n_data_dims eq 4) then begin
            ;; mask_dims[0:2] must equal data_dims
            if (n_mask_dims ne 3) then begin
                void = dialog_message('Data Dimensions does not equal Mask Dimensions.', /error, /center)
                ptr_free, nif.voxel_data
                return
            endif else begin
                dim_diff = total(abs(data_dims[0:2] - mask_dims))
                if (dim_diff ne 0) then begin
                    void = dialog_message('Data Dimensions does not equal Mask Dimensions.', /error, /center)
                    ptr_free, nif.voxel_data                    
                    return
                endif
            endelse
        endif
        
        pbar = obj_new('progressbar', title='Applying Filter', text='Applying Background Mask', /fast_loop, /nocancel)
        pbar->Start
        for a = 0, adim-1 do begin
        
            (*project.dataarray[project.ci].state1)[*,*,*,a] *= *nif.voxel_data
            
            if (a mod 10 eq 0) then pbar->Update, float(a)/adim * 100.0
            
        endfor
        pbar->Destroy
        
        project.procpramarray[project.ci].state_2 = 0
        ptr_free, nif.voxel_data 
        
    endif else begin
    
        void = dialog_message('NIFTI file could not be read.', /error, /center)
        
    endelse

end


; Subroutine name: mas_remove_background_filter
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Alter p_data acording the the mask data

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

pro mas_remove_background_filter, p_data
    COMPILE_OPT IDL2
    HEAP_GC
    ;bring in the appropriate global variables
    COMMON scan_data
    CI = project.ci

    ;IF the remove background active flag is not set and
    ;the pointer to the masks is not valid then skip this step
    IF project.procPramArray[ci].remove_background_active EQ 0 THEN RETURN

    slice_axis =  project.procPramArray[ci].slice_axis
    sort_dir    =  project.procPramArray[ci].sort_dir
    key_adim  =  project.procPramArray[ci].key_adim
    fdim_start    =  project.procPramArray[ci].fdim_start
    pdim_start  =  project.procPramArray[ci].pdim_start
    sdim_start  =  project.procPramArray[ci].sdim_start

    dims = size(*p_data, /dimension)
    print, dims

    progressbar = Obj_New('progressbar', Color='red', Text='Background Filter',/NOCANCEL)
    progressbar -> Start

    CASE (size(*p_data))[0] OF
       2:BEGIN
         if slice_axis eq 0 then begin
                ;slice along the read_phase

            p_image = ptr_new((*project.dataArray[CI].state1)[*,*,sdim_start,key_adim])

          (*p_data)[*,*] = (*p_data)[*,*]  * $
              FLOAT(*(generate_background_mask( p_image )))

            end else if slice_axis eq 1 then begin
                ;slice along the read_slice

            p_image = PTR_NEW(REFORM((*project.dataArray[CI].state1)[*,pdim_start,*,key_adim]))


          (*p_data)[*,*] = (*p_data)[*,*]  * $
              FLOAT(*(generate_background_mask( p_image )))

            end else if slice_axis eq 2 then begin
                ;slice along the phase_slice
          p_image = PTR_NEW(REFORM((*project.dataArray[CI].state1)[fdim_start,*,*,key_adim]))


          (*p_data)[*,*] = (*p_data)[*,*]  * $
              FLOAT(*(generate_background_mask( p_image )))
         end

       END

       3:BEGIN

        ;which dimension do they want to traverse
            if sort_dir eq 0 then begin
                ;sort the data by the sdim
                if slice_axis eq 0 then begin
                    ;slice along the read_phase
                  FOR ii=0, dims[2]-1 DO BEGIN
                    p_image = ptr_new((*project.dataArray[CI].state1)[*,*,ii,key_adim])

                 (*p_data)[*,*,ii] = (*p_data)[*,*,ii]  * $
                   FLOAT(*(generate_background_mask( p_image )))


              ENDFOR

                end else if slice_axis eq 1 then begin
                    ;slice along the read_slice
                    print,'creting a mask for read_slice'
              FOR ii=0, dims[2]-1 DO BEGIN

                    p_image = PTR_NEW(REFORM((*project.dataArray[CI].state1)[*,ii,*,key_adim]))

                 mask = FLOAT(*(generate_background_mask( p_image )))

                 (*p_data)[*,*,ii] = (*p_data)[*,*,ii]  * mask

              ENDFOR

                end else if slice_axis eq 2 then begin
                    ;slice along the phase_slice
              FOR ii=0, dims[2]-1 DO BEGIN

                    p_image = PTR_NEW(REFORM((*project.dataArray[CI].state1)[ii,*,*,key_adim]))

                 mask = FLOAT(*(generate_background_mask( p_image )))

                 (*p_data)[*,*,ii] = (*p_data)[*,*,ii]  * mask

              ENDFOR
                end


            end else begin
                ;sort the data by the adim
          ;for traversing the adim the mask stays the same for each slice.

                if slice_axis eq 0 then begin
                    ;slice along the read_phase

                 p_image = ptr_new((*project.dataArray[CI].state1)[*,*,sdim_start,key_adim])
                 mask = FLOAT(*(generate_background_mask( p_image )))

              FOR ii=0, dims[2]-1 DO BEGIN

                 (*p_data)[*,*,ii] = (*p_data)[*,*,ii]  * mask

                   progressBar -> Update, (float(ii)/float(dims[2]-1))*100.0

              ENDFOR

                end else if slice_axis eq 1 then begin
                    ;slice along the read_slice
              p_image = ptr_new(REFORM((*project.dataArray[CI].state1)[*,pdim_start,*,key_adim]))
                 mask = FLOAT(*(generate_background_mask( p_image )))

              FOR ii=0, dims[2]-1 DO BEGIN

                 (*p_data)[*,*,ii] = (*p_data)[*,*,ii]  * mask

              ENDFOR

                end else if slice_axis eq 2 then begin
                    ;slice along the phase_slice
              p_image = ptr_new(REFORM((*project.dataArray[CI].state1)[fdim_start,*,*,key_adim]))
                 mask = FLOAT(*(generate_background_mask( p_image )))

              FOR ii=0, dims[2]-1 DO BEGIN

                 (*p_data)[*,*,ii] = (*p_data)[*,*,ii]  * mask

              ENDFOR
                end

            end


       END
       4:BEGIN
         PRINT,'Processing for 4d data'

         FOR ss=0, dims[2]-1 DO BEGIN
          FOR aa=0, dims[3]-1 DO $
              (*p_data)[*,*,ss,aa] = (*p_data)[*,*,ss,aa]  * $
                 FLOAT((*project.procPramArray[ci].remove_background_masks)[*,*,ss])

         END
       END

       ELSE: RETURN

    ENDCASE



    progressbar -> Destroy

END


; Subroutine name: mas_background_histo
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This procedure will make a histogram for each image in the slice range and the
; adim range.

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

pro mas_background_histo, slice_range, adim_range
    COMPILE_OPT IDL2
    HEAP_GC
    ;bring in the appropriate global variables
    COMMON scan_data
    CI = project.ci

    mas_load_state_1

    ;IF the range to display is across slice or adim then make a square
    ;view grid
    ;this is IF they want one slice then view grid is 1 by 1
    IF slice_range[0] EQ slice_range[1] AND adim_range[0] EQ adim_range[1] THEN BEGIN
       ;PRINT, 'single slice'
       VIEW_GRID = [1,1]

    ENDIF ELSE IF adim_range[0] EQ adim_range[1] THEN BEGIN
       ;PRINT, 'across slices'
       side_length = CEIL(SQRT(slice_range[1]-slice_range[0]+1))
       VIEW_GRID = [side_length, side_length]

    ENDIF ELSE IF slice_range[0] EQ slice_range[1] THEN BEGIN
       ;PRINT, 'across adim'
       side_length = CEIL(SQRT(adim_range[1]- adim_range[0]+1))
       VIEW_GRID = [side_length, side_length]

    ENDIF ELSE BEGIN
       ;PRINT, 'aross all slice range and adim'
       VIEW_GRID = [slice_range[1]-slice_range[0]+1,adim_range[1]- adim_range[0]+1 ]
    END

    progressbar = Obj_New('progressbar', Color='red', Text='Creating Histograph Display',/NOCANCEL)
    progressbar -> Start

    FOR aa = adim_range[0], adim_range[1] DO BEGIN

       progressBar -> Update, (float(aa)/float(adim_range[1]))*100.0

        FOR ss = slice_range[0],slice_range[1]  DO BEGIN

           name = 'S='+strtrim(ss+1,2)+' A='+strtrim(aa+1,2)
           ;PRINT, name

           ;IF this is the first plot then set up the plot window
           IF AA EQ adim_range[0] AND SS EQ slice_range[0] THEN $
             iplot, histogram((*project.dataArray[ci].state1)[*,*,ss,aa]), $
              /HISTOGRAM ,VIEW_GRID=VIEW_GRID, XTITLE=name $
           ELSE $
             iplot, histogram((*project.dataArray[ci].state1)[*,*,ss,aa]), $
             /HISTOGRAM, VIEW_GRID=VIEW_GRID, XTITLE=name, /VIEW_NEXT
        ENDFOR
    ENDFOR

    progressbar -> Destroy

END


; Subroutine name: mas_background_activate
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Turn of or off the sensitivity of the widgets based on arguments passed in.

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

pro mas_background_activate, sensitive, wWidget
    COMMON scan_data
    CI = project.ci


    WIDGET_control, Widget_Info(wWidget, FIND_BY_UNAME='BACKGROUND_THRESH_LOWER_BOUND'), SENSITIVE=sensitive

    WIDGET_control, Widget_Info(wWidget, FIND_BY_UNAME='BACKGROUND_THRESH_UPPER_BOUND'), SENSITIVE=sensitive



    IF project.imndArray[CI].sdim EQ 1 THEN BEGIN
       WIDGET_control, Widget_Info(wWidget, FIND_BY_UNAME='BACKGROUND_KEY_SLICE'), SENSITIVE=0


    ENDIF ELSE BEGIN

       WIDGET_control, Widget_Info(wWidget, FIND_BY_UNAME='BACKGROUND_KEY_SLICE'), SENSITIVE=sensitive

    END

    IF project.imndArray[CI].adim EQ 1 THEN $
       WIDGET_control, Widget_Info(wWidget, FIND_BY_UNAME='BACKGROUND_KEY_ADIM'), SENSITIVE=0 $
    ELSE $
       WIDGET_control, Widget_Info(wWidget, FIND_BY_UNAME='BACKGROUND_KEY_ADIM'), SENSITIVE=sensitive

    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='BACKGROUND_HISTO_CURRENT'), SENSITIVE=sensitive

    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='BACKGROUND_HISTO_SLICE'), SENSITIVE=sensitive

    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='BACKGROUND_HISTO_ARRAY'), SENSITIVE=sensitive

    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='BACKGROUND_HISTO_ALL'), SENSITIVE=sensitive

    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='BACKGROUND_PREVIEW'), SENSITIVE=sensitive

END


; Subroutine name: mas_background_filter_redraw
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This procedure will refresh the values in the background filter tool from memory.

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

pro mas_background_filter_redraw

    COMPILE_OPT IDL2
    ;;bring in the appropriate global variables
    COMMON scan_data
    COMMON common_widgets
    CI = project.ci
    
    IF N_ELEMENTS(REMOVE_BACKGROUND_BASE) NE 1 THEN RETURN

    IF WIDGET_INFO( REMOVE_BACKGROUND_BASE,  /VALID_ID ) THEN BEGIN
        
        wWidget =  REMOVE_BACKGROUND_BASE
        
        WIDGET_CONTROL, WIDGET_INFO(wWidget, FIND_BY_UNAME='BACKGROUND_ACTIVATE'), $
          SET_BUTTON=project.procPramArray[ci].remove_background_active
        
        IF project.procPramArray[ci].remove_background_active EQ 0 THEN BEGIN
            
            mas_background_activate, 0, wWidget
            
            WIDGET_CONTROL, WIDGET_INFO(wWidget, FIND_BY_UNAME='BACKGROUND_THRESH_DRAW'), $
              GET_VALUE=drawID, SENSITIVE = 0, $
              DRAW_XSIZE = 1, $
              DRAW_YSIZE = 1
            WSET , drawID
            
            tv, [128]
            
            WIDGET_CONTROL, WIDGET_INFO(wWidget, FIND_BY_UNAME='BACKGROUND_EROSION_DILATION_DRAW'), $
              GET_VALUE=drawID, SENSITIVE = 0, $
              DRAW_XSIZE = 1, $
              DRAW_YSIZE = 1
            WSET , drawID
            
            tv, [128]
            
        ENDIF ELSE BEGIN
            
            mas_background_activate, 1, wWidget
            key_adim  = project.procPramArray[ci].key_adim
            key_slice = project.procPramArray[ci].key_slice
            threshold = project.procPramArray[ci].threshold
            
            ;;set the min labels to look right.
            slice_range_min = min((*project.dataArray[CI].state1), max = slice_range_max )
            
            sz_data = size((*project.dataArray[ci].state1))
            
            WIDGET_CONTROL, WIDGET_INFO(wWidget, FIND_BY_UNAME='BACKGROUND_THRESH_LOWER_BOUND'), $
              set_value = [threshold[0], slice_range_min, threshold[1]]
            
            WIDGET_CONTROL, WIDGET_INFO(wWidget, FIND_BY_UNAME='BACKGROUND_THRESH_UPPER_BOUND'), $
              set_value = [threshold[1], threshold[0], slice_range_max]
            
            
            WIDGET_CONTROL, WIDGET_INFO(wWidget, FIND_BY_UNAME='BACKGROUND_KEY_SLICE'), $
              set_value = key_slice, $
              SET_SLIDER_MAX=sz_data[3]-1
            
            
            WIDGET_CONTROL, WIDGET_INFO(wWidget, FIND_BY_UNAME='BACKGROUND_KEY_ADIM'), $
              set_value = key_adim, $
              SET_SLIDER_MAX=project.imndArray[CI].adim-1
            
            IF PTR_VALID(project.procPramArray[ci].remove_background_thres_image) THEN BEGIN
                ;;is there an image to draw to the draw widget.
                
                
                WIDGET_CONTROL, WIDGET_INFO(wWidget, FIND_BY_UNAME='BACKGROUND_THRESH_DRAW'),$
                  GET_VALUE=drawID
                
                ;old_geom = widget_info(drawID, /geometry)
                ;if (old_geom.draw_xsize ne sz_data[1] and old_geom.draw_ysize ne sz_data[2]) then begin
                ;    widget_control, drawID, DRAW_XSIZE = sz_data[1], DRAW_YSIZE = sz_data[2]
                ;endif
                
                WSET , drawID
                
                TV, BYTSCL( *(project.procPramArray[ci].remove_background_thres_image),max=1, min=0)
                
                WIDGET_CONTROL, WIDGET_INFO(wWidget, FIND_BY_UNAME='BACKGROUND_EROSION_DILATION_DRAW'),$
                  GET_VALUE=drawID

                ;old_geom = widget_info(drawID, /geometry)
                ;if (old_geom.draw_xsize ne sz_data[1] and old_geom.draw_ysize ne sz_data[2]) then begin
                ;    print, 'resize'
                ;    widget_control, drawID, DRAW_XSIZE = sz_data[1], DRAW_YSIZE = sz_data[2]
                ;endif
                
                WSET , drawID
                
                p_image = ptr_new( (*project.dataArray[ci].state1) [*,*, key_slice, key_adim] )
                
                p_mask = erosion_dilation( project.procPramArray[ci].remove_background_thres_image, 3)
                
                
                p_mask = AUTO_SELECT_REGION( p_image, p_mask)
                
                
                p_mask = dilation_erosion( p_mask, 3)
                
                
                tv, BYTSCL( *p_mask , max=1, min=0)
                
            END ELSE BEGIN
                
                WIDGET_CONTROL, WIDGET_INFO(wWidget, FIND_BY_UNAME='BACKGROUND_THRESH_DRAW'), GET_VALUE=drawID
                
                
            END
            
            
        ENDELSE
        
    ENDIF
    
END
    
    
; Subroutine name: mas_background_filter_event
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting


pro mas_background_filter_event, event

    ;error handling.
    catch, error_status
    if (error_status ne 0) then begin
       help, calls= trace_back

       dummy = DIALOG_MESSAGE(['Please send a copy of this error message to the developer of this program.', $
         !ERROR_STATE.msg,trace_back], /ERROR, $
           TITLE='Error in WID_BASE_MAIN_event')
       RETURN
    endif

    COMPILE_OPT IDL2
    HEAP_GC
    ;bring in the appropriate global variables
    COMMON scan_data
    CI = project.ci


    wWidget =  Event.top
    name = widget_info(event.id, /uname)

    CASE name of
        
        'BACKGROUND_ACTIVATE': BEGIN
            project.procPramArray[ci].remove_background_active = event.select
            project.procPramArray[ci].state_2 = 0
            mas_load_state_2
        END
        
        'BACKGROUND_THRESH_LOWER_BOUND': BEGIN
            WIDGET_CONTROL, event.id, get_value = sTemp
            project.procPramArray[ci].threshold[0] = float(sTemp)
            project.procPramArray[ci].state_2 = 0
            update_threshold_display
        END
        
        'BACKGROUND_THRESH_UPPER_BOUND': BEGIN
            WIDGET_CONTROL, event.id, get_value = sTemp
            project.procPramArray[ci].threshold[1] = float(sTemp)
            project.procPramArray[ci].state_2 = 0
            update_threshold_display
        END
        
        'BACKGROUND_KEY_SLICE': BEGIN
            project.procPramArray[ci].key_slice = event.value > 0
            update_threshold_display
        END
        
        'BACKGROUND_KEY_ADIM': BEGIN
            project.procPramArray[ci].key_adim= event.value > 0
            project.procPramArray[ci].state_2 = 0
            update_threshold_display
        END
        
        'BACKGROUND_HISTO_CURRENT': BEGIN
            key_slice = project.procPramArray[ci].key_slice
            key_adim = project.procPramArray[ci].key_adim
            mas_background_histo, [key_slice,key_slice], [key_adim,key_adim]
        END
        
        'BACKGROUND_HISTO_SLICE': BEGIN
            ;;PRINT, 'create histogram for all slices in the slice range'
            key_adim = project.procPramArray[ci].key_adim
            mas_background_histo, [0,project.imndArray[ci].sdim-1], [key_adim,key_adim]
        END
        
        'BACKGROUND_HISTO_ARRAY': BEGIN
            PRINT, 'create histogram for all adims at the current slice'
            key_slice = project.procPramArray[ci].key_slice
            adim = project.imndArray[ci].adim-1
            mas_background_histo, [key_slice,key_slice], [0,adim]
        END
        
        'BACKGROUND_HISTO_ALL': BEGIN
            ;;PRINT, 'create histogram for all adims AND all slices'
            adim = project.imndArray[ci].adim-1
            mas_background_histo, [0,project.imndArray[ci].sdim-1], [0,adim]
        END
        
        'BACKGROUND_PREVIEW': BEGIN
            mas_load_state_1
            project.procPramArray[ci].key_slice = 0
            sz_data = size((*project.dataArray[ci].state1))
            FOR ii = 0, sz_data[3] -1 DO BEGIN
                
                project.procPramArray[ci].key_slice = ii
                update_threshold_display
                mas_background_filter_redraw
               
                wait, .0625
            ENDFOR
        END
        
    ENDCASE
    
    mas_background_filter_redraw
    
END


; Subroutine name: mas_background_filter
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

pro mas_background_filter, event
    COMPILE_OPT IDL2
    ;bring in the appropriate global variables
    COMMON scan_data
    COMMON common_widgets

    CI = project.ci
    mas_load_state_1

    proj_size = size(*project.dataarray[project.ci].state1)
    
    ;;If the ADT ROI window is already on the screen then don't make another one.
    IF N_ELEMENTS(REMOVE_BACKGROUND_BASE) EQ 1 THEN BEGIN
        IF  WIDGET_INFO( REMOVE_BACKGROUND_BASE,  /VALID_ID ) eq 1 THEN return
    ENDIF

    title = string ('MAS Background Filter')
    REMOVE_BACKGROUND_BASE = widget_base(TITLE=title, $
                                         UVALUE = 'REMOVE_BACKGROUND_BASE', $
                                         XOFFSET=420 ,YOFFSET=0  ,COLUMN=1)

    ;;Activate
    base = widget_base(REMOVE_BACKGROUND_BASE, /NONEXCLUSIVE)
    void = widget_button(base, value='ACTIVATE',UNAME = 'BACKGROUND_ACTIVATE')


    ;;histogram
    base = widget_base(REMOVE_BACKGROUND_BASE, row=2, FRAME=1)
    void = widget_label(base, value='Show Histogram(s)', /ALIGN_CENTER)
    base2 = widget_base(base, row=1)
    void = widget_button(base2, value='Key Image',UNAME = 'BACKGROUND_HISTO_CURRENT')
    void = widget_button(base2, value='Slice',UNAME = 'BACKGROUND_HISTO_SLICE')
    void = widget_button(base2, value='Array',UNAME = 'BACKGROUND_HISTO_ARRAY')
    void = widget_button(base2, value='All ',UNAME = 'BACKGROUND_HISTO_ALL')

    ;;key image
    base = widget_base(REMOVE_BACKGROUND_BASE, row=2, FRAME=1)
    void = widget_label(base, value='Key Image', /ALIGN_CENTER)
    base2 = widget_base(base, row=1)

    IF project.imndArray[CI].sdim EQ 1 THEN BEGIN
        void = widget_slider(base2, value=1, MINIMUM = 0, MAXIMUM=1, UNAME ='BACKGROUND_KEY_SLICE', title='Slice')
    ENDIF ELSE BEGIN
        void = widget_slider(base2, value=project.procPramArray[ci].key_slice, $
                             MINIMUM = 0, MAXIMUM=project.imndArray[CI].sdim-1, $
                             UNAME ='BACKGROUND_KEY_SLICE', title='Slice')
    ENDELSE

    IF project.imndArray[CI].adim EQ 1 THEN BEGIN
        void = widget_slider(base2, value=1, $
                             MINIMUM = 0, MAXIMUM=1, SENSITIVE=0, $
                             UNAME ='BACKGROUND_KEY_ADIM', title='Adim' )
    ENDIF ELSE BEGIN
        void = widget_slider(base2, value=project.procPramArray[ci].key_adim, $
                             MINIMUM = 0, MAXIMUM=project.imndArray[CI].adim-1, $
                             UNAME ='BACKGROUND_KEY_ADIM', title='Adim')
    ENDELSE

    ;;threshold
    base = widget_base(REMOVE_BACKGROUND_BASE, row=3, FRAME=1)
    base1 = widget_base(base, row=1)
    
    thr_low = CW_FSLIDER(base1 , UNAME='BACKGROUND_THRESH_LOWER_BOUND' $
                         ,MINIMUM=0.0 $
                         ,MAXIMUM=1.0  $
                         ,VALUE=0  $
                         ,TITLE ='Low Threshold', /EDIT )
    
    thr_high = CW_FSLIDER(base1 , UNAME='BACKGROUND_THRESH_UPPER_BOUND' $
                          ,MINIMUM=0.0 $
                          ,MAXIMUM=1.0  $
                          ,VALUE=0  $
                          ,TITLE ='Hi Threshold', /EDIT )
    
;    thr_base = widget_base(REMOVE_BACKGROUND_BASE, /row, /align_center)

;    thr_draw_base = widget_base(thr_base, /column, /align_center)
    thr_draw = WIDGET_DRAW(base, XSIZE= proj_size[1], $;project.imndArray[CI].fdim, $
                           YSIZE=proj_size[2], $ ;project.imndArray[CI].pdim, $
                           UNAME='BACKGROUND_THRESH_DRAW'  )
;    lbl = widget_label(thr_draw_base, value="Threshold Mask", /align_center)

;    er_dl_draw_base = widget_base(thr_base, /column)
    er_dl_draw = WIDGET_DRAW(base, XSIZE= proj_size[1], $;project.imndArray[CI].fdim, $
                             YSIZE= proj_size[2], $; project.imndArray[CI].pdim, $
                             UNAME='BACKGROUND_EROSION_DILATION_DRAW'  )
;    lbl = widget_label(er_dl_draw_base, value="Final Mask", /align_center)

    ;;generate and preview buttons
    base = widget_base(REMOVE_BACKGROUND_BASE, COLUMN=3)
    btn_preview = widget_button(base, value='Preview',UNAME = 'BACKGROUND_PREVIEW')
    
    state = ptr_new({ thr_low:thr_low, $
                      thr_high:thr_high, $
                      thr_draw:thr_draw, $
                      er_dl_draw:er_dl_draw }, /no_copy)
                      
    WIDGET_CONTROL, REMOVE_BACKGROUND_BASE, set_uvalue=state
    WIDGET_CONTROL, REMOVE_BACKGROUND_BASE, /realize
    
    mas_background_filter_redraw
    
    xmanager, 'mas_background_filter',REMOVE_BACKGROUND_BASE,/no_block, GROUP_LEADER=WID_BASE_MAIN

END


; Subroutine name: mas_background_filter_file
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting


PRO mas_background_filter_file

END
