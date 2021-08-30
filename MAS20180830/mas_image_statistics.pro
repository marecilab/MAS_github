;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

function reformat_stat_display_string, display_string_array

    header = display_string_array[0]
    
    hdr_arr = strsplit(header, string(09B), /extract)
    nfields = n_elements(hdr_arr)
    nrows   = n_elements(display_string_array)
    
    data_matrix = strarr(nfields, nrows)
    for r = 0, nrows-1 do begin
        data_matrix[*,r] = strtrim(strsplit(display_string_array[r], string(09B), /extract), 2)
    endfor

    ret = strarr(nrows)
    for f = 0, nfields-1 do begin
        lens = strlen(data_matrix[f,*])
        maxlen = max(lens)
        diffs = maxlen - lens
        for r = 0, nrows-1 do begin
            if (diffs[r] gt 0) then begin
                if (r eq 0) then begin ;; pad the header the same as the data rows for now
                    data_matrix[f,r] = data_matrix[f,r]+strjoin(replicate(' ', diffs[r]))
                endif else begin
                    data_matrix[f,r] = data_matrix[f,r]+strjoin(replicate(' ', diffs[r]))
;;                    data_matrix[f,r] = strjoin(replicate(' ', diffs[r]))+data_matrix[f,r]
;;  uncomment previous line to "pre-pad" the fields insteda of post-padding them.
                endelse
            endif
        endfor
    endfor
    
    for r = 0, nrows-1 do begin
        ret[r] = strjoin(data_matrix[*,r], string(09B))
    endfor
    
    return, ret

end

pro display_stats_event, event

    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        void = dialog_message(['The file could not be written', $
            'Please make sure that you have write permissions.'], $
            /error, /center)
        return
    endif
    
    uname = widget_info(event.id, /uname)
    if (uname ne 'btn_save') then return
    
    widget_control, event.top, get_uvalue=title
    
    textarea = widget_info(event.top, find_by_uname='display_area')
    widget_control, textarea, get_value=text
    
    dest = dialog_pickfile(file=title, default_extension='txt')
    
    if (dest ne '') then begin
        openw, lun, dest, /get_lun
        printf, lun, text
        close, lun
        free_lun, lun
    endif

    catch, /cancel
    
end


; Subroutine name: calculate_actual_string_length
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; The function takes in an input string and calculates the string length
; in number of characters displayed on the screen with reference to tab's in
; particular. tab's are set to equal 4 spaces

; Editing Information:
    ;Edited by HS 2006/10/04
    ;Fix spelling mistakes and commenting.

function calculate_actual_string_length, string_to_measure
    void = where(byte(string_to_measure) eq 09,count)
    return, strlen(string_to_measure)+count*5
end


; Subroutine name: display_stats
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; This procedure will take in a string array and display that string array on the screen
; so that it can be copied into another program , and data will be tab seperated.
; title will the be the title of the widget main

; Editing Information:
    ;Edited by HS 2006/10/04
    ;Fix spelling mistakes and commenting.

pro display_stats, display_string_array, title, reformat_columns=reformat_columns

    if (keyword_set(reformat_columns)) then begin
        display_string_array = reformat_stat_display_string(display_string_array)
    endif
    ;in order to display the text properly we need to resize the text widget accordingly
    sz_string_array = size(display_string_array)

    ;setting the y direction on the widget is easy now all we have to do is
    ; set the x size this will be the length of the largest string
    x_size = max(strlen(display_string_array))


    for ii=0, sz_string_array[1]-1 do begin
       x_size2 = calculate_actual_string_length(display_string_array[ii])
       if x_size2 gt x_size then x_size=x_size2
    end

;    print, 'print stats x_size=',x_size
;    print,'size of display_string_array',sz_string_array

    base = widget_base(/COLUMN, TITLE=title, mbar=mbar, uvalue=title)

    btn_file = widget_button(mbar, value='File', /menu)
    btn_save = widget_button(btn_file, uname='btn_save', value='Save...')

    if sz_string_array[0] gt 1 then begin
       y_Size = sz_string_array[2]+1
    end else begin
       y_Size = sz_string_array[1]+1
    end

    ;there may be alot to display and i don't want it to exceed the size of the visable area on the screen.
    screen_size = get_screen_size()
    wid_uname = 'display_area'
    
    if x_size gt 200 and y_Size gt 70 then begin
       ;if both are too big then cut the size down and let them scroll the text widget.
       w=widget_text(base, /SCROLL, value = display_string_array, uname=wid_uname, $
                     SCR_XSIZE =screen_size[0]-100, SCR_YSIZE =screen_size[1]-100 )

    end else if x_size gt 200 then begin

       ;if only x is larger then let them scroll in the x dir and set y approatly.
       w=widget_text(base,/SCROLL,value = display_string_array, uname=wid_uname, $
                     SCR_XSIZE =screen_size[0]-100, YSIZE=y_Size )

    end else if y_Size gt 70 then begin

       ;if only y is larger then let them scorll in the y dir and set x approatly.
       w=widget_text(base,/SCROLL,value = display_string_array, uname=wid_uname, $
                     XSIZE=x_size, SCR_YSIZE =screen_size[1]-100 )

    end else begin
       ;both sizes are fine.
       w=widget_text(base,/SCROLL,value = display_string_array, uname=wid_uname, $
                     XSIZE=x_size, YSIZE=y_Size )

    end

    widget_control, base, /realize

    xmanager, 'display_stats', base,/no_block

end

; Subroutine name: do_image_statistics
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04
    ;Fix spelling mistakes and commenting.
    ;Edited by MK 2015/01/11
    ; Added SNR to the statistics
    
pro do_image_statistics, ev
    COMPILE_OPT IDL2
    ;bring in the appropriate global variables
    COMMON scan_data
    CI = project.ci

    sdimStart = project.procPramArray[CI].sdim_Start
    adimStart = project.procPramArray[CI].adim_Start
    sort_dir = project.procpramarray[ci].sort_dir

    ;copy the image locally. this image will have the roi drawn on it
    ;process the image the way the user has set up.
    
    if (project.procpramarray[ci].no_transform_roi eq 1) then begin

        ;make sure the image is in state1
        mas_load_state_2
    
        image = *project.dataArray[project.ci].state2
    
    endif else begin
    
        mas_load_state_2, native_orientation=1
        
        image = *project.dataarray[project.ci].state2
        
        ;case project.procpramarray[ci].slice_axis of
        ;    0: image = (*project.dataarray[project.ci].state1)[*,*,sdimStart,adimStart]
        ;    1: image = (*project.dataarray[project.ci].state1)[*,sdimStart,*,adimStart]
        ;    2: image = (*project.dataarray[project.ci].state1)[sdimStart,*,*,adimStart]
        ;endcase
        
        image = temporary(reform(image))
    
    endelse
    
    sz_image = size(image)

;    if sz_image[0] ge 3 then begin
;       void = dialog_message(['Multi Slice display is not supported yet.' $
;                              ,'Please use single display when doing image statistics'],/error, /center )
;       return
;   end

    if project.roi.ci eq 0 then begin
        update_status_bar,'Please select a different item in the ROI list'
        return
    end
    
    ;appropriate to make sure that the there is an roi to draw on.
    if not ptr_valid(project.roi.pROIs[project.roi.ci]) then begin
        void =  dialog_message('Please either make an ROI or change selection',/error, /center)
        return
    end

    ;get the data object from the roi program.
    ;regions = *project.roi.pROIs[project.roi.ci]->get(/all)
    regions = (*project.roi.pROIs[project.roi.ci])

    ;appropriate to make sure the object regions is valid
    if not(total(obj_valid(regions))) then begin
       ;update_status_bar, 'No region returned from ROI tool'
       return
    end

    if (size(regions))[0] eq 0 then begin
        numRoi = 1
    endif else begin
        numRoi =  (size(regions))[1]
    endelse

    ;make an array to hold all the names of the arrays
    name_array =  strarr(numRoi)

    ; make a pointer array to hold all the masks
    mask = ptrarr(numRoi)
    
    img_dims = [ (size(image))[1], (size(image))[2] ]

    ;mas_roi_transform, regions, img_dims, /current
    for temp=0 , numRoi-1 do begin
        mask[temp] = ptr_new( regions[temp]->ComputeMask(dimensions=img_dims , MASK_RULE=2) )
        regions[temp] -> GETPROPERTY, name= name
        name_array[temp] = name
    end
    ;mas_roi_transform, regions, img_dims, /native

    
    nslices = sz_image[0] eq 3 ? sz_image[3] : 1
    if (nslices eq 1) then begin
        loop_start = (sort_dir eq 0) ? sdimStart : adimStart
        loop_end   = loop_start
    endif else begin
        loop_start = 0
        loop_end = nslices-1
    endelse
    
    header_printed = 0
    
; ================== Part that Ali Added ========================
    
    for sl = 0, nslices-1 do begin
    
        slice_image = reform(image[*,*,sl])
        stats_string = strarr(numRoi)
        if (header_printed eq 0) then begin
            final_output = ((sort_dir eq 0) ? 'Slice#': 'Array#')+string([09B])+ $
                           ((sort_dir eq 0) ? 'Array#': 'Slice#')+string([09B])+ $
                            'ROI#'+string([09B])+ $
                            'Signal'+string([09B])+ $
                            'Mean'+string([09B])+ $
                            'Stddev'+string([09B])+ $
                            'SNR'   +string([09B])+ $
                            'Sum Of Squares'+string([09B])+ $
                            'Min'+string([09B])+ $
                            'Max'+string([09B])+ $
                            'Var'+string([09B])+ $
                            'Count'
            header_printed = 1
        endif 
        
        maskCounter = 0
        for maskCounter= 0 , numRoi-1 do begin
        
            IMAGE_STATISTICS  , slice_image, mask=*(mask[maskCounter])  $
                , MEAN=image_mean $
                , STDDEV=image_stddev $
                , SUM_OF_SQUARES=image_sum_of_squares $
                , MINIMUM=image_min $
                , MAXIMUM=image_max $
                , VARIANCE=image_var $
                , Count=image_count
                
            project.procPramArray[project.ci].image_var = image_var
            
            ;Calculating the Signal.
            Image_Signal=(sqrt(image_sum_of_squares/image_count-2*image_var))
            Image_SNR = Image_Signal/image_stddev
            stats_string[maskCounter] = $
                string(sl+loop_start,format='(I03)') + string([09B]) + $
                string((sort_dir eq 0) ? adimStart : sdimStart ,format='(I03)') + string([09B]) + $
                name_array[maskCounter]  + string([09B]) + $
                string(Image_Signal) + string([09B]) + $
                string(image_mean)   + string([09B]) + $
                string(image_stddev) + string([09B]) + $
                string(Image_SNR) + string([09B]) + $
                string(image_sum_of_squares) + string([09B]) + $
                string(image_min) + string([09B]) + $
                string(image_max) + string([09B]) + $
                string(image_var) + string([09B]) + $
                string(image_count)
                
        endfor
        
        final_output = (n_elements(final_output) eq 0) ? stats_string : [ final_output, stats_string ]
        
    endfor

    display_stats, final_output, 'Image Statistics' ;;stats_string, 'Image Statistics'
   
end


; Subroutine name: plot_histogram
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; This procedure will plot the histogram for the mask specified
; the default number of bins is 100

; Editing Information:
    ;Edited by HS 2006/10/04
    ;Fix spelling mistakes and commenting.
    ;Edited by MK 2014/11/30
    ; Bug: Image array indexed using circular brackets which IDL interprets as the image function.
    ; Solution : Circular brackets replaced with square brackets

pro plot_histogram, image, mask, title, roi_name


    IMAGE_STATISTICS, image, mask=mask, MINIMUM=image_min, MAXIMUM=image_max
    nbins = 100

    ;we have to calculate the independent var from image max to image mean in steps of
    ;binsize
    BINSIZE = (image_max - image_min) / (NBINS - 1)
    indep_var = fltarr(nbins)
    for ii=0, nbins-1 do indep_var[ii] = ii*BINSIZE+image_min


    histoplot = HISTOGRAM(image[where(mask)],nbins=nbins)

    ;name the window
    windowName = 'B Value Vs mean least squares fit ROI'


    IPLOT, indep_var, histoplot,  /HISTOGRAM, $
       title = title, $
       IDENTIFIER=iPlot_indentifer, $
       NAME= roi_name, $
       XTitle = 'Density for '+roi_name, $
       YTitle = 'Occurrences'

end


; Subroutine name: do_image_histogram
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:


; Editing Information:
    ;Edited by HS 2006/10/04
    ;Fix spelling mistakes and commenting.
    ;Changed name from hisotgram to histogram.

pro do_image_histogram, ev
    COMPILE_OPT IDL2
    ;bring in the appropriate global variables
    COMMON scan_data
    CI = project.ci

    sdimStart = project.procPramArray[CI].sdim_Start
    adimStart = project.procPramArray[CI].adim_Start

    ;make sure the image is in state1
    mas_load_state_2

    if (project.procpramarray[ci].no_transform_roi eq 1) then begin
        ;; copy the image locally. this image will have the roi drawn on it
        ;; process the image the way the user has set up.
    
        image = *project.dataArray[project.ci].state2
    
    endif else begin
    
        case project.procpramarray[ci].slice_axis of
            0: image = (*project.dataarray[project.ci].state1)[*,*,sdimStart,adimStart]
            1: image = (*project.dataarray[project.ci].state1)[*,sdimStart,*,adimStart]
            2: image = (*project.dataarray[project.ci].state1)[sdimStart,*,*,adimStart]
            endcase
            image = temporary(reform(image))
    
    endelse

    sz_image = size(image)

    if sz_image[0] ge 3 then begin
       void = dialog_message(['multi Slice display is not supported yet.' $
                    ,'Please use single display when doing image statistics'],/error, /center )
       return
    end

    if project.roi.ci eq 0 then begin
        update_status_bar,'Please select a different item in the ROI list'
        return
    end

    ;get the data object from the roi program.
;    widget_control, ev.id, get_uvalue=og
;    regions = og->get(/all)
    regions = (*project.roi.pROIs[project.roi.ci])


    ;appropriate to make sure the object regions is valid
    if not(total(obj_valid(regions))) then begin
       ;update_status_bar, 'No region returned from ROI tool'
       return
    end

    ;print, 'size of regions', size(regions)

    if     (size(regions))[0] eq 0 then numRoi = 1 $
    else numRoi =  (size(regions))[1]

    ;print, 'numRoi=',numRoi

    ;we have to make an array to hold the names of the roi's
    name_array = strarr(numRoi)

    ; make a pointer array to hold all the masks
    mask = ptrarr(numRoi)

    for temp=0 , numRoi-1 do begin
       mask[temp] = ptr_new(regions[temp] -> ComputeMask( dimensions =[ (size(image))[1], (size(image))[2] ] , MASK_RULE=2))
       regions[temp] -> GETPROPERTY, name= name
       name_array[temp] = name
    end

    maskCounter = 0
    for maskCounter= 0 , numRoi-1 do begin

       windowName = 'Histogram For '+name_array[maskCounter] ;strtrim(maskCounter,2)

       plot_histogram, image , *(mask[maskCounter]) , windowName, name_array[maskCounter]


    end

end


; Subroutine name: mas_image_statistics
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:


; Editing Information:
    ;Edited by HS 2006/10/04
    ;Fix spelling mistakes and commenting.

pro mas_image_statistics
    COMPILE_OPT IDL2

    ;bring in the appropriate global variables
    COMMON scan_data
    CI = project.ci

    sdimStart = project.procPramArray[CI].sdim_Start
    adimStart = project.procPramArray[CI].adim_Start

    ;make sure the image is in state1
    mas_load_state_2

    ;copy the image locally. this image will have the roi drawn on it
    ;process the image the way the user has set up.
    image = *project.dataArray[project.ci].state2
    mas_windowing , image
    ;mas_make_sheet_display, image

    sz_image = size(image)

    if sz_image[0] ge 3 then begin
       void = dialog_message(['multi Slice display is not supported yet.' $
                    ,'Please use single display when doing image statistics'],/error, /center )
       return
    end


    tlb=widget_base(/column)

;    print, 'project.roi_valid',project.roi_valid
;    help, where(obj_isa(obj_valid(),'idlgrroi'))
;    if project.roi_valid eq 1 then begin
;        ;mas_xroi, image , outgroup=og,apptlb=tlb, group=tlb, title='slice:'+strtrim(sdimStart+1,2), $
;            regions_in=*project.roi
;
;    end else begin
            ;mas_xroi, image , outgroup=og,apptlb=tlb, group=tlb, title='slice:'+strtrim(sdimStart+1,2)
;    end


    b=widget_button(tlb,value='Print Stats', event_pro='do_image_statistics', uvalue=og)
    b=widget_button(tlb,value='Histogram', event_pro='do_image_histogram', uvalue=og)


end
