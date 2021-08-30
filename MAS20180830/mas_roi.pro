;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

;+
; :Description:
;    For an open scan with ROIs it returns the 
;    number of regions contained in an ROI set.
;
;
;
;
;
; :Author: wtriplett
;-
function mas_roi_get_number_of_rois

    compile_opt idl2
    common scan_data
    
    ci = project.ci
    
    if (ptr_valid(project.roi.prois[project.roi.ci])) then begin

        if project.roi.ci ne 0 then begin
            oroi = (*project.roi.prois[project.roi.ci])
            return, n_elements(oroi)
        endif
    
    endif
    
    return, 0L
    
end

;+
; :Description:
;    For an open scan, this function returns the name or an array
;    of names for the regions in the currently selected ROI set.
;
;
;
; :Keywords:
;    region_num - optional, specify which region number to get the name of.
;
; :Author: wtriplett
;-
function mas_roi_get_current_name, region_num=region_num

    compile_opt idl2
    common scan_data
    
    ci = project.ci
    
    if (ptr_valid(project.roi.prois[project.roi.ci])) then begin

        if project.roi.ci ne 0 then begin
            oroi = (*project.roi.prois[project.roi.ci])
            if (n_elements(region_num) eq 0) then begin
                names = strarr(n_elements(oroi))
                for r = 0, n_elements(oroi)-1 do begin
                    oroi[r]->GetProperty, name=roi_name
                    names[r] = roi_name
                 endfor
                 return, names
            endif else begin
                oroi[region_num]->GetProperty, name=roi_name
                return, roi_name
            endelse
        endif
    
    endif
    
    return, ''
    

end

;+
; :Description:
;    Returns a mask with the regions of interest active, that can be
;    applied to the currently selected imagery.
;
; :Params:
;    dim - the dimension of the mask. If not provided, they will be
;          computed from the data dimensions and the currently selected
;          slice axis
;
; :Keywords:
;    region_num - optionally specify a single region to receive the mask of in
;                 an ROI set with multiple regions
;    crop       - this returns the minimum dimensions needed to contain the region.
;    no_transform - set this keyword to have the function ignore MAS's current
;                   zoom and rotate options, and return a mask that can be applied
;                   directly to the data.
;
; :Author: wtriplett
;-
function mas_roi_get_current_mask, dim, $
                                   region_num=region_num, $
                                   crop=crop, $
                                   no_transform=no_transform

    compile_opt idl2
    common scan_data
    
    ci = project.ci
    
    if (n_elements(dim) eq 0) then begin
        case project.procpramarray[ci].slice_axis of
            0: begin
                dim = [ project.imndarray[ci].fdim, project.imndarray[ci].pdim ]
                dim *= [ project.procpramarray[ci].freq_interp, project.procpramarray[ci].phase_interp ]
            end
            1: begin
                dim = [ project.imndarray[ci].fdim, project.imndarray[ci].sdim ]
                dim *= [ project.procpramarray[ci].freq_interp, project.procpramarray[ci].slice_interp ]
            end
            2: begin
                dim = [ project.imndarray[ci].pdim, project.imndarray[ci].sdim ]
                dim *= [ project.procpramarray[ci].phase_interp, project.procpramarray[ci].slice_interp ]
            end
        endcase
    endif
    
    xdim = dim[0]
    ydim = dim[1]

    if (ptr_valid(project.roi.prois[project.roi.ci])) then begin

        mask = ptr_new(bytarr(xdim, ydim), /no_copy)

        if project.roi.ci ne 0 then begin
            oroi = (*project.roi.prois[project.roi.ci])
            if (n_elements(region_num) ne 0) then begin
                roi_start = region_num
                roi_end   = region_num
            endif else begin
                roi_start = 0
                roi_end   = n_elements(oroi)-1
            endelse
            for roi = roi_start, roi_end do begin
                
                if not keyword_set(no_transform) then mas_roi_transform, oRoi, dim, /current

                *mask = *mask OR oroi[roi]->computeMask(mask_rule=2, $
                                                        dimensions=[xdim, ydim])

                if not keyword_set(no_transform) then mas_roi_transform, oRoi, dim, /native
            endfor
        endif

    endif else begin
        return, ptr_new()
    endelse

    sz = size(*mask, /dimensions)
    
    crop = [0, 0, sz[0]-1, sz[1]-1]

    yrange = total((*mask), 1)
    xrange = total((*mask), 2)
    
    tmp = where(yrange ne 0, tmp_ct)
    if (tmp_ct ne 0) then begin
        max_y = max(tmp, min=min_y)
        crop[1] = min_y-1 > 0
        crop[3] = max_y+2 < (sz[1]-1)
    endif

    tmp = where(xrange ne 0, tmp_ct)
    if (tmp_ct ne 0) then begin
        max_x = max(tmp, min=min_x)
        crop[0] = min_x-1 > 0
        crop[2] = max_x+2 < (sz[0]-1)
    endif
    
    return, mask

end

pro mas_roi_transform, oroi, img_dims, native=native, current=current

    common scan_data
    ci = project.ci
    ;;return

    if (keyword_set(native) and keyword_set(current)) then return
    if (not keyword_set(native) and not keyword_set(current)) then return

    ctr = img_dims/2.0

    voxel_dim = [ project.imndarray[ci].f_voxsz/project.procpramarray[ci].freq_interp , $
                  project.imndarray[ci].p_voxsz/project.procpramarray[ci].phase_interp, $
                  project.imndarray[ci].s_voxsz/project.procpramarray[ci].slice_interp ]
    voxel_dim /= min(voxel_dim)

    case project.procpramarray[ci].rotate_direction of
        0: rot_ctr = img_dims/2.0
        1: rot_ctr = [min(img_dims), min(img_dims)]/2.0
        2: rot_ctr = img_dims/2.0
        3: rot_ctr = [max(img_dims), max(img_dims)]/2.0
        else:
    endcase

    ;; the order of the transformations matters
    if (keyword_set(native)) then begin
        ;; "Back to native (relative to state1) orientation)
        case project.procpramarray[ci].rotate_direction of
            1: inv_rotate = -90
            2: inv_rotate = 180
            3: inv_rotate = 90
            else: inv_rotate = 0
        endcase
        
        case project.procpramarray[ci].flip_direction of
            1: reverse_dir = 1
            2: reverse_dir = 2
            else: reverse_dir = 0
        endcase
        
        x_reduce = 1.0/project.procpramarray[ci].x_zoom
        y_reduce = 1.0/project.procpramarray[ci].y_zoom
        
        for i=0, n_elements(oroi)-1 do begin
            
            if (not obj_valid(oroi[i])) then continue
            
            if (reverse_dir eq 1) then begin
                oroi[i]->rotate, [0,1,0], 180, center=[ctr[0], ctr[1]]
            endif else if (reverse_dir eq 2) then begin
                oroi[i]->rotate, [1,0,0], 180, center=[ctr[0], ctr[1]]
            endif

            oroi[i]->rotate, [0,0,1], inv_rotate, center=[rot_ctr[0], rot_ctr[1]]
            
            oroi[i]->scale, [x_reduce, y_reduce]
            
        endfor
        
    endif else begin
        ;; From state1 orientation to current
        case project.procpramarray[ci].rotate_direction of
            1: for_rotate = 90
            2: for_rotate = 180
            3: for_rotate = -90
            else: for_rotate = 0
        endcase
            
        case project.procpramarray[ci].flip_direction of
            1: reverse_dir = 1
            2: reverse_dir = 2
            else: reverse_dir = 0
        endcase
        
        x_reduce = project.procpramarray[ci].x_zoom
        y_reduce = project.procpramarray[ci].y_zoom
        
        for i=0, n_elements(oroi)-1 do begin
            
            if (not obj_valid(oroi[i])) then continue

            oroi[i]->scale, [x_reduce, y_reduce]
            
            oroi[i]->rotate, [0,0,1], for_rotate, center=[rot_ctr[0], rot_ctr[1]]
            
            if (reverse_dir eq 1) then begin
                oroi[i]->rotate, [0,1,0], 180, center=[ctr[0], ctr[1]]
            endif else if (reverse_dir eq 2) then begin
                oroi[i]->rotate, [1,0,0], 180, center=[ctr[0], ctr[1]]
            endif
            
        endfor

    endelse

end
; Subroutine name: mas_roi_mask
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This procedure will take in a pointer to a data object and apply the maks on top of it.

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

PRO mas_roi_mask, p_data, crop=crop, no_transform=no_transform

    common scan_data
    CI = project.ci

    mask_index = project.roi.mask
    IF mask_index LE 0 THEN RETURN

    sz_p_data = size(*p_data)

    no_transform = keyword_set(no_transform) ? 1 : 0

    pmask = mas_roi_get_current_mask([sz_p_data[1], sz_p_data[2]], $
                                    no_transform=no_transform)

    if (not ptr_valid(pmask)) then return

    print, max(*pmask)
    CASE sz_p_data[0] OF

        2: (*p_data) *= FLOAT(*pmask)
        3: FOR ii=0,  sz_p_data[3]-1 DO (*p_data)[*,*,ii] *= FLOAT(*pmask)
        ELSE: PRINT, 'CASE NOT HANDLED'
        
    ENDCASE
    
    ptr_free, pmask
    
END


; Subroutine name: mas_roi_get_image
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This function will return a pointer to the image to draw the roi on.
; It will keep ADT images in mind and the selections made on that control panel also.

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

function mas_roi_get_image
    common scan_data
    CI = project.ci

    ;populate the correct image to draw an roi on.
    if project.scan_Open_Flag eq 0 then begin
        ;if there is no scan opened then create an image and let them draw on that.

        pimage = ptr_new(bytscl(dist(128)))

    end else begin
        ;let the user draw on there own selected image.

        ;if the image is an ADT image let them draw the roi on their selected image.
        if (project.procpramarray[ci].adt_proccess_flag ne 0 or $
            project.imndArray[ci].image_type eq 3 or $
            project.imndArray[ci].image_type eq 11 or $
            project.imndarray[ci].image_type eq 16) then begin

            ;don't do the fitting if the want to draw the image on a diffusion weighted image.
            if project.procPramArray[ci].adt_display_type ne 4 then begin

                ;load the ADT image
                do_adt_regress_mono

                ;make the image
                case project.procpramArray[ci].slice_axis of
                   0: slice = project.procPramArray[ci].sdim_start
                   1: slice = project.procPramArray[ci].pdim_start
                   2: slice = project.procPramArray[ci].fdim_start
                endcase
                image = adt_make_image(slice)

                ;if project.procpramArray[ci].slice_axis eq 0 then $
                ;    image = adt_make_image (project.procPramArray[ci].sdim_start)
                ;if project.procpramArray[ci].slice_axis eq 1 then $
                ;    image = adt_make_image (project.procPramArray[ci].pdim_start)
                ;if project.procpramArray[ci].slice_axis eq 2 then $
                ;    image = adt_make_image (project.procPramArray[ci].fdim_start)

                ;;rotate, flip and zoom the image

                pimage = ptr_new(image)
                mas_rotate_flip, pimage
                mas_zoom, pimage

            end else begin;; The image is ADT, just that diffusion weighted is chosen

                ;load the image in to the structure and prepare for processing
                mas_load_state_2

                ;if the image has 3 dimensions then it must be multi slice data
                ;and drawing an roi on mulitple slice data is not allowed.
                ;although i don't see why not it just easier this way.
                if (size(*project.dataArray[ci].state2))[0] gt 2 then begin
                    void = dialog_message(['Multiple slice data sets are not permited', $
                                           'to have ROIs on them yet'],/error, /center)
                    return,-1
                end

                image = (*project.dataArray[ci].state2)[*,*]
                mas_windowing, image
                pimage = ptr_new(image,/no_copy)
            end

        end else begin ;; The image is _NOT_ an ADT image at all
            ; basically just repeat the above block of code

            ;load the image in to the structure and prepare for processing
            mas_load_state_2

            ;if the image has 3 dimensions then it must be multi slice data
            ;and drawing an roi on multiple slice data is not allowed.
            ;although i don't see why not it just easier this way.
            if (size(*project.dataArray[ci].state2))[0] gt 2 then begin
               void = dialog_message(['Multiple slice data sets are not permited', $
                                      'to have ROIs on them yet'],/error, /center)
                return,-1
            end

            image = (*project.dataArray[ci].state2)[*,*]
            mas_windowing, image
            pimage = ptr_new(image,/no_copy)
        end
     end

    return, pimage
end


; Subroutine name: mas_roi_add
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This procedure will add a pointer to an roi object to the roi structure.
; It will add the roi data to the end of the structure.

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_roi_add, p_oRoi
    COMMON scan_data

     ;we have to make sure that there is enough room in the structure for the new roi.
     if project.roi.ni ge project.roi.max_num_roi then begin
        void = dialog_message(['Max number of ROI has been reached.', $
                               'Please remove an ROI before adding another.'], $
                              /error, /center)
        return
     end

    ;we need to check to make sure they drawed an roi.
    ;if they drew nothing then exit quietly.
    ;check to make sure the object regions is valid
    if not(total(obj_valid(*p_oRoi))) then begin
       update_status_bar, 'No region returned from ROI tool'
       return
    end
    ;if they get past this point then clean up the status bar.
    update_status_bar, ''


    ;now that we are sure that we have room for the current ROI we can add it to the structure.
    project.roi.pROIs[project.roi.ni] = p_oRoi


    ;add the name to the next index on the roi selction screen.
    *project.roi.pdisplay_names = [ [*project.roi.pdisplay_names],$
                                    ['ROI set#'+strtrim(string(project.roi.ni),2)] ]
    ;*project.roi.pdisplay_names += [['ROI set#'+strtrim(string(project.roi.ni),2)]]

    ;move the current index to the one just added.
    project.roi.ci  = project.roi.ni
    ;increment the next index
    project.roi.ni ++

end


; Subroutine name: mas_new_roi
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Create a new roi.

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_new_roi, event
    COMMON scan_data
    CI = project.ci

    ;get the image to draw the roi on.
    pImage = mas_roi_get_image ()

    if not(ptr_valid(pImage)) then return

    img_dims = size(*pimage, /dimensions)

    ;block to make sure that no other xroi pros can be running at the same time
    if project.roi.xroi_running eq 1 then begin
        update_status_bar,'Close current ROI program before opening another'
        return
    end
    ;project.roi.xroi_running = 1
    xroi, *pimage , regions_out=oRoi, /block, /modal ,group = event.top
    project.roi.xroi_running = 0
    
    if (project.procpramarray[ci].no_transform_roi eq 0) then begin
        mas_roi_transform, oRoi, img_dims, /native
    endif
    
    ;add the roi
    mas_roi_add, ptr_new(oRoi)    
end


; Subroutine name: mas_roi_edit
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Edit an existing roi

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_roi_edit, roi_index, event
    COMMON scan_data
    CI = project.ci


    ;get the image to draw the roi on.
    pImage =  mas_roi_get_image()
    
    if (not ptr_valid(pImage)) then return

    img_dims = (size(*pImage, /dimensions))[0:1]

    ;gather the current roi object
    oRoi =  *project.roi.pROIs[roi_index]

    ;lock out the roi so it cant be called more than once.
    ;project.roi.xroi_running = 1

    if (project.procpramarray[ci].no_transform_roi eq 1) then begin
    
        xroi, *pImage, regions_in = oRoi, regions_out=oRoi, /block, /modal ,group = event.top
        project.roi.xroi_running = 0
        
    endif else begin
    
        mas_roi_transform, oRoi, img_dims, /current

        xroi, *pImage, regions_in = oRoi, regions_out=oRoi, /block, /modal ,group = event.top
        project.roi.xroi_running = 0

        mas_roi_transform, oRoi, img_dims, /native
    
    endelse
    
    ;we need to check to make sure they drawed an roi.
    ;if they drew nothing then exit quietly.
    ;check to make sure the object regions is valid
    if not(total(obj_valid(oRoi))) then begin
       update_status_bar, 'No region returned from ROI tool'
       return
    end

    project.roi.pROIs[roi_index] = ptr_new(oRoi)

end


; Subroutine name: mas_roi_save
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Save the current index to a file

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_roi_save, event
    common scan_data
    roiCI = project.roi.ci

    ;you can't choose the new roi
    if roiCI eq 0 then return

    file = (Dialog_PickFile( Title = 'Save Current ROI', path=project.current_path))

    if file eq '' then return

    roi = *project.roi.pROIs[roiCI]
    save, roi, FILENAME=file+'.sav'


end


; Subroutine name: mas_roi_load
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This procedure will load a roi from a file and store it locally at the end of the roi pointer

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_roi_load, event
    COMMON scan_data


    roi_file = dialog_pickfile( PATH=project.current_Path)

    if roi_file eq '' then return

    RESTORE, roi_file , RESTORED_OBJECTS = ORoi, /RELAXED_STRUCTURE_ASSIGNMENT

    mas_roi_add, ptr_new(ORoi,/no_copy)

    WIDGET_CONTROL, Widget_Info(event.top, FIND_BY_UNAME='w_roi_list')$
        ,set_value = *project.roi.pdisplay_names, SET_LIST_SELECT=project.roi.ci


end


; Subroutine name: mas_roi_delete
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Delete the current roi and shift all the rest into its place.

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_roi_delete, event

    COMMON scan_data
    roiCI = project.roi.ci
    roiNI = project.roi.ni
    wWidget = event.top
    ;you are not allowed to remove the New Roi from the list.
    if roiCI eq 0 then return

    ;since we are removing one we won't have to do this one
    roiNI --

    ;start at the current index and go until the next index
    for ii = roiCI, roiNI do begin

        ;this pointer gets the next one in line
        project.roi.pROIs[ii] = project.roi.pROIs[ii +1 ]

        if not (ii ge roiNI) then $
            (*project.roi.pdisplay_names)[0,ii]=(*project.roi.pdisplay_names)[0,ii +1 ]

    end

    ;printing out some values to make sure that the data is being deleted.
    ;for cc=0, 10 do  help,project.roi.pROIs[cc]

    temp = (*project.roi.pdisplay_names)[0,0:ii-2]

    project.roi.pdisplay_names = ptr_new(temp)

    project.roi.ci --
    project.roi.ni --

    if project.roi.ci eq 0 then begin
        WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_new_edit'), set_value ='New'
        WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_delete'), SENSITIVE=0
        WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_save'), SENSITIVE=0
        WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_name'), set_value='', SENSITIVE=0

    end

    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_list'), $
        set_value = *project.roi.pdisplay_names, SET_LIST_SELECT=project.roi.ci

end


; Subroutine name: mas_roi_event
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This is the event call back.

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_roi_event, event
    COMMON scan_data

    catch, error_status
    if (error_status ne 0) then begin
        catch, /cancel
        help, calls= trace_back
        
        dummy = DIALOG_MESSAGE(['Please send a copy of this error message to the developer of this program.', $
                                !ERROR_STATE.msg,trace_back], /ERROR, $
                               TITLE='Error in WID_BASE_MAIN_event')
       RETURN
    endif


    CI = project.ci
    wWidget =  Event.top

    ;if xroi is running then don't let them open another one.
    ;this is a lock to keep race conditions from happening.
    if project.roi.xroi_running eq 1 then begin
        update_status_bar,'ROI editor alreay running. Please close before editing another.'
        return
    end
    ;clean up status bar.
    update_status_bar,''

    case Event.id of

        Widget_Info(wWidget, FIND_BY_UNAME='w_roi_list'): begin

            ;if they just click on the list change the name of the new/edit button.
            if event.clicks eq 1 then begin

                if event.index eq 0 then begin
                ;if they click on the first item in the list

                    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_new_edit'), set_value ='New'
                    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_delete'), SENSITIVE=0
                    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_save'), SENSITIVE=0
                    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_name'), set_value='', SENSITIVE=0

                    ;change the index to the current roi selcted.
                    project.roi.ci = event.index

                end else begin
                ;if they click on any other item in the list.
                    ;change the label on the button to respresent what is the current course of action.
                    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_new_edit'), set_value ='Edit'
                    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_delete'), SENSITIVE=1
                    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_save'), SENSITIVE=1
                    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_name'), $
                        set_value=(*project.roi.pdisplay_names)[event.index], $
                        SENSITIVE=1

                    ;change the index to the current roi selected.
                    project.roi.ci = event.index

                end

            end


            ;if the user double clicks then let them edit or draw a new roi.
            if event.clicks eq 2 then begin
                ;did they select the new ROI if they did then let them draw a new roi
                if event.index eq 0 then begin

                    ;change the name of the roi button.
                    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_new_edit'), set_value ='New '

                    ;let the user create a new roi.
                    mas_new_roi, event

                    ;now we update the selection list widget.
                    WIDGET_CONTROL, Event.id, set_value = *project.roi.pdisplay_names, SET_LIST_SELECT=project.roi.ci
              WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_new_edit'), set_value ='Edit'
              WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_delete'), SENSITIVE=1
                    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_save'), SENSITIVE=1
                    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_name'), $
                        set_value=(*project.roi.pdisplay_names)[project.roi.ci], $
                        SENSITIVE=1


                ;they double clicked and want to edit an existing roi.
                end else begin

                    WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_new_edit'), set_value ='Edit'
                    ;change the index to the current roi selected.
                    project.roi.ci = event.index
                    ;allow them to edit the current index.
                    mas_roi_edit, event.index, event
                    if (project.procPramArray[CI].flow_ecc_stat eq 1) then project.procPramArray[CI].flow_proccess_flag = 0
                end

            end



        end
        Widget_Info(wWidget, FIND_BY_UNAME='w_roi_name'): begin

            ;for this event we will get the name from the text box and rewrite the name in the roi list box.
            ;extract the new name from the text box
            WIDGET_CONTROL, event.id, GET_VALUE=name
            ;store that name to the structure
            (*project.roi.pdisplay_names)[project.roi.ci] = name
            ;now we update the selection list widget.
            WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_list'), $
                set_value = *project.roi.pdisplay_names, SET_LIST_SELECT=project.roi.ci

        end

        Widget_Info(wWidget, FIND_BY_UNAME='w_roi_new_edit'): begin

            ;if they select the first item in the scan list let them draw a new roi
            if project.roi.ci eq 0 then begin

                ;let the user create a new roi.
                mas_new_roi, event

                ;now we update the selection list widget.
                WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_list')$
                    , set_value = *project.roi.pdisplay_names, SET_LIST_SELECT=project.roi.ci

                WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_new_edit'), set_value ='Edit'
          WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_delete'), SENSITIVE=1
                WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_save'), SENSITIVE=1
                WIDGET_CONTROL, Widget_Info(wWidget, FIND_BY_UNAME='w_roi_name'), $
                    set_value=(*project.roi.pdisplay_names)[project.roi.ci], $
                    SENSITIVE=1


            end else begin
            ;if they select any other item in the list allow the usere to edit the roi.
                mas_roi_edit, project.roi.ci, event                
            end


        end
        Widget_Info(wWidget, FIND_BY_UNAME='w_roi_load'): begin
            mas_roi_load, event

        end
        Widget_Info(wWidget, FIND_BY_UNAME='w_roi_save'): begin
            mas_roi_save, event

        end
        Widget_Info(wWidget, FIND_BY_UNAME='w_roi_delete'): begin
            mas_roi_delete, event

        end
       Widget_Info(wWidget, FIND_BY_UNAME='w_roi_mask'): begin
            project.roi.mask = event.index
         project.procPramArray[ci].state_2 = 0

        end

       Widget_Info(wWidget, FIND_BY_UNAME='w_roi_xform'): begin
         project.procPramArray[ci].no_transform_roi = event.select
        end

      Widget_Info(wWidget, FIND_BY_UNAME='roi_graph_toggle'): begin
         project.procPramArray[ci].curvefit_graph = event.select
        end
    
    endcase
    mas_redraw_GUI

end


pro mas_roi
end
