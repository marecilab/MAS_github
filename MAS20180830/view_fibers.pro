;; $Id$
;; Copyright 2008 University of Florida. All Rights Reserved

;*****************************************************************************************************
;
; NAME:
;   make_fibers, filename
;
; PURPOSE:
;   read in a tracts.trs data file to make a fiber array
;
; ARGUMENTS:
; filename - the file to read in. Usually 'tracts.trs'
;
; MODIFICATION HISTORY:
;  20080421 - Created BT
;
;*****************************************************************************************************


;; Subroutine name: make_fibers
;; Created by: BT 20080421
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine: read tracts.trs file and create fiber array
;;
;; Editing Information:
;;
pro make_fibers, filename, fiber_data=fiber_data, transform=transform, use_new_method=use_new_method
;;function make_fibers, filename

    if not file_test(filename) then begin
        junk = dialog_message('Cannot locate: '+filename, /center, /error)
        return
    endif

    restore, filename

    filedir   = file_dirname(filename)
    tractfile = file_basename(filename)

    if (stregex(tractfile, '^(.+)\.[Tt][Rr][Ss]$') ne -1) then begin
        junk = stregex(tractfile, '^(.+)\.[Tt][Rr][Ss]$', /extract, /subexpr)
        roifile = junk[1]+'.ROI'
    endif else begin
        roifile = ''
    endelse

    if (n_elements(branch_levels) eq 0) then begin
        branch_levels = [0]
    endif
    
    if (not keyword_set(transform)) then begin
        transform = diag_matrix(replicate(1.0, 4))
    endif
    
    if (not keyword_set(use_new_method)) then begin
    ;; a single fiber track is stored in thie file as:
    ;; mid_point->end_point->mid_point->start_point. we
    ;; have to reverse the 2nd half and place it in front
    ;; of the first half
        
        nfibers = n_elements(start_points-1)
        fibers = ptrarr(nfibers-1)

        for j = 1L, nfibers-1 do begin
            
            beg = start_points[j-1]
            fin = start_points[j] - 1;; 'end' is reserved word.
            mid = mid_points[j-1]
            
            tmp = [ [ reverse(rrx[beg+1:mid-1],1), rrx[mid:fin] ], $
                    [ reverse(rry[beg+1:mid-1],1), rry[mid:fin] ], $
                    [ reverse(rrz[beg+1:mid-1],1), rrz[mid:fin] ] ]
            
            fibers[j-1] = ptr_new(transpose(tmp), /no_copy)
            
        endfor
    
    endif else begin
    ;; This is the routine for the new fiber tracking method, in which the
    ;; tracts are automatically flipped before stored.

        nfibers = n_elements(start_points)
        fibers = ptrarr(nfibers-1)
        
        sp_ind = 0L
        for j = 1L, nfibers-1 do begin
           if (sp_ind eq nfibers) then break
           tmp = [ [rrx[ start_points[sp_ind]:start_points[sp_ind+1]] ] , $
                   [rry[ start_points[sp_ind]:start_points[sp_ind+1]] ] , $
                   [rrz[ start_points[sp_ind]:start_points[sp_ind+1]] ] ]
           sp_ind += 2
           fibers[j-1] = ptr_new(transpose(tmp), /no_copy)
        endfor
    
    endelse

    ;; pass back some other information about the
    ;; tract job
    fdata = { fibers:    temporary(fibers), $
              seed_points: temporary(mid_points), $
              transform: transform, $
              branch_levels: branch_levels, $
              matrix:    matrix, $
              nfibers:   nfibers, $
              filedir:   filedir, $
              tractfile: tractfile, $
              roifile:   roifile }

    fiber_data = ptr_new(fdata, /no_copy)

end

;; Subroutine name: rebuild_fiber_cache, state
;; Created by: BT 20080421
;; Calling Information:
;;   state - the struct containing the state information. see below
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine: read in the raw fiber data array and cache
;;   information about roi containment, length, etc.
;;
;; Editing Information:
;;
pro rebuild_fiber_cache, state

    progressbar = OBJ_NEW('progressbar', Color='red', $
                          Title='Building fiber cache', $
                          /NOCANCEL)
    progressbar->SetProperty, TEXT='Building data cache...'

    progressbar -> Start
    
    fiber_data = (*state).fib_rawdata
    
    oroi_model = (*state).model->getByName('ROI_model')

    branch_levels = (*state).fib_branch_levels
    print, branch_levels

    curr_branch_lvl = 0
    n_branch_levels = n_elements(branch_levels)
    transform = invert((*state).transform)
    
    size_reduce = 1
    roi_work = bytarr([(*state).matrix/size_reduce, (*state).roi_count])
    fib_arr  = 0
    for r = 0, (*state).roi_count-1 do begin
        roi_name = strcompress('ROI_'+string(r), /remove_all)
        roi = oroi_model->getByName(roi_name)
        roi->getProperty, DATA=roi_pts
        for ro = 0L, n_elements(roi_pts)/3 - 1 do begin
            pt = floor(roi_pts[*,ro])/size_reduce
            roi_work[pt[0],pt[1],pt[2], r] = 1B
        endfor
    endfor
   
    ;; check roi containment and length of fibers
    for j = 0L, n_elements(fiber_data)-1 do begin
        
        progressbar->Update, (float(j)/n_elements(fiber_data)) * 100.0

        ; free the current cache
        if (ptr_valid( (*state).fib_cache[j] )) then begin
            ptr_free, (*state).fib_cache[j]
        endif

        fib_c = ulong64(0)
        
        curr_fiber = fiber_data[j]
        
        if (n_branch_levels gt 1 and j eq branch_levels[curr_branch_lvl]) then begin
            curr_branch_lvl = (curr_branch_lvl + 1) < (n_branch_levels-1)
        endif

        if not ptr_valid(curr_fiber) then continue

        ;; measured in number of coordinates.
        fsize = n_elements(*curr_fiber)/3
        
        for st = 0L, fsize-1 do begin
            pt = ((floor((*curr_fiber)[*,st]) > [0,0,0]) < ((*state).matrix-1))/size_reduce
            for ro = 0, (*state).roi_count-1 do begin
                if (roi_work[pt[0],pt[1],pt[2],ro] ne 0) then begin
                    fib_c = fib_c OR (ulong64(2)^ro)
                    continue
                endif
            endfor
        endfor
        
        if ( (*state).roi_count eq 0 ) then begin
            fib_c = 1
        endif

        if fib_c eq 0 then begin
            fib_c = 1
        endif

        if (curr_branch_lvl ne 0) then fib_c = 1

        (*state).fib_cache[j] = ptr_new({roi_hits:fib_c, $
                                         mid_point: (*state).fib_midpoints[j], $
                                         branch_level:curr_branch_lvl, $
                                         length: fsize}, /no_copy)

    endfor

    progressbar->Destroy

end

;; Subroutine name: rebuild_fiber_model, state
;; Created by: BT 20080421
;; Calling Information:
;;   state - the struct containing the state information. see below
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine: rebuilds fiber graphics model. usually
;;  in response to a GUI event.
;;
;; Editing Information:
;;
pro rebuild_fiber_model, state

    compile_opt idl2
    fiber_data = (*state).fib_rawdata

    colors = ['Khaki',      'Red',       'Green',        'Lime Green', 'Blue', 'Maroon', 'Cyan',  $
              'Magenta',    'Orange',    'Yellow',       'Orchid',     'Violet', $
              'Deep Pink',  'Tomato',    'Dark Green',   'Olive',$
              'Purple',     'Goldenrod', 'Saddle Brown', 'Coral',$
              'Brown',      'Navy',      'Light Cyan',   'Lavender',   'Green Yellow',$
              'Lawn Green', 'Honeydew']
    
    
    omodel_fibs_par = (*state).model->getByName('FIB_parent_model')
    
    obj_destroy, omodel_fibs_par->getByName('FIB_child_model')
    
    omodel_fibs_chl = obj_new('idlgrmodel', name='FIB_child_model', lighting=0)
    
    ;; set up fibers model
    fibs        = objarr((*state).fib_nfibers)
    n_discarded = 0
    n_fibs      = 0
    n_rois      = (*state).roi_count

    oroi_model = (*state).model->getByName('ROI_model')
    fstep = (*state).fib_fstep
    
    terminate = 0

    for j = 0, n_elements(fiber_data)-1, fstep do begin

        fib_c = ulong64(0)
    
        curr_fiber = fiber_data[j]

        ;; this should never happen
        if not ptr_valid(curr_fiber) then continue

        ;; check the length
        fsize = (*( (*state).fib_cache[j] )).length
        if fsize lt (*state).fib_len_threshold then continue

        ;; check whether to draw this fiber based on which
        ;; ROIs the user wants to see
        fib_c = (*( (*state).fib_cache[j] )).roi_hits
        mask = (*state).roi_mask and fib_c

        branch_lvl = (*( (*state).fib_cache[j] )).branch_level
        ;;print, branch_lvl

        mask = (*state).roi_mask and fib_c
        ;; Fiber passes through ANY ROIs the user wants?
        ;if ((*state).roi_mask and fib_c) eq 0 then continue
        ;; Fiber passes through ALL ROIs the user wants?
        ;if ((*state).roi_constrain) and (mask ne (*state).roi_mask) then continue

        if ((*state).roi_constrain) then begin
        
            if (mask ne (*state).roi_mask) then continue
            if (mask ne fib_c) then continue
            
        endif else begin
            
            if (mask eq 0) then continue
        
        endelse
        
        ;; choose a color
        ;;tmp = (fib_c + 3*branch_lvl) mod n_elements(colors)
        tmp = (fib_c) mod n_elements(colors)
        
        fiber_color = FSC_Color(colors[tmp], /triple)
        ;;fiber_color = fiber_color / (branch_lvl+1)

        pivot = (*state).fib_midpoints[j]
        fiblen = n_elements(*curr_fiber)/3
        
        low = (pivot - pivot * float((*state).fib_max_len)/100.0) > 0
        hi  = (pivot + (fiblen - pivot) * float((*state).fib_max_len)/100.0) < fiblen
        fiblen_new = (hi-low) > 1

        if (1 && fiblen_new lt fiblen) then begin
            fib_range = indgen(fiblen_new) + low
            fib_range = (fib_range > 0) < (fiblen - 1)
            fiber = (*curr_fiber)[*,[fib_range]]
        endif else begin
            fiber = *curr_fiber
        endelse
        
        fibs[n_fibs++] = obj_new('idlgrpolyline', fiber, $
                                 name='fiber_'+strcompress(string(j), /remove_all), $
                                 color=transpose(fiber_color), $
                                 alpha_channel=(*state).fib_alpha,           $
                                 depth_test_disable=0,      $
                                 thick=(*state).fib_thick)
    endfor

    ;; collapse array
    if (n_fibs gt 0) then begin
        ofibs = fibs[0:n_fibs-1]
        fibs  = 0
        omodel_fibs_chl->add, ofibs
    endif

    omodel_fibs_par->add, omodel_fibs_chl

end

;; Subroutine name: view_fibers
;;
;; Created by: BT 20080421
;; Calling Information:
;;   fiber_data    - output from make_fibers
;;   thick         - initial tract thickness (def: 1)
;;   alpha         - initial alpha for tracts(def: .75)
;;   show_axis     - whether or not to show  (def: 1)
;;   len_threshold - tract length threshold  (def: 20)
;;   fsteps        - fibers to skip (def 0)
;;   in_mas        - whether this is being called from within mas or
;;                   from the idl terminal
;;supplemental_omodel - a precreated omodel containing something that
;;                      should be displayed along with the fibers
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine: creates the visualization environment and
;; displays it
;;
;; Editing Information:
;;
pro view_fibers, fiber_data,         $ 
                   thick=thick,        $
                   alpha=alpha,        $
                   show_axis=show_axis,$
                   len_threshold=len_threshold, $
                   fstep=fstep, $
                   in_mas=in_mas, $
                   project_dataset=project_dataset, $
                   supplemental_omodel=supplemental_omodel
          
    common scan_data, project
    forward_function get_dir_slash

    if (keyword_set(in_mas)) then begin
        in_mas = 1
        rotate_Direction = project.procPramArray[project.ci].rotate_Direction
        flip_Direction   = project.procPramArray[project.ci].flip_Direction
        interp_factors   = double([ project.procPramArray[project.ci].freq_Interp, $
                                     project.procPramArray[project.ci].phase_Interp, $
                                     project.procPramArray[project.ci].slice_Interp ] )
    endif else begin
        in_mas = 0
        rotate_Direction = 0
        flip_Direction   = 0
        interp_factors   = [ 1D,1D,1D ]
    endelse
    
     ;; in mas, we have the imagery available
    if (in_mas) then begin
        image_datapool = project.dataarray[project.ci].frac_ani
    
        if (keyword_set(project_dataset)) then begin
            (*fiber_data).transform = diag_matrix(1D/[interp_factors, 1D]) # (*fiber_data).transform
            (*fiber_data).matrix = (size(*project.dataarray[project_dataset].state1, /dimensions))[0:2]
            image_datapool = project.dataarray[project_dataset].state1
        endif
    endif
    
    matrix   = (*fiber_data).matrix
    roi_file = (*fiber_data).filedir+(get_dir_slash())[1]+(*fiber_data).roifile

    ;; Sanity check for ROI file and read if available.
    looping = 0
    num_rois = 0
    ROI_FILE_TEST:
    roi_file_info = file_info(roi_file, /noexpand_path)

    if (not roi_file_info.exists) then begin
       ;; file doesn't exist
       yesno = 'No'
       if (looping eq 0) then begin
          ;; ask only once. if "yes" execution will start from
          ;; ROI_FILE_TEST with the newly user-selected filename.
          yesno = dialog_message(['ROI file could not be found:', $
                                  roi_file, $
                                  'Do you want to locate it?'], $
                                 /center, /question)
       endif

       if (yesno eq 'Yes') then begin
          roi_file = dialog_pickfile(path=project.current_path)
          looping = 1
          goto, ROI_FILE_TEST
       endif else begin
          ;; we need to do something here. We should find a way to
          ;; nullify the ROI display feature and just display the
          ;; streamlines in one color.
          void = dialog_message(['No valid ROI file could be loaded.', $
                                 'ROI coloring and information will be disabled.'], $
                                /error, /center)
       endelse
       
    endif else if (not roi_file_info.read) then begin
       ;; file exists but can't be read
       void = dialog_message(['ROI file is not readable:', $
                              roi_file, $
                              'ROI coloring and information will be diabled.'], $
                             /error, /center)
    endif else begin
       restore, roi_file
       num_rois = n_elements(roi_objects)
    endelse

    if not keyword_set(thick) then begin
        thick = 1.0
    endif

    if not keyword_set(alpha) then begin
        alpha = 0.75
    endif

    if not keyword_set(len_threshold) then begin
        len_threshold = 20
    endif

    if not keyword_set(fstep) then begin
        fstep = 1
    endif

    ;; Set up the various models
    omodel          = obj_new('idlgrmodel', name='MAIN_model', lighting=0)

    omodel_fibs_par = obj_new('idlgrmodel', name='FIB_parent_model', depth_test_disable=2, lighting=0, $
                               transform=(*fiber_data).transform)
    omodel_fibs_chl = obj_new('idlgrmodel', name='FIB_child_model', depth_test_disable=2, lighting=0)

    omodel_PS       = obj_new('idlgrmodel', name="PS_image", depth_test_disable=2, lighting=0)
    omodel_FP       = obj_new('idlgrmodel', name="FP_image", depth_test_disable=2, lighting=0)
    omodel_FS       = obj_new('idlgrmodel', name="FS_image", depth_test_disable=2, lighting=0)
    omodel_im       = obj_new('idlgrmodel', name='IMG_model', lighting=0)

    omodel_axes     = obj_new('idlgrmodel', name='AXES_model', lighting=0)

    oroi_model      = obj_new('idlgrmodel', name='ROI_model', depth_test_disable=2, lighting=1, $
                               transform=(*fiber_data).transform)
    
    ;; order is important
    omodel->add, omodel_fibs_par, position=0
    omodel->add, oroi_model,      position=1
    omodel->add, omodel_im,       position=2    
    omodel->add, omodel_axes,     position=3
    
    omodel_fibs_par->add, omodel_fibs_chl
    omodel_im->add, omodel_PS
    omodel_im->add, omodel_FP
    omodel_im->add, omodel_FS
    
    omodel_FS->rotate, [1,0,0], 90
    omodel_PS->rotate, [1,0,0], 90
    omodel_PS->rotate, [0,0,1], 90
    
    if keyword_set(supplemental_omodel) and obj_valid(supplemental_omodel) then begin
        omodel->add, supplemental_omodel, position=4
    endif

    ;; rotate the image model as necessary
    case rotate_direction of
        0: deg = 0
        1: begin
            deg = 90
            omodel_im->translate, -matrix[1]/2, -matrix[0]/2, 0
            omodel_im->rotate, [0,0,1], deg
            omodel_im->translate, matrix[0]/2, matrix[1]/2, 0
         end

        2: begin
            deg = 180
            omodel_im->translate, -matrix[0]/2, -matrix[1]/2, 0
            omodel_im->rotate, [0,0,1], deg
            omodel_im->translate, matrix[0]/2, matrix[1]/2, 0
        end

        3: begin
            deg = 270
            omodel_im->translate, -matrix[1]/2, -matrix[0]/2, 0
            omodel_im->rotate, [0,0,1], deg
            omodel_im->translate, matrix[0]/2, matrix[1]/2, 0
        end

    endcase

;;  model             - the main idlgrmodel object. contains
;;                      everything else
;;  tlb               - the tlb for the main viewport. returned by
;;                      xobjview
;;  matrix            - the dimensions (from mas) of the environment
;;  roi_on            - 1 ~ user want to see ROIs drawn, 0 ~ otherwise
;;  roi_count         - number of ROIs 
;;  roi_constrain     - 1 ~ user wants to see tracts running through
;;                      every visiable ROI, 0 ~ any ROI
;;  roi_identify      - 1 ~ user wants to see ROI numbers near ROIs, 0
;;                      ~ otherwise
;;  roi_mask          - n-bit vector (n = roi_count), where bit i=1 if
;;                      user wants to see ROI #i, zero otherwise. 
;;                      0<=i<=(roi_count-1)
;;  in_mas            - we're running from mas
;;  active_FP         - active freq-phase slice number
;;  active_FS         - active freq-slice slice number
;;  active_PS         - active phase-slice slice number
;;  image_alpha       - alpha percent for image display
;;  image_contrast    - contrast for image display
;;  image_datapool    - pointer to the data array that contains the
;;                      image data. don't free this pointer!
;;  fibers_front      - show fibers in front of imagery, otherwise
;;                      show imagery in front of fibers
;;  fib_rawdata       - raw fiber data (from make_fibers)
;;  fib_cache         - cache pointer
;;  fib_nfibers       - number of fibers
;;  fib_thick         - fiber thickness
;;  fib_alpha         - fiber alpha
;;  fib_fstep         - display increment
;;  fib_len_threshold - minimum length

    state = ptr_new({ model:omodel, $
                      tlb:0, $
                      matrix:matrix, $
                      transform: (*fiber_data).transform, $
                      roi_on:1,  $
                      roi_count:0, $
                      roi_constrain:0, $
                      roi_identify:1, $
                      roi_mask:0l, $
                      in_mas:in_mas, $
                      active_FP:-1, $
                      active_FS:-1, $
                      active_PS:-1, $
                      min_intensity: long(0), $
                      max_intensity: long(0), $
                      image_alpha: .75, $
                      image_contrast: 0, $
                      image_datapool: ptr_valid(image_datapool) ? image_datapool : ptr_new(), $
                      fibers_front:0, $
                      fib_rawdata:(*fiber_data).fibers, $
                      fib_midpoints: (*fiber_data).seed_points, $
                      fib_cache: ptrarr((*fiber_data).nfibers),$
                      fib_nfibers:(*fiber_data).nfibers, $
                      fib_branch_levels:(*fiber_data).branch_levels, $
                      fib_thick: thick, $
                      fib_alpha: alpha, $
                      fib_fstep: fstep, $
                      fib_max_len: 100, $
                      fib_len_threshold: len_threshold }, /no_copy)

            
    ptr_free, fiber_data

    ;; set up axis
    if keyword_set(show_axis) then begin
    
        oaxfnt = obj_new('idlgrfont', size=4)
        oxaxis = obj_new('idlgraxis', 0, range=[0,matrix[0]], color=[200,200,200])
        oyaxis = obj_new('idlgraxis', 1, range=[0,matrix[1]], color=[200,200,200])
        ozaxis = obj_new('idlgraxis', 2, range=[0,matrix[2]], color=[200,200,200])
        
        oxaxis->getProperty, ticktext=tmptxt
        tmptxt->setProperty, font=oaxfnt
        oyaxis->getProperty, ticktext=tmptxt
        tmptxt->setProperty, font=oaxfnt
        ozaxis->getProperty, ticktext=tmptxt
        tmptxt->setProperty, font=oaxfnt

        omodel_axes->add, oxaxis
        omodel_axes->add, oyaxis
        omodel_axes->add, ozaxis

    endif

    ;; set up ROI model
    if (num_rois gt -1) then begin

    colors = ['Khaki',      'Red',       'Green',        'Lime Green', 'Blue', 'Maroon', 'Cyan',  $
              'Magenta',    'Orange',    'Yellow',       'Orchid',     'Violet', $
              'Deep Pink',  'Tomato',    'Dark Green',   'Olive',$
              'Purple',     'Goldenrod', 'Saddle Brown', 'Coral',$
              'Brown',      'Navy',      'Light Cyan',   'Lavender',   'Green Yellow',$
              'Lawn Green', 'Honeydew']

        (*state).roi_count = ++num_rois

        orois = objarr(num_rois)
        oroi_lbls = objarr(num_rois)
        
        oroi_lbl_model  = obj_new('idlgrmodel', depth_test_disable=2, $
                                  name='ROI_LBL_model', lighting=0, hide=0)
        ofnt_lbl = obj_new('idlgrfont', size=9)

        orois[0] = obj_new('idlgrroi', $
                           name='ROI_0', $
                           depth_test_disable=2, $
                           [ [0,0,0] ], $
                           style=0, $
                           thick=4, $
                           type=2, $
                           color=[0,0,0], $
                           alpha_channel=0.0)
        oroi_lbls[0] = obj_new('idlgrtext', $
                               'B', $
                               locations=[0,0,0], $
                               /onglass,    $
                               font=ofnt_lbl,$
                               color=[255,255,255])
        
        vol = bytarr(matrix)
        
        for i = 1, num_rois-1 do begin 

            tmp = roi_objects[i-1]

            if (ptr_valid(tmp)) then begin
                
               ;; orois[i] = obj_new('idlgrroigroup', name=strcompress('ROI_'+string(i),/remove_all))
                
               ; orois[i] = obj_new('idlgrroi', $
               ;                    round(*tmp), $ ;;(*tmp)[0,*], (*tmp)[1,*], (*tmp)[2,*], $
               ;                    name=strcompress('ROI_'+string(i),/remove_all), $
               ;                    depth_test_disable=2, $
               ;                    style=0, $
               ;                    thick=4, $
               ;                    type=2, $
               ;                    color=[100, 100, 100], $
               ;                    alpha_channel=.6)

               ; junk = orois[i]->computeGeometry(centroid=c)
                c = reform((*tmp)[*,0])
                
                for f = 0L, n_elements(*tmp)/3 - 1 do begin
                    pt = floor((*tmp)[*,f])
                    vol[pt[0],pt[1],pt[2]] = 1B
                endfor
                roi_color = transpose(FSC_Color(colors[(i+2) mod n_elements(colors)], /triple))
                shade_volume, vol, 0, v, p, /low
                orois[i] = obj_new('idlgrpolygon', v, polygons=p, $
                                    color=roi_color, alpha=0.4, $
                                    name=strcompress('ROI_'+string(i),/remove_all))
                vol[*] = 0B
                oroi_lbls[i] = obj_new('idlgrtext', $
                                       string(i), $
                                       locations=c, $
                                       /onglass,    $
                                       font=ofnt_lbl,$
                                       color=[255,255,255])
            endif

            ptr_free, tmp

        endfor

        oroi_lbl_model->add, oroi_lbls
        oroi_model->add, [orois, oroi_lbl_model]

        ;; start with all rois visible
        (*state).roi_mask = long(2^num_rois - 1)

        if (n_elements(ptr_seeds) ne 0) then begin
           ptr_free, ptr_seeds[0]
           ptr_free, ptr_seeds
        endif

    endif

    (*state).max_intensity = max(*image_datapool, min=min_int)
    (*state).min_intensity = min_int
    
    ;; build up the stuff
    rebuild_fiber_cache, state
    rebuild_fiber_model, state

    xobjview, omodel, background=[0,0,0], tlb=tlb, xsize=800, ysize=600, renderer=1
    
    (*state).tlb = tlb

    ;; this blocks, so we cleanup after the user closes window
    view_fibers_GUI, state=state

    ;; cleanup
;    obj_destroy, omodel
;    if keyword_set(show_axis) then obj_destroy, oaxfnt
;    ptr_free, (*state).fib_rawdata
;    ptr_free, (*state).fib_cache
;    ptr_free, state

end
