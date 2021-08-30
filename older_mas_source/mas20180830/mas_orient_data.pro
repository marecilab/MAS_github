;; $Id$
;;

;; The purpose of these routines is to take data originating from
;; where ever (PAR/REC, fid, DICOM, etc.) and determine the anatomical
;; orientation.
;;
;; After state1 has been loaded, the data is in a 3 or 4 dimensional
;; array, depening on the type of scan. The first three indices
;; [i,j,k] represent a coordinate system that is relative to an
;; anatomical coordinate system [x,y,z] by some possibly unknown
;; transformation. 
;;
;; In the above case, the x-, y-, and z-axes may coorespond to
;; some permutation of (Left-Right, Anterior-Posterior,
;; Inferior-Superior). 
;;
;; If the transformation can be deduced from the parameters provided
;; by the file, then our job here is to provide a way to map from
;; anatomical coordinates (x,y,z) to in-data indices (i,j,k), and to
;; provide a consistent way to label image displays.
;;
;; If the transformation is unknown, then we should be able to allow
;; the user to orient the image based on their knowledge of the
;; acquisition process.
;;
;; If neither of the above are possible, then no adjustments should be
;; made, and the data will be displayed in the (i,j,k) coordinate system.
;;

;; I_axis_min = L ( => R); def: Imin ( => Imax)
;; J_axis_min = I ( => S); def: Jmin ( => Jmax)
;; K_axis_min = A ( => P); def: Kmin ( => Kmax)

;; project.procpramarray[project.ci].orientation_Imin
;; project.procpramarray[project.ci].orientation_Jmin
;; project.procpramarray[project.ci].orientation_Kmin

;+
; :Description:
;    This is a "constant" function that returns the
;    axes vectors. They are ordered to correspond to
;    the numeric code.
;
; :Author: btt
;-
function mas_orient_get_axes_vecs

    axes_vecs = [ [  0, -1,  0 ], $  ;; A
                  [  1,  0,  0 ], $  ;; L -
                  [  0,  0,  1 ], $  ;; I
                  [  0,  1,  0 ], $  ;; P -
                  [ -1,  0,  0 ], $  ;; R
                  [  0,  0, -1 ] ]   ;; S -
                  
   return, transpose(axes_vecs)

end

;+
; :Description:
;    This is a "constant" function that returns
;    the axes letters corresponding to the orientation
;    numeric code.
;
; :Author: btt
;-
function mas_orient_get_axes_ltrs

    axes_ltrs = ['P', 'R', 'S', 'A', 'L', 'I']
    
    return, axes_ltrs

end

;;; Not used
;function mas_orient_get_anat_indices, $
;   transverse=transverse, $
;   coronal=coronal,       $
;   sagittal=sagittal
;
;    common scan_data
;
;    ci = project.ci
;    
;    ijk = [project.procpramarray[ci].orientation_Imin, $
;           project.procpramarray[ci].orientation_jmin, $
;           project.procpramarray[ci].orientation_Kmin]
;    
;    if (keyword_set(coronal)) then begin
;        temp = ijk eq 1 or ijk eq 4 or ijk eq 2 or ijk eq 5
;     endif else if (keyword_set(transverse)) then begin
;        temp = ijk eq 1 or ijk eq 4 or ijk eq 0 or ijk eq 3
;     endif else if (keyword_set(sagittal)) then begin
;        temp = ijk eq 2 or ijk eq 5 or ijk eq 0 or ijk eq 3
;     endif else begin
;        ;; at least on direction must be specified
;     endelse
;     
;     indices = where(temp ne 0)
;
;     return, indices
;     
;end

;+
; :Description:
;    Constant function mapping the proper name for
;    an axis to its numeric code..
;    
;    Given an axis as a keyword, function will return
;    the numeric code which can be used elsewhere.
;
; :Keywords: One of these must be specified.
;    anterior
;    superior
;    posterior
;    left
;    right
;    inferior
;
; :Author: btt
;-
function mas_orient_get_id, anterior=anterior, $
                            superior=superior, $
                            posterior=posterior, $
                            left=left, $
                            right=right, $
                            inferior=inferior


  if keyword_set(anterior)  then return, 0
  if keyword_set(left)      then return, 1
  if keyword_set(inferior)  then return, 2
  if keyword_set(posterior) then return, 3
  if keyword_set(right)     then return, 4
  if keyword_set(superior)  then return, 5

  return, -1
  
end

;+
; :Description:
;    Given one numeric code for an axis, this function
;    will return the numeric code for its opposite.
;    
;    That is, given the code for "left", function will return
;    code for "right", etc.
;    
; :Params:
;    direction - the orientation code for which to obtain
;                the opposite
;
; :Author: btt
;-
function mas_orient_get_opposite, direction

    return, (direction + 3) mod 6

end

;+
; :Description:
;    Given a numeric code, this function will return
;    the axis letter or, optionally, the entire word.
;
; :Params:
;    id - numeric code to look up
;
; :Keywords:
;    longname - set this keyword to get the whole word
;               instead of the first letter
;
; :Author: btt
;-
function mas_orient_get_label_from_id, id, longname=longname

   ;; for opposite label, id_opp = (id + 3) mod 6
    case id of
       0: return, keyword_set(longname) ? 'Anterior'  : 'A'
       1: return, keyword_set(longname) ? 'Left'      : 'L'
       2: return, keyword_set(longname) ? 'Inferior'  : 'I'
       3: return, keyword_set(longname) ? 'Posterior' : 'P'
       4: return, keyword_set(longname) ? 'Right'     : 'R'
       5: return, keyword_set(longname) ? 'Superior'  : 'S'
    endcase

end

function mas_orient_get_axes_array, orientation, labels=labels

    ;; 0 => Anterior  => [  0  1  0 ]
    ;; 1 => Left      => [ -1  0  0 ]
    ;; 2 => Inferior  => [  0  0  1 ]
    ;; 3 => Posterior => [  0 -1  0 ]
    ;; 4 => Right     => [  1  0  0 ]
    ;; 5 => Superior  => [  0  0 -1 ]
    
    axes = [0,1,2,3,4,5]

    if (arg_present(labels)) then begin
        labels = strarr(n_elements(axes))
        for i = 0, n_elements(axes)-1 do begin
            labels[i] = mas_orient_get_label_from_id(axes[i])
        endfor
    endif
    
    return, axes

 end

;+
; :Description:
;    Given a 3x3 orientation matrix, such as from a DICOM file
;    or a NIFTI sform, this function generates the internal
;    numeric codes for the orientation or the standard letter
;    codes for the orientation.
;    
;    For example, given
;     (R->L)   (A->P)  (S->I)
;      -1       0       0
;       0       1       0
;       0       0       1
;    
;    This function returns: 4 3 2 (see above), or L A S
;                           
; :Params:
;    matrix - a 3x3 orientation matrix 
;
; :Keywords:
;    lettercode - set this keyword to return the letters 
;                 rather than the numeric code.
;
; :Author: btt
;-
function mas_orient_from_matrix, matrix, lettercode=lettercode

    axes_vecs = mas_orient_get_axes_vecs()
    axes_ltrs = mas_orient_get_axes_ltrs()


    orient = intarr(3)

    for j = 0, 2 do begin
       for i = 0, 5 do begin
          tot = total(abs(axes_vecs[i,*] - transpose(matrix[j,*])))
          if (tot eq 0.0) then begin
             orient[j] = i
             break
          endif
       endfor
    endfor
    
    if (keyword_set(lettercode)) then begin
        return, [axes_ltrs[orient[0]], axes_ltrs[orient[1]], axes_ltrs[orient[2]]]
    endif else begin
        return, orient
    endelse
end

;+
; :Description:
;    Given a three element array containing the orientation
;    numeric code, or a three character string containing
;    the axes code ('LAS', 'RPI', etc). this function will
;    return a 3x3 matrix corresponding to that code.
;
; :Params:
;    ijkmin - the three element array containing the numeric codes
;             if this is not specified, then the current scan's
;             codes will be used.
;
; :Keywords:
;    lettercode - set this keyword to a string containing the
;                 orientation code to get a matrix for. EG:
;                 
;     IDL> print, mas_orient_get_matrix(lettercode='RAS')
;           1       0       0
;           0       1       0
;           0       0       1             
;
; :Author: btt
;-
function mas_orient_get_matrix, ijkmin, lettercode=lettercode

    common scan_data

    axes_vecs = mas_orient_get_axes_vecs()
    axes_ltrs = mas_orient_get_axes_ltrs()
    
;    axes_vecs = transpose(axes_vecs)

    if (keyword_set(lettercode)) then begin
        ijkmin = bytarr(3)
        tmp = string(transpose(byte(lettercode)))
        for i = 0, n_elements(tmp)-1 do begin
            ijkmin[i] = where(tmp[i] eq axes_ltrs, count)
            if (count eq 0) then begin
                print, "Orientation couldn't be computed!"
            endif
        endfor
    endif
    if (n_elements(ijkmin) eq 0) then begin
       ci = project.ci
       ijkmin = [ project.procpramarray[ci].orientation_Imin, $
                  project.procpramarray[ci].orientation_Jmin, $
                  project.procpramarray[ci].orientation_Kmin]
       ;print, ijkmin
    endif

    
    return,  [ axes_vecs[ijkmin[0],*], $
               axes_vecs[ijkmin[1],*], $
               axes_vecs[ijkmin[2],*] ]
end


;+
; :Description:
;    This function is used by MAS to generate single letter
;    labels for image axes.
;    
;    It returns a four-element array containing the labels
;    for the image plane in this order:
;      TOP, RIGHT, BOTTOM, LEFT
;
;    Set the keyword ORIENT_STRING to a variable that will
;    contain the orientation name of the current slice.
;    
;    If the keyword APPLY_FLIP_ROTATE is set, the results
;    will automatically be adjusted to match any display
;    rotating or flipping set in the mas display panel.
;    
; :Params:
;    slice_axis - the axis through which the data is sliced:
;                 0 = "freq/phase"
;                 1 = "freq/slice"
;                 2 = "phase/slice"
;
; :Keywords:
;    apply_flip_rotate - set this to have the results adjusted 
;                        for rotation / flipping
;    orient_string - set this to a variable to retrieve the proper
;                    name for the image's orientation
;
; :Author: btt
;-
function mas_orient_get_axis_labels, slice_axis, $
                                     apply_flip_rotate=apply_flip_rotate, $
                                     orient_string=orient_string

    common scan_data

    ci = project.ci

    orient_string = 'UNKNOWN'
    
    ;; if we have data that is not oriented, then just return empty strings
    if project.procpramarray[ci].have_orientation eq 0 then return, strarr(4)

    if (n_elements(slice_axis) eq 0) then begin

       slice_axis = project.procpramarray[ci].slice_axis

    end else if slice_axis lt 0 or slice_axis gt 2 then begin

       ;; error - slice axis out of range

    endif

    Imin = project.procpramarray[project.ci].orientation_Imin
    Imax = mas_orient_get_opposite(Imin)
    Jmin = project.procpramarray[project.ci].orientation_Jmin
    Jmax = mas_orient_get_opposite(Jmin)
    Kmin = project.procpramarray[project.ci].orientation_Kmin
    Kmax = mas_orient_get_opposite(Kmin)

    ;; proceeds clockwise
    case slice_axis of 

       0: begin ;; "Freq/Phase" 
          
          axis = [Jmax, Imax, Jmin, Imin]
          
       end

       1: begin ;; "Freq/Slice"

          axis = [Kmax, Imax, Kmin, Imin]

       end

       2: begin

          axis = [Kmax, Jmax, Kmin, Jmin]

       end
       
    endcase
       
    axis_str = strarr(4)
    for i = 0,3 do begin
       axis_str[i] = mas_orient_get_label_from_id(axis[i])
    endfor
    axis = temporary(axis_str)
    
    ;; determine the name of the view. For example, sagittal views
    ;; contain sup, inf, ant, post.
    sorted_axis = strjoin(axis[sort(axis)], '')
    case sorted_axis of 
        'AIPS' : orient_string = 'SAGITTAL'
        'ILRS' : orient_string = 'CORONAL'
        'ALPR' : orient_string = 'TRANSVERSE'
        else   : orient_string = 'UNKNOWN'
    endcase

    if (not keyword_set(apply_flip_rotate)) then begin
       
       return, axis

    endif

    rotate_dir = project.procpramarray[ci].rotate_direction

    case rotate_dir of
       
       0 : begin ;; No rotation
       end
       
       1: begin ;; 90 degrees
          axis = shift(axis, -1)
       end

       2: begin ;; 180 degrees
          axis = shift(axis, 2)
       end

       3: begin ;; 270 degrees
          axis = shift(axis, -3)
       end

    endcase

    flip_dir = project.procpramarray[ci].flip_direction

    case flip_dir of

       0: begin ;; No Flip
       end

       1: begin ;; Horizontal
          axis = [axis[0], axis[3], axis[2], axis[1]]
       end

       2: begin ;; Vertical
          axis = [axis[2], axis[1], axis[0], axis[3]]
       end
       
    endcase


    return, axis

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;+
; :Description:
;    Event handler for selecting a "Slice Axis" in the "Attach Orientation"
;    window (Display menu).
;
; :Params:
;    event - event generated by selection widget
;   
; :Author: btt
;-
pro mas_orient_slice_axis_event, event

    common scan_data

    ci = project.ci

    widget_control, event.top, get_uvalue=state
    
    ;; we want to adjust the display of the slice for apparent voxel size
    data_dim = size(*project.dataarray[ci].state1, /dimensions)
    voxel_dim = [ project.imndarray[ci].f_voxsz / project.procpramarray[ci].freq_interp, $
                  project.imndarray[ci].p_voxsz / project.procpramarray[ci].phase_interp, $
                  project.imndarray[ci].s_voxsz / project.procpramarray[ci].slice_interp ]
    voxel_dim/=min(voxel_dim)

    case event.index of

       0: begin
          (*state).slice_axis = 0
          dim = [ data_dim[0] * voxel_dim[0], data_dim[1] * voxel_dim[1] ]
          slice_dim = data_dim[2]
          xmin = (*state).Imin
          ymin = (*state).Jmin
       end
       1: begin
          (*state).slice_axis = 1
          dim = [ data_dim[0] * voxel_dim[0], data_dim[2] * voxel_dim[2] ]
          slice_dim = data_dim[1]
          xmin = (*state).Imin
          ymin = (*state).Kmin
       end
       2: begin
          (*state).slice_axis = 2
          dim = [ data_dim[1] * voxel_dim[1], data_dim[2] * voxel_dim[2] ]
          slice_dim = data_dim[0]
          xmin = (*state).Jmin
          ymin = (*state).Kmin
       end

    endcase

    slice = fix(slice_dim/2.) > 0 ;; choose the center slice by default

    widget_control, (*state).image_window, scr_xsize=dim[0], scr_ysize=dim[1]
    widget_control, (*state).image_slider, set_slider_min=0, set_slider_max=slice_dim-1, set_value=slice
    ;; generating this event manually sets the proper image int he viewer.
    mas_orient_image_select_event, { id:(*state).image_slider, top:event.top, handler:0L, value:slice }

    orients = mas_orient_get_axes_array(event.index, labels=labels)
    
    widget_control, (*state).wid_TOP,    set_value=labels, set_droplist_select=mas_orient_get_opposite(ymin)
    widget_control, (*state).wid_RIGHT,  set_value=labels, set_droplist_select=mas_orient_get_opposite(xmin)
    widget_control, (*state).wid_BOTTOM, set_value=labels, set_droplist_select=ymin
    widget_control, (*state).wid_LEFT,   set_value=labels, set_droplist_select=xmin
    
end

;+
; :Description:
;    Event handler that gets and displays the slice data 
;    in the window (Display->Attach Orientation).
;
; :Params:
;    event
;
; :Author: btt
;-
pro mas_orient_image_select_event, event

    common scan_data
    ci = project.ci

    widget_control, event.top, get_uvalue=state
    widget_control, event.id, get_value=val
   
    data_ptr = project.dataarray[project.ci].state1
    voxel_dim = [ project.imndarray[ci].f_voxsz / project.procpramarray[ci].freq_interp, $
                  project.imndarray[ci].p_voxsz / project.procpramarray[ci].phase_interp, $
                  project.imndarray[ci].s_voxsz / project.procpramarray[ci].slice_interp ]
    voxel_dim/=min(voxel_dim)

    data_dim = size(*data_ptr, /dimensions)

    val = fix(val) > 0
    case (*state).slice_axis of
       0: img = congrid(reform((*data_ptr)[*,*,val,0]), $
                        data_dim[0]*voxel_dim[0], data_dim[1]*voxel_dim[1],$
                        interp=0)
       1: img = congrid(reform((*data_ptr)[*,val,*,0]), $
                        data_dim[0]*voxel_dim[0], data_dim[2]*voxel_dim[2],$
                        interp=0)
       2: img = congrid(reform((*data_ptr)[val,*,*,0]), $
                        data_dim[1]*voxel_dim[1], data_dim[2]*voxel_dim[2],$
                        interp=0)
    endcase

    image = bytscl(reform(img))
    dim = size(image, /dimensions)

    widget_control, (*state).image_window, get_value=win

    wset, win
    tv, image

end

;+
; :Description:
;    Event handler for various buttons.
;
; :Params:
;    event
;
; :Author: btt
;-
pro mas_orient_button_event, event

    common scan_data

    ci = project.ci
    
    widget_control, event.top, get_uvalue=state

    name = widget_info(event.id, /uname)

    case name of 
       'btn_apply': begin
       
          testmatrix = mas_orient_get_matrix([(*state).Imin,(*state).Jmin,(*state).Kmin])
          if (determ(testmatrix, /check) eq 0) then begin
            void = dialog_message(["Invalid orientation", $
                                   "Please make sure that all axes have been oriented",$
                                   "and that there are no duplicate axis orientations."],$
                                    /error, /center) 
            return
          endif
          project.procpramarray[ci].orientation_Imin = (*state).Imin
          project.procpramarray[ci].orientation_Jmin = (*state).Jmin
          project.procpramarray[ci].orientation_Kmin = (*state).Kmin
          project.procpramarray[ci].have_orientation = 1
          orient_file = project.imndarray[ci].file_path+'.orientation.dat'
          write_ok = 1
          if (file_test(orient_file)) then begin
             if (not file_test(orient_file, /write)) then begin
                write_ok = 0
             endif
          endif else if (not file_test(file_dirname(orient_file), /write)) then begin
             write_ok = 0
          endif
          if (write_ok eq 1) then begin
             openw, lun, orient_file, /get_lun
             printf, lun, [ project.procpramarray[ci].orientation_Imin, $
                            project.procpramarray[ci].orientation_Jmin, $
                            project.procpramarray[ci].orientation_Kmin ]
             close, lun
             free_lun, lun
          endif else begin
             junk = dialog_message('Could not write orientation file.', /error, /center)
          endelse

       end

       'btn_cancel': begin

       end
       
       else: print, 'mas_orient_button_event: recd event from unknown widget: '+name
    endcase

    ptr_free, state

    widget_control, event.top, /destroy
    
end

;+
; :Description:
;    Event handler for events generated by changing the drop list
;    to choose the orientation for a single axis. This automatically
;    updates the opposite axis.
;
; :Params:
;    event
;
;
;
; :Author: btt
;-
pro mas_orient_data_event, event

    widget_control, event.top, get_uvalue=state

    ref = mas_orient_get_axes_array((*state).orientation)

    TOP    = widget_info((*state).wid_TOP, /droplist_select) 
    RIGHT  = widget_info((*state).wid_RIGHT, /droplist_select)
    BOTTOM = widget_info((*state).wid_BOTTOM, /droplist_select)
    LEFT   = widget_info((*state).wid_LEFT, /droplist_select)

    case event.id of

       (*state).wid_TOP: begin 
          index = event.index
          other = RIGHT

          TOP    = index
          BOTTOM = ( index + 3 ) mod 6

          if (other eq index or other eq (index+3) mod 6) then begin
             LEFT  = ( index + 1 ) mod 6
             RIGHT = ( index + 4 ) mod 6
          endif

       end

       (*state).wid_LEFT: begin
          index = event.index
          other = TOP

          LEFT  = index
          RIGHT = ( index + 3 ) mod 6

          if (other eq index or other eq (index+3) mod 6) then begin
             TOP    = ( index + 1 ) mod 6
             BOTTOM = ( index + 4 ) mod 6
          endif

       end

       else: return

    endcase

    (*state).axes = [ref[TOP],   $
                     ref[RIGHT], $
                     ref[BOTTOM],$
                     ref[LEFT]]

    case (*state).slice_axis of

       0: begin
          (*state).Imin = ref[LEFT]
          (*state).Jmin = ref[BOTTOM]
       end
       1: begin
          (*state).Imin = ref[LEFT]
          (*state).Kmin = ref[BOTTOM]
       end
       2: begin
          (*state).Jmin = ref[LEFT]
          (*state).Kmin = ref[BOTTOM]
       end

    endcase

    widget_control, (*state).wid_TOP, set_droplist_select = TOP
    widget_control, (*state).wid_RIGHT, set_droplist_select = RIGHT
    widget_control, (*state).wid_BOTTOM, set_droplist_select = BOTTOM
    widget_control, (*state).wid_LEFT, set_droplist_select = LEFT

 end


;+
; :Description:
;    Main GUI setup for the data orientation setting. Found 
;    in Display->Attach Orientation.
;
; :Author: btt
;-
pro mas_orient_data

    common scan_data
    common common_widgets

    ci = project.ci
    if (not ptr_valid(project.dataarray[ci].state1)) then return

    data_dims = size(*project.dataarray[ci].state1, /dimensions)

    base = widget_base(title='Orient Data', /column, group_leader=WID_BASE_MAIN, /modal)

    base_2col = widget_base(base, row=2)
    
    base_leftcol = widget_base(base_2col, /column)

    base_rightcol = widget_base(base_2col, /column)

    initial_orientation = 1

    Imin = project.procpramarray[ci].orientation_Imin
    Imax = mas_orient_get_opposite(Imin)
    Jmin = project.procpramarray[ci].orientation_Jmin
    Jmax = mas_orient_get_opposite(Jmin)
    Kmin = project.procpramarray[ci].orientation_Kmin
    Kmax = mas_orient_get_opposite(Kmin)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    b = widget_base(base_leftcol, /frame, /row)
    l = widget_label(b, value='Slice Axis:')
    axis_select = widget_droplist(b, value=['Freq/Phase', 'Freq/Slice', 'Phase/Slice'], $
                                  uname='axis_select', event_pro='mas_orient_slice_axis_event')

    b = widget_base(base_leftcol, /frame, /row)
    l = widget_label(b, value="Display Slice:")
    slice_slider = widget_slider(b, min=0, max=data_dims[2]-1, value=(fix(data_dims[2]/2.)-1 > 0), $
                                 uname='slice_slider', event_pro='mas_orient_image_select_event')

    b = widget_base(base_leftcol, /frame, /column, /align_center)
    
    orients = mas_orient_get_axes_array(initial_orientation, labels=labels)

    TOP_select = widget_droplist(b, value=labels, /align_center, $
                                 uname='TOP_select')

    widget_control, TOP_select, set_droplist_select=Jmax

    bb = widget_base(b, /row, /align_center)
    LEFT_select = widget_droplist(bb, value=labels, /align_center, $
                                  uname='LEFT_select')

    widget_control, LEFT_select, set_droplist_select=Imin
    
    drawwin = widget_draw(bb, scr_xsize=data_dims[0], scr_ysize=data_dims[1], /align_center, $
                          graphics_level=1)

    RIGHT_select = widget_droplist(bb, value=labels, /align_center, $
                                   uname='RIGHT_select', sensitive=0)

    widget_control, RIGHT_select, set_droplist_select=Imax

    BOTTOM_select = widget_droplist(b, value=labels, /align_center, $
                                    uname='BOTTOM_select', sensitive=0)

    widget_control, BOTTOM_select, set_droplist_select=Jmin
    
    btn_base = widget_base(base, /row, /align_center)
    btn_apply  = widget_button(btn_base, value="Apply Orientation", uname='btn_apply', $
                              event_pro='mas_orient_button_event')
    btn_cancel = widget_button(btn_base, value="Discard Changes", uname='btn_cancel', $
                              event_pro='mas_orient_button_event')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    state = ptr_new({ slice_axis: 0, $
                      current_slice: 50, $ ;slice, $
                      orientation: initial_orientation, $
                      Imin: Imin, $
                      Jmin: Jmin, $
                      Kmin: Kmin, $
                      axes:orients, $
                      image_window: drawwin, $
                      image_slider: slice_slider, $
                      wid_TOP: TOP_select, $
                      wid_RIGHT: RIGHT_select, $
                      wid_BOTTOM: BOTTOM_select, $
                      wid_LEFT: LEFT_select}) ;, $
    

    widget_control, base, set_uvalue=state

    widget_control, base, /realize

    mas_orient_image_select_event, { id:(*state).image_slider, top:base, handler:0L, value:(*state).current_slice }

    xmanager, 'mas_orient_data',  base, /no_block

end

pro mas_orient_compute_transform, from=from, to=to, lettercode=lettercode

    


end
