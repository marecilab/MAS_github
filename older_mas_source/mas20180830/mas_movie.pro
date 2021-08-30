;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

;*****************************************************************************************************
;
; NAME:
;   mas_movie_se
;
; PURPOSE:
;
;
; ARGUMENTS:
;
;
; MODIFICATION HISTORY:
;
;
;*****************************************************************************************************
pro mas_movie_se
    COMPILE_OPT IDL2

    COMMON scan_data
    ci = project.ci

    sz_image = size(*project.dataArray[ci].signal_Enhancement)
    x_zoom = project.procpramArray[ci].x_zoom
    y_zoom = project.procpramarray[ci].y_zoom
    rotate_Direction = project.procPramArray[project.ci].rotate_Direction


    if rotate_Direction eq 0 or rotate_Direction eq 2 then begin
        sx = x_zoom * sz_image[1]
        sy = y_zoom * sz_image[2]
    end else begin
        sx = x_zoom * sz_image[2]
        sy = y_zoom * sz_image[1]

    end


    if sz_image[0] eq 2 then update_status_bar , 'only 1 frame in movie'
    if sz_image[0] eq 3 then begin


       base = WIDGET_BASE(TITLE = 'Animation Widget')
       animate = CW_ANIMATE(base, sx, sy, sz_image[3], /NO_KILL )
       WIDGET_CONTROL, /REALIZE, base
       FOR ii=0,sz_image[3]-1 DO begin
         p_data = ptr_new((*project.dataArray[ci].signal_Enhancement)[*,*,ii])
         mas_zoom, p_data
         mas_rotate_flip, p_data
         CW_ANIMATE_LOAD, animate, FRAME=ii, IMAGE = (*p_data)
       end

       update_status_bar , ' '

       ;WIDGET_CONTROL, /REALIZE, base

       CW_ANIMATE_GETP, animate, pixmap_vect
       CW_ANIMATE_RUN, animate, 25
       XMANAGER, 'CW_ANIMATE Demo', base, EVENT_HANDLER = 'EHANDLER', /NO_BLOCK


    end



end


;*****************************************************************************************************
;
; NAME:
;   mas_movie
;
; PURPOSE:
;
;
; ARGUMENTS:
;
;
; MODIFICATION HISTORY:
;   2008-03-11 - added /NO_BLOCK to XMANAGER call. 
;
;
;*****************************************************************************************************
pro mas_movie
    heap_gc

    COMPILE_OPT IDL2
    COMMON scan_data
    ci = project.ci

    mas_load_state_2

    ;max_intensity = max(*(project.dataArray[ci].state2))
    ;min_intensity = min(*(project.dataArray[ci].state2))
    min_intensity = project.procPramArray[ci].min_display_thresh
    max_intensity = project.procPramArray[ci].max_display_thresh

    image = 0
    sz_image = size(*project.dataArray[ci].state2)

    if sz_image[0] eq 2 then update_status_bar , 'only 1 frame in movie'
    if sz_image[0] eq 3 then begin

       base = WIDGET_BASE(TITLE = 'Animation Widget', column=1)
       animate = CW_ANIMATE(base, sz_image[1], sz_image[2], sz_image[3], /NO_KILL )
       junk = widget_text(base, value='', /editable)
    
; Creating a progressbar
progressbar = Obj_New('progressbar', Color='red', Text='Loading Frames',/NOCANCEL)
progressbar -> Start
progressbartemp = float(1.0/sz_image[3])
progressbar -> Update, progressbartemp*100

       ;WIDGET_CONTROL, /REALIZE, base
       FOR ii=0,sz_image[3]-1 DO begin
         image = (*project.dataArray[ci].state2)[*,*,ii]
         mas_windowing,image ,max_intensity, min_intensity

         CW_ANIMATE_LOAD, animate, FRAME=ii, IMAGE = image

        progressbar_load = progressbartemp*ii
         progressbar -> Update, progressbar_load*100
       end
       update_status_bar , ' '

		progressbar -> Destroy
		WIDGET_CONTROL, /REALIZE, base

       CW_ANIMATE_GETP, animate, pixmap_vect
       ; 25 % run speed
       CW_ANIMATE_RUN, animate, 25
       XMANAGER, 'CW_ANIMATE Demo', base, EVENT_HANDLER = 'EHANDLER', /NO_BLOCK


    end

end

; Subroutine name: mas_multi_movie
; Created by: Garrett Astary
; Calling Information: Called by mas_multi_movie_gui_event when user selects Play button or Export button

; Bugs or Important Comments to Developers:


; Purpose of subroutine: Creates multi slice movie and then either a)loads the movie using CW animate for display
;                        on the screen or b)exports the movie as a series of .tif files where each .tif is a frame in the movie 

; Editing Information:

pro mas_multi_movie, params
    heap_gc

    COMPILE_OPT IDL2
    COMMON scan_data
    ci = project.ci
    
    ;Get chosen slice information from info array
    info = params.info
    ;We need to check to see how many slices the user wants to veiw
    chosen = where(info ge 0)
    num_slices = (size(chosen))[1] 
    chosen_slices = info[chosen]
    
    ;In order to get the array dim for multiple slices we will have to change the sdim_start selection such that
    ;mas_windowing organizes the appropriate slice, we want to store the original sdim_start and change it back when 
    ;we finish making the movie
    sdim_start = project.procPramArray[ci].sdim_start
    x_fov = project.imndArray[project.ci].p_fov
    y_fov = project.imndArray[project.ci].f_fov
    
    if num_slices lt 1 then return
    ;Multi Slice Movie should not be sensitive if adim = 0, included a check here in case
    if project.imndarray[ci].adim eq 0 then return
    ;If the user wants to create a multi slice movie, we know they want the scan oriented in the Array dim
    ;and multiple slices. We temporarily set these flags and reverse them at the end of the .pro
    sort_dir = project.procPramArray[CI].sort_dir
    project.procPramArray[CI].sort_dir = 1
    
    single_Multi_flag = project.procPramArray[project.ci].single_Multi_flag
    project.procPramArray[project.ci].single_Multi_flag = 1
    
    ;The interp and zoom factors are used when creating the .tif files
    ;We want to ensure these are correct regardless of the user rotating the image
    rotate_direction = project.procPramArray[project.ci].rotate_direction
    if rotate_direction eq 0 or 2 then begin
      x_interp = project.procPramArray[project.ci].freq_interp
      x_zoom = project.procPramArray[project.ci].x_zoom
      y_interp = project.procPramArray[project.ci].phase_interp
      y_zoom = project.procPramArray[project.ci].y_zoom
    endif else begin
      x_interp = project.procPramArray[project.ci].phase_interp
      x_zoom = project.procPramArray[project.ci].y_zoom
      y_interp = project.procPramArray[project.ci].freq_interp
      y_zoom = project.procPramArray[project.ci].x_zoom
    endelse
    
    ;import the first chosen slice, get dimensions from this slice
    ;dimensions will be constant for each slice
    mas_load_state_1
    project.procPramArray[ci].sdim_start = chosen_slices[0]
    temp = mas_subset(project.dataArray[ci].state1)
    mas_zoom, temp
    mas_rotate_flip, temp
    sz_image = size(*temp)
    xsize = sz_image[1]
    ysize = sz_image[2]
    zsize = sz_image[3]
    slice1 = (*temp)
    
    ;Array that will contain the read,phase,array information for all chosen slices
    stored_slices = dblarr(xsize,ysize,zsize,num_slices)
    stored_slices[*,*,*,0] = slice1
    
    ;If the number of chosen slices is greater than 1, than import them all into the stored_slices array now
    if num_slices gt 1 then begin
      for i=1, num_slices-1 do begin
        mas_load_state_1
        project.procPramArray[ci].sdim_start = chosen_slices[i]
        temp = mas_subset(project.dataArray[ci].state1)
        mas_zoom, temp
        mas_rotate_flip,temp
        new_slice = (*temp)
        stored_slices[*,*,*,i] = new_slice
      endfor
    endif
    
    ;The next set of if statements will determine how many slices the user wants to show and then create a composite image
    ;based on this selection. The xfactor and yfactor variables approriately size the animation frame by increasing the xdim 
    ;or ydim depending on how the composite layout is constructed
    if num_slices eq 1 then begin
      composite = dblarr(xsize, ysize, zsize)
      composite[*,*,*] = stored_slices[*,*,*,0]
      xfactor = 1
      yfactor = 1
    endif else if num_slices eq 2 then begin
      
      composite = dblarr(2*xsize, ysize, zsize)
      
      for i=0, zsize-1 do begin
      composite[0:xsize-1,*, i] = stored_slices[*,*,i,0]
      composite[xsize:2*xsize-1,*, i] = stored_slices[*,*,i,1]
      endfor
      
      xfactor = 2
      yfactor = 1
    endif else if num_slices eq 3 then begin
    
      composite = dblarr(3*xsize, ysize, zsize)
    
      for i=0, zsize-1 do begin
        composite[0:xsize-1,*, i] = stored_slices[*,*,i,0]
        composite[xsize:2*xsize-1,*, i] = stored_slices[*,*,i,1]
        composite[2*xsize:3*xsize-1,*, i] = stored_slices[*,*,i,2]
      endfor
      
      xfactor = 3
      yfactor = 1
    endif else if num_slices eq 4 then begin
      
      composite = dblarr(2*xsize, 2*ysize, zsize)
      for i=0, zsize-1 do begin
        ;interestingly enough, the CW_ANIMATE function loads frames such that left-most voxel is 0 and
        ;the bottom-most voxel is zero. Thus, to load a slice in the top right quadrant for the movie we want to load
        ;it into the bottom right quadrant of the array
        composite[0:xsize-1,ysize:2*ysize-1, i] = reform(stored_slices[*,*,i,0])
        composite[xsize:2*xsize-1,ysize:2*ysize-1, i] = reform(stored_slices[*,*,i,1])
        composite[0:xsize-1,0:ysize-1, i] = reform(stored_slices[*,*,i,2])
        composite[xsize:2*xsize-1,0:ysize-1, i] = reform(stored_slices[*,*,i,3])
      endfor
      
      xfactor = 2
      yfactor = 2
    endif else if num_slices eq 5 then begin
      
      composite = dblarr(3*xsize, 2*ysize, zsize)
      
      for i=0, zsize-1 do begin
        composite[0:xsize-1,ysize:2*ysize-1, i] = stored_slices[*,*,i,0]
        composite[xsize:2*xsize-1,ysize:2*ysize-1, i] = stored_slices[*,*,i,1]
        composite[2*xsize:3*xsize-1,ysize:2*ysize-1, i] = stored_slices[*,*,i,2]
        composite[0:xsize-1,0:ysize-1, i] = stored_slices[*,*,i,3]
        composite[xsize:2*xsize-1,0:ysize-1, i] = stored_slices[*,*,i,4]
      endfor
      xfactor =3
      yfactor = 2
    endif else begin
      composite = dblarr(3*xsize, 2*ysize, zsize)
      
      for i=0, zsize-1 do begin
        composite[0:xsize-1,ysize:2*ysize-1, i] = stored_slices[*,*,i,0]
        composite[xsize:2*xsize-1,ysize:2*ysize-1, i] = stored_slices[*,*,i,1]
        composite[2*xsize:3*xsize-1,ysize:2*ysize-1, i] = stored_slices[*,*,i,2]
        composite[0:xsize-1,0:ysize-1, i] = stored_slices[*,*,i,3]
        composite[xsize:2*xsize-1,0:ysize-1, i] = stored_slices[*,*,i,4]
        composite[2*xsize:3*xsize-1,0:ysize-1, i] = stored_slices[*,*,i,5]
      endfor
      xfactor =3
      yfactor = 2
    endelse
    
    
    
    
    ;byte scaling the time-stack of images
    max_intensity = max(composite)
    min_intensity = min(composite)
    
    
    ;Check to see if user wants to view the movie
    if params.export_flag eq 0 then begin
      base = WIDGET_BASE(TITLE = 'Animation Widget', column=1)
      animate = CW_ANIMATE(base, xfactor*sz_image[1], yfactor*sz_image[2], sz_image[3], /NO_KILL )
      junk = widget_text(base, value='', /editable)
  
      ; Creating a progressbar
      progressbar = Obj_New('progressbar', Color='red', Text='Loading Frames',/NOCANCEL)
      progressbar -> Start
      progressbartemp = float(1.0/sz_image[3])
      progressbar -> Update, progressbartemp*100
  
  
        
         FOR ii=0,sz_image[3]-1 DO begin
           image = composite[*,*,ii]
           mas_windowing,image, .75*max_intensity, min_intensity 
           CW_ANIMATE_LOAD, animate, FRAME=ii, IMAGE = image
  
          progressbar_load = progressbartemp*ii
           progressbar -> Update, progressbar_load*100
         end
         update_status_bar , ' '
  
      progressbar -> Destroy
      WIDGET_CONTROL, /REALIZE, base
  
         CW_ANIMATE_GETP, animate, pixmap_vect
         ; 25 % run speed
         CW_ANIMATE_RUN, animate, 25
         XMANAGER, 'CW_ANIMATE Demo', base, EVENT_HANDLER = 'EHANDLER', /NO_BLOCK
       endif else begin
         ;user selects directory to store .tif files, they as x.tif where x is consecutive numbers starting at 1
         fpath = dialog_pickfile(/directory, path=project.current_path)
         XRESOL = xfactor*xsize/(xfactor*x_interp*x_zoom*x_fov * (1/2.54) )
         YRESOL = yfactor*ysize/(yfactor*y_interp*y_zoom*y_fov * (1/2.54) )
         composite = reverse(composite,2)
         for ii=0,sz_image[3]-1 do begin
              image = composite[*,*,ii]
              mas_windowing, image, .75*max_intensity, min_intensity 
              tifname = fpath+strcompress(string(ii+1),/remove_all)+'.tif'
              print, tifname
              WRITE_TIFF, tifname, image, XRESOL=XRESOL ,YRESOL=YRESOL , orientation = 1
              
         endfor
       endelse
;reset the initial user settings and redraw MAS to reflect these
project.procPramArray[ci].sdim_start = sdim_start
project.procPramArray[CI].sort_dir = sort_dir
project.procPramArray[project.ci].single_Multi_flag = single_Multi_flag
mas_redraw_GUI
end

; Subroutine name: mas_multi_movie_gui_event
; Created by: Garrett Astary
; Calling Information: Event handler for mas_multi_movie_gui

; Bugs or Important Comments to Developers:


; Purpose of subroutine: Handles events from mas_multi_movie_gui, stores slices selections in info array located in the params
;                        structure. 

; Editing Information:

pro mas_multi_movie_gui_event, event

 ;bring in global variables
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    
    ;grab parameter structure that is stored in the uvalue of the mas_multi_movie_gui
    widget_control, event.top, get_uvalue=params, /no_copy
    info = params.info
    CASE event.id OF

      ; input for R1 value of contrast agent, value stored in info structure
      msm_slice_1_droplist: BEGIN
        slice1 = widget_info(msm_slice_1_droplist, /DROPLIST_SELECT) - 1
        info[0] = slice1
      END
      
       msm_slice_2_droplist: BEGIN
        slice2 = widget_info(msm_slice_2_droplist, /DROPLIST_SELECT) - 1
        info[1] = slice2
      END
      
       msm_slice_3_droplist: BEGIN
       slice3 = widget_info(msm_slice_3_droplist, /DROPLIST_SELECT) - 1
        info[2] = slice3
      END
      
       msm_slice_4_droplist: BEGIN
       slice4 = widget_info(msm_slice_4_droplist, /DROPLIST_SELECT) - 1
       info[3] = slice4
      END
      
       msm_slice_5_droplist: BEGIN
       slice5 = widget_info(msm_slice_5_droplist, /DROPLIST_SELECT) - 1
       info[4] = slice5
      END
      
       msm_slice_6_droplist: BEGIN
       slice6 = widget_info(msm_slice_6_droplist, /DROPLIST_SELECT) - 1
       info[5] = slice6
      END
      
      msm_play_button: BEGIN
      ;get information from each dropdown, if user did not change a dropdown menu the info array will not be updated
      ;we need to grab the droplist selection from each box and store it in info
      slice1 = widget_info(msm_slice_1_droplist, /DROPLIST_SELECT) - 1
      info[0] = slice1
      slice2 = widget_info(msm_slice_2_droplist, /DROPLIST_SELECT) - 1
      info[1] = slice2
      slice3 = widget_info(msm_slice_3_droplist, /DROPLIST_SELECT) - 1
      info[2] = slice3
      slice4 = widget_info(msm_slice_4_droplist, /DROPLIST_SELECT) - 1
      info[3] = slice4
      slice5 = widget_info(msm_slice_5_droplist, /DROPLIST_SELECT) - 1
      info[4] = slice5
      slice6 = widget_info(msm_slice_6_droplist, /DROPLIST_SELECT) - 1
      info[5] = slice6
      test = where(info ge 0)
      params.info = info
      ;if the user hits play without selecting any slices, we just return
      if (size(test))[0] eq 0 then begin
        widget_control, event.top, set_uvalue=params  
        return
      endif
      ;user wants to play movie
      params.export_flag = 0
      mas_multi_movie, params
      END
      
      msm_export_button: BEGIN
      slice1 = widget_info(msm_slice_1_droplist, /DROPLIST_SELECT) - 1
      info[0] = slice1
      slice2 = widget_info(msm_slice_2_droplist, /DROPLIST_SELECT) - 1
      info[1] = slice2
      slice3 = widget_info(msm_slice_3_droplist, /DROPLIST_SELECT) - 1
      info[2] = slice3
      slice4 = widget_info(msm_slice_4_droplist, /DROPLIST_SELECT) - 1
      info[3] = slice4
      slice5 = widget_info(msm_slice_5_droplist, /DROPLIST_SELECT) - 1
      info[4] = slice5
      slice6 = widget_info(msm_slice_6_droplist, /DROPLIST_SELECT) - 1
      info[5] = slice6
      test = where(info ge 0)
      params.info = info
      if (size(test))[0] eq 0 then begin
        widget_control, event.top, set_uvalue=params  
        return
      endif
      params.export_flag = 1
      mas_multi_movie, params
      END
      
    ENDCASE
;after updating the info array after a user selection we want to place it back into the uvalue of the mas_multi_movie_gui
params.info = info    
widget_control, event.top, set_uvalue=params    
end

; Subroutine name: mas_multi_movie_gui
; Created by: Garrett Astary
; Calling Information: Called from mas.pro when user selects Multi Slice Movie (widget button uname = W_MENU_multi_slice_movie)  
;                      from Display dropdown

; Bugs or Important Comments to Developers:


; Purpose of subroutine: Creates widget base and widgets that allow user to choose which slices (up to 6 total) to add
;                        to the multislice movie

; Editing Information:
pro mas_multi_movie_gui
    HEAP_GC
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    CI = project.CI
    sdim = project.imndArray[CI].sdim
    
    msm_window_base = widget_base(column=1, TITLE='Multi Slice Movie', xoffset=420, /base_align_center)
    msm_label = widget_label(msm_window_base, value='Select Slices to add to Movie:')
    
    msm_slice_base = widget_base(msm_window_base, column = 2)
    
    msm_slice1_base = widget_base(msm_slice_base, column = 2)
    msm_slice_1_label = widget_label(msm_slice1_base, value='Slice 1:')
    ;creates array [0,....,sdim-1] to be included in droplist
    slices = indgen(sdim)
    slices = strcompress(string(slices),/remove_all)
    ;droplist now looks like [None, 0, ....,sdim-1]
    choices = ['None', slices]
    msm_slice_1_droplist =  widget_droplist(msm_slice1_base, value=choices, uname = 'msm_slice_1_droplist')
    
    msm_slice2_base = widget_base(msm_slice_base, column = 2)
    msm_slice_2_label = widget_label(msm_slice2_base, value='Slice 2:')
    msm_slice_2_droplist =  widget_droplist(msm_slice2_base, value=choices, uname = 'msm_slice_2_droplist')
    
    msm_slice3_base = widget_base(msm_slice_base, column = 2)
    msm_slice_3_label = widget_label(msm_slice3_base, value='Slice 3:')
    msm_slice_3_droplist =  widget_droplist(msm_slice3_base, value=choices, uname = 'msm_slice_3_droplist')
    
    msm_slice4_base = widget_base(msm_slice_base, column = 2)
    msm_slice_4_label = widget_label(msm_slice4_base, value='Slice 4:')
    msm_slice_4_droplist =  widget_droplist(msm_slice4_base, value=choices, uname = 'msm_slice_4_droplist')
    
    msm_slice5_base = widget_base(msm_slice_base, column = 2)
    msm_slice_5_label = widget_label(msm_slice5_base, value='Slice 5:')
    msm_slice_5_droplist =  widget_droplist(msm_slice5_base, value=choices, uname = 'msm_slice_5_droplist')
    
    msm_slice6_base = widget_base(msm_slice_base, column = 2)
    msm_slice_6_label = widget_label(msm_slice6_base, value='Slice 6:')
    msm_slice_6_droplist =  widget_droplist(msm_slice6_base, value=choices, uname = 'msm_slice_6_droplist')
    
    ;After choosing desired slices, user can either view the movie or export as a series of .tif files for creating a movie
    ;in another application (Quicktime Pro, etc.)
    msm_action_base = widget_base(msm_window_base, column = 2)
    msm_play_button = widget_button(msm_action_base, value = 'Play')
    msm_export_button = widget_button(msm_action_base, value = 'Export to .tif')

    widget_control, msm_window_base, /realize
    ;Parameter structure
    ;params = {info:    6 element vector array containing the selections from the 6 dropdown menus
    ;                   Each element is equal to the (dropdown selection - 1)
    ;                   Ex:  If 'None' is chosen by dropdown 4 then info[3] = -1
    ;                        If  4 is chosen for dropdown 6 then info[5] = 3
    ;          export_flag: 0 - not exporting, user wants to view movie
    ;                       1 - will export movie as a series of consecutively numbered .tif files              
    params = {info:intarr(6), export_flag:0}
    
    widget_control, msm_window_base, set_uvalue=params, /no_copy

    xmanager, 'mas_multi_movie_gui', msm_window_base

end