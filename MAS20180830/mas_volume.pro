;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved


; Subroutine name: mas_slicer_3d
; Created by
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_slicer_3d
    ;print, 'slicer 3d'
    COMMON scan_data
    CI = project.ci

    ;if the data is not arranged in a 3d format then make the data in a 3d format
    if project.procPramArray[CI].single_Multi_Flag  eq  0 then begin
       project.procPramArray[CI].single_Multi_Flag = 1
       project.procPramArray[CI].state_2 = 0
    end

    ;slicer needs a 3d data set.
    mas_load_state_2
    mas_redraw_gui

    ;; Save the current color table because
    ;; slicer3d is going to change it
    curr_ct = 0
    tvlct, curr_ct, /get

    if (size(*project.dataArray[CI].state2))[0] eq 3 then  $
       SLICER3, project.dataArray[CI].state2 $
    else update_status_bar,'Slicer 3D unable to process data'

    ;; restore previously saved color table
    tvlct, curr_ct

end


; Subroutine name: mas_ivolume
; Created by
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Launches the 3D IDL slicer with the loaded data.

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_ivolume

    COMMON scan_data
    CI = project.ci

    ;if the data is not arranged in a 3d format then make the data in a 3d format
    if project.procPramArray[CI].single_Multi_Flag  eq  0 then begin
       project.procPramArray[CI].single_Multi_Flag = 1
       project.procPramArray[CI].state_2 = 0
    end

    ;slicer needs a 3d data set.
    mas_load_state_2
    mas_redraw_gui

    print, size(project.dataArray[CI].state2)


    ivolume, *project.dataArray[CI].state2

;    if (size(*project.dataArray[CI].state2))[0] eq 3 then  $
;
;    else update_status_bar,'Slicer 3D unable to process data'


end

; Subroutine name: mas_volume_slicer
; Created by: Bill T.
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Displays a 3D dataset in an enhanced
;     slicer view

; Editing Information:

pro mas_volume_slicer

    common scan_data

    ci = project.ci

    interp_f = project.procpramarray[ci].freq_interp
    interp_p = project.procpramarray[ci].phase_interp
    interp_s = project.procpramarray[ci].slice_interp

    x_zoom = project.procpramarray[ci].x_zoom
    y_zoom = project.procpramarray[ci].y_zoom

    voxdims = [project.imndarray[ci].f_voxsz/interp_f/x_zoom,$
               project.imndarray[ci].p_voxsz/interp_p/y_zoom,$
               project.imndarray[ci].s_voxsz/interp_s]
    voxdims /= min(voxdims)
    voxsz_f = voxdims[0]
    voxsz_p = voxdims[1]
    voxsz_s = voxdims[2]

    if (project.imndarray[ci].adim gt 1) then begin

        if (project.procpramarray[ci].single_multi_flag eq 1) then begin

            mas_load_state_2
            mas_3d_slicer, project.dataarray[ci].state2, f_scale=voxsz_f, p_scale=voxsz_p, s_scale=voxsz_s

        endif else begin

            ;; we need to get the volume data into state_2
            project.procpramarray.single_multi_flag = 1
            project.procpramarray[ci].state_2 = 0

            mas_load_state_2
            mas_3d_slicer, project.dataarray[ci].state2, f_scale=voxsz_f, p_scale=voxsz_p, s_scale=voxsz_s

            project.procpramarray.single_multi_flag = 0
            project.procpramarray[ci].state_2 = 0
            mas_redraw_GUI

        endelse

    endif else begin
        
        mas_3d_slicer, project.dataarray[ci].state1, f_scale=voxsz_f, p_scale=voxsz_p, s_scale=voxsz_s

    endelse

end
; Subroutine name: mas_measure_volume_event
; Created by
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting


PRO mas_measure_volume_event, Event
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets

    CI = project.CI
    case Event.id of

       ;store the values returned from the slider widgets
        vol_start: project.procPramArray[CI].vol_start_slice = event.value
        vol_stop : project.procPramArray[CI].vol_stop_slice = event.value

       ;this will handle the pressing of the go button on the volumetrics main widget.
        vol_draw: begin
         ;if the start gt stop dont let the user start drawing roi to measure volume.
         if  project.procPramArray[CI].vol_start_slice gt project.procPramArray[CI].vol_stop_slice then begin
          void = dialog_message(['Start slice is greater than stop slice, please correct.','Aborting volumetrics processing'], /error, /center)
          return
         end
         n_pixels = 0
         mas_load_state_1

         for ii=project.procPramArray[CI].vol_start_slice-1 , project.procPramArray[CI].vol_stop_slice-1 do begin

          ;prep the image for displaying in the roi tool
          image = reform((*project.dataArray[ci].state1)[*,*,ii,0])
          mas_windowing, image
          pImage = ptr_new( image, /no_copy)
            mas_zoom, pImage

          ;let the user draw an roi
          xroi, *pImage , regions_out=regions, /block, /modal, group = event.top

          ;we need to chk to make sure they drawed an roi.
             ;if they drew nothing then skip to the next slice quietly.
             ;chk to make sure the object regions is valid
             if not(total(obj_valid(regions))) then begin
                continue
             end else begin

              if  (size(regions))[0] eq 0 then numRoi = 1 $
              else numRoi =  (size(regions))[1]

              ; make a pointer array to hold all the masks
              mask = ptrarr(numRoi)

              ;extract the masks


              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;what if the roi's  over lap each other
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              for temp=0 , numRoi-1 do begin
                 ;extract the mask where the inside the roi
                  mask[temp] = ptr_new(regions[temp] -> ComputeMask( dimensions =[ (size((*pimage)))[1], (size((*pimage)))[2] ] , MASK_RULE=2))

                 ;count the number of pixels in the mask
                  IMAGE_STATISTICS , *pimage, mask=*(mask[temp]) , count=roi_pixels, /LABELED

                 ;cont then number of pixels and add them to the total num of pixels.
                 n_pixels +=   roi_pixels[(size(roi_pixels))[1]-1]
                 print, n_pixels

              endfor

          end

         endfor

          ;how much volume is one pixel.
          ;volume = [(xZoom*xInterpolation*Ffov/Fdim) * (yZoom* yInterpolation*Pfov/pDim) * (zInterpolation*Sfov/sDim) ] *num_pixels

          x_zoom = float(project.procPramArray[CI].x_zoom)
          freq_interp= float(project.procPramArray[CI].freq_interp)
          fdim = float(project.imndarray[ci].fdim)
          f_fov = float(project.imndarray[ci].f_fov)

          xdim = f_fov/( fdim * x_zoom * freq_interp )
          ;print, 'xdim',xdim


          y_zoom= float(project.procPramArray[CI].y_zoom)
          Phase_interp= float(project.procPramArray[CI].Phase_interp)
          pdim = float(project.imndarray[ci].pdim)
          p_fov = float(project.imndarray[ci].p_fov)

          ydim = p_fov/( pdim * y_zoom * Phase_interp)
          ;print, 'ydim',ydim

          Slice_interp= float(project.procPramArray[CI].Slice_interp)
          sdim = float(project.imndarray[ci].sdim)
          s_fov = float(project.imndarray[ci].s_fov)

          zdim = s_fov/( sdim * slice_interp )
          ;print, 'zdim',zdim


          pixel_volume = xdim * ydim * zdim
          ;print, 'pixel_volume',pixel_volume

          roi_volume = pixel_volume * n_pixels

          ;print, 'Volume = '+strtrim(roi_volume,2)
          stemp = strarr(2)
          stemp[0] = 'Volume'
          stemp[1] = strtrim(roi_volume,2)
          display_stats, stemp, 'Volumetrics'


       end
    end
end


; Subroutine name: mas_measure_volume
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This procedure will create the base for all volume measuring widgets will be attached.

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_measure_volume
    HEAP_GC
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    CI = project.CI

    if project.scan_Open_Flag eq 0 then begin
       void = dialog_message(['Please open a scan before using Volumetrics tool.','Aborting Volumetrics tool.'], /error, /center)
       return
    end

    IF N_ELEMENTS(vol_window_base) EQ 1 THEN $
    IF     widget_info( vol_window_base,  /VALID_ID ) eq 1 THEN return

    title = string ('Volumetrics measurement: Scan #', strtrim(CI,2))
    vol_window_base = widget_base(TITLE=title, SCR_XSIZE=200 , XOFFSET=420 ,YOFFSET=0, COLUMN=1  )

    vol_start = widget_slider(vol_window_base, MINIMUM=1, FRAME=0, MAXIMUM=project.imndarray[ci].sdim , $
                            TITLE = 'Start'  , $
                            VALUE= project.procPramArray[CI].vol_start_slice )

    vol_stop = widget_slider(vol_window_base, MINIMUM=1, FRAME=0, MAXIMUM=project.imndarray[ci].sdim , $
                            TITLE = 'Stop' ,  SENSITIVE=SENSITIVE , $
                            VALUE= project.procPramArray[CI].vol_stop_slice )

    vol_draw = widget_button(vol_window_base,value='Draw' ,SENSITIVE=SENSITIVE )

    widget_control, vol_window_base, /realize

    xmanager, 'mas_measure_volume',vol_window_base,/no_block, GROUP_LEADER=WID_BASE_MAIN

end


pro mas_volume
end
