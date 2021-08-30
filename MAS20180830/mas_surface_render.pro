; Subroutine name: mas_surface_render
; Created by: Kevin M.
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Displays a 3D dataset in an enhanced 
;     surface view by calling the mas_3d_surface function

; Editing Information:

pro mas_surface_render

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
            mas_3d_surface, project.dataarray[ci].state2, f_scale=voxsz_f, p_scale=voxsz_p, s_scale=voxsz_s

        endif else begin

            ;; we need to get the volume data into state_2
            project.procpramarray.single_multi_flag = 1
            project.procpramarray[ci].state_2 = 0

            mas_load_state_2
            mas_3d_surface, project.dataarray[ci].state2, f_scale=voxsz_f, p_scale=voxsz_p, s_scale=voxsz_s

            project.procpramarray.single_multi_flag = 0
            project.procpramarray[ci].state_2 = 0
            mas_redraw_GUI

        endelse

    endif else begin
        
        mas_3d_surface, project.dataarray[ci].state1, f_scale=voxsz_f, p_scale=voxsz_p, s_scale=voxsz_s

    endelse

end
