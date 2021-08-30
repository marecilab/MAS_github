;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

; Subroutine name: mas_zoom
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_zoom, p_data, scale_voxels=scale_voxels
    COMPILE_OPT IDL2
    COMMON scan_data

    ci = project.ci

    if (keyword_set(scale_voxels)) then begin

       slice_axis = project.procPramarray[ci].slice_axis
       
       voxel_dim = [ project.imndarray[ci].f_voxsz/project.procpramarray[ci].freq_interp , $
                     project.imndarray[ci].p_voxsz/project.procpramarray[ci].phase_interp, $
                     project.imndarray[ci].s_voxsz/project.procpramarray[ci].slice_interp ]
       voxel_dim /= min(voxel_dim)
       
       case slice_axis of
          
          0: begin
             x_voxdim = voxel_dim[0]
             y_voxdim = voxel_dim[1]
          end
          
          1: begin
             x_voxdim = voxel_dim[0]
             y_voxdim = voxel_dim[2]
          end
          
          2: begin
             x_voxdim = voxel_dim[1]
             y_voxdim = voxel_dim[2]
          end
          
       endcase
    endif else begin
       x_voxdim = 1.
       y_voxdim = 1.
    endelse

    ;; bring in the data that is necessary.
    x_zoom = project.procPramArray[ci].x_zoom * x_voxdim
    y_zoom = project.procPramArray[ci].y_zoom * y_voxdim

    sz_scan       = size(*p_data)

    new_sx = sz_scan[1] * x_zoom
    new_sy = sz_scan[2] * y_zoom

    update_status_bar, 'Zooming'

    if sz_scan[0] eq 2 then begin
        ;;(*p_DATA) = REBIN((*p_DATA), new_sx, new_sy, /SAMPLE )
        (*p_DATA) = CONGRID((*p_DATA), new_sx, new_sy, INTERP=0)
    end
    if sz_scan[0] eq 3 then begin
       ;;(*p_DATA) = REBIN((*p_DATA), new_sx, new_sy, sz_scan[3],
       ;;/SAMPLE )
       tmp = fltarr(new_sx, new_sy, sz_scan[3])
       for sl = 0, sz_scan[3]-1 do begin
          tmp[*,*,sl] = CONGRID(reform((*p_DATA)[*,*,sl]), new_sx, new_sy, INTERP=0)
       endfor
       ptr_free, p_data
       p_data = ptr_new(tmp, /no_copy)
    end

    update_status_bar, ''


end
