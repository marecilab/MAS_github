

;; $Id$
;;
; Subroutine name: adt_show_gradient_cleanup
; Created by:
; Calling Information:
;
;   tlb: top-level base
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; 
;    Cleans up pointers and objects associated with the gradient direction viewer
;
; Editing Information:

pro adt_show_gradient_cleanup, tlb

    widget_control, tlb, get_uvalue=uvalue

    obj_destroy, (*uvalue).otxtmodel
    obj_destroy, (*uvalue).overts
    obj_destroy, (*uvalue).osurf_hi
    if (obj_valid((*uvalue).osurf_lo)) then obj_destroy, (*uvalue).osurf_lo
    obj_destroy, (*uvalue).pts
    ptr_free, uvalue

end

; Subroutine name: adt_show_gradient_angle_event
; Created by:
; Calling Information:
;
;   ev: the button press event
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; 
;    Handles click of "Hide/Show Angles" button on the 3D viewer top 
;    bar.
;
; Editing Information:

pro adt_show_gradient_angle_event, ev

    widget_control, ev.id, get_uvalue=uvalue
    (*uvalue).otxtmodel->getProperty, hide=isHidden
    (*uvalue).otxtmodel->setProperty, hide=(1 - isHidden)

    if (isHidden eq 0) then begin
       widget_control, ev.id, set_value='Show Angles'
    endif else begin
       widget_control, ev.id, set_value='Hide Angles'
    endelse

    xobjview, refresh=ev.top

end

; Subroutine name: adt_show_gradient_profile
; Created by:
; Calling Information:
;
;   GRADIENTS - user-supplied gradients in substitute for the
;               current mas scan's gradients
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:
  ; Edited by mkulam 
  ; 1) Estimate for threshold for high and low b-values changed from mu(b) to mu(b)-sigma(b) 
  ;    treating the low b-values as outliers in the normal b-value distribution
  ; 2) Included the effect of selecting array data for the fit
  
pro adt_show_gradient_profile, gradients=gradients

    common scan_data

    if (keyword_set(gradients)) then begin

       theta = ptr_new(gradients[1,*])
       phi   = ptr_new(gradients[0,*])

    endif else begin

       theta = project.imndarray[project.ci].angle_theta
       phi   = project.imndarray[project.ci].angle_phi

       if not (ptr_valid(theta) and ptr_valid(phi)) then begin
          junk = dialog_message(['Gradient profile can only be shown for diffusion scans.'], $
                                /error, /center)
          return
       endif

    endelse

    bvals = project.imndarray[project.ci].bval_array
    
    if ptr_valid(project.procpramArray[project.ci].array_select) then begin       
      bvals = ptr_new((*bvals)[*project.procpramArray[project.ci].array_select])
      theta = ptr_new((*theta)[*project.procpramArray[project.ci].array_select])
      phi = ptr_new((*phi)[*project.procpramArray[project.ci].array_select])   
    endif  
            
    bthres = mean(*bvals)-stddev(*bvals)

    high_b = where(*bvals gt bthres, count_high)
    low_b  = where(*bvals le bthres, count_low)
    total_b = count_high + count_low
    if (count_high le 0) then begin
        junk = dialog_message(['No non-zero bvalues found. Aborting.'], $
                              /error, /center)
        return
    endif else if (count_high eq 1) then begin
        junk = dialog_message(['Diffusion scan must have more than', $
                               'gradient direction. Aborting.'], $
                                /error, /center)
        return
    endif

    gr           = fltarr(3, n_elements(*bvals))
    gr[0,high_b] = (*phi)[high_b]
    gr[1,high_b] = 90.0 - (*theta)[high_b]
    gr[2,high_b] = 5

    gr[0,low_b]  = (*phi)[low_b]
    gr[1,low_b]  = 90.0 - (*theta)[low_b]
    gr[2,low_b]  = 3.5

    gr_actual           = fltarr(3, n_elements(*bvals))
    gr_actual[0,high_b] = (*phi)[high_b]
    gr_actual[1,high_b] =  (*theta)[high_b]
    gr_actual[2,high_b] = 1
    gr_actual[0,low_b]  = (*phi)[low_b]
    gr_actual[1,low_b]  = (*theta)[low_b]
    gr_actual[2,low_b]  = 1

    gr_rect = cv_coord(from_sphere=gr, /to_rect, /degrees)
    gr_rect_act = cv_coord(from_sphere=gr_actual, /to_rect, /degrees)

    catch, qhull_error1
    if (qhull_error1 ne 0) then begin
        catch, /cancel
        void = dialog_message(['Unable to create visualization.', $
                               'There must be at least three unique gradient directions.'], $
                               /error, /center)
        return
    endif
    
    oldverts = gr_rect[*,high_b]
    qhull, oldverts, tetrahedra, /delaunay
    newconn=tetra_surface(oldverts, tetrahedra)
    osurf_hi = obj_new('IDLgrPolygon', oldverts, polygons=newconn, style=1, Color=[255,0,0])
    ;;osurf_hi = obj_new('IDLgrPolygon', oldverts, alpha_channel=0.75, polygons=newconn, style=2, Color=[128,0,0])
    
    catch, /cancel
    
    catch, qhull_error
    if (qhull_error ne 0) then begin
        catch, /cancel
        goto, SKIP_QHULL
    endif
    
    if (count_low ge 2) then begin
        oldverts = gr_rect[*,low_b]
        qhull, oldverts, tetrahedra, /delaunay
        newconn=tetra_surface(oldverts, tetrahedra)
        osurf_low = obj_new('IDLgrPolygon', oldverts, polygons=newconn, style=1, Color=[0,0,255])
        ;;osurf_low = obj_new('IDLgrPolygon', oldverts, polygons=newconn, style=2, Color=[0,0,128])
    endif

    catch, /cancel
    SKIP_QHULL:
    
    overts = objarr(total_b)

    pts = obj_new('idlgrpolygon', gr_rect, thick=4, style=0, depth_test_disable=2)
    otxt_model = obj_new('idlgrmodel')
    ofnt = obj_new('idlgrfont', size=8)

    for i = 0, total_b-1 do begin
        overts[i] = obj_new('idlgrpolyline', [ [0,0,0], [gr_rect[*,i]] ], thick=1, $
                          name=strcompress('('+string(gr[0,i])+','+string(gr[1,i])+')', /remove_all), $
                          linestyle=0, color=[90,90,90], alpha_channel=0.2, depth_test_disable=2)

        otxt_model->add, obj_new('idlgrtext', location=[gr_rect[*,i]+0.05], color=[0,0,0], /onglass, font=ofnt,$
                                 alpha_channel=0.6, $
                                 strcompress(string(gr_actual[1,i], format='(3F0.3)')+','+$
                                             string(gr_actual[0,i], format='(3F0.3)'), /remove_all))

    endfor

    omodel = obj_new('idlgrmodel', depth_test_disable=2)

    if (count_low ge 2 and obj_valid(osurf_low)) then omodel->add, osurf_low
    omodel->add, osurf_hi
    
    omodel->add, overts
    omodel->add, pts
    omodel->add, otxt_model

    xobjview, omodel, xsize=600, ysize=600, tlb=tlb

    tb = widget_info(tlb, find_by_uname="xobjview:_resetbase")

    if (widget_info(tb, /valid_id)) then begin

       btn = widget_button(tb, $
                           value='Hide Angles', $
                           uname='ANGLE_BTN', $
                           event_pro='adt_show_gradient_angle_event', $
                           uvalue=ptr_new({ otxtmodel: otxt_model, $
                                            overts: overts, $
                                            osurf_hi: osurf_hi, $
                                            osurf_lo: obj_valid(osurf_low) ? osurf_low : 0, $
                                            pts: pts }) )
       xmanager, 'adt_show_gradient_angle', btn, cleanup='adt_show_gradient_cleanup', /no_block

    endif

end


; Subroutine name: get_orient_rgb_axis
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
;
; This function will map the RGB orientation maps to the color scheme
; from Sinisa Pajevic and Carlo Pierpaoli
;
; Editing Information:

FUNCTION get_orient_rgb_axis, tx_matrix=tx_matrix, slice_axis=slice_axis

    COMMON scan_data
    forward_function mas_orient_get_axis_labels
    
    CI = project.ci

    ident_4    = diag_matrix(replicate(1.0D,4))
    tx_matrix  = ident_4; dblarr(4,4)

    if (arg_present(slice_axis)) then begin
        slice_axis = (fix(slice_axis) > 0) < 2
    endif else begin
        slice_axis = project.procpramArray[CI].slice_axis
    endelse

    orientation = mas_orient_get_axis_labels(0, orient_string=orient_str)

    ;; if set to Arbitrary or Oblique then assume Transverse Colors
    ;CASE strupcase(project.imndArray[CI].orientation[0]) OF
    case orient_str of
        'TRANSVERSE': begin;; OK
            axis_orient = [0,1,2]
            ;; orientation[0] indicates the image-top anatomical orientation
            if (orientation[0] eq 'L' or orientation[0] eq 'R') then begin
                ;; y-axis is left-right
                axis_orient[0] = 1 & axis_orient[1] = 0
            endif
            case slice_axis of
                0: begin
                end
                1: begin
                    tx_matrix[*,0] = ident_4[0, *]
                    tx_matrix[*,1] = ident_4[2, *]
                    tx_matrix[*,2] = ident_4[1, *]
                    axis_orient = [axis_orient[0], axis_orient[2], axis_orient[1]]
                end
                2: begin
                    tx_matrix[*,0] = ident_4[1, *]
                    tx_matrix[*,1] = ident_4[2, *]
                    tx_matrix[*,2] = ident_4[0, *]
                    axis_orient = [axis_orient[1], axis_orient[2], axis_orient[0]]
                end
            endcase
        end

        'CORONAL'   : begin
            axis_orient = [0,2,1] ;; [0,1,2]
            ;; orientation[0] indicates the image-top anatomical orientation
            if (orientation[0] eq 'L' or orientation[0] eq 'R') then begin
                ;; y-axis is left-right
                axis_orient[0] = 2 & axis_orient[1] = 0
            endif            
            case slice_axis of
                0: begin
                end
                1: begin
                    axis_orient = [axis_orient[0], axis_orient[2], axis_orient[1]]
                    tx_matrix[*,0] = ident_4[0, *]
                    tx_matrix[*,1] = ident_4[2, *]
                    tx_matrix[*,2] = ident_4[1, *]
                end
                2: begin
                    axis_orient = [axis_orient[1], axis_orient[2], axis_orient[0]]
                    tx_matrix[*,0] = ident_4[1, *]
                    tx_matrix[*,1] = ident_4[2, *]
                    tx_matrix[*,2] = ident_4[0, *]
                end
            endcase
        end

        'SAGITTAL'  : begin
            axis_orient = [1, 2, 0]; [0,1,2] ;;;;;;[1,2,0] ; [2,1,0]; [1,0,2]; ([2,0,1] orig)
            ;; orientation[0] indicates the image-top anatomical orientation
            if (orientation[0] eq 'A' or orientation[0] eq 'P') then begin
                ;; y-axis is anterior-posterior
                axis_orient[0] = 2 & axis_orient[1] = 1
            endif
            case slice_axis of
                0: begin
                end
                1: begin
                    axis_orient = [axis_orient[0], axis_orient[2], axis_orient[1]]
                    tx_matrix[*,0] =  ident_4[0, *]
                    tx_matrix[*,1] =  ident_4[2, *] ;; -
                    tx_matrix[*,2] =  ident_4[1, *]
                end
                2: begin
                    axis_orient = [axis_orient[1], axis_orient[2], axis_orient[0]]
                    tx_matrix[*,0] =  ident_4[1, *]
                    tx_matrix[*,1] =  ident_4[2, *] ;; -
                    tx_matrix[*,2] =  ident_4[0, *]
                end
            endcase
        end

        ELSE: begin
            axis_orient = [0,1,2]
            case slice_axis of
                0: ;; noop
                1: begin
                    axis_orient = [axis_orient[0], axis_orient[2], axis_orient[1]]
                    tx_matrix[*,0] = ident_4[axis_orient[0], *]
                    tx_matrix[*,1] = ident_4[axis_orient[1], *]
                    tx_matrix[*,2] = ident_4[axis_orient[2], *]
                end
                2: begin
                    axis_orient = [axis_orient[1], axis_orient[2], axis_orient[0]]
;;                    axis_orient = [axis_orient[1], axis_orient[0], axis_orient[2]]
                    tx_matrix[*,0] = ident_4[axis_orient[0], *]
                    tx_matrix[*,1] = ident_4[axis_orient[1], *]
                    tx_matrix[*,2] = ident_4[axis_orient[2], *]
                end
            endcase
        end

    ENDCASE

    tx_matrix[3,3] = 1.0D
;    print, 'tx_matrix:'
;    print, tx_matrix

    RETURN, axis_orient
END

function adt_get_orientation_color, x, y, sl

    common scan_data

    ci = project.ci
    slice_axis = project.procpramarray[ci].slice_axis

    glyph_color = float([0,0,0])

    if project.procpramarray[ci].adt_proccess_flag eq 0 then begin
        return, glyph_color
    endif

    axis_orient = get_orient_rgb_axis(slice_axis=0)

    case slice_axis of
        0: begin
            fa    = reform((*project.dataarray[ci].frac_ani)[x,y,sl])
            evecs = transpose(reform((*project.dataarray[ci].eign_vec)[x,y,sl,*,*]))
            ev_axis = [ 0, 1, 2 ]
        end
        1: begin
            fa    = reform((*project.dataarray[ci].frac_ani)[x,sl,y])
            evecs = transpose(reform((*project.dataarray[ci].eign_vec)[x,sl,y,*,*]))
            ev_axis = [ 0, 2, 1 ]
        end
        2: begin
            fa    = reform((*project.dataarray[ci].frac_ani)[sl,x,y])
            evecs = transpose(reform((*project.dataarray[ci].eign_vec)[sl,x,y,*,*]))
            ev_axis = [ 1, 2, 0 ]
        end
    endcase

    glyph_color[axis_orient[0]] = abs([fa * evecs[0,ev_axis[0]]])
    glyph_color[axis_orient[1]] = abs([fa * evecs[0,ev_axis[1]]])
    glyph_color[axis_orient[2]] = abs([fa * evecs[0,ev_axis[2]]])

;    glyph_color = bytscl(abs([fa * evecs[0,ev_axis[0]], $
;                              fa * evecs[0,ev_axis[1]], $
;                              fa * evecs[0,ev_axis[2]]]))
;    glyph_color = glyph_color[axis_orient]
    
    return, bytscl(glyph_color)
end


; Subroutine name: orient
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
;
; This function will return the orientation maps to be used during display routines
;
; Editing Information:
    ;Edited by HS 2006/10/05
    ;Fix spelling mistakes and commenting.

function orient, fa, ev, fdim, pdim, ix, iy, slice_axis

    ir=fdim
    ip=pdim
    im=fltarr(ir, ip, 3)

	; orient DOES NOT needs to be rotated depending on the slicing selection
	axisid = get_orient_rgb_axis(slice_axis=0)

    case slice_axis of 
        0: ev_axis = [ 0, 1, 2 ]
        1: ev_axis = [ 0, 2, 1 ]
        2: ev_axis = [ 1, 2, 0 ]
    end
    
    ;im[*,*,0]=CONGRID(fa*abs(ev[*,*,axisid[0]]), ir, ip, CUBIC=-0.5)
    ;im[*,*,1]=CONGRID(fa*abs(ev[*,*,axisid[1]]), ir, ip, CUBIC=-0.5)
    ;im[*,*,2]=CONGRID(fa*abs(ev[*,*,axisid[2]]), ir, ip, CUBIC=-0.5)

    im[*,*,axisid[0]]=CONGRID(fa*abs(ev[*,*,ev_axis[0]]), ir, ip, CUBIC=-0.5)
    im[*,*,axisid[1]]=CONGRID(fa*abs(ev[*,*,ev_axis[1]]), ir, ip, CUBIC=-0.5)
    im[*,*,axisid[2]]=CONGRID(fa*abs(ev[*,*,ev_axis[2]]), ir, ip, CUBIC=-0.5)

    bim=bytscl(im)
    return, bim
end


; Subroutine name: color_mapper
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This function returns the 24-bit (3 dimensional) pixel values
; of a given grayscale slice according to the given color table.
; Editing Information:
    ;Edited by HS 2006/10/05
    ;Fix spelling mistakes and commenting.

FUNCTION color_mapper, im, fdim, pdim, table

  	color_im=BYTARR(fdim,pdim,3)
  	LOADCT, table
  	TVLCT, r, g, b, /get

	pix_idx = im[INDGEN(fdim*pdim, /LONG)]

	img_r = r[pix_idx]
	img_b = b[pix_idx]
	img_g = g[pix_idx]

	color_im = INTARR(fdim,pdim,3)
	color_im[*,*,0] = REFORM(TEMPORARY(img_r),fdim,pdim)
	color_im[*,*,1] = REFORM(TEMPORARY(img_g),fdim,pdim)
	color_im[*,*,2] = REFORM(TEMPORARY(img_b),fdim,pdim)

	RETURN, color_im
END


; draws the orientation circle
function circle, fdim, pdim, void, slice_axis=slice_axis

    dx=floor(fdim/2)
    dy=floor(pdim/2)
    im=fltarr(fdim, pdim, 3)
    r=ceil(1./3.*min([fdim, pdim]))
    
    ;; colorball needs to be rotated depending on the slicing selection according to Pierpaoli
    if (n_elements(slice_axis) ne 0) then begin
        axisid = get_orient_rgb_axis(slice_axis=slice_axis)
    endif else begin
        axisid = get_orient_rgb_axis()
    endelse
    
    for j=0, r-1 do begin
       for i=-r+1, r-1 do begin

          if i^2+j^2 le r^2 then begin
             x=float(i)/float(r)
             y=float(j)/float(r)
             im[dx+i,dy+j,axisid[2]]=sqrt(1.-x^2-y^2)
             im[dx+i,dy+j,axisid[1]]=y
             im[dx+i,dy+j,axisid[0]]=abs(x)
             im[dx+i,dy-j,axisid[2]]=sqrt(1.-x^2-y^2)
             im[dx+i,dy-j,axisid[1]]=y
             im[dx+i,dy-j,axisid[0]]=abs(x)
          endif

       endfor
    endfor
    
    return, bytscl(im)
end


; Subroutine name: adt_rotate_flip_display
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
;
; This routine will take in an ADT display sheet and flip and/or rotate each element individually.
;
; Editing Information:
    ;Edited by HS 2006/10/05
    ;Fix spelling mistakes and commenting.


pro adt_rotate_flip_display, sheet
    COMMON scan_data
    CI = project.ci

    HEAP_GC

    sz = size(*sheet)

    sz_x = sz[1]/3
    sz_y = sz[2]/4
    new_sz_x = sz_x
    new_sz_y = sz_y
    ;if you have to rotate the image by 90 or 270 then swap the x and y directions.
    if project.procPramArray[ci].rotate_Direction eq 1 or $
        project.procPramArray[ci].rotate_Direction eq 3 then begin
;        print, 'swaping x and y sizes'
        new_sz_y = sz[1]/3
        new_sz_x = sz[2]/4

    end

    ;create a new sheet to store the sheet after it is rotated.
    new_sheet = ptr_new(bytarr(new_sz_x * 3 , new_sz_y * 4, 3))

    for ii=0, 2 do begin
       for jj=0, 3 do begin

         temp =  ptr_new((*sheet)[ii*sz_x:(ii+1)*sz_x-1,jj*sz_y:(jj+1)*sz_y-1,*])
         mas_rotate_flip ,temp
         (*new_sheet)[ii*new_sz_x:(ii+1)*new_sz_x-1,jj*new_sz_y:(jj+1)*new_sz_y-1,*] = *temp

       end
    end
    sheet = new_sheet
end


; Subroutine name: adt_zoom_display
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
;
; adt_zoom_display takes in an ADT display sheet and zoom's it to the users specified
; value in MAS
;
; Editing Information:
    ;Edited by HS 2006/10/05
    ;Fix spelling mistakes and commenting.

pro adt_zoom_display, sheet
    HEAP_GC
    COMMON scan_data
    CI = project.ci

    ;zoom the image
    x_zoom    = project.procPramArray[ci].x_zoom
    y_zoom    = project.procPramArray[ci].y_zoom

    sz_sheet = size(*sheet)
    *sheet = rebin (*sheet, x_zoom*sz_sheet[1], y_zoom*sz_sheet[2],sz_sheet[3], /SAMPLE  )

end


; Subroutine name: adt_make_sheet
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; adt_make_sheet is used to return a sheet of images in color.
;
; Editing Information:
    ;Edited by HS 2006/10/05
    ;Fix spelling mistakes and commenting.

function adt_make_sheet , slice_number
    ;do garbage collection
    HEAP_GC
    COMMON scan_data
    CI = project.ci

    ;bring in variables
    scale_aver_d_max_by = project.procPramArray[CI].avg_d_max
    scale_aver_d_min_by = project.procPramArray[CI].avg_d_min
    scale_dia_max_by = project.procPramArray[CI].dia_max
    scale_dia_min_by = project.procPramArray[CI].dia_min
    scale_off_dia_max_by = project.procPramArray[CI].off_dia
    scale_fa_max_by = project.procPramArray[CI].fa_max
    slice_axis = project.procpramArray[ci].slice_axis

    sz_adt = size(*project.dataArray[CI].adt)


    ;if the whole data set has not been fitted then just grap the important slice.
    if (size(*project.dataArray[CI].adt))[0] eq 3 then begin
        ;if the image has not already been zoomed then the slice_number is incorrect
        ;so take the slice_number selection and divide that by Slice_zoom

         adt = (*project.dataArray[CI].adt)[*,*,*]
         fa = (*project.dataArray[CI].frac_Ani)[*,*]
         evecs=reform((*project.dataArray[CI].eign_Vec)[*,*,*,0])
         aver_d = (*project.dataArray[CI].Avg_Dif)[*,*]


        ;make a sheet to hold all the individual images
        sheet = bytarr(3*sz_adt[1],4*sz_adt[2],3)

        ;nice variables with low number of characters
        sx = sz_adt[1]
        sy = sz_adt[2]

        ;byte scale all images that are to be displayed
        byt_aver_d = Bytscl(aver_d, max=scale_aver_d_max_by, min=scale_aver_d_min_by)
        byt_S0 = Bytscl(reform(ADT[*,*,0]))
        byt_Dxx = Bytscl(ADT[*,*,1]>0.,max=scale_dia_max_by, min=scale_dia_min_by)
        byt_Dyy = Bytscl(ADT[*,*,2]>0.,max=scale_dia_max_by, min=scale_dia_min_by)
        byt_Dzz = Bytscl(ADT[*,*,3]>0.,max=scale_dia_max_by, min=scale_dia_min_by)
        byt_Dxy = Bytscl(ADT[*,*,4],max=scale_off_dia_max_by, min=-scale_off_dia_max_by)
        byt_Dxz = Bytscl(ADT[*,*,5],max=scale_off_dia_max_by, min=-scale_off_dia_max_by)
        byt_Dyz = Bytscl(ADT[*,*,6],max=scale_off_dia_max_by, min=-scale_off_dia_max_by)
        byt_fa = Bytscl(fa, max=scale_fa_max_by)

        ;Now collate the images into a single sheet
        for color=0,2 do begin

          sheet[0     :sx-1   ,3*sy  :4*sy-1, color] = byt_Dxx
          sheet[sx    :2*sx-1 ,3*sy  :4*sy-1, color] = byt_Dxy
          sheet[2*sx  :3*sx-1 ,3*sy  :4*sy-1, color] = byt_Dxz
          sheet[0     :sx-1   ,2*sy  :3*sy-1, color] = byt_S0
          sheet[sx    :2*sx-1 ,2*sy  :3*sy-1, color] = byt_Dyy
          sheet[2*sx  :3*sx-1 ,2*sy  :3*sy-1, color] = byt_Dyz
          sheet[0     :sx-1   ,sy    :2*sy-1, color] = byt_aver_d
          sheet[2*sx  :3*sx-1 ,sy    :2*sy-1, color] = byt_Dzz

        endfor

          sheet[0     :sx-1   ,0     :sy-1,*] = color_mapper(byt_fa , sx, sy, 0)
          temp_1 = orient(fa, evecs[*,*,*], sx, sy, sx, sy, slice_axis)
          sheet[sx    :2*sx-1 ,0     :sy-1,*] = Bytscl(temp_1, max=scale_fa_max_by*255.0)
          sheet[2*sx  :3*sx-1 ,0     :sy-1,*]=circle(sx, sy, slice_axis)

        return, sheet
    end

    if (size(*project.dataArray[CI].adt))[0] eq 4 then begin
        ;since the whole data set has been fitted we can extract different slices
        ;that is what we are doing here.

        if slice_axis eq 0 then begin
            ;if the image has not already been zoomed then the slice_number is incorrect
            ;so take the slice_number selection and divide that by Slice_zoom


            adt = reform((*project.dataArray[CI].adt)[*,*,slice_number,*])
            fa = reform((*project.dataArray[CI].frac_Ani)[*,*,slice_number])
            evecs=reform((*project.dataArray[CI].eign_Vec)[*,*,slice_number,*,0])
            aver_d = reform((*project.dataArray[CI].Avg_Dif)[*,*,slice_number])

            ;set the size of a single freq. phase image
            sx = sz_adt[1]
            sy = sz_adt[2]

        end else if slice_axis eq 1 then begin


            adt = reform((*project.dataArray[CI].adt)[*,slice_number,*,*])
            fa = reform((*project.dataArray[CI].frac_Ani)[*,slice_number,*])
            evecs=reform((*project.dataArray[CI].eign_Vec)[*,slice_number,*,*,0])
            aver_d = reform((*project.dataArray[CI].Avg_Dif)[*,slice_number,*])

            ;set the size of a single freq. slice image
            sx = sz_adt[1]
            sy = sz_adt[3]

        end else if slice_axis eq 2 then begin


            adt = reform((*project.dataArray[CI].adt)[slice_number,*,*,*])
            fa = reform((*project.dataArray[CI].frac_Ani)[slice_number,*,*])
            evecs=reform((*project.dataArray[CI].eign_Vec)[slice_number,*,*,*,0])
            aver_d = reform((*project.dataArray[CI].Avg_Dif)[slice_number,*,*])

            ;set the size of a single phase slice image
            sx = sz_adt[2]
            sy = sz_adt[3]

        end

        ;make a sheet to hold all the individual images
        sheet = bytarr(3*sx,4*sy,3)



        ;byte scale all images that are to be displayed
        byt_aver_d = Bytscl(aver_d, max=scale_aver_d_max_by, min=scale_aver_d_min_by)
        byt_S0 = Bytscl(reform(ADT[*,*,0]))
        byt_Dxx = Bytscl(ADT[*,*,1]>0.,max=scale_dia_max_by, min=scale_dia_min_by)
        byt_Dyy = Bytscl(ADT[*,*,2]>0.,max=scale_dia_max_by, min=scale_dia_min_by)
        byt_Dzz = Bytscl(ADT[*,*,3]>0.,max=scale_dia_max_by, min=scale_dia_min_by)
        byt_Dxy = Bytscl(ADT[*,*,4],max=scale_off_dia_max_by, min=-scale_off_dia_max_by)
        byt_Dxz = Bytscl(ADT[*,*,5],max=scale_off_dia_max_by, min=-scale_off_dia_max_by)
        byt_Dyz = Bytscl(ADT[*,*,6],max=scale_off_dia_max_by, min=-scale_off_dia_max_by)
        byt_fa = Bytscl(fa, max=scale_fa_max_by)

        ;Now correlate the images into a single sheet
        for color=0,2 do begin

          sheet[0     :sx-1   ,3*sy  :4*sy-1, color] = byt_Dxx
          sheet[sx    :2*sx-1 ,3*sy  :4*sy-1, color] = byt_Dxy
          sheet[2*sx  :3*sx-1 ,3*sy  :4*sy-1, color] = byt_Dxz
          sheet[0     :sx-1   ,2*sy  :3*sy-1, color] = byt_S0
          sheet[sx    :2*sx-1 ,2*sy  :3*sy-1, color] = byt_Dyy
          sheet[2*sx  :3*sx-1 ,2*sy  :3*sy-1, color] = byt_Dyz
          sheet[0     :sx-1   ,sy    :2*sy-1, color] = byt_aver_d
          sheet[2*sx  :3*sx-1 ,sy    :2*sy-1, color] = byt_Dzz

        endfor

          sheet[0     :sx-1   ,0     :sy-1,*] = color_mapper(byt_fa , sx, sy, 0)
          temp_1 = orient(fa, evecs[*,*,*], sx, sy, sx, sy, slice_axis)
          sheet[sx    :2*sx-1 ,0     :sy-1,*] = Bytscl(temp_1, max=scale_fa_max_by*255.0)
          sheet[2*sx  :3*sx-1 ,0     :sy-1,*]=circle(sx, sy, slice_axis)

        return, sheet

    end
    void = dialog_message('Error in creating a sheet!', /center, /error)
    return, -1
end


; Subroutine name: adt_make_image
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
;
;Instead of making a sheet of images make a single image to display
;
; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.


function adt_make_image , slice_number, float_flag
    HEAP_GC
    COMMON scan_data
    CI = project.ci

    ;if they didnt define float flag then define it for them
    if N_ELEMENTS(float_flag) eq 0 then float_flag = 0

    ;bring in the common variables
    adt_display_type =     project.procPramArray[CI].adt_display_type
    adt_proccess_flag = project.procPramArray[CI].adt_proccess_flag
    ADT_Start =    project.procPramArray[CI].ADT_Start
    adt_display_multi =    project.procPramArray[CI].adt_display_multi
    scale_aver_d_max_by = project.procPramArray[CI].avg_d_max
    scale_aver_d_min_by = project.procPramArray[CI].avg_d_min
    scale_dia_max_by = project.procPramArray[CI].dia_max
    scale_dia_min_by = project.procPramArray[CI].dia_min
    scale_off_dia_max_by = project.procPramArray[CI].off_dia
    scale_fa_max_by = project.procPramArray[CI].fa_max
    slice_axis = project.procpramArray[ci].slice_axis

    ;; the data must already be fitted, and should be up to someone
    ;; else up the call stack to have fitted it.
    if adt_proccess_flag eq 0 then return, -1

    if slice_axis ne 0 and adt_proccess_flag eq 0 then $
        void = dialog_message(['Since whole data set has not been fitted...',$
                               'Assuming freq. phase direction to slice data'])

    sz_adt = size(*project.dataArray[project.CI].adt)

    ;what type of image does the user want?
    if adt_display_type eq 0 then begin
       ;user want tensor image

       ;grab the specified slice of interest and type of image
       ;if the user wants a color trace movie they they are going to get it.
       ;has the whole data set been processed
        if adt_proccess_flag eq 0 then begin
           if ADT_Start eq 7 then $
             image = reform((*project.dataArray[project.CI].adt)[*,*,1:3]) $
           else $
             image = reform((*project.dataArray[project.CI].adt)[*,*,ADT_Start])
        end else begin
            ;the whole data set has been processed
            ;since this has been done we can extract slices in different planes
            CASE ADT_Start OF
                7: begin
                    if slice_axis eq 0 then $
                        image = reform((*project.dataArray[project.CI].adt)[*,*,slice_number,1:3]) $
                    else if slice_axis eq 1 then $
                        image = reform((*project.dataArray[project.CI].adt)[*,slice_number,*,1:3]) $
                    else if slice_axis eq 2 then $
                        image = reform((*project.dataArray[project.CI].adt)[slice_number,*,*,1:3])
                END
                8: BEGIN
                 if slice_axis eq 0 then begin
                        fa = reform((*project.dataArray[project.CI].frac_Ani)[*,*,slice_number])
                        evecs=reform((*project.dataArray[CI].eign_Vec)[*,*,slice_number,*,0])
                    end else if slice_axis eq 1 then begin
                        fa = reform((*project.dataArray[project.CI].frac_Ani)[*,slice_number,*])
                        evecs=reform((*project.dataArray[CI].eign_Vec)[*,slice_number,*,*,0])
                    end else if slice_axis eq 2 then begin
                        fa = reform((*project.dataArray[project.CI].frac_Ani)[slice_number,*,*])
                        evecs=reform((*project.dataArray[CI].eign_Vec)[slice_number,*,*,*,0])
                    end
                end
                9: BEGIN
                    image = circle(200,200,0)
                end

                else: begin
                    if slice_axis eq 0 then $
                        image = reform((*project.dataArray[project.CI].adt)[*,*,slice_number,ADT_Start]) $
                    else if slice_axis eq 1 then $
                        image = reform((*project.dataArray[project.CI].adt)[*,slice_number,*,ADT_Start]) $
                    else if slice_axis eq 2 then $
                        image = reform((*project.dataArray[project.CI].adt)[slice_number,*,*,ADT_Start])
          END
            ENDCASE
        end

        sz_image = size(image)

        ;if the users just want a floating point image then let them quit now.
        if float_flag eq 1 then return, image

       ;scale the image so that it looks correct
       ; if ADT_Start eq 0 then image = color_mapper( Bytscl(image), sz_image[1], sz_image[2], 8)
       if ADT_Start eq 0 then image = color_mapper( Bytscl(image), sz_image[1], sz_image[2], 0)
       if ADT_Start eq 1 then image = Bytscl(image, min=scale_dia_min_by ,max=scale_dia_max_by)
       if ADT_Start eq 2 then image = Bytscl(image, min=scale_dia_min_by ,max=scale_dia_max_by)
       if ADT_Start eq 3 then image = Bytscl(image, min=scale_dia_min_by ,max=scale_dia_max_by)
       if ADT_Start eq 4 then image = Bytscl(image, min=-scale_off_dia_max_by ,max=scale_off_dia_max_by)
       if ADT_Start eq 5 then image = Bytscl(image, min=-scale_off_dia_max_by ,max=scale_off_dia_max_by)
       if ADT_Start eq 6 then image = Bytscl(image, min=-scale_off_dia_max_by ,max=scale_off_dia_max_by)
       if ADT_Start eq 7 then begin
           byt_Dxx = Bytscl(image[*,*,0]>0.,max=scale_dia_max_by, min=scale_dia_min_by)
           byt_Dyy = Bytscl(image[*,*,1]>0.,max=scale_dia_max_by, min=scale_dia_min_by)
           byt_Dzz = Bytscl(image[*,*,2]>0.,max=scale_dia_max_by, min=scale_dia_min_by)
           image[*,*,0] = bytscl(byt_Dzz)
           image[*,*,1] = bytscl(byt_Dyy)
           image[*,*,2] = bytscl(byt_Dxx)

       end
       if ADT_Start eq 8 then begin

			CASE slice_axis OF
				0: BEGIN
					sx = sz_adt[1]
					sy = sz_adt[2]
				END
				1: BEGIN
					sx = sz_adt[1]
					sy = sz_adt[3]
				END
				2: BEGIN
					sx = sz_adt[2]
					sy = sz_adt[3]
				END
			ENDCASE

        	image = orient(fa, evecs , sx, sy, sx, sy, slice_axis)
            image = Bytscl(image, max=scale_fa_max_by*255.0)
       end

       return, image

    end

    if adt_display_type eq 2 then begin
        ;user wants average diffusivity image
        ;dont let the user change slice directions if they have not processed the whole image.
        if sz_adt[0] eq 3 then begin
            aver_d = (*project.dataArray[CI].Avg_Dif)[*,*]
            if float_flag eq 1 then return, aver_d
            byt_aver_d = Bytscl(aver_d, max=scale_aver_d_max_by, min=scale_aver_d_min_by)
			return, byt_aver_d
        end
        if sz_adt[0] eq 4 then begin

            ;what plane do you want to slice the image by
            if slice_axis eq 0 then begin

                aver_d = reform((*project.dataArray[CI].Avg_Dif)[*,*,slice_number])
                if float_flag eq 1 then return, aver_d
                byt_aver_d = Bytscl(aver_d, max=scale_aver_d_max_by, min=scale_aver_d_min_by)
				return, byt_aver_d

            end else if slice_axis eq 1 then begin
                aver_d = reform((*project.dataArray[CI].Avg_Dif)[*,slice_number,*])
          if float_flag eq 1 then return, aver_d
                byt_aver_d = Bytscl(aver_d, max=scale_aver_d_max_by, min=scale_aver_d_min_by)
				return, byt_aver_d

            end else if slice_axis eq 2 then begin
                aver_d = reform((*project.dataArray[CI].Avg_Dif)[slice_number,*,*])
          if float_flag eq 1 then return, aver_d
                byt_aver_d = Bytscl(aver_d, max=scale_aver_d_max_by, min=scale_aver_d_min_by)
				return, byt_aver_d
            end
        end
    end

    if adt_display_type eq 3 then begin
        ;user wants fractional anisotropy image
        if sz_adt[0] eq 3 then begin
            fa = (*project.dataArray[CI].frac_Ani)[*,*]
            if float_flag eq 1 then return, fa
            byt_fa = Bytscl(fa, max=scale_fa_max_by)
			return, byt_fa
        end

        if sz_adt[0] eq 4 then begin
            ;what plane do you want to slice the image by
            if slice_axis eq 0 then begin
                fa = reform((*project.dataArray[CI].frac_Ani)[*,*,slice_number])
                if float_flag eq 1 then return, fa
                byt_fa = Bytscl(fa, max=scale_fa_max_by)
				return, byt_fa

            end else if slice_axis eq 1 then begin
                fa = reform((*project.dataArray[CI].frac_Ani)[*,slice_number,*])
          if float_flag eq 1 then return, fa
                byt_fa = Bytscl(fa, max=scale_fa_max_by)
				return, byt_fa

            end else if slice_axis eq 2 then begin
                fa = reform((*project.dataArray[CI].frac_Ani)[slice_number,*,*])
          if float_flag eq 1 then return, fa
                byt_fa = Bytscl(fa, max=scale_fa_max_by)
				return, byt_fa
            end
        end
    end

    if adt_display_type eq 4 then begin

        adim_start = project.procpramarray[project.ci].adim_start
        case slice_axis of
            0: image = reform((*project.dataArray[CI].state1)[*,*,slice_number,adim_start])
            1: image = reform((*project.dataArray[CI].state1)[*,slice_number,*,adim_start])
            2: image = reform((*project.dataArray[CI].state1)[slice_number,*,*,adim_start])
            else:
        endcase

        ;; it is expected that this is flipped/rotated as the it used
        ;; to duplicate state2 data. Unforunately.
        ptmp = ptr_new(image)
        mas_rotate_flip, ptmp
        image = *ptmp
        ptr_free, ptmp

        if float_flag eq 1 then return, image

        mas_windowing, image

        return, image
    endif

;;    Replaced by above block, BT 20080606
;     if adt_display_type eq -99 then begin ;; was 4
;        mas_load_state_2

;        if (size((*project.dataArray[CI].state2)))[0] eq 2 then $
;          image = reform((*project.dataArray[project.CI].state2)[*,*]) $
;        else $
;          image = reform((*project.dataArray[project.CI].state2)[*,*,slice_number])

;        ;if there is a need for a floating point image then return that else scale the image.
;        if float_flag eq 1 then return, image

;        mas_windowing , image

;        return, image
;     end

    IF adt_display_type EQ 5 THEN BEGIN
        ;;user want 2D, lattice anisotropy image

        result = project.procPramArray[CI].latt_cur_slice
        is3D = project.procPramArray[CI].use_3D_lattice

        if (result ne slice_number) then begin
            adt_calc_lattice_anisotrophy, slice_number, 0
            project.procPramArray[CI].use_3D_lattice = 0
            project.procPramArray[CI].latt_cur_slice = slice_number

        endif else begin

            if (project.procPramArray[CI].use_3D_lattice eq 1) then begin

                adt_calc_lattice_anisotrophy, slice_number, 0
                project.procPramArray[CI].use_3D_lattice = 0
                project.procPramArray[CI].latt_cur_slice = slice_number

            endif

        endelse

        IF sz_adt[0] EQ 3 THEN BEGIN
            la = (*project.dataArray[CI].latt_Ani)[*,*]
        ENDIF ELSE IF sz_adt[0] EQ 4 THEN BEGIN
            ;; la calculation is freq_phase, vs phase_slice vs freq_slice dependent
            ;; only stores for current slice axis, the la matrix will always have the
            ;; current slice axis as the third dimension
            la = REFORM((*project.dataArray[CI].latt_Ani)[*,*,slice_number])
        ENDIF

        IF float_flag EQ 1 THEN RETURN, la
        byt_la = BYTSCL(la)
        RETURN, byt_la
    END

    IF adt_display_type EQ 6 THEN BEGIN
        ;;user want 3D lattice anisotropy image

        result = project.procPramArray[CI].latt_cur_slice
        is3D = project.procPramArray[CI].use_3D_lattice

        if (result ne slice_number) then begin

            adt_calc_lattice_anisotrophy, slice_number, 1

            project.procPramArray[CI].use_3D_lattice = 1
            project.procPramArray[CI].latt_cur_slice = slice_number

        endif else begin

            if (project.procPramArray[CI].use_3D_lattice eq 0) then begin

                adt_calc_lattice_anisotrophy, slice_number, 1

                project.procPramArray[CI].use_3D_lattice = 1
                project.procPramArray[CI].latt_cur_slice = slice_number

            endif

        endelse

        IF sz_adt[0] EQ 3 THEN BEGIN
            la = (*project.dataArray[CI].latt_Ani)[*,*]
        ENDIF ELSE IF sz_adt[0] EQ 4 THEN BEGIN
            ;; la calculation is freq_phase, vs phase_slice vs freq_slice dependent
            ;; only stores for current slice axis, the la matrix will always have the
            ;; current slice axis as the third dimension
            la = REFORM((*project.dataArray[CI].latt_Ani)[*,*,slice_number])
        ENDIF

        IF float_flag EQ 1 THEN RETURN, la
        byt_la = BYTSCL(la)
        RETURN, byt_la
    END

    ;;if they are trying to use the procedure for what it isnt meant for return -1
    void = dialog_message(['Error in adt_make_image','Try using the other display button'], /error, /center)
    return, -1
end


; Subroutine name: adt_movie
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

PRO adt_movie
    HEAP_GC
    COMMON scan_data
    CI = project.ci

    do_adt_regress_mono

    if project.procpramArray[CI].adt_proccess_flag eq 0 then return

    adt_display_type =     project.procPramArray[CI].adt_display_type
    adt_proccess_flag = project.procPramArray[CI].adt_proccess_flag
    ADT_Start =    project.procPramArray[CI].ADT_Start
    adt_display_multi =    project.procPramArray[CI].adt_display_multi

    x_zoom    = project.procPramArray[ci].x_zoom
    y_zoom    = project.procPramArray[ci].y_zoom


    slice_axis = project.procpramArray[ci].slice_axis

    scale_aver_d_max_by = project.procPramArray[CI].avg_d_max
    scale_aver_d_min_by = project.procPramArray[CI].avg_d_min
    scale_dia_max_by = project.procPramArray[CI].dia_max
    scale_dia_min_by = project.procPramArray[CI].dia_min
    scale_off_dia_max_by = project.procPramArray[CI].off_dia
    scale_fa_max_by = project.procPramArray[CI].fa_max

    sz_adt     =   size(*project.dataArray[project.CI].adt)
    ;print, 'sz_adt',sz_adt
    update_status_bar , 'Loading Movie'

    if adt_display_type eq 4 then begin
       update_status_bar, 'Please change display selection.'
       return
    end
    if sz_adt[0] lt 4 then begin
       ;print,(size(*project.dataArray[project.CI].adt))[0]
       update_status_bar,'Whole data set not fited. Please change Tensor Processing'
       return
    end


    if adt_display_multi eq 0 then begin
       ;the user wants a single image movie

       ;now measure the the size of one of the images and
       ;assume that all the images will be the same size.
       p_image = ptr_new(adt_make_image(1))
       mas_rotate_flip ,p_image
       mas_zoom, p_image
       sx = (size(*p_image))[1]
       sy = (size(*p_image))[2]

        if slice_axis eq 0 then begin
            sz = sz_adt[3]

        end else if slice_axis eq 1 then begin
            sz = sz_adt[2]

        end else if slice_axis eq 2 then begin
            sz = sz_adt[1]
        end

	progressbar = Obj_New('progressbar', Color='red', Text='Loading Frames',/NOCANCEL)
	progressbar -> Start
	progressbartemp = float(1.0/sz)
	progressbar -> Update, progressbartemp*100


        base = WIDGET_BASE(TITLE = project.imndarray[ci].display_name);'Animation Widget')

        ;XINTERANIMATE, SET=[sx, sy, sz], /SHOWLOAD, /BLOCK ,MPEG_QUALITY=80
        animate = CW_ANIMATE(base, sx, sy, sz, /NO_KILL )
;        WIDGET_CONTROL, /REALIZE, base


        ; Load the images into XINTERANIMATE:
        for ii=0, sz-1 do begin
            p_image = ptr_new(adt_make_image(ii))

            mas_rotate_flip ,p_image

            ;zoom the image
            mas_zoom, p_image

            CW_ANIMATE_LOAD, animate, FRAME=ii, IMAGE = *p_image

 			progressbar_load = progressbartemp*ii
         	progressbar -> Update, progressbar_load*100
        end

		progressbar -> Destroy

        ;Play the animation:
		WIDGET_CONTROL, /REALIZE, base
        CW_ANIMATE_GETP, animate, pixmap_vect
        CW_ANIMATE_RUN, animate, 25
        XMANAGER, 'CW_ANIMATE Demo', base, EVENT_HANDLER = 'EHANDLER', /NO_BLOCK



    end
    if adt_display_multi eq 1 then begin
       ;the user wants a sheet movie

      ;now measure the the size of one of the images and
       ;assume that all the images will be the same size.
       p_sheet = ptr_new(adt_make_sheet(0))
       adt_rotate_flip_display, p_sheet
       adt_zoom_display, p_sheet

       sx = (size(*p_sheet))[1]
       sy = (size(*p_sheet))[2]

        ;setting the right dimensions depending up fit_zoom_flag and what
        ;direction to traverse the data.
        if slice_axis eq 0 then begin

            sz = sz_adt[3]

        end else if slice_axis eq 1 then begin

            sz = sz_adt[2]

        end else if slice_axis eq 2 then begin

            sz = sz_adt[1]
        end

; Creating a progressbar
progressbar = Obj_New('progressbar', Color='red', Text='Loading Frames',/NOCANCEL)
progressbar -> Start
progressbartemp = float(1.0/sz)
progressbar -> Update, progressbartemp*100


        base = WIDGET_BASE(TITLE = project.imndarray[ci].display_name); 'Animation Widget')

       ;XINTERANIMATE, SET=[3*sx, 4*sy, sz], /SHOWLOAD, /BLOCK ,MPEG_QUALITY=80
        animate = CW_ANIMATE(base, sx, sy, sz, /NO_KILL )
;        WIDGET_CONTROL, /REALIZE, base


       ; Load the images into XINTERANIMATE:


        for ii=0, sz-1 do begin
            p_sheet = ptr_new(adt_make_sheet(ii))

            adt_rotate_flip_display, p_sheet

            ;zoom the image
            adt_zoom_display, p_sheet


            ;XINTERANIMATE, FRAME = ii , IMAGE = sheet
            CW_ANIMATE_LOAD, animate, FRAME=ii, IMAGE = *p_sheet

            loadct, 0,  /SILENT

            progressbartemp += progressbartemp
         	progressbar -> Update, progressbartemp*100

       end

progressbar -> Destroy

       ; Play the animation:
       ;XINTERANIMATE, /KEEP_PIXMAPS

       WIDGET_CONTROL, /REALIZE, base

        CW_ANIMATE_GETP, animate, pixmap_vect
        CW_ANIMATE_RUN, animate, 25
        XMANAGER, 'CW_ANIMATE Demo', base, EVENT_HANDLER = 'EHANDLER', /NO_BLOCK

    end

END


; Subroutine name: ADT_MPEG_WIDGET_EVENT
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

pro ADT_MPEG_WIDGET_EVENT, event
    COMPILE_OPT IDL2
    HEAP_GC
    ;bring in the appropriate global variables
    COMMON scan_data
    CI = project.ci


    wWidget =  Event.top

    case Event.id of


       Widget_Info(wWidget, FIND_BY_UNAME='ADT_MPEG_FILE_BUTTON'):begin
         ;the user clicked the button to select a file
         ;have the users select a file
         file_name = dialog_pickfile(PATH=project.current_Path, /WRITE, FILTER='*.mpg')

         ;check to see if the users selected nothing
         if file_name eq '' then return

         ;check to see if they put a .mpg in the file name.if they didnt then add it
         if strpos(file_name, '.mpg') lt 0  then begin
          file_name +='.mpg'
         end

         ;store the path for later.
         project.current_Path = file_name

         ;display the file name in the text widget
         widget_control, Widget_Info(wWidget, FIND_BY_UNAME='ADT_MPEG_FILE_NAME'), $
          SET_VALUE = file_name
       end

       Widget_Info(wWidget, FIND_BY_UNAME='ADT_MPEG_CREATE'):begin
         ;create the mpeg file

         ;first gather all the info from the mpeg widgets
         ;i should do some error checking to make sure that the values are actually number and not strings
         ;for rep_frame and frame_gap
         WIDGET_control, Widget_Info(wWidget, FIND_BY_UNAME='ADT_MPEG_FILE_NAME'), get_value = file_name
         if file_name eq 'Pick a file' then begin
          void = dialog_message( ['Please select a file name'])
                 return
         end

         WIDGET_control, Widget_Info(wWidget, FIND_BY_UNAME='ADT_MPEG_REP_FRAME'), get_value = sTemp
          rep_frame = (fix(sTemp))[0]

         WIDGET_control, Widget_Info(wWidget, FIND_BY_UNAME='ADT_MPEG_FRAME_GAP'), get_value = sTemp
          frame_gap = (fix(sTemp))[0]

         WIDGET_control, Widget_Info(wWidget, FIND_BY_UNAME='ADT_MPEG_QUALITY'), get_value = sTemp
          quality = (fix(sTemp))[0]

         if quality lt 0 or quality gt 100 then begin
          void = dialog_message( ['Please change quality from 0 to 100'])
                 return
         end

         ;check to see if frame_gap is the withen the correct range.
         if frame_gap le 0 then begin
          void = dialog_message( ['The Frame Gap can not be less than 1' $
                              ,'Please change Frame Gap'])
              return

         end

         ;check to see if frame_gap is the withen the correct range.
         if rep_frame le 0 then begin
          void = dialog_message( ['The Frame Repetions can not be less than 1' $
                              ,'Please change Frame Gap'])
              return

         end


           do_adt_regress_mono

           adt_display_type =     project.procPramArray[CI].adt_display_type
           adt_proccess_flag = project.procPramArray[CI].adt_proccess_flag
           ADT_Start =    project.procPramArray[CI].ADT_Start
           adt_display_multi =    project.procPramArray[CI].adt_display_multi

           x_zoom    = project.procPramArray[ci].x_zoom
           y_zoom    = project.procPramArray[ci].y_zoom


           slice_axis = project.procpramArray[ci].slice_axis

           scale_aver_d_max_by = project.procPramArray[CI].avg_d_max
           scale_aver_d_min_by = project.procPramArray[CI].avg_d_min
           scale_dia_max_by = project.procPramArray[CI].dia_max
           scale_dia_min_by = project.procPramArray[CI].dia_min
           scale_off_dia_max_by = project.procPramArray[CI].off_dia
           scale_fa_max_by = project.procPramArray[CI].fa_max

           sz_adt     =   size(*project.dataArray[project.CI].adt)
         ;print, sz_adt


           if adt_display_type eq 4 then begin
              void = dialog_message( [   'Diffusion-Weighted mpeg creation is handled under the export menu' $
                              ,'Please change display selection'])
              return
           end
           if sz_adt[0] lt 4 then begin
              ;print,(size(*project.dataArray[project.CI].adt))[0]
              void = dialog_message( ['Whole data set not fited.',' Please change Tensor Processing'])
              return
           end


           if slice_axis eq 0 then begin
               sx = sz_adt[1]*x_zoom
               sy = sz_adt[2]*y_zoom
               sz = sz_adt[3]

           end else if slice_axis eq 1 then begin
               sx = sz_adt[1]*x_zoom
               sy = sz_adt[3]*y_zoom
               sz = sz_adt[2]

           end else if slice_axis eq 2 then begin
               sx = sz_adt[2]*x_zoom
               sy = sz_adt[3]*y_zoom
               sz = sz_adt[1]
           end



         if adt_display_multi eq 0 then begin
              ;the user wants a single image movie
                 progressbar = Obj_New('progressbar', Color='red', Text='Loading Frames',/NOCANCEL)
             progressbar -> Start

          rgb_image = bytarr(3,sx,sy)

          ; Open an MPEG sequence:
             mpegID = MPEG_OPEN( [sx,sy] , FILENAME=file_name, quality=quality , IFRAME_GAP=frame_gap)

          ; Add the frames:
             counter=0
             for ii=0,(sz-1) do begin

                  image = ptr_new(adt_make_image(ii))

                  mas_rotate_flip ,image

                   mas_zoom, image

              rgb_image[0,*,*] = (*image)[*,*,0]
              rgb_image[1,*,*] = (*image)[*,*,1]
              rgb_image[2,*,*] = (*image)[*,*,2]

              progressBar -> Update, (float(ii)/float(sz-1))*100.0


              for jj=0, rep_frame-1 do begin

                 MPEG_PUT, mpegID, IMAGE=rgb_image , FRAME=counter, /ORDER;;, /COLOR
                 counter++

                  end

             end

          progressbar -> Destroy
          progressbar = Obj_New('progressbar', Color='red', Text=string('Saving:', file_name),/NOCANCEL)
          progressbar -> Start

                                ; Save the MPEG sequence in the file myMovie.mpg:
          MPEG_SAVE, mpegID, FILENAME=mpegFileLocation

                                ; Close the MPEG sequence:
          MPEG_CLOSE, mpegID

          progressbar -> Destroy

      end
         if adt_display_multi eq 1 then begin
              ;the user wants a sheet movie

          progressbar = Obj_New('progressbar', Color='red', Text='Loading Frames',/NOCANCEL)
             progressbar -> Start

          rgb_image = bytarr(3, 3*sx, 4*sy)

          ; Open an MPEG sequence:
             mpegID = MPEG_OPEN( [3*sx, 4*sy] , FILENAME=file_name, quality=quality , IFRAME_GAP=frame_gap)

          ; Add the frames:
             counter=0
             for ii=0,(sz-1) do begin

                  sheet = ptr_new(adt_make_sheet(ii))

                   adt_rotate_flip_display, sheet

                   ;zoom the image
                   adt_zoom_display, sheet

              rgb_image[0,*,*] = (*sheet)[*,*,0]
              rgb_image[1,*,*] = (*sheet)[*,*,1]
              rgb_image[2,*,*] = (*sheet)[*,*,2]

              progressBar -> Update, (float(ii)/float(sz-1))*100.0

              for jj=0, rep_frame-1 do begin

                 MPEG_PUT, mpegID, IMAGE=rgb_image , FRAME=counter, /ORDER
                 counter++

                  end

             end

          progressbar -> Destroy
          progressbar = Obj_New('progressbar', Color='red', Text='Writing MPEG File...',/NOCANCEL)
             progressbar -> Start


             ; Save the MPEG sequence in the file myMovie.mpg:
             MPEG_SAVE, mpegID, FILENAME=mpegFileLocation

             ; Close the MPEG sequence:
             MPEG_CLOSE, mpegID

          progressbar -> Destroy


           end

       end


    else:
       ;If there is an event that i dont capture then do nothing.
       ;the 3 text boxes i dont care about.
    endcase

end


; Subroutine name: ADT_MPEG_WIDGET
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This will create the widget interface to control the movie creation process.

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

PRO ADT_MPEG_WIDGET

    COMPILE_OPT IDL2
    HEAP_GC
    ;bring in the appropriate global variables
    COMMON scan_data
    COMMON common_widgets
    CI = project.ci



    ;If the ADT ROI window is already on the screen then dont make another one.
    IF N_ELEMENTS(ADT_MPEG_WINDOW_BASE) EQ 1 THEN $
        IF     widget_info( ADT_MPEG_WINDOW_BASE,  /VALID_ID ) eq 1 THEN return

    ;if the current scan does not have a good b-matrix then make every thing insensitive
    SENSITIVE =  ptr_valid (project.imndArray[CI].b_matrix)

    title = string ('ADT MPEG Creation')
    ADT_MPEG_WINDOW_BASE = widget_base(TITLE=title, $
            UVALUE = 'ADT_MPEG_WINDOW_BASE', $
            XOFFSET=420 ,YOFFSET=637  ,ROW=5)

    void = widget_button(ADT_MPEG_WINDOW_BASE, value='Choose file name',UNAME = 'ADT_MPEG_FILE_BUTTON',XSIZE = 100)
    void = widget_text(ADT_MPEG_WINDOW_BASE, value='Pick a file' ,XSIZE = 40, /EDITABLE ,UNAME ='ADT_MPEG_FILE_NAME')



    void = widget_label(ADT_MPEG_WINDOW_BASE, value='Repeat Frames (10)',XSIZE = 100, /ALIGN_RIGHT)
    void = widget_text(ADT_MPEG_WINDOW_BASE, value='10' ,XSIZE = 40, /EDITABLE ,UNAME ='ADT_MPEG_REP_FRAME')

    void = widget_label(ADT_MPEG_WINDOW_BASE, value='Frame Gap (2)',XSIZE = 100, /ALIGN_RIGHT)
    void = widget_text(ADT_MPEG_WINDOW_BASE, value='2' ,XSIZE = 40, /EDITABLE ,UNAME ='ADT_MPEG_FRAME_GAP')

    void = widget_label(ADT_MPEG_WINDOW_BASE, value='Quality (80)',XSIZE = 100, /ALIGN_RIGHT)
    void = widget_text(ADT_MPEG_WINDOW_BASE, value='80' ,XSIZE = 40, /EDITABLE ,UNAME ='ADT_MPEG_QUALITY')

    void = widget_button(ADT_MPEG_WINDOW_BASE, value='Create MPEG',UNAME = 'ADT_MPEG_CREATE',XSIZE = 100)

    widget_control, ADT_MPEG_WINDOW_BASE, /realize

    xmanager, 'ADT_MPEG_WIDGET',ADT_MPEG_WINDOW_BASE,/no_block, GROUP_LEADER=ADT_WINDOW_BASE


END

; Subroutine name: display_adt
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.
    ;Edited by CD 2007/02/20
    ;if 2D image, display using the enhance display window from mas_display

pro _draw_scroll_h, ev
    widget_control, ev.top, get_uvalue=draw
    widget_control, (*draw), xsize=ev.x, ysize=ev.y
    widget_control, ev.top, xsize=ev.x, ysize=ev.y
end

pro _draw_scroll_c, tlb
    widget_control, tlb, get_uvalue=draw, /no_copy
    ptr_free, draw
end

; Subroutine name: display_adt
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; Widget to display ADT results
;
; Editing Information:
; Edited by Magdoom 2015/10/20
; Fixed bug with multi display

pro display_adt

    COMMON scan_data

    print, 'Inside display_adt'
    HEAP_GC
    do_adt_regress_mono

    CI = project.CI

    ;; likely that adt processing had an error or was cancelled.
    if (project.procpramarray[ci].adt_proccess_flag eq 0) then return

    sdim_start = project.procPramArray[CI].sdim_start
    pdim_start = project.procPramArray[CI].pdim_start
    fdim_start = project.procPramArray[CI].fdim_start

    ffov = project.imndArray[ci].f_fov
    pfov = project.imndArray[ci].p_fov
    sfov = project.imndArray[ci].s_fov

    if (project.procpramarray[ci].down_sample eq 1) then begin
        Freq_interp  = 1.0/project.procPramArray[ci].Freq_interp
        Phase_interp = 1.0/project.procPramArray[ci].Phase_interp
        Slice_interp = 1.0/project.procPramArray[ci].Slice_interp
    endif else begin
        Freq_interp  = project.procPramArray[ci].Freq_interp
        Phase_interp = project.procPramArray[ci].Phase_interp
        Slice_interp = project.procPramArray[ci].Slice_interp
    endelse

    slice_axis = project.procpramArray[ci].slice_axis

    title = 'Diffusion Tensor & Parameters Slice#'+strtrim(sdim_start,2)

    ;check to see what type of display the users wants either a sheet of multiple images
    ;or a display of a single image
    IF project.procPramArray[CI].adt_display_multi EQ 0 THEN BEGIN

        ;;if the user wants to display all the slices for a scan
        ;;this button connects the multi button from the main display
        ;;so the user can get all the slices for the image they want.

        ;; note that ADT_Start == 9 means Color Ball only, so disply/save
        ;; a single image
       IF project.procPramArray[CI].single_Multi_flag EQ 0 or project.procPramArray[CI].ADT_start eq 9 THEN BEGIN

            ;we have to make sure that we are selecting the right slice and that
            ;fov is initilized correctly
            IF slice_axis EQ 0 THEN BEGIN
                slice_number = sdim_start
                x_fov = ffov*Freq_interp
                y_fov = pfov*Phase_interp

            ENDIF ELSE IF slice_axis EQ 1 THEN BEGIN
                slice_number = pdim_start
                x_fov = ffov*Freq_interp
                y_fov = sfov*Slice_interp

            ENDIF ELSE IF slice_axis EQ 2 THEN BEGIN
                slice_number = fdim_start
                x_fov = pfov*Phase_interp
                y_fov = sfov*Slice_interp

            ENDIF

           p_display = ptr_new(adt_make_image(slice_number))
           p_raw     = ptr_new(adt_make_image(slice_number,1))
           disp_float_labels = 1

           ;has the image already been zoomed? if they are trying to display a diffusion
           ; weighted image which comes from state2 then there is no need to zoom the image b/c it has
           ; already been done.

           IF project.procPramArray[CI].adt_display_type ne 4 then begin
               mas_zoom, p_display
               mas_zoom, p_raw

			; *************************************
			; HS - 20061111
			; Now we have to rotate every image except the diffusion weighted imaged
			; *************************************
               mas_rotate_flip, p_display
               mas_rotate_flip, p_raw

           ENDIF

           ;; CD - This will eventually call Display_DrawRes_Create below

       END ELSE BEGIN
            ;the user want all slices from a specified data set.
           ;print, 'The user wants a slice movie'
           print, 'The user wants a multi-slice image'
           sz = size(*project.dataArray[project.CI].adt)
           ;we have to make sure that we are selecting the right slice and that
           ;fov is initilized correctly
           IF slice_axis EQ 0 THEN BEGIN
               slice_number = sdim_start
               x_fov = ffov*Freq_interp
               y_fov = pfov*Phase_interp
               ;;how many slices in the x and y direction to display?
               totalFrames = sz[3]
                
           END ELSE IF slice_axis EQ 1 THEN BEGIN
               slice_number = pdim_start
               x_fov = ffov*Freq_interp
               y_fov = sfov*Slice_interp
               ;;how many slices in the x and y direction to display?
               totalFrames = sz[2]
           END ELSE IF slice_axis EQ 2 THEN BEGIN
               slice_number = fdim_start
               x_fov = pfov*Phase_interp
               y_fov = sfov*Slice_interp
               ;;how many slices in the x and y direction to display?
               totalFrames = sz[1]
           END
           
           ;; this gets overwritten later on.
           xSize = ceil(sqrt(totalFrames))
           ySize = ceil(float(totalFrames)/float(xSize))

           ;; gets overwritten later -- this is just to get the size
           ;;                           of a single image
           p_display = ptr_new(adt_make_image(0))

           ;;has the image already been zoomed? if they are trying to display a diffusion
           ;; weighted image which comes from state2 then there is no need to zoom the image b/c it has
           ;; already been done.

           IF project.procPramArray[CI].adt_display_type ne 4 then begin
               mas_zoom, p_display

			; *************************************
			; HS - 20061111
			; Now we have to rotate every image except the diffusion weighted imaged
			; *************************************

               mas_rotate_flip, p_display
           ENDIF

           sz_image = size(*p_display)
           ;; now that we have the size, we can free this ptr
           ptr_free, p_display

		; ========= HAVE TO CHANGE RIGHT HERE TO INCREASE THE WINDOW SIZE
		;             if sz_image[0] eq 3 then begin
		;    ;nr            Window, /FREE, Xsize=xSize *sz_image[1], Ysize=ySize *sz_image[2],  $
		;    ;                TITLE=title
		;             end else begin
		;    ;nr            Window, /FREE, Xsize=xSize *sz_image[1], Ysize=ySize *sz_image[2],  $
		;    ;                TITLE=title
		;             end

		; HS - 20061108
		; Create a window that can fit all the images being displayed.
		; I will first get the dimensions of displaying half the images
		; in the xdirection and the remaining on the next row.
		; If this turns out to be more than the x axis screen dimension
		; I will instead divide the images in 3 columns.
           screen_dimensions = GET_SCREEN_SIZE(RESOLUTION=resolution)

;;           xSize = sz_image[1]*ceil(totalFrames/2)
;;           ySize = sz_image[2]*fix(totalFrames/(totalFrames/2))
;;
;;           if xSize ge screen_dimensions[0] then begin
;;               xSize = sz_image[1]*ceil(totalFrames/3)
;;               ySize = sz_image[2]*fix(totalFrames/(totalFrames/3))
;;           endif

;;;;;;;;;;;;;;;
           ;; # of images that will fit in x direction
           xsize = floor(float(screen_dimensions[0])*0.9/sz_image[1])
           ysize = ceil(float(totalframes)/float(xsize))
           ;; y-size of the scroll viewport
           ysc_size = min([screen_dimensions[1]-150, ysize*sz_image[2]])

           base = widget_base(/tlb_size_events, title='ADT Multiple Slice Display')

	   dw = widget_draw(base, $
                            x_scroll_size=xsize*sz_image[1]+3, y_scroll_size=ysc_size, $
                            xsize=xsize*sz_image[1], ysize=ysize*sz_image[2], $
                            /scroll, retain=2)

           draw = ptr_new(dw, /no_copy)
;**           widget_control, base, /realize
;**           widget_control, base, set_uvalue=draw
;**           widget_control, *draw, get_value=drawwin
;**           XMANAGER, 'draw_scroll', base, /NO_BLOCK, $
;**             event_handler='_draw_scroll_h', cleanup='_draw_scroll_c'

;**           XMANAGER, 'draw_scroll_cnt', *draw, /NO_BLOCK
;**           wset, drawwin
;;;;;;;;;;;;;;;;

           if (sz_image[0] eq 3) then begin
               im = bytarr(xsize*sz_image[1], ysize*sz_image[2], 3)
           endif else begin
               im = bytarr(xsize*sz_image[1], ysize*sz_image[2])
           endelse

           for ii=0,totalFrames-1 do begin
               p_display = ptr_new(adt_make_image(ii))

               if project.procPramArray[CI].adt_display_type ne 4 then begin
                   mas_zoom, p_display
                   mas_rotate_flip, p_display
               endif

               xpos = (ii mod xsize)  * sz_image[1]
               ypos = floor(ii/xsize) * sz_image[2]
               if (sz_image[0] eq 3) then begin
                   im[xpos:xpos+sz_image[1]-1, ypos:ypos+sz_image[2]-1,*] = reverse((*p_display)[*,*,*], 2)
               endif else begin
                   im[xpos:xpos+sz_image[1]-1, ypos:ypos+sz_image[2]-1] = reverse((*p_display)[*,*], 2)
               endelse

           end

           if (sz_image[0] eq 3) then begin
               tmp = bytarr(3, xsize*sz_image[1], ysize*sz_image[2])
               tmp[0,*,*] = im[*,*,0]
               tmp[1,*,*] = im[*,*,1]
               tmp[2,*,*] = im[*,*,2]
               im = tmp
               im = temporary(reverse(im, 3))
               orient_circle = circle(100,100,1)
               pcircle = ptr_new(orient_circle, /no_copy)
               mas_rotate_flip, pcircle
           endif else begin
               im = temporary(reverse(im, 2))
               pcircle = ptr_new()
           endelse

;**           im = temporary(reverse(im, 2))
;**           if sz_image[0] eq 3 then tv, im, true=3 else tv, im

           mas_display_multi, im, tab_title='ADT Multislice', /standalone, fov_units=0, $
             extra_thumbnail = (ptr_valid(pcircle)) ? temporary(*pcircle) : 0

           loadct, 0,  /SILENT

           return

			; CD - this code does not call Display_DrawRes_Create below
			; returns above, thus do not create p_raw (floating point data)
			; to be passed to the display routine
       END

   END ELSE BEGIN
        ;sheet display
        ;we have to make sure that we are selecting the right slice and that
        ;fov is initilized correctly
       IF slice_axis EQ 0 THEN BEGIN
           slice_number = sdim_start
           x_fov = ffov
           y_fov = pfov
       END ELSE IF slice_axis EQ 1 THEN BEGIN
           slice_number = pdim_start
           x_fov = ffov
           y_fov = sfov
       END ELSE IF slice_axis EQ 2 THEN BEGIN
           slice_number = fdim_start
           x_fov = pfov
           y_fov = sfov
       END

       p_display = ptr_new(adt_make_sheet(slice_number))
       disp_float_labels = 0

		; HS 20061113
		; There is no need to rotate things anymore
       adt_rotate_flip_display, p_display
       adt_zoom_display, p_display

   ENDELSE

   sz_sheet = size((*p_display))

    ;does the user want to diplay using iImage or window
   IF project.iImage_flag EQ 1 THEN BEGIN

       iImage, (*p_display) $
         ,IMAGE_DIMENSIONS=[x_fov,y_fov] $
         ,xtitle='cm' $
         ,ytitle='cm' $
         ,title = display_title

   END ELSE BEGIN
		;CD - use enhanced window in mas_display
		;need to find out what floating point to display

       if ( (size(*p_display))[0] eq 3 ) then begin
           orient_circle = circle(100,100,1)
           pcircle = ptr_new(orient_circle, /no_copy)
           mas_rotate_flip, pcircle
       endif else begin
           pcircle = ptr_new()
       endelse

       IF disp_float_labels EQ 1 THEN BEGIN
;           Display_DrawResS_Create, *p_display, title, ffov, pfov, FPVALS=*p_raw
           mas_display_multi, *p_display, $
             tab_title=title, $
             fov_x=ffov, fov_y=pfov, fov_units=1, $
             fp_vals=*p_raw, /standalone, $
             extra_thumbnail = (ptr_valid(pcircle)) ? temporary(*pcircle) : 0
       ENDIF ELSE BEGIN
;           Display_DrawResS_Create, *p_display, title, ffov, pfov
           mas_display_multi, *p_display, $
             tab_title=title, $
             fov_x=ffov, fov_y=pfov, fov_units=1, $
             /standalone, $
             extra_thumbnail = (ptr_valid(pcircle)) ? temporary(*pcircle) : 0
       ENDELSE
   END

   ptr_free, p_display

END


; Subroutine name: adt_image_statistics_event
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.


PRO adt_image_statistics_event, event

    COMMON scan_data
    COMMON common_widgets

    ;error handling.

    CI = project.ci

    case Event.id of

      bS0: begin
     project.procPramArray[CI].adt_img_stat_flg_array[0] = event.select
    end
       bXX: begin
     project.procPramArray[CI].adt_img_stat_flg_array[1] = event.select
    end
    bYY: begin
     project.procPramArray[CI].adt_img_stat_flg_array[2] = event.select
    end
    bZZ: begin
     project.procPramArray[CI].adt_img_stat_flg_array[3] = event.select
    end
    bXY: begin
     project.procPramArray[CI].adt_img_stat_flg_array[4] = event.select
    end
    bXZ: begin
     project.procPramArray[CI].adt_img_stat_flg_array[5] = event.select
    end
    bYZ: begin
     project.procPramArray[CI].adt_img_stat_flg_array[6] = event.select
    end
    bAD: begin
     project.procPramArray[CI].adt_img_stat_flg_array[7] = event.select
    end
       bFA: begin
     project.procPramArray[CI].adt_img_stat_flg_array[8] = event.select
    end
    bE1: begin
     project.procPramArray[CI].adt_img_stat_flg_array[9] = event.select
    end
    bE2: begin
     project.procPramArray[CI].adt_img_stat_flg_array[10] = event.select
    end
    bE3: begin
     project.procPramArray[CI].adt_img_stat_flg_array[11]  = event.select
    end
    endcase
end


; Subroutine name: adt_do_image_histogram
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.


pro adt_do_image_histogram, ev
    COMPILE_OPT IDL2
    ;bring in the appropriate global variables
    COMMON scan_data
    CI = project.ci

    sdim_Start = project.procPramArray[CI].sdim_Start
    pdim_start = project.procPramArray[CI].pdim_start
    fdim_start = project.procPramArray[CI].fdim_start

    adt_display_type =     project.procPramArray[CI].adt_display_type
    ADT_Start =    project.procPramArray[CI].ADT_Start
    EigenVec_Start =   project.procPramArray[CI].EigenVec_Start
    img_stat_flgs = project.procPramArray[CI].adt_img_stat_flg_array

    ;check to make sure that images have been fitted
    do_adt_regress_mono

    ;copy the image locally. this image will have the roi drawn on it
    ;process the image the way the user has set up.

    ;bring in all the images for 1 sheet of ADT display
    sz_adt =size(*project.dataArray[CI].adt)
    sz_frac_Ani = size(*project.dataArray[CI].frac_Ani )
    sz_Avg_Dif =  size(*project.dataArray[CI].Avg_Dif)
    sz_eign_Vec =  size(*project.dataArray[CI].eign_Vec)

    if sz_adt[0] eq 3 then begin
       ;single slice passed in


       flt_ADT = fltarr(sz_adt[1],sz_adt[2],12)

       flt_ADT[*,*,0:6] = (*project.dataArray[CI].adt)[*,*,*]
       flt_ADT[*,*,7]   = (*project.dataArray[CI].Avg_Dif)[*,*]
       flt_ADT[*,*,8]   = (*project.dataArray[CI].frac_Ani )[*,*]
       flt_ADT[*,*,9:11]= (*project.dataArray[CI].eign_val )[*,*,*]

    end
    if sz_adt[0] eq 4 then begin


       ;How did they slice the data
        ;freq. phase
        if project.procpramArray[ci].slice_axis eq 0 then begin
            flt_ADT = fltarr(sz_adt[1],sz_adt[2],12)

            flt_ADT[*,*,0:6] = (*project.dataArray[CI].adt)[*,*,sdim_start,*]
            flt_ADT[*,*,7]   = (*project.dataArray[CI].Avg_Dif)[*,*,sdim_start]
            flt_ADT[*,*,8]   = (*project.dataArray[CI].frac_Ani )[*,*,sdim_start]
            flt_ADT[*,*,9:11]= reform((*project.dataArray[CI].eign_val )[*,*,sdim_start,*])

        end

        ;freq. slice
        if project.procpramArray[ci].slice_axis eq 1 then begin
            flt_ADT = fltarr(sz_adt[1],sz_adt[3],12)

            flt_ADT[*,*,0:6] = reform((*project.dataArray[CI].adt)[*,pdim_start,*,*])
            flt_ADT[*,*,7]   = reform((*project.dataArray[CI].Avg_Dif)[*,pdim_start,*])
            flt_ADT[*,*,8]   = reform((*project.dataArray[CI].frac_Ani )[*,pdim_start,*])
            flt_ADT[*,*,9:11]= reform((*project.dataArray[CI].eign_val )[*,pdim_start,*,*])

        end

        ;phase slice
        if project.procpramArray[ci].slice_axis eq 2 then begin
            flt_ADT = fltarr(sz_adt[2],sz_adt[3],12)

            flt_ADT[*,*,0:6] = reform((*project.dataArray[CI].adt)[fdim_start,*,*,*])
            flt_ADT[*,*,7]   = reform((*project.dataArray[CI].Avg_Dif)[fdim_start,*,*])
            flt_ADT[*,*,8]   = reform((*project.dataArray[CI].frac_Ani )[fdim_start,*,*])
            flt_ADT[*,*,9:11]= reform((*project.dataArray[CI].eign_val )[fdim_start,*,*,*])

        end

    end

    ;zoom the image
    x_zoom = 1;project.procPramArray[ci].x_zoom
    y_zoom = 1;project.procPramArray[ci].y_zoom


    ;rotate, flip and zoom the image
    ;pflt_ADT = ptr_new(flt_ADT,/no_copy)
    ;mas_rotate_flip, pflt_ADT
    ;mas_zoom, pflt_ADT
    ;flt_adt = *pflt_ADT
    sz_flt_ADT = size(flt_adt)



    ;get the data object from the roi program.
;    widget_control, ev.id, get_uvalue=og
;    regions = og->get(/all)
    ;check to make sure that the there is an roi to draw on.
    if not ptr_valid(project.roi.pROIs[project.roi.ci]) then begin
        void =  dialog_message('Please either make an ROI or change selection',/error)
        return
    end
    regions = (*project.roi.pROIs[project.roi.ci])


    ;check to make sure the object regions is valid
    if not(total(obj_valid(regions))) then begin
       ;update_status_bar, 'No region returned from ROI tool'
       return
    end

    ;print, 'size of regions', size(regions)

    if     (size(regions))[0] eq 0 then numRoi = 1 $
    else numRoi =  (size(regions))[1]

    ;we need an array to hold the roi names
    name_array = strarr(numRoi)

    ; make a pointer array to hold all the masks
    mask = ptrarr(numRoi)

    for temp=0 , numRoi-1 do begin
       mask[temp] = ptr_new(regions[temp] -> ComputeMask( dimensions =[ sz_flt_ADT[1], sz_flt_ADT[2] ] , MASK_RULE=2))
       regions[temp] -> GETPROPERTY, name= name
       name_array[temp] = name
    end


    for maskCounter= 0 , numRoi-1 do begin
       for ii = 0 , 11 do begin


         if img_stat_flgs[ii] eq 1 then begin

          if ii eq 0 then begin ;SO
              windowName = 'S0 Histogram for '+name_array[maskCounter]

              plot_histogram, flt_ADT[*,*,ii], *(mask[maskCounter]), windowName, name_array[maskCounter]
          end
          if ii eq 1 then begin ;XX
              windowName = 'XX Histogram for '+name_array[maskCounter]

              plot_histogram, flt_ADT[*,*,ii], *(mask[maskCounter]), windowName, name_array[maskCounter]

          end
          if ii eq 2 then begin ;YY
              windowName = 'YY Histogram for '+name_array[maskCounter]

              plot_histogram, flt_ADT[*,*,ii], *(mask[maskCounter]), windowName, name_array[maskCounter]

          end
          if ii eq 3 then begin ;ZZ
              windowName = 'ZZ Histogram for '+name_array[maskCounter]

              plot_histogram, flt_ADT[*,*,ii], *(mask[maskCounter]), windowName, name_array[maskCounter]

          end
          if ii eq 4 then begin ;XY
              windowName = 'XY Histogram for '+name_array[maskCounter]

              plot_histogram, flt_ADT[*,*,ii], *(mask[maskCounter]), windowName, name_array[maskCounter]

          end
          if ii eq 5 then begin ;XZ
              windowName = 'XZ Histogram for '+name_array[maskCounter]

              plot_histogram, flt_ADT[*,*,ii], *(mask[maskCounter]), windowName, name_array[maskCounter]

          end
          if ii eq 6 then begin ;YZ
              windowName = 'YZ Histogram for '+name_array[maskCounter]

              plot_histogram, flt_ADT[*,*,ii], *(mask[maskCounter]), windowName, name_array[maskCounter]

          end
          if ii eq 7 then begin ;AD
              windowName = 'AD Histogram for '+name_array[maskCounter]

              plot_histogram, flt_ADT[*,*,ii], *(mask[maskCounter]), windowName, name_array[maskCounter]

          end
          if ii eq 8 then begin ;FA
              windowName = 'FA Histogram for '+name_array[maskCounter]

              plot_histogram, flt_ADT[*,*,ii], *(mask[maskCounter]), windowName, name_array[maskCounter]

          end
          if ii eq 9 then begin ;E1
              windowName = 'E1 Histogram for '+name_array[maskCounter]

              plot_histogram, flt_ADT[*,*,ii], *(mask[maskCounter]), windowName, name_array[maskCounter]

          end
          if ii eq 10 then begin ;E2
              windowName = 'E2 Histogram for '+name_array[maskCounter]

              plot_histogram, flt_ADT[*,*,ii], *(mask[maskCounter]), windowName, name_array[maskCounter]

          end
          if ii eq 11 then begin ;E3
              windowName = 'E3 Histogram for '+name_array[maskCounter]

              plot_histogram, flt_ADT[*,*,ii], *(mask[maskCounter]), windowName, name_array[maskCounter]

          end
         endif
       endfor

    endfor

end


; Subroutine name: count_repetions
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This function will count the number of repetions in an array and
; then return values and their occurences.
; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

function count_repetions, array

    ;copy the array locally
    rep_array = array

    length =  (size(array))[1]
    if length lt 1 then return, -1

    ;make an array to hold the values of the repetions
    reps = fltarr(length)
    ;make an array to hold couter for how many times this value occured.
    count= lonarr(length)


    done = 0
    ii = 0

    while not done do begin

       ;store the value to remove
       reps[ii] = rep_array[0]

       ;print, 'check_value',reps[ii]
       ;print, 'size', size(rep_array)


       location = where( rep_array ne reps[ii] , count_temp, NCOMPLEMENT= n_complement )

       ;store how many were removed.
       count[ii] = n_complement

       ;print, 'number of reps =', count[ii]

       if count_temp eq 0 then begin
         done = 1
         break
       end

       rep_array = rep_array[location]
       ii++

    end

    ;print, reps
    ;print, 'num pixels',total(count)

    return, [[reps[0:ii]],[count[0:ii]]]
end


; Subroutine name: mask_indices
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This function will locate the duplicated mask array indices and remove then and return the
; pure array indices.
; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

function mask_indices , image, mask , x_zoom, y_zoom
        ind = ARRAY_INDICES(image , mask )
        ;ind = float(ind)
        ind[0,*] /= float(x_zoom)
        ind[1,*] /= float(y_zoom)
       ;now remove the duplicates from the indices array x and y component.

    ii=0
    done = 0

    ;make an array to hold all the values.
    x_indices = intarr((size(ind))[2])
    y_indices = intarr((size(ind))[2])
    x_ind = ind[0,*]
    y_ind = ind[1,*]


    while not(done) do begin
         ;print, 'looking for ',ind[0,0],ind[1,0]

         x_indices[ii] = x_ind[0]
         y_indices[ii] = y_ind[0]

         ;find all the duplicates for the first entry
         location = where(x_indices[ii] eq x_ind and y_indices[ii] eq y_ind, COMPLEMENT= comp, count )


       ; print, 'size comp',

         ;if there were none found then we are at the end of the array.
         if (size( comp))[0] eq 0 then begin
          done = 1
          break
         end
         ;remove the duplicates
         x_ind = x_ind[comp]
         y_ind = y_ind[comp]


         ii++


    end

;   print, 'x', x_indices[0:ii]
;   print, 'y', y_indices[0:ii]
    return, [[x_indices[0:ii]],[y_indices[0:ii]]]
end


; Subroutine name: add_string_to_array
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This procedure will take in an array and change its size by 1 then add in the next string,'str'
; This is designed for a 1d string array.
; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.


pro add_string_to_array, array, str
    sz_arr = (size(array))[1]
    new_arr = strarr(sz_arr+1)
    new_arr = [[array],[str]]
    array = new_arr
    ;print, array
end

; Subroutine name: adt_subset_selected_slice
; Created by: Bill Triplett 20080523
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
;   grab all the ADT information for a slice

; Editing Information:

function adt_subset_selected_slice, want_slice=want_slice, native=native

    compile_opt idl2
    common scan_data

    ci = project.ci

    params = ptr_new(project.procpramarray[ci])

    sdim_Start       = (*params).sdim_Start
    pdim_start       = (*params).pdim_start
    fdim_start       = (*params).fdim_start
    adt_display_type = (*params).adt_display_type
    ADT_Start        = (*params).ADT_Start
    EigenVec_Start   = (*params).EigenVec_Start
    img_stat_flgs    = (*params).adt_img_stat_flg_array
    x_zoom           = (*params).x_zoom
    y_zoom           = (*params).y_zoom


    ;check to make sure that images have been fitted
    do_adt_regress_mono

    ;copy the image locally. this image will have the roi drawn on it
    ;process the image the way the user has set up.

    ;bring in all the images for 1 sheet of ADT display
    sz_adt      = size(*project.dataArray[CI].adt)
    sz_frac_Ani = size(*project.dataArray[CI].frac_Ani )
    sz_Avg_Dif  = size(*project.dataArray[CI].Avg_Dif)
    sz_eign_Val = size(*project.dataArray[CI].eign_Val)

    if sz_adt[0] eq 3 then begin
       ;single slice passed in
       flt_ADT = fltarr(sz_adt[1],sz_adt[2],12)

       flt_ADT[*,*,0:6]   = (*project.dataArray[CI].adt)[*,*,*]
       flt_ADT[*,*,7]   = (*project.dataArray[CI].Avg_Dif)[*,*]
       flt_ADT[*,*,8]   = (*project.dataArray[CI].frac_Ani )[*,*]
       flt_ADT[*,*,9:11]= (*project.dataArray[CI].eign_val )[*,*,*]

    end
    if sz_adt[0] eq 4 then begin

        ;How did they slice the data
        ;freq. phase
        if project.procpramArray[ci].slice_axis eq 0 then begin
            flt_ADT = fltarr(sz_adt[1],sz_adt[2],12)
            if (keyword_set(want_slice)) then sdim_start = want_slice
            flt_ADT[*,*,0:6] = (*project.dataArray[CI].adt)[*,*,sdim_start,*]
            flt_ADT[*,*,7]   = (*project.dataArray[CI].Avg_Dif)[*,*,sdim_start]
            flt_ADT[*,*,8]   = (*project.dataArray[CI].frac_Ani )[*,*,sdim_start]
            flt_ADT[*,*,9:11]= reform((*project.dataArray[CI].eign_val )[*,*,sdim_start,*])

        end

        ;freq. slice
        if project.procpramArray[ci].slice_axis eq 1 then begin
            flt_ADT = fltarr(sz_adt[1],sz_adt[3],12)
            if (keyword_set(want_slice)) then pdim_start = want_slice
            flt_ADT[*,*,0:6] = reform((*project.dataArray[CI].adt)[*,pdim_start,*,*])
            flt_ADT[*,*,7]   = reform((*project.dataArray[CI].Avg_Dif)[*,pdim_start,*])
            flt_ADT[*,*,8]   = reform((*project.dataArray[CI].frac_Ani )[*,pdim_start,*])
            flt_ADT[*,*,9:11]= reform((*project.dataArray[CI].eign_val )[*,pdim_start,*,*])

        end

        ;phase slice
        if project.procpramArray[ci].slice_axis eq 2 then begin
            flt_ADT = fltarr(sz_adt[2],sz_adt[3],12)
            if (keyword_set(want_slice)) then fdim_start = want_slice
            flt_ADT[*,*,0:6] = reform((*project.dataArray[CI].adt)[fdim_start,*,*,*])
            flt_ADT[*,*,7]   = reform((*project.dataArray[CI].Avg_Dif)[fdim_start,*,*])
            flt_ADT[*,*,8]   = reform((*project.dataArray[CI].frac_Ani )[fdim_start,*,*])
            flt_ADT[*,*,9:11]= reform((*project.dataArray[CI].eign_val )[fdim_start,*,*,*])

        end

    end

    ;rotate, flip and zoom the image
    pflt_ADT = ptr_new(flt_ADT,/no_copy)
    if not keyword_set(native) then begin
        mas_rotate_flip, pflt_ADT
        mas_zoom, pflt_ADT
    endif

    return, pflt_ADT
end

pro adt_export_roi_data, event

    compile_opt idl2
    common scan_data
    ci = project.ci

    pflt_adt = adt_subset_selected_slice(native=(1-project.procpramarray[ci].no_transform_roi))
    sz_flt_adt = size(*pflt_adt)
    img_stat_flgs = project.procPramArray[CI].adt_img_stat_flg_array
;;MASKCHANGE
    x_zoom = 1;project.procPramArray[CI].x_zoom
    y_zoom = 1;project.procPramArray[CI].y_zoom

    curr_slice = project.procPramArray[CI].sdim_start
    curr_slice = (curr_slice lt 10) ? '0'+string(curr_slice) : string(curr_slice)
    temp   = strsplit(project.imndArray[ci].file_path, (get_dir_slash())[0], /extract)

    ;; these are the only two image_types that should be calling this procedure
    case project.imndarray[ci].image_type of

        11: file_name = temp[n_elements(temp)-1]
        3 : file_name = temp[n_elements(temp)-2]
        else: file_name = project.imndArray[ci].file_path

    endcase

    ;;check to make sure that the there is an roi to draw on.
    if not ptr_valid(project.roi.pROIs[project.roi.ci]) then begin
        void =  dialog_message('Please either make an ROI or change selection',/error)
        return
    end
    regions = (*project.roi.pROIs[project.roi.ci])

    ;; check to make sure the object regions is valid
    if not(total(obj_valid(regions))) then begin
        update_status_bar, 'No region returned from ROI tool'
        return
    end

    if (size(regions))[0] eq 0 then begin
        numRoi = 1
    endif else begin
        numRoi = (size(regions))[1]
    endelse

    ;make an array to hold the name of the roi's
    name_array = strarr(numRoi)

    ; make a pointer array to hold all the masks
    mask = ptrarr(numRoi)
    dim  = [ sz_flt_ADT[1], sz_flt_ADT[2] ]
    for temp=0 , numRoi-1 do begin
;;        mask[temp] = mas_roi_get_current_mask(dim, region_num=temp)
        mask[temp] = ptr_new(regions[temp]->computeMask(dimensions=dim,$
                                                        MASK_RULE=2) )
        regions[temp]->getProperty, name=name
        name_array[temp] = name
    end

    ;; header labels
    lblarr = ['S0', 'XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ', 'AD', 'FA', 'E1', 'E2', 'E3']

    ;; set up the header
    header = ['ROI', 'x', 'y', lblarr[where(img_stat_flgs ne 0)]]
    str_data = strjoin(header, string(9B))


    for maskCounter=0,numRoi-1 do begin

        maskPixels = where(*(mask[maskCounter]) gt 0, count)
        if (count eq 0) then continue ;; mask contains no pixels

        array_indices_of_mask = mask_indices((*pflt_adt)[*,*,0], $
                                             maskPixels, $
                                             x_zoom, y_zoom)

        n_coords = (size(array_indices_of_mask))[1]

        str = strarr(3 + (size(where(img_stat_flgs ne 0)))[1])

        ;; build the string array line-by-line
        for coord = 0,n_coords-1 do begin

            x = array_indices_of_mask[coord,0]
            y = array_indices_of_mask[coord,1]

            str[0] = name_array[maskCounter]
            str[1] = string(x)
            str[2] = string(y)

            str_mark = 0
            for ii=0, n_elements(lblarr)-1 do begin
                if img_stat_flgs[ii] eq 0 then continue
                str[str_mark + 3] = string((*pflt_adt)[x,y,ii])
                str_mark++
            endfor

            str_data = [str_data,  strjoin(strcompress(str, /remove_all), string(9B))]

        endfor

    endfor

    display_stats, transpose(str_data), strcompress(file_name+'_'+curr_slice, /remove_all)

end

; Subroutine name: adt_do_image_statistics
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.
    ;Replaced with new version 2009/06/28


pro adt_do_image_statistics, ev

    COMPILE_OPT IDL2
    ;bring in the appropriate global variables
    COMMON scan_data
    CI = project.ci

    sdim_Start = project.procPramArray[CI].sdim_Start
    pdim_start = project.procPramArray[ci].pdim_start
    fdim_start = project.procPramArray[ci].fdim_start
    sdim = project.iMNDArray[CI].sdim
    pdim = project.IMNDArray[ci].pdim
    fdim = project.IMNDArray[ci].fdim
    adt_display_type = project.procPramArray[CI].adt_display_type
    ADT_Start = project.procPramArray[CI].ADT_Start
    EigenVec_Start = project.procPramArray[CI].EigenVec_Start
    img_stat_flgs = project.procPramArray[CI].adt_img_stat_flg_array
    x_zoom = 1;project.procPramArray[CI].x_zoom
    y_zoom = 1;project.procPramArray[CI].y_zoom

    adt_param_names = ['S0','XX','YY','ZZ','XY','XZ','YZ','AD','FA','E1','E2','E3']
    
    ;before we process any info we need to see if anybody checked an image they want the
    ;roi maped to.
    selected_params = where (img_stat_flgs eq 1,n_selected_images)
    if n_selected_images eq 0 then begin
       void = dialog_message('Please select an image to map the ROI onto')
       return
    end

    ;;this section tell whether or not the user pressed the detailed
    ;;stats button and if so give them detailed info.
    wWidget =  ev.top
    detailed_flag = 0
    if ev.id eq Widget_Info(wWidget, FIND_BY_UNAME='Print Stats Detailed') then begin
       detailed_flag = 1
    end

    ;;check to make sure that images have been fitted
    do_adt_regress_mono

    ;;this is just to check to make sure they now they are fitting data that is not
    ;;in the direction they are thinking.
    ;;make the image
    case project.procpramarray[ci].slice_axis of
       0: begin
          void = adt_make_image (sdim_start)
          if (project.procpramarray[ci].single_multi_flag) then begin
             sl_start = 0
             sl_end = sdim-1
          endif else begin
             sl_start = sdim_start
             sl_end = sdim_start
          endelse
          sl_selected = sdim_start
       end
       1: begin
          void = adt_make_image (pdim_start)
          if (project.procpramarray[ci].single_multi_flag) then begin
             sl_start = 0
             sl_end = pdim-1
          endif else begin
             sl_start = pdim_start
             sl_end = pdim_start
          endelse
          sl_selected = pdim_start
       end
       2: begin
          void = adt_make_image (fdim_start)
          if (project.procpramarray[ci].single_multi_flag) then begin
             sl_start = 0
             sl_end = fdim-1
          endif else begin
             sl_start = fdim_start
             sl_end = fdim_start
          endelse
          sl_selected = fdim_start
       end
    endcase

    ;copy the image locally. this image will have the roi drawn on it
    ;process the image the way the user has set up.

    ;bring in all the images for 1 sheet of ADT display
    sz_adt      = size(*project.dataArray[CI].adt)
    sz_frac_Ani = size(*project.dataArray[CI].frac_Ani )
    sz_Avg_Dif  = size(*project.dataArray[CI].Avg_Dif)
    sz_eign_Val = size(*project.dataArray[CI].eign_Val)

    pflt_ADT   = adt_subset_selected_slice(native=(1-project.procpramarray[ci].no_transform_roi))
    flt_ADT    = *pflt_ADT
    ptr_free, pflt_ADT
    sz_flt_adt = size(flt_ADT)

    ;get the data object from the roi program.
;    widget_control, ev.id, get_uvalue=og
;    regions = og->get(/all)

    ;check to make sure that the there is an roi to draw on.
    if not ptr_valid(project.roi.pROIs[project.roi.ci]) then begin
        void =  dialog_message('Please either make an ROI or change selection',/error)
        return
    end
    regions = (*project.roi.pROIs[project.roi.ci])

    ;check to make sure the object regions is valid
    if not(total(obj_valid(regions))) then begin
       update_status_bar, 'No region returned from ROI tool'
       return
    end

    ;print, 'size of regions', size(regions)

    if     (size(regions))[0] eq 0 then numRoi = 1 $
    else numRoi =  (size(regions))[1]

    ;make an array to hold the name of the roi's
    roi_name_array = strarr(numRoi)
    param_name_array = adt_param_names[selected_params]
    
    ; make a pointer array to hold all the masks
    mask = ptrarr(numRoi)
    dim = [ sz_flt_ADT[1], sz_flt_ADT[2] ]

    for temp=0 , numRoi-1 do begin
        if (not obj_valid(regions[temp])) then continue
;;        mask[temp] = mas_roi_get_current_mask(dim, region_num=temp)
        mask[temp] = ptr_new(regions[temp]->ComputeMask(dimensions=dim, $
                                                        MASK_RULE=2))

        regions[temp] -> GETPROPERTY, name= name
        roi_name_array[temp] = name
    end

    ;setting the size of the display string array
    ;DSAC =  display string array counter to keep place where the string
    ;belongs in the array. after a new line is added the counter should increase by 1
    ;DDSAC =  detailed display string array counter to keep place where the string
    ;belong in the array, after a new line is added the coutner should increase by 1
    SL = sdim_start

    ;; [ col = [mean, std, sum, sum sq, min, max, var], row = adt parameter * numRoi ]
    stat_array = fltarr(8, n_selected_images, numRoi, SL_end-SL_start+1)
       

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    sl_counter = 0
    for SL = SL_Start, SL_end do begin

       pflt_ADT   = adt_subset_selected_slice(native=(1-project.procpramarray[ci].no_transform_roi), want_slice=SL)
       flt_ADT    = *pflt_ADT
       ptr_free, pflt_ADT

       ;; loop over each ADT parameter the user selected
       param_ct = 0
       for ii = 0, n_elements(img_stat_flgs)-1 do begin
       
           if (img_stat_flgs[ii] eq 0) then continue
           
           ;;now loop throught all possible images and the mask counter
           for maskCounter = 0, numRoi-1 do begin

              maskPixels = where( *(mask[maskCounter]) gt 0, count)
              if count eq 0 then continue

                IMAGE_STATISTICS, flt_ADT[*,*,ii], mask=*(mask[maskCounter])  $
                               , MEAN=image_mean $
                               , STDDEV=image_stddev $
                               , DATA_SUM=image_sum $
                               , SUM_OF_SQUARES=image_sum_of_squares $
                               , MINIMUM=image_min $
                               , MAXIMUM=image_max $
                               , VARIANCE=image_var

                stat_array[*, param_ct, maskCounter, sl_counter] = [ image_mean, image_stddev, image_sum, image_sum_of_squares, image_min, image_max, image_var, count]
                
            endfor
            
            param_ct++
            
       endfor
       
       sl_counter++
       
    endfor
    
    print_format = '(E0, E0, E0, E0, E0, E0, E0, E0)'
    htab         = string(09B)
    crlf         = string(13B)+string(10B)
    header       = 'Slice'+htab+'ADT Parameter'+htab+'Region Name'+htab+ $
                    'Mean'+htab+'Std. Dev.'+htab+'Sum'+htab+'Sum of Squares'+htab+ $
                    'Min'+htab+'Max'+htab+'Variance'+htab+'Count'
                    
    output_lines = ['File: '+project.imndarray[ci].file_path, '']
    line         = 0

    ;; Param => Method => ROI
    for slice = 0, sl_end-sl_start do begin
        for param = 0, n_elements(param_name_array)-1 do begin
            output_lines = [output_lines, header]
            for roi = 0, n_elements(roi_name_array)-1 do begin
                pref = string(slice+sl_start, format='(I0)')+htab+param_name_array[param]+htab+roi_name_array[roi]
                stat = strjoin(string(stat_array[*, param, roi, slice], format='(E0)'), htab)
                output_lines = [ output_lines, pref+htab+stat ]
            endfor
            output_lines = [output_lines, '']
         endfor
         output_lines = [output_lines, '']
     endfor
    
    display_stats, output_lines, 'ADT Statistics'
    ;;print, stat_array, format=print_format

end

; Subroutine name: adt_do_image_statistics
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.
    ;Replaced by BT 2009/07/28

;pro adt_do_image_statistics_old, ev
;    COMPILE_OPT IDL2
;    ;bring in the appropriate global variables
;    COMMON scan_data
;    CI = project.ci
;
;    sdim_Start = project.procPramArray[CI].sdim_Start
;    pdim_start = project.procPramArray[ci].pdim_start
;    fdim_start = project.procPramArray[ci].fdim_start
;    sdim = project.iMNDArray[CI].sdim
;    pdim = project.IMNDArray[ci].pdim
;    fdim = project.IMNDArray[ci].fdim
;    adt_display_type = project.procPramArray[CI].adt_display_type
;    ADT_Start = project.procPramArray[CI].ADT_Start
;    EigenVec_Start = project.procPramArray[CI].EigenVec_Start
;    img_stat_flgs = project.procPramArray[CI].adt_img_stat_flg_array
;    x_zoom = 1;project.procPramArray[CI].x_zoom
;    y_zoom = 1;project.procPramArray[CI].y_zoom
;
;    ;before we process any info we need to see if anybody checked an image they want the
;    ;roi maped to.
;    void = where (img_stat_flgs eq 1,n_selected_images)
;    if n_selected_images eq 0 then begin
;       void = dialog_message('Please select an image to map the ROI onto')
;       return
;    end
;
;    ;;this section tell whether or not the user pressed the detailed
;    ;;stats button and if so give them detailed info.
;    wWidget =  ev.top
;    detailed_flag = 0
;    if ev.id eq Widget_Info(wWidget, FIND_BY_UNAME='Print Stats Detailed') then begin
;       detailed_flag = 1
;    end
;
;    ;;check to make sure that images have been fitted
;    do_adt_regress_mono
;
;    ;;this is just to check to make sure they now they are fitting data that is not
;    ;;in the direction they are thinking.
;    ;;make the image
;    case project.procpramarray[ci].slice_axis of
;       0: begin
;          void = adt_make_image (sdim_start)
;          if (project.procpramarray[ci].single_multi_flag) then begin
;             sl_start = 0
;             sl_end = sdim-1
;          endif else begin
;             sl_start = sdim_start
;             sl_end = sdim_start
;          endelse
;          sl_selected = sdim_start
;       end
;       1: begin
;          void = adt_make_image (pdim_start)
;          if (project.procpramarray[ci].single_multi_flag) then begin
;             sl_start = 0
;             sl_end = pdim-1
;          endif else begin
;             sl_start = pdim_start
;             sl_end = pdim_start
;          endelse
;          sl_selected = pdim_start
;       end
;       2: begin
;          void = adt_make_image (fdim_start)
;          if (project.procpramarray[ci].single_multi_flag) then begin
;             sl_start = 0
;             sl_end = fdim-1
;          endif else begin
;             sl_start = fdim_start
;             sl_end = fdim_start
;          endelse
;          sl_selected = fdim_start
;       end
;    endcase
;
;    ;copy the image locally. this image will have the roi drawn on it
;    ;process the image the way the user has set up.
;
;    ;bring in all the images for 1 sheet of ADT display
;    sz_adt      = size(*project.dataArray[CI].adt)
;    sz_frac_Ani = size(*project.dataArray[CI].frac_Ani )
;    sz_Avg_Dif  = size(*project.dataArray[CI].Avg_Dif)
;    sz_eign_Val = size(*project.dataArray[CI].eign_Val)
;
;    pflt_ADT   = adt_subset_selected_slice(native=(1-project.procpramarray[ci].no_transform_roi))
;    flt_ADT    = *pflt_ADT
;    ptr_free, pflt_ADT
;    sz_flt_adt = size(flt_ADT)
;
;    ;get the data object from the roi program.
;;    widget_control, ev.id, get_uvalue=og
;;    regions = og->get(/all)
;
;    ;check to make sure that the there is an roi to draw on.
;    if not ptr_valid(project.roi.pROIs[project.roi.ci]) then begin
;        void =  dialog_message('Please either make an ROI or change selection',/error)
;        return
;    end
;    regions = (*project.roi.pROIs[project.roi.ci])
;
;    ;check to make sure the object regions is valid
;    if not(total(obj_valid(regions))) then begin
;       update_status_bar, 'No region returned from ROI tool'
;       return
;    end
;
;    ;print, 'size of regions', size(regions)
;
;    if     (size(regions))[0] eq 0 then numRoi = 1 $
;    else numRoi =  (size(regions))[1]
;
;    ;make an array to hold the name of the roi's
;    name_array = strarr(numRoi)
;
;    ; make a pointer array to hold all the masks
;    mask = ptrarr(numRoi)
;    dim = [ sz_flt_ADT[1], sz_flt_ADT[2] ]
;
;    for temp=0 , numRoi-1 do begin
;        if (not obj_valid(regions[temp])) then continue
;;;        mask[temp] = mas_roi_get_current_mask(dim, region_num=temp)
;        mask[temp] = ptr_new(regions[temp]->ComputeMask(dimensions=dim, $
;                                                        MASK_RULE=2))
;
;        regions[temp] -> GETPROPERTY, name= name
;        name_array[temp] = name
;    end
;
;    ;setting the size of the display string array
;    ;DSAC =  display string array counter to keep place where the string
;    ;belongs in the array. after a new line is added the counter should increase by 1
;    ;DDSAC =  detailed display string array counter to keep place where the string
;    ;belong in the array, after a new line is added the coutner should increase by 1
;    SL = sdim_start
;
;    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;    for SL = SL_Start, SL_end do begin
;
;       pflt_ADT   = adt_subset_selected_slice(native=(1-project.procpramarray[ci].no_transform_roi), want_slice=SL)
;       flt_ADT    = *pflt_ADT
;       ptr_free, pflt_ADT
;
;       DSAC = 0
;       DDSAC = 0
;
;       display_string_array = strarr(1+numRoi*(n_selected_images+2))
;       detailed_string_array = strarr(1)
;
;       display_string_array[DSAC] = 'Image'+string([09B])$
;                                    +'mean'+string([09B])+'stddev'+string([09B])$
;                                    +'sum'+string([09B])+'sumOfSquares'+string([09B])$
;                                    +'min'+string([09B])+'max'+string([09B])+'var'
;       DSAC = DSAC+1
;
;       ;;now loop throught all possible images and the mask counter
;       maskCounter = 0
;       for maskCounter= 0 , numRoi-1 do begin
;          display_string_array[DSAC] = ''
;          DSAC = DSAC+1
;          display_string_array[DSAC] = name_array[maskCounter]+string([09B])+ $
;                                       'Slice '+strcompress(string(SL), /remove_all)+string([09B])+$
;                                       project.imndarray[ci].file_path
;          DSAC = DSAC+1
;
;          add_string_to_array, detailed_string_array, name_array[maskCounter]+string([09B])+$
;                               'Slice '+strcompress(string(SL), /remove_all)+string([09B])+$
;                               project.imndarray[ci].file_path
;
;          ;;this is a flag to tell whether or not the location array has been calculated for the
;          ;;current region.
;          generate_location_once = 0
;
;
;          for ii = 0 , 11 do begin
;
;             if img_stat_flgs[ii] eq 0 then continue
;             maskPixels = where( *(mask[maskCounter]) gt 0, count)
;             if count eq 0 then continue
;
;             IMAGE_STATISTICS, flt_ADT[*,*,ii], mask=*(mask[maskCounter])  $
;                               , MEAN=image_mean $
;                               , STDDEV=image_stddev $
;                               , DATA_SUM=image_sum $
;                               , SUM_OF_SQUARES=image_sum_of_squares $
;                               , MINIMUM=image_min $
;                               , MAXIMUM=image_max $
;                               , VARIANCE=image_var
;
;             temp = flt_ADT[*,*,ii]
;             values_and_count = count_repetions(temp[maskPixels])
;
;
;             ;;kyle padgett wants to divede the number of count by (x zoom and y zoom)
;             ;;print,values_and_count[*,1]
;             values_and_count[*,1] /= float(x_zoom*y_zoom)
;             ;;print,values_and_count[*,1]
;
;
;             ;;now kyle padgett wants to look at the location of where these pixels are.
;             ;;but only print it once per roi.
;             if generate_location_once eq 0 then begin
;                array_indices_of_mask = mask_indices (temp, maskPixels, x_zoom, y_zoom)
;                x_loc = strtrim(array_indices_of_mask[*,0],2)
;                x_loc = 'x pos'+STRING(9B)+STRJOIN(x_loc,STRING(9B))
;                add_string_to_array, detailed_string_array, x_loc
;                y_loc = strtrim(array_indices_of_mask[*,1],2)
;                y_loc = 'y pos'+STRING(9B)+STRJOIN(y_loc,STRING(9B))
;                add_string_to_array, detailed_string_array, y_loc
;
;                generate_location_once = 1
;             end
;
;             if ii eq 0 then begin
;                display_string_array[DSAC] = string( FORMAT='(%"S0\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e")' $
;                                                     ,image_mean , image_stddev , image_sum , image_sum_of_squares, image_min, image_max, image_var )
;                DSAC = DSAC+1
;
;                add_string_to_array, detailed_string_array, 'S0'
;             end
;             if ii eq 1 then begin
;                display_string_array[DSAC] = string( FORMAT='(%"XX\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e")' $
;                                                     ,  image_mean , image_stddev , image_sum , image_sum_of_squares, image_min, image_max, image_var )
;                DSAC = DSAC+1
;
;                add_string_to_array, detailed_string_array, 'XX'
;
;             end
;             if ii eq 2 then begin
;                display_string_array[DSAC] = string( FORMAT='(%"YY\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e")' $
;                                                     ,  image_mean , image_stddev , image_sum , image_sum_of_squares, image_min, image_max, image_var )
;                DSAC = DSAC+1
;
;                add_string_to_array, detailed_string_array, 'YY'
;             end
;             if ii eq 3 then begin
;                display_string_array[DSAC] = string( FORMAT='(%"ZZ\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e")' $
;                                                     ,  image_mean , image_stddev , image_sum , image_sum_of_squares, image_min, image_max, image_var )
;                DSAC = DSAC+1
;
;                add_string_to_array, detailed_string_array, 'ZZ'
;             end
;             if ii eq 4 then begin
;                display_string_array[DSAC] = string( FORMAT='(%"XY\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e")' $
;                                                     ,  image_mean , image_stddev , image_sum , image_sum_of_squares, image_min, image_max, image_var )
;                DSAC = DSAC+1
;
;                add_string_to_array, detailed_string_array, 'XY'
;             end
;             if ii eq 5 then begin
;                display_string_array[DSAC] = string( FORMAT='(%"XZ\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e")' $
;                                                     ,  image_mean , image_stddev , image_sum , image_sum_of_squares, image_min, image_max, image_var )
;                DSAC = DSAC+1
;
;                add_string_to_array, detailed_string_array, 'XZ'
;             end
;             if ii eq 6 then begin
;                display_string_array[DSAC] = string( FORMAT='(%"YZ\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e")' $
;                                                     ,  image_mean , image_stddev , image_sum , image_sum_of_squares, image_min, image_max, image_var )
;                DSAC = DSAC+1
;
;                add_string_to_array, detailed_string_array, 'YZ'
;             end
;             if ii eq 7 then begin
;               display_string_array[DSAC] = string( FORMAT='(%"AD\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e")' $
;                                                    ,  image_mean , image_stddev , image_sum , image_sum_of_squares, image_min, image_max, image_var )
;               DSAC = DSAC+1
;
;               add_string_to_array, detailed_string_array, 'AD'
;            end
;             if ii eq 8 then begin
;                display_string_array[DSAC] = string( FORMAT='(%"FA\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e")' $
;                                                     ,  image_mean , image_stddev , image_sum , image_sum_of_squares, image_min, image_max, image_var )
;                DSAC = DSAC+1
;
;                add_string_to_array, detailed_string_array, 'FA'
;             end
;             if ii eq 9 then begin
;                display_string_array[DSAC] = string( FORMAT='(%"E1\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e")' $
;                                                     ,  image_mean , image_stddev , image_sum , image_sum_of_squares, image_min, image_max, image_var )
;                DSAC = DSAC+1
;
;                add_string_to_array, detailed_string_array, 'E1'
;             end
;             if ii eq 10 then begin
;                display_string_array[DSAC] = string( FORMAT='(%"E2\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e")' $
;                                                     ,  image_mean , image_stddev , image_sum , image_sum_of_squares, image_min, image_max, image_var )
;                DSAC = DSAC+1
;
;                add_string_to_array, detailed_string_array, 'E2'
;             end
;             if ii eq 11 then begin
;                display_string_array[DSAC] = string( FORMAT='(%"E3\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e")' $
;                                                     ,  image_mean , image_stddev , image_sum , image_sum_of_squares, image_min, image_max, image_var )
;                DSAC = DSAC+1
;
;                add_string_to_array, detailed_string_array, 'E3'
;             end
;
;             ;;trim off the extra space for the values and add in tabs
;             values = strtrim(values_and_count[*,0],2)
;             values = 'values'+STRING(9B)+STRJOIN(values,STRING(9B))
;             add_string_to_array, detailed_string_array, values
;
;             ;;trim off the extra space for the count and add in tabs
;             count = strtrim(values_and_count[*,1],2)
;             count = 'Count'+STRING(9B)+STRJOIN(count,STRING(9B))
;             add_string_to_array, detailed_string_array, count
;
;          endfor
;
;       endfor
;
;       IF (n_elements(display_string_array_big) eq 0) then begin
;          display_string_array_big = display_string_array
;       ENDIF ELSE BEGIN
;          display_string_array_big = [ display_string_array_big, display_string_array ]
;       ENDELSE
;
;       IF (SL eq sl_selected) then begin
;          detailed_string_array_big = detailed_string_array
;       ENDIF
;
;    ENDFOR
;
;    display_stats, display_string_array_BIG, 'ADT Image Stats'
;
;    if detailed_flag then begin
;       display_stats, detailed_string_array_big, 'ADT Image Stats'
;    end
;
;   HEAP_GC
;
;end

; Subroutine name: adt_ROI
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This procedure will make a display for the roi mapping selection.

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

pro adt_ROI
    COMPILE_OPT IDL2
    HEAP_GC
    ;bring in the appropriate global variables
    COMMON scan_data
    COMMON common_widgets
    CI = project.ci

    ;If the ADT ROI window is already on the screen then dont make another one.
    IF N_ELEMENTS(ADT_ROI_WINDOW_BASE) EQ 1 THEN $
        IF     widget_info( ADT_ROI_WINDOW_BASE,  /VALID_ID ) eq 1 THEN return

    ;if the current scan does not have a good b-matrix then make every thing insensitive
    SENSITIVE =  ptr_valid (project.imndArray[CI].b_matrix)

    title = string ('ADT ROI Mapping Selection Scan #', strtrim(CI,2))
    ADT_ROI_WINDOW_BASE = widget_base(TITLE=title, $
            UVALUE = 'ADT_ROI_WINDOW_BASE', $
            XOFFSET=420 ,YOFFSET=637  )

    tlb=widget_base(ADT_ROI_WINDOW_BASE, /column)

    b=widget_button(tlb,value='Print Stats', event_pro='adt_do_image_statistics',uname='Print Stats', uvalue=og)
    b=widget_button(tlb,value='Print Stats Detailed', event_pro='adt_do_image_statistics',uname='Print Stats Detailed', uvalue=og, sensitive=0)
    b=widget_button(tlb,value='Export Selected Data', event_pro='adt_export_roi_data', uvalue=og)
    b=widget_button(tlb,value='Plot Histogram', event_pro='adt_do_image_histogram', uvalue=og)

    tlb2=widget_base(tlb, /row)

    ex_base_col1 = Widget_Base(tlb2, UNAME='ex_base_col1' ,COLUMN=1 ,/NONEXCLUSIVE)
    ex_base_col2 = Widget_Base(tlb2, UNAME='ex_base_col2' ,COLUMN=1 ,/NONEXCLUSIVE)
    ex_base_col3 = Widget_Base(tlb2, UNAME='ex_base_col3' ,COLUMN=1 ,/NONEXCLUSIVE)
    ex_base_col4 = Widget_Base(tlb2, UNAME='ex_base_col4' ,COLUMN=1 ,/NONEXCLUSIVE)



    bXX = Widget_Button(ex_base_col1, UNAME='XX' ,VALUE='XX', event_pro='adt_image_statistics_event')

    bXY = Widget_Button(ex_base_col1, UNAME='XY' ,VALUE='XY', event_pro='adt_image_statistics_event')
    bE1 = Widget_Button(ex_base_col1, UNAME='E1' ,VALUE='E1', event_pro='adt_image_statistics_event')


    bYY = Widget_Button(ex_base_col2, UNAME='YY' ,VALUE='YY', event_pro='adt_image_statistics_event')
    bXZ = Widget_Button(ex_base_col2, UNAME='XZ' ,VALUE='XZ', event_pro='adt_image_statistics_event')
    bE2 = Widget_Button(ex_base_col2, UNAME='E2' ,VALUE='E2', event_pro='adt_image_statistics_event')


    bZZ = Widget_Button(ex_base_col3, UNAME='ZZ' ,VALUE='ZZ', event_pro='adt_image_statistics_event')
    bYZ = Widget_Button(ex_base_col3, UNAME='YZ' ,VALUE='YZ', event_pro='adt_image_statistics_event')
    bE3 = Widget_Button(ex_base_col3, UNAME='E3' ,VALUE='E3', event_pro='adt_image_statistics_event')

    bS0 = Widget_Button(ex_base_col4, UNAME='S0' ,VALUE='S0', event_pro='adt_image_statistics_event')
    bFA = Widget_Button(ex_base_col4, UNAME='FA' ,VALUE='FA', event_pro='adt_image_statistics_event')
    bAD = Widget_Button(ex_base_col4, UNAME='AD' ,VALUE='AD', event_pro='adt_image_statistics_event')

    if project.procPramArray[CI].adt_img_stat_flg_array[0] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bS0
    if project.procPramArray[CI].adt_img_stat_flg_array[1] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bXX
    if project.procPramArray[CI].adt_img_stat_flg_array[2] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bYY
    if project.procPramArray[CI].adt_img_stat_flg_array[3] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bZZ
    if project.procPramArray[CI].adt_img_stat_flg_array[4] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bXY
    if project.procPramArray[CI].adt_img_stat_flg_array[5] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bXZ
    if project.procPramArray[CI].adt_img_stat_flg_array[6] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bYZ
    if project.procPramArray[CI].adt_img_stat_flg_array[7] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bAD
    if project.procPramArray[CI].adt_img_stat_flg_array[8] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bFA
    if project.procPramArray[CI].adt_img_stat_flg_array[9] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bE1
    if project.procPramArray[CI].adt_img_stat_flg_array[10] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bE2
    if project.procPramArray[CI].adt_img_stat_flg_array[11] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bE3


    widget_control, ADT_ROI_WINDOW_BASE, /realize

    xmanager, 'adt_image_statistics_event',ADT_ROI_WINDOW_BASE,/no_block, GROUP_LEADER=ADT_WINDOW_BASE

end


; Subroutine name: adt_save_tif
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

    ;Edited by BT 2008/02/20 to use
    ;display routines to generate image data so that
    ;rotations and flips in main mas gui are honored

pro adt_save_tif, file=file, data_only=data_only, data_rcpt=data_rcpt; 'file' not implemented
    HEAP_GC
    COMPILE_OPT IDL2
    COMMON scan_data
    CI = project.ci
    sdim = project.imndArray[CI].sdim
    do_adt_regress_mono

    tiff_directory = DIALOG_PICKFILE( Title = 'Select the directory and filename for the Tiff image', $
                                      path=project.current_path, $
                                      /WRITE)

    if (tiff_directory eq '') then return ; user didn't cancel the save window

    pos = stregex(tiff_directory, '\.[Tt][Ii][Ff]?$', length=len)
    if (pos gt 0) then begin
        tiff_directory = strmid(tiff_directory, 0, pos)
    endif

    orient_filename = tiff_directory + '_orientation'

    sdim_start = project.procPramArray[CI].sdim_start
    pdim_start = project.procPramArray[CI].pdim_start
    fdim_start = project.procPramArray[CI].fdim_start

    Freq_interp  = project.procPramArray[ci].Freq_interp
    Phase_interp = project.procPramArray[ci].Phase_interp
    Slice_interp = project.procPramArray[ci].Slice_interp

    ffov = project.imndArray[ci].f_fov
    pfov = project.imndArray[ci].p_fov
    sfov = project.imndArray[ci].s_fov

    slice_axis = project.procpramArray[ci].slice_axis
    ADT_Start = project.procPramArray[CI].ADT_Start

    ;;  just save the color ball by itself
    if (ADT_start eq 9) then begin
        image = circle(200,200,0)
        pimage = ptr_new(bytarr(3,200,200))
        (*pimage)[0,*,*] = image[*,*,0]
        (*pimage)[1,*,*] = image[*,*,1]
        (*pimage)[2,*,*] = image[*,*,2]
        mas_rotate_flip, pimage
        filename = tiff_directory+'.tif'

        WRITE_TIFF, filename, *pimage
        ptr_free, pimage
        return
    endif

    adt_size = size(*project.dataArray[ci].adt)

    IF slice_axis EQ 0 THEN BEGIN
        slice_number = sdim_start
        num_slices = adt_size[3]
        x_fov = ffov*Freq_interp
        y_fov = pfov*Phase_interp

    ENDIF ELSE IF slice_axis EQ 1 THEN BEGIN
        slice_number = pdim_start
        num_slices = adt_size[2]
        x_fov = ffov*Freq_interp
        y_fov = sfov*Slice_interp

    ENDIF ELSE IF slice_axis EQ 2 THEN BEGIN
        slice_number = fdim_start
        num_slices = adt_size[1]
        x_fov = pfov*Phase_interp
        y_fov = sfov*Slice_interp
    ENDIF

    ;; check to see what type of display the users wants either a sheet of multiple images
    ;; or a display of a single image
    if project.procPramArray[CI].adt_display_multi eq 1 then begin

        ;; user wants all DTI parameters on one sheet
        ;; get the dti sheet and orient it
        img_tmp = ptr_new(adt_make_sheet(slice_number))
        adt_rotate_flip_display, img_tmp
        adt_zoom_display, img_tmp

        ;; the image array must be rotated for the tiff writer
        *img_tmp = temporary(reverse(*img_tmp, 2))

        img_sz = size(*img_tmp)
        img_data = bytarr(img_sz[3], img_sz[1], img_sz[2])

        file_name = $
          tiff_directory + "_" + $
          strcompress(string(slice_number+1), /REMOVE_ALL) + ".tif"

        ;; WRITE_TIFF wants [CHANNEL, X_PIXELS, Y_PIXELS]
        img_data[0,*,*] = (*img_tmp)[*,*,0]
        img_data[1,*,*] = (*img_tmp)[*,*,1]
        img_data[2,*,*] = (*img_tmp)[*,*,2]
        ptr_free, img_tmp

        ;; compute the tiff resolution from the fov
        x_res = (img_sz[1]/3 * 2.54) / x_fov
        y_res = (img_sz[2]/4 * 2.54) / y_fov

        WRITE_TIFF, file_name, img_data, XRESOL=x_res, YRESOL=y_res

    endif else begin

        ;; if the user wants to display all the slices for a scan
        ;; this button connects the multi button from the main display
        ;; so the user can get all the slices for the image they want.
        if project.procPramArray[CI].single_Multi_flag eq 0 then begin

            sl_start = slice_number
            sl_end = slice_number
            orient_filename += strcompress('_'+string(slice_number+1)+'.tif', /remove_all)

        endif else begin

            sl_start = 0
            sl_end = num_slices-1
            orient_filename += strcompress('_1-'+string(num_slices)+'.tif', /remove_all)

        endelse

        for slice = sl_start, sl_end do begin

            img_tmp = ptr_new(adt_make_image(slice))

            if project.procPramArray[CI].adt_display_type ne 4 then begin
                mas_zoom, img_tmp
                mas_rotate_flip, img_tmp
            endif

            file_name = $
              tiff_directory + "_" + $
              strcompress(string(slice+1), /REMOVE_ALL) + ".tif"

            img_sz = size(*img_tmp)

            if (img_sz[0] eq 2) then begin
                ;; image is greyscale
                img_data = bytarr(img_sz[1], img_sz[2])
                img_data = reverse(*img_tmp,2)
            endif else begin
                ;; image is rgb
                img_data = bytarr(img_sz[3], img_sz[1], img_sz[2])
                img_data[0,*,*] = reverse((*img_tmp)[*,*,0],2)
                img_data[1,*,*] = reverse((*img_tmp)[*,*,1],2)
                img_data[2,*,*] = reverse((*img_tmp)[*,*,2],2)
            endelse

            x_res = (img_sz[1] * 2.54) / x_fov
            y_res = (img_sz[2] * 2.54) / y_fov

            WRITE_TIFF, file_name, img_data, XRESOL=x_res, YRESOL=y_res

            if (ptr_valid(img_tmp)) then ptr_free, img_tmp
        endfor

        ;; Write the color ball, properly oriented
;        slash = get_dir_slash()
;        slashpos = STRPOS(tiff_directory, slash[1], /REVERSE_SEARCH)
;        split_directory = STRMID(tiff_directory,0, slashpos+1)
;        rotate_Direction = project.procPramArray[project.ci].rotate_Direction
;        case rotate_Direction of
;
;            0: begin
;                orient_filename = split_directory+'OrientSphere_0deg.tif'
;            end
;
;            1: begin
;                orient_filename = split_directory+'OrientSphere_90deg.tif'
;            end
;
;            2: begin
;                orient_filename = split_directory+'OrientSphere_180deg.tif'
;            end
;
;            3:begin
;                orient_filename = split_directory+'OrientSphere_270deg.tif'
;            end
;        end

;        if (FILE_TEST(orient_filename) eq 0) then begin

        if (project.procPramArray[CI].ADT_Start eq 8) then begin
            orient = ptr_new(circle(adt_size[1], adt_size[2], slice_axis))
            mas_rotate_flip, orient
            img_data = bytarr(3, img_sz[1], img_sz[2])
            img_data[0,*,*] = reverse((*orient)[*,*,0],2)
            img_data[1,*,*] = reverse((*orient)[*,*,1],2)
            img_data[2,*,*] = reverse((*orient)[*,*,2],2)
            ptr_free, orient

            WRITE_TIFF, orient_filename, img_data, XRESOL=x_res, YRESOL=y_res
        end

    endelse

    ;project.procPramArray[CI].adt_display_multi = adt_display_multi_old

end

; Subroutine name: adt_save_flt
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

pro adt_save_flt
    HEAP_GC
    COMPILE_OPT IDL2
    COMMON scan_data
    CI = project.ci

    do_adt_regress_mono

    ;save the data to flt files

   answer = 'Yes'
   if project.procpramArray[project.CI].adt_proccess_flag eq 0 then $
     answer = dialog_message('The whole ADT data set was not processed. Do you wish to save .flt files?', /QUESTION $
                                 ,TITLE='ADT .dwi Save')
   if answer eq 'Yes' then begin
     flt_directory = DIALOG_PICKFILE( Title = 'Set directory and choose file name for *.dwi files' $
             ,path=project.current_path)

     if flt_directory ne '' then begin


      update_status_bar,'saving flt 1/7'

      pImage = ptr_new(*project.dataArray[project.CI].adt)
      mas_rotate_flip, pImage

      sz_alp = size( *project.dataArray[project.CI].adt)
      OpenW, lun0, flt_directory+'.adt', /SWAP_IF_LITTLE_ENDIAN, /get_lun
      WriteU, lun0, sz_alp, *pImage
      close, lun0
      Free_LUN, lun0


;          ; Write the error result
;          update_status_bar,'saving flt 2/7'
;          sz = size(sigma)
;          OpenW, lun0, flt_directory+'sigma.flt', /SWAP_IF_LITTLE_ENDIAN, /get_lun
;          WriteU, lun0, sz, sigma
;        close, lun0
;          Free_LUN, lun0


      ;Calculate other "Basser" parameters
      ; (see Pierpaoli and Basser, MRM 36, 893-906)
      update_status_bar,'saving flt 3/7'

      Avg_Dif = ptr_new(*project.dataArray[project.CI].Avg_Dif)
      mas_rotate_flip, Avg_Dif
      sz_av_D = size(*Avg_Dif)

      OpenW, lun0, flt_directory+'.avdif', /SWAP_IF_LITTLE_ENDIAN, /get_lun
      WriteU, lun0, sz_av_D, *Avg_Dif
      close, lun0
      Free_LUN, lun0

      eign_val = ptr_new(*project.dataArray[project.CI].eign_val)
      mas_rotate_flip, eign_val
      sz_evals = size(*eign_val)

      update_status_bar,'saving flt 4/7'
      OpenW, lun0, flt_directory+'.eval', /SWAP_IF_LITTLE_ENDIAN, /get_lun
      WriteU, lun0, sz_evals, *eign_val
      close, lun0
      Free_LUN, lun0

      eign_Vec = ptr_new(*project.dataArray[project.CI].eign_Vec)
      mas_rotate_flip, eign_Vec
      sz_evecs = size(eign_Vec)

      update_status_bar,'saving flt 5/7'
      OpenW, lun0, flt_directory+'.evec', /SWAP_IF_LITTLE_ENDIAN, /get_lun
      WriteU, lun0, sz_evecs, *eign_Vec
      close, lun0
      Free_LUN, lun0


      frac_Ani = ptr_new(*project.dataArray[project.CI].frac_Ani)
      mas_rotate_flip, frac_Ani

      update_status_bar,'saving flt 6/7'
      OpenW, lun0, flt_directory+'.fa', /SWAP_IF_LITTLE_ENDIAN, /get_lun
      WriteU, lun0, sz_av_D, *frac_Ani
      close, lun0
      Free_LUN, lun0



      s0 = ptr_new((*project.dataArray[project.CI].adt)[*,*,*,0])
      mas_rotate_flip, s0

      update_status_bar,'saving flt 7/7'
      OpenW, lun0, flt_directory+'.s0', /SWAP_IF_LITTLE_ENDIAN, /get_lun
      WriteU, lun0, sz_av_D, *s0
      close, lun0
      Free_Lun, lun0


    write_geo_info_file, flt_directory, sz_alp

    ;flt_directory = DIALOG_PICKFILE( Title = 'Choose directory for *.flt files' $
             ;,path=flt_directory ,/directory )

    ;if flt_directory ne '' then begin

          OpenW, lun0, flt_directory+'Fract_Anisotropy.flt', /SWAP_IF_LITTLE_ENDIAN, /get_lun
          WriteU, lun0, sz_av_D, *frac_Ani
          close, lun0
          Free_Lun, lun0

          OpenW, lun0, flt_directory+'Eigenvectors.flt', /SWAP_IF_LITTLE_ENDIAN, /get_lun
          WriteU, lun0, size(*project.dataArray[project.CI].eign_Vec) $
                   ,*project.dataArray[project.CI].eign_Vec
          close, lun0
          Free_Lun, lun0

          OpenW, lun0, flt_directory+'Aver_D.flt', /SWAP_IF_LITTLE_ENDIAN, /get_lun
          WriteU, lun0, sz_av_D, *Avg_Dif
          close, lun0
          Free_Lun, lun0

          OpenW, lun0, flt_directory+'S0.flt', /SWAP_IF_LITTLE_ENDIAN, /get_lun
          WriteU, lun0, sz_av_D, *s0
          close, lun0
          Free_Lun, lun0


      ;end



      update_status_bar,''
     end
   end
end

;;; Function is never called by any routine in mas
; FUNCTION linear_regress, X, Y, n_obs, n_coeffs, CONST=const, SIGMA=sigma

; 	ON_ERROR,2

; 	xmean  = TOTAL(X,2)/n_obs
; 	ymean  = TOTAL(Y)/n_obs

; 	xx = X - REBIN(xmean,n_coeffs,n_obs)

; 	sigmax = SQRT( TOTAL(xx^2,2) / (n_obs-1) )
; 	sigmay = SQRT(TOTAL((Y - ymean)^2)/(n_obs-1))

; 	array  = INVERT((xx # TRANSPOSE(xx)) / ((n_obs-1)*sigmax # sigmax))

; 	correlation = (TEMPORARY(xx) # (Y - ymean)) / (sigmax * sigmay * (n_obs-1))

; 	coeff = (correlation # array)*(sigmay/sigmax)
; 	const = ymean - TOTAL(coeff*xmean)

; 	yfit  = (coeff # X) + const
; 	freen = n_coeffs-n_obs-1 > 1

; 	chisq  = TOTAL((Y - yfit)^2)
; 	varnce = chisq/freen

; 	sigma = SQRT(array[LINDGEN(n_coeffs)*(n_coeffs+1)]*varnce/((n_obs-1)*sigmax^2))

; 	RETURN, coeff
; END


; Subroutine name: do_adt_regress_mono
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This is the subroutine that loads the data from state1 and

pro do_adt_regress_mono

    common scan_data

    if (project.procpramarray[project.ci].adt_regress_method eq 1) then begin
        do_adt_regress_standard
    endif else begin
        do_adt_regress_mow
    endelse

end

pro do_adt_regress_standard_singlevoxel, ijk, eval=eval, evec=evec, fa=fa, ad=ad, fit=fit, bval_threshold=bval_threshold

    common scan_data, project
    
    ci = project.ci
    sig = project.dataarray[ci].state1
    
    bvals = *project.imndarray[ci].bval_array
    
    if (n_elements(bval_threshold) ne 0) then begin
         use_bvals = where(bvals lt bval_threshold, nb)
         if (nb eq 0) then begin
            message, 'No b values eligible using the present threshold.'
         endif
    endif else begin
        use_bvals = lindgen(n_elements(bvals))
    endelse
    
    bmatrix = (*project.imndArray[CI].b_matrix)[*,use_bvals]
    bvals = bvals[use_bvals]
    
    if (n_elements(ijk) eq 0) then begin
        i = project.procpramarray[ci].fdim_start
        j = project.procpramarray[ci].pdim_start
        k = project.procpramarray[ci].sdim_start
    endif else begin
        i = ijk[0]
        j = ijk[1]
        k = ijk[2]
    endelse
    
    ;number of files
    n     = n_elements(bvals)
    
    ;rank of the tensor
    n_bvals = n_elements(bvals)
    
    ; redefine parameters here since easier to follow regression
    X      = temporary(-bmatrix)
    
    n_obs    = n_bvals
    n_coeffs = 6
    
    xmean  = TOTAL(X,2)/n_obs
    xx     = X - REBIN(xmean,n_coeffs,n_obs)
    sigmax = SQRT( TOTAL(xx^2,2) / (n_obs-1) )
    varx   = sigmax^2
    array  = INVERT(MATRIX_MULTIPLY(xx, xx, /BTRANSPOSE) / ((n_obs-1)*sigmax # sigmax))
    
    freen  = n_coeffs-n_obs-1 > 1
    
    y_iter  = REFORM(ALOG( (*sig)[i, j, k, 0:n_bvals-1] ))
    y_iter = y_iter[use_bvals]
    
    void = where( FINITE(y_iter) ne 1, n_nan )
    if n_nan gt 0 then return
    
    ymean  = TOTAL(y_iter)/n_obs
    yy_iter = y_iter - ymean
    
    sigmay = SQRT(TOTAL(yy_iter^2)/(n_obs-1))
    correlation = (xx # yy_iter) / (sigmax * sigmay * (n_obs-1))
    
    coeff = (correlation # array)*(sigmay/sigmax)
    const = ymean - TOTAL(coeff*xmean)
    if (total(finite(coeff, /nan)) ne 0) then return
    
    ADT = reform([ [exp(const)], [coeff] ])
    
    A_matrix = [ [ ADT[1], ADT[4], ADT[5] ], $
                 [ ADt[4], ADT[2], ADT[6] ], $
                 [ ADT[5], ADT[6], ADT[3] ] ]
    
    ev = 1 ; initialize eigenvectors (see IDL function EIGENQL)
    eigenval = EIGENQL(A_matrix, EIGENVECTORS = ev, /overwrite)
    eval = eigenval
    evec = ev
    
    avg_eign = TOTAL(eigenval) / 3.
    avg_diff = abs(avg_eign) ;; the above is equal to the previous assignment
    frac_anis = ( (sqrt(3./2.) * sqrt(total((eigenval - avg_eign)^2)/total(eigenval^2))) ) < 1.

    fa = frac_anis
    ad = avg_diff
    fit = adt
    
end

; Subroutine name: do_adt_regress_standard
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This is the subroutine that loads the data from state1 and computes
; the diffusion tensor using standard linear regression

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.
    ;Edited by CD 2007/06/28
    ;normalize data prior to regression in all cases
    ;Edited by MK 2014/11/11
    ;Added feature which allows one to choose the diffusion directions to fit
pro do_adt_regress_standard, save_residuals=save_residuals

    COMPILE_OPT IDL2
    COMMON common_widgets
    COMMON scan_data
    CI = project.ci

    loadct, 0, /silent
   
    if project.procpramArray[CI].adt_proccess_flag eq 0 then begin

        ;***********************        
        if (ptr_valid(project.dataArray[ci].adt)) then ptr_free, project.dataArray[ci].adt
        if (ptr_valid(project.dataArray[ci].eign_val)) then ptr_free, project.dataArray[ci].eign_val
        if (ptr_valid(project.dataArray[ci].eign_Vec)) then ptr_free, project.dataArray[ci].eign_Vec
        if (ptr_valid(project.dataArray[ci].frac_Ani)) then ptr_free, project.dataArray[ci].frac_Ani
        if (ptr_valid(project.dataArray[ci].Avg_Dif)) then ptr_free, project.dataArray[ci].Avg_Dif
        if (ptr_valid(project.dataarray[ci].fit_residuals)) then ptr_free, project.dataarray[ci].fit_residuals
        if (ptr_valid(project.dataArray[ci].fibers_adt)) then ptr_free, project.dataArray[ci].fibers_adt
        if (ptr_valid(project.dataArray[ci].fibers_eign_val)) then ptr_free, project.dataArray[ci].fibers_eign_val
        if (ptr_valid(project.dataArray[ci].fibers_eign_Vec)) then ptr_free, project.dataArray[ci].fibers_eign_Vec
        if (ptr_valid(project.dataArray[ci].fibers_frac_Ani)) then ptr_free, project.dataArray[ci].fibers_frac_Ani
        if (ptr_valid(project.dataArray[ci].fibers_Avg_Dif)) then ptr_free, project.dataArray[ci].fibers_Avg_Dif

        heap_gc

        if (not ptr_valid(project.imndarray[ci].bval_array)) then begin
            junk = dialog_message(['No B values exist for this scan.', $
                                   'ADT cannot be processed with a valid', $
                                   'B matrix and B values'], /error, /center)
            return
        endif

    	;check to see what regression type the user wants to use then call
    	; the appropriate one if not -1 which is the default.
        if project.procPramArray[CI].regression_type ne -1 then begin
            mas_do_smooth
            return
        endif
                
           mas_load_state_1
           data = *project.dataarray[ci].state1
            
           progressbar = Obj_New('progressbar', Color='red', Text='ADT Fitting (Standard Method)', /fast_loop)
           progressbar -> Start
        
           sz_sig = size(data)

           sdim_start = project.procpramArray[ci].sdim_start
           adim_start = project.procpramArray[ci].adim_start

            fdim = sz_sig[1]
            pdim = sz_sig[2]
            sdim = sz_sig[3]
            adim = sz_sig[4]

            ;now we need to generate a mask for our data set.
            ;since we have multiple slices we have to generate it from the beginning
            ;and loop through then freq then phase then slice. that way we dont
            ;have to recreate them unnecessarily.

            mask_background_flag =  project.procPramArray[CI].remove_background_active

            ;; note that if we need to mask the ADT fit with either an
            ;; ROI or a background mask, we will have to duplicate the
            ;; state1 data, at least _for now_. It doesn't seem a
            ;; common practice to use ROI or background mask, though.
            IF (mask_background_flag or project.roi.mask) then begin
                sig = ptr_new(data)
                IF mask_background_flag EQ 1 THEN BEGIN
                    print, "Masking background..."
                    mask = bytarr(sz_sig[1], sz_sig[2])
                    key_adim = project.procPramArray[CI].key_adim

                    FOR ss=0, sz_sig[3]-1 DO BEGIN
                        p_image_mask = PTR_NEW(REFORM(data[*,*,ss,key_adim]))
                        mask[*,*]= *(generate_background_mask( p_image_mask))
                        FOR AA =0, sz_sig[4]-1 DO BEGIN
                            (*sig)[*,*,ss,aa] *= mask
                        END
                    ENDFOR

                ENDIF

                ;;use the roi mask
                IF project.roi.mask GE 1 THEN begin
                    print, "Masking with ROI..."
                    FOR ii=0, sz_sig[3]-1 DO BEGIN
                        FOR jj=0, sz_sig[4]-1 DO BEGIN
                            p_image = ptr_new((*sig)[*,*,ii,jj])
                            mas_roi_mask, p_image
                            (*sig)[*,*,ii] = *p_image
                        ENDFOR
                    ENDFOR

                ENDIF

            ENDIF ELSE BEGIN
                ;; Since we're not going to modify the state1 data, we
                ;; can just use the state1 pointer explicitly to save memory.
                sig = ptr_new(data)
            ENDELSE

            bvals = *project.imndArray[CI].b_matrix           
            ;number of files
            n     = adim
            ;rank of the tensor
            n_bvals = project.imndArray[CI].n_bvals
            
            if ptr_valid(project.procpramArray[CI].array_select) then begin
              sig = ptr_new((*sig)[*,*,*,*project.procpramArray[CI].array_select])
              bvals = (bvals)[*,*project.procpramArray[CI].array_select]
              n = (size(*project.procpramArray[CI].array_select,/Dimension))[0]
              n_bvals = n
            endif  
                        
            k     = sdim_start
            user_thr = project.procPramArray[CI].adtThreshold
            thr = max(*sig) * user_thr

            ; Make arrays for the alp, ADT and sigma matrices as well as 'Basser' result
            ln_sig_1d = Make_Array(1, adim)
            fit_result = Make_Array(1, 7)
            ADT = Make_Array(fdim, pdim, sdim, 7, value = 0.0)
            sigma = Make_Array(fdim, pdim, sdim, 6, value = 0.0)
            A_matrix = Make_array(3, 3)
            evals = Make_array(fdim, pdim, sdim, 3)
            evecs = Make_array(fdim, pdim, sdim, 3, 3)
            avg_diff = Make_array(fdim, pdim, sdim)
            frac_anis = Make_Array(fdim,pdim,sdim)

            ; the following code was taken from the regress.pro function that comes
            ; as a part of IDL.  This code involves the calculation of the mean, stddev
            ; and variance of both the X and Y terms in the linear regression Y = mX + b.
            ; since X is equal to the -bmatrix (variable = -bvals) in all of the regressions
            ; in the loop that repeats for each fdim,pdim,sdim the calcuations of the mean,'
            ; stddev, and variance of the X term can be moved outside the loop and only
            ; calculate once, this results in 20% more efficient code.  In addition, the
            ; overhead in the regress.pro function from IDL that allows for a weighted
            ; regression can be removed, this results in 10% more efficient code
            ; for details on the actual algorithm for solving the regression, see regress.pro
            ; in the IDL library
            ; redefine parameters here since easier to follow regression
            X      = temporary(-bvals)

            ;; commented by BT. This gets assigned on a slice-by-slice
            ;; basis in the inner for loop below
;            Y      = temporary(alog(sig))

            n_obs    = n_bvals
            n_coeffs = 6

            xmean  = TOTAL(X,2)/n_obs
            xx     = X - REBIN(xmean,n_coeffs,n_obs)
            sigmax = SQRT( TOTAL(xx^2,2) / (n_obs-1) )
            varx   = sigmax^2
            array  = INVERT(MATRIX_MULTIPLY(xx, xx, /BTRANSPOSE) / ((n_obs-1)*sigmax # sigmax))

            freen  = n_coeffs-n_obs-1 > 1
            ln_thr = alog(thr)

            start_time = systime(1)

            ;;refit = make_array(fdim,pdim,sdim,adim, value=0.0)

            FOR i = 0, fdim-1 DO BEGIN
                FOR j = 0, pdim-1 DO BEGIN
                    FOR k = 0, sdim-1 DO BEGIN

                        if (*sig)[i, j, k, adim_start] le thr then continue

                        y_iter  = REFORM(ALOG( (*sig)[i, j, k, 0:n_bvals-1] ))

                        void = where( FINITE(y_iter) ne 1, n_nan )
                        if n_nan gt 0 then  continue

                        ymean  = TOTAL(y_iter)/n_obs
                        yy_iter = y_iter - ymean

                        sigmay = SQRT(TOTAL(yy_iter^2)/(n_obs-1))
                        correlation = (xx # yy_iter) / (sigmax * sigmay * (n_obs-1))

                        coeff = (correlation # array)*(sigmay/sigmax)
                        const = ymean - TOTAL(coeff*xmean)
                        if (total(finite(coeff, /nan)) ne 0) then continue
                        ADT[i, j, k, 0] = exp(const)
                        ADT[i, j, k, 1:6] = coeff

                        A_matrix = [[ADT[i,j,k,1], ADT[i,j,k,4], ADT[i,j,k,5]], $
                                    [ADt[i,j,k,4], ADT[i,j,k,2], ADT[i,j,k,6]], $
                                    [ADT[i,j,k,5], ADT[i,j,k,6], ADT[i,j,k,3]]]

                        ev = 1 ; initialize eigenvectors (see IDL function EIGENQL)
                        eigenval = EIGENQL(A_matrix, EIGENVECTORS = ev, /overwrite)
                        evals[i, j, k, *] = eigenval
                        evecs[i, j, k, *, *] = ev

                        ;;; Avg. diffusion may be calc'd from the avg eigenvalues. invariant
                        ;;; scalar measures
                        ;; see Basser and Pierpaoli, "Microstructural features measured using diffusion
                        ;; tensor imaging", J. Magn. Reson. B 111, 209-219 (1996)
                        ;; the change below is equivalent mathmatically to [sqrt(1.-I2/(I1^2-2.*I2)) < 1.]

                        avg_eign = TOTAL(eigenval) / 3.
                        avg_diff[i,j,k] = abs(avg_eign) ;; the above is equal to the previous assignment
                        frac_anis[i,j,k] = ( (sqrt(3./2.) * sqrt(total((eigenval - avg_eign)^2)/total(eigenval^2))) ) < 1.

                    ENDFOR

                ENDFOR

                progressBar -> Update, (float(i)/float(fdim-1))*100.0

                if (progressBar->checkcancel()) then begin
                    progressbar -> Destroy
                    project.procpramArray[project.CI].adt_proccess_flag = 0
                    heap_gc
                    return
                endif

            ENDFOR

            project.procpramArray[project.CI].adt_proccess_flag = 1
            progressbar -> Destroy            

;****************************
; ***********************************************************
; =============== PART THAT SAVES THE PROJECT ADT DATA ======
; ***********************************************************

        update_status_bar, 'Saving structures to memory'

        project.dataArray[project.CI].adt = ptr_new(ADT,/no_copy)
        project.dataArray[project.CI].eign_val = ptr_new(evals,/no_copy)
        project.dataArray[project.CI].eign_Vec = ptr_new(evecs,/no_copy)
        project.dataArray[project.CI].frac_Ani = ptr_new(frac_anis,/no_copy)
        project.dataArray[project.CI].Avg_Dif =  ptr_new(avg_diff,/no_copy)
        
        update_status_bar, ''

    endif

    ; HS - 20061031
    ; This is to force the gui to recalculate the sensitivity of the Fiber Track Mapping
    ; button in the menu.
    mas_redraw_gui

   ;this is a way for me to clean up my errors b/c i didn't write much of this code. i just
   ;scrubbed it.
    heap_gc

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  Experimenal, not used at the moment.  
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro do_adt_regress_mow

    COMPILE_OPT IDL2
    COMMON common_widgets
    COMMON scan_data
    CI = project.ci

    if project.procpramArray[CI].state_1 eq 0 then begin
        project.procpramArray[CI].adt_proccess_flag = 0
    endif

    if project.procpramArray[CI].adt_proccess_flag eq 0 then begin

        mas_load_state_1

        if (ptr_valid(project.dataArray[ci].adt)) then ptr_free, project.dataArray[ci].adt
        if (ptr_valid(project.dataArray[ci].eign_val)) then ptr_free, project.dataArray[ci].eign_val
        if (ptr_valid(project.dataArray[ci].eign_Vec)) then ptr_free, project.dataArray[ci].eign_Vec
        if (ptr_valid(project.dataArray[ci].frac_Ani)) then ptr_free, project.dataArray[ci].frac_Ani
        if (ptr_valid(project.dataArray[ci].Avg_Dif)) then ptr_free, project.dataArray[ci].Avg_Dif

        if (ptr_valid(project.dataArray[ci].fibers_adt)) then ptr_free, project.dataArray[ci].fibers_adt
        if (ptr_valid(project.dataArray[ci].fibers_eign_val)) then ptr_free, project.dataArray[ci].fibers_eign_val
        if (ptr_valid(project.dataArray[ci].fibers_eign_Vec)) then ptr_free, project.dataArray[ci].fibers_eign_Vec
        if (ptr_valid(project.dataArray[ci].fibers_frac_Ani)) then ptr_free, project.dataArray[ci].fibers_frac_Ani
        if (ptr_valid(project.dataArray[ci].fibers_Avg_Dif)) then ptr_free, project.dataArray[ci].fibers_Avg_Dif

        heap_gc

        if (not ptr_valid(project.imndarray[ci].bval_array)) then begin
            junk = dialog_message(['No B values exist for this scan.', $
                                   'ADT cannot be processed with a valid', $
                                   'B matrix and B values'], /error, /center)
            return
        endif

        progressbar = Obj_New('progressbar', Color='red', Text='ADT Fitting (MOW Method)', /fast_loop)
        progressbar -> Start

        sz_sig = size(*project.dataArray[CI].state1)

        adim_start = project.procpramArray[ci].adim_start

        fdim = sz_sig[1]
        pdim = sz_sig[2]
        sdim = sz_sig[3]
        adim = sz_sig[4]

        sig = project.dataArray[CI].state1
        bvals = *project.imndArray[CI].b_matrix
        ;;number of files
        n     = adim

        ;;rank of the tensor
        n_bvals = project.imndArray[CI].n_bvals

        user_thr = project.procPramArray[CI].adtThreshold
        max_sig = max(*sig)
        thr = max_sig * user_thr

        ;; Make arrays for the alp, ADT and sigma matrices as well as 'Basser' result
        ln_sig_1d = Make_Array(1, adim)
        ADT = Make_Array(fdim, pdim, sdim, 7, value = 0.0)
        A_matrix = Make_array(3, 3)
        evals = Make_array(fdim, pdim, sdim, 3)
        evecs = Make_array(fdim, pdim, sdim, 3, 3)
        avg_diff = Make_array(fdim, pdim, sdim)
        frac_anis = Make_Array(fdim,pdim,sdim)
        n_obs    = n_bvals
        n_coeffs = 6

        start_time = systime(1)

        a_mat = fltarr(7,n_bvals)
        a_mat[1:6,*] = bvals * 1e-0
        id = replicate(-1.0, n_obs)

        FOR i = 0, fdim-1 DO BEGIN
            FOR j = 0, pdim-1 DO BEGIN
                FOR k = 0, sdim-1 DO BEGIN

                    if (*sig)[i, j, k, adim_start] le thr then continue

;;                    y_iter  = REFORM(ALOG( (*sig)[i, j, k, 0:n_bvals-1] ))
                    y_iter  = REFORM((*sig)[i, j, k, 0:n_bvals-1] )
                    y_iter  = 1.0/sqrt(y_iter)

                    void = where( FINITE(y_iter) ne 1, n_nan )
                    if n_nan gt 0 then  continue

                    a_mat[0,*] = y_iter

                    temp = invert(matrix_multiply(a_mat, a_mat, /btranspose)) # a_mat # id

                    ADT[i,j,k,0] = temp[0]^2
                    ADT[i,j,k,1:6] = 2*temp[1:6]

                    A_matrix = [[ADT[i,j,k,1], ADT[i,j,k,4], ADT[i,j,k,5]], $
                                [ADT[i,j,k,4], ADT[i,j,k,2], ADT[i,j,k,6]], $
                                [ADT[i,j,k,5], ADT[i,j,k,6], ADT[i,j,k,3]]]

                    ev = 1 ; initialize eigenvectors (see IDL function EIGENQL)
                    eigenval = EIGENQL(A_matrix, EIGENVECTORS = ev, /overwrite)
                    evals[i, j, k, *] = eigenval
                    evecs[i, j, k, *, *] = ev
                    avg_eign = TOTAL(eigenval) / 3.
                    avg_diff[i,j,k] = abs(avg_eign) ;; the above is equal to the previous assignment
                    frac_anis[i,j,k] = ( (sqrt(3./2.) * sqrt(total((eigenval - avg_eign)^2)/total(eigenval^2))) ) < 1.

                    progressBar -> Update, (float(i)/float(fdim-1))*100.0

                endfor
            endfor

            if (progressBar->checkcancel()) then begin
                progressbar -> Destroy
                project.procpramArray[project.CI].adt_proccess_flag = 0
                heap_gc
                return
            endif

        endfor
        project.procpramArray[project.CI].adt_proccess_flag = 1
        progressbar -> Destroy
        print, "Time: "+string(systime(1) - start_time)
        update_status_bar, 'Saving structures to memory'

; These will be the "backup" data which will not be edited by any routine in MAS.
        project.dataArray[project.CI].adt = ptr_new(ADT)
        project.dataArray[project.CI].eign_val = ptr_new(evals)
        project.dataArray[project.CI].eign_Vec = ptr_new(evecs)
        project.dataArray[project.CI].frac_Ani = ptr_new(frac_anis)
        project.dataArray[project.CI].Avg_Dif = ptr_new(avg_diff)

; These are the ones that will be edited in MAS and Fibers
        project.dataArray[project.CI].fibers_adt = ptr_new(ADT,/no_copy)
        project.dataArray[project.CI].fibers_eign_val = ptr_new(evals,/no_copy)
        project.dataArray[project.CI].fibers_eign_Vec = ptr_new(evecs,/no_copy)
        project.dataArray[project.CI].fibers_frac_Ani = ptr_new(frac_anis,/no_copy)
        project.dataArray[project.CI].fibers_Avg_Dif = ptr_new(avg_diff,/no_copy)

        update_status_bar, ''

    endif

    ; HS - 20061031
    ; This is to force the gui to recalculate the sensitivity of the Fiber Track Mapping
    ; button in the menu.
    mas_redraw_gui

   ;this is a way for me to clean up my errors b/c i didn't write much of this code. i just
   ;scrubbed it.
    heap_gc

end

; Subroutine name: adt_calc_lattice_anisotrophy
; Created by: CD 9/10/07
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; executes lattice anistrophy calculation defined
; in Pajevic and Pierpaoli MRM 42:526-540 (1999)

; Editing Information:

PRO adt_calc_lattice_anisotrophy, slice_number, LA_3D_FLAG

    HEAP_GC
    COMMON scan_data
    CI = project.ci

	slice_axis = project.procpramArray[CI].slice_axis

	do_adt_regress_mono

	  ADT   = *project.dataArray[CI].adt        ;(r by p by s by 7)
    EVALS = *project.dataArray[CI].eign_val   ;(r by p by s by 3)
    EVECS = *project.dataArray[CI].eign_Vec   ;(r by p by s by 3 by 3) (one ev per row)

	CASE slice_axis OF
		0: BEGIN
			xres = project.imndArray[CI].f_fov
			yres = project.imndArray[CI].p_fov
			zres = project.imndArray[CI].s_fov
		END

		1: BEGIN
			ADT   = TRANSPOSE(ADT,[0,2,1,3])
			EVALS = TRANSPOSE(EVALS,[0,2,1,3])
			EVECS = TRANSPOSE(EVECS,[0,2,1,3,4])

			xres = project.imndArray[CI].f_fov
			yres = project.imndArray[CI].s_fov
			zres = project.imndArray[CI].p_fov
		END

		2: BEGIN
			ADT   = TRANSPOSE(ADT,[1,2,0,3])
			EVALS = TRANSPOSE(EVALS,[1,2,0,3])
			EVECS = TRANSPOSE(EVECS,[1,2,0,3,4])

			xres = project.imndArray[CI].p_fov
			yres = project.imndArray[CI].s_fov
			zres = project.imndArray[CI].f_fov
		END
	ENDCASE

	sz_sig = SIZE(ADT)
	xdim = sz_sig[1]
	ydim = sz_sig[2]
	zdim = sz_sig[3]

	xres = xres / xdim
	yres = yres / ydim
	zres = zres / zdim

	; if slice_axis switched since last lattice aniso calculation, have to reprocess
	; if the same see if that slice has already been calculated
	;IF slice_axis NE project.procpramArray[CI].latt_slice_index THEN BEGIN
		project.procpramArray[CI].latt_slice_status = PTR_NEW(INTARR(zdim))
		(*project.procpramArray[CI].latt_slice_status)[slice_number] = 1

		project.procpramArray[CI].latt_slice_index = slice_axis

		LA = DBLARR(xdim,ydim,zdim) ; stored in fdim,pdim,sdim format regardless of slice selection
	;ENDIF ELSE BEGIN
	;	IF (*project.procpramArray[CI].latt_slice_status)[slice_number] EQ 1 THEN RETURN

	;	LA = *project.dataArray[CI].latt_Ani
	;ENDElSE

	slice_start = slice_number
	slice_end = slice_start

	num_rot_elements = 9

	egvc_rot_one = [0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8]
	egvl_rot_one = [0,0,0,1,1,1,2,2,2,0,0,0,1,1,1,2,2,2,0,0,0,1,1,1,2,2,2]

	egvc_rot_two = [0,1,2,3,4,5,6,7,8,3,4,5,6,7,8,0,1,2,6,7,8,0,1,2,3,4,5]
	egvl_rot_two = [0,0,0,1,1,1,2,2,2,1,1,1,2,2,2,0,0,0,2,2,2,0,0,0,1,1,1]

	idx_thr = [0,1,2]

	sz_rot = (SIZE(egvc_rot_one))[1]

	FOR k=slice_start, slice_end DO BEGIN

		IF LA_3D_FLAG EQ 1 AND slice_number NE 0 AND slice_number NE zdim-1 THEN BEGIN
			progressbar = Obj_New('progressbar', Color='red', Text='3D LA Calculation Slice '+STRTRIM(STRING(k+1),2),/NOCANCEL)
		ENDIF ELSE BEGIN
			progressbar = Obj_New('progressbar', Color='red', Text='2D LA Calculation Slice '+STRTRIM(STRING(k+1),2),/NOCANCEL)
		ENDELSE

    	progressbar -> Start

		lattice_n = FLTARR(3+2*(xdim-1),3+2*(ydim-1))

		; define num_neighbors
		num_3D_neighbors = 18
		lattice_3D  = FLTARR(xdim, ydim, num_3D_neighbors)
		lattice_wgt = FLTARR(xdim, ydim, num_3D_neighbors)

		FOR j=0, ydim-2 DO BEGIN
			FOR i=0, xdim-2 DO BEGIN

				; calculate d_dp for each voxel with itself
				egvc_slf_one = [(EVECS(i,j,k,*,*))[egvc_rot_one],(EVECS(i+1,j,k,*,*))[egvc_rot_one],(EVECS(i,j+1,k,*,*))[egvc_rot_one],(EVECS(i+1,j+1,k,*,*))[egvc_rot_one]]
				egvl_slf_one = [(EVALS(i,j,k,*))[egvl_rot_one],(EVALS(i+1,j,k,*))[egvl_rot_one],(EVALS(i,j+1,k,*))[egvl_rot_one],(EVALS(i+1,j+1,k,*))[egvl_rot_one]]

				egvc_slf_two = [(EVECS(i,j,k,*,*))[egvc_rot_two],(EVECS(i+1,j,k,*,*))[egvc_rot_two],(EVECS(i,j+1,k,*,*))[egvc_rot_two],(EVECS(i+1,j+1,k,*,*))[egvc_rot_two]]
				egvl_slf_two = [(EVALS(i,j,k,*))[egvl_rot_two],(EVALS(i+1,j,k,*))[egvl_rot_two],(EVALS(i,j+1,k,*))[egvl_rot_two],(EVALS(i+1,j+1,k,*))[egvl_rot_two]]

				d_dp_slf = egvc_slf_one * SQRT(egvl_slf_one) * egvc_slf_two * SQRT(egvl_slf_two)

				; neighbor i,j (middle), reference voxel
				d_dp_ref_slf = 0
				FOR idx=0, 8 DO d_dp_ref_slf = d_dp_ref_slf + TOTAL(d_dp_slf[(idx_thr+(idx*3))])^2

				; neighbor i+1,j (right)
				d_dp_r_slf = 0
				FOR idx=9, 17 DO d_dp_r_slf = d_dp_r_slf + TOTAL(d_dp_slf[(idx_thr+(idx*3))])^2

				; neighbor i,j+1 (top)
				d_dp_t_slf = 0
				FOR idx=18, 26 DO d_dp_t_slf = d_dp_t_slf + TOTAL(d_dp_slf[(idx_thr+(idx*3))])^2

				; neighbor i+1,j+1 (top right)
				d_dp_tr_slf = 0
				FOR idx=27, 35 DO d_dp_tr_slf = d_dp_tr_slf + TOTAL(d_dp_slf[(idx_thr+(idx*3))])^2

				; calculate d_dp for each voxel with reference voxel
				egvc_one = [(EVECS(i,j,k,*,*))[egvc_rot_one],(EVECS(i,j,k,*,*))[egvc_rot_one],(EVECS(i,j,k,*,*))[egvc_rot_one]]
				egvl_one = [(EVALS(i,j,k,*))[egvl_rot_one],(EVALS(i,j,k,*))[egvl_rot_one],(EVALS(i,j,k,*))[egvl_rot_one]]

 				egvc_two = [(EVECS(i+1,j,k,*,*))[egvc_rot_two],(EVECS(i,j+1,k,*,*))[egvc_rot_two],(EVECS(i+1,j+1,k,*,*))[egvc_rot_two]]
 				egvl_two = [(EVALS(i+1,j,k,*))[egvl_rot_two],(EVALS(i,j+1,k,*))[egvl_rot_two],(EVALS(i+1,j+1,k,*))[egvl_rot_two]]

				trace_ref = ADT[i,j,k,1] + ADT[i,j,k,2] + ADT[i,j,k,3]
				trace_r   = ADT[i+1,j,k,1] + ADT[i+1,j,k,2] + ADT[i+1,j,k,3]
				trace_t   = ADT[i,j+1,k,1] + ADT[i,j+1,k,2] + ADT[i,j+1,k,3]
				trace_tr  = ADT[i+1,j+1,k,1] + ADT[i+1,j+1,k,2] + ADT[i+1,j+1,k,3]

				d_dp_iter = egvc_one * SQRT(egvl_one) * egvc_two * SQRT(egvl_two)

				; neighbor i+1,j (right)
				d_dp_r = 0.0
				FOR idx=0, 8 DO d_dp_r = d_dp_r + TOTAL(d_dp_iter[(idx_thr+(idx*3))])^2
				d_dp_ital_r = d_dp_r - (1.0/3 * trace_ref * trace_r)

				; neighbor i, j+1 (top)
				d_dp_t = 0.0
				FOR idx=9, 17 DO d_dp_t = d_dp_t + TOTAL(d_dp_iter[(idx_thr+(idx*3))])^2
				d_dp_ital_t = d_dp_t - (1.0/3 * trace_ref * trace_t)

				; neighbor i+1, j+1 (top right)
				d_dp_tr = 0.0
				FOR idx=18, 26 DO d_dp_tr = d_dp_tr + TOTAL(d_dp_iter[(idx_thr+(idx*3))])^2
				d_dp_ital_tr = d_dp_tr - (1.0/3 * trace_ref * trace_tr)

				lattice_n[2*i+2,2*j+1] = SQRT((3*d_dp_ital_r)/(8*d_dp_r)) + (3.0/4)*(d_dp_ital_r/(SQRT(d_dp_ref_slf*d_dp_r_slf)))
				lattice_n[2*i+1,2*j+2] = SQRT((3*d_dp_ital_t)/(8*d_dp_t)) + (3.0/4)*(d_dp_ital_t/(SQRT(d_dp_ref_slf*d_dp_t_slf)))
				lattice_n[2*i+2,2*j+2] = SQRT((3*d_dp_ital_tr)/(8*d_dp_tr)) + (3.0/4)*(d_dp_ital_tr/(SQRT(d_dp_ref_slf*d_dp_tr_slf)))

				; if we want 3D Lattice Anistropy, calculate it, since only calculating one slice at time
				; to maintain reasonable timeframe, there is no overlap between calculations for the 18
				; neighbors in other slices

				IF LA_3D_FLAG EQ 1 AND $
					slice_number NE 0 AND slice_number NE zdim-1 AND $
					i NE 0 AND j NE 0 $
				THEN BEGIN

					; -----------------------------
					; BEGIN 3D NEIGHBOR DEFINITIONS
					;------------------------------

					egvc_rot = [[egvc_rot_one],[egvc_rot_two],[egvc_rot_one],[egvc_rot_two]]
					egvl_rot = [[egvl_rot_one],[egvl_rot_two],[egvl_rot_one],[egvl_rot_two]]

					k_slice  = [k+1, k+1, k-1, k-1]
					i_slice  = [i,i,i,i+1,i+1,i+1,i-1,i-1,i-1]
					j_slice  = [j,j+1,j-1,j,j+1,j-1,j,j+1,j-1]

					; ---------------------------
					; END 3D NEIGHBOR DEFINITIONS
					;----------------------------

					sz_k   = (SIZE(k_slice))[1]
					sz_ij  = (SIZE(i_slice))[1]

					egvc_slf = FLTARR(sz_rot*sz_ij,sz_k)
					egvl_slf = FLTARR(sz_rot*sz_ij,sz_k)

					egvc_vxl = FLTARR(sz_rot*sz_ij,sz_k)
					egvl_vxl = FLTARR(sz_rot*sz_ij,sz_k)

					trace_mat = FLTARR(num_3D_neighbors)

					FOR kidx = 0, sz_k-1 DO BEGIN
						egvc_iter = egvc_rot[*,kidx]
						egvl_iter = egvl_rot[*,kidx]
						k_iter    = k_slice[kidx]

						FOR ijidx = 0, sz_ij-1 DO BEGIN
							egvc_slf[ijidx*sz_rot:(ijidx+1)*sz_rot-1 ,kidx] = (EVECS[i_slice[ijidx],j_slice[ijidx],k_iter,*,*])[egvc_iter]
							egvl_slf[ijidx*sz_rot:(ijidx+1)*sz_rot-1 ,kidx] = (EVALS[i_slice[ijidx],j_slice[ijidx],k_iter,*])[egvl_iter]

							IF kidx MOD 2 EQ 0 THEN BEGIN
								egvc_vxl[ijidx*sz_rot:(ijidx+1)*sz_rot-1 ,kidx] = (EVECS[i,j,k,*,*])[egvc_iter]
								egvl_vxl[ijidx*sz_rot:(ijidx+1)*sz_rot-1 ,kidx] = (EVALS[i,j,k,*])[egvl_iter]

								trace_mat[(kidx/2)*sz_k+ijidx] = ADT[i_slice[ijidx],j_slice[ijidx],k_iter,1] + $
									ADT[i_slice[ijidx],j_slice[ijidx],k_iter,2] + ADT[i_slice[ijidx],j_slice[ijidx],k_iter,3]

							ENDIF ELSE BEGIN
								egvc_vxl[ijidx*sz_rot:(ijidx+1)*sz_rot-1 ,kidx] = (EVECS[i_slice[ijidx],j_slice[ijidx],k_iter,*,*])[egvc_iter]
								egvl_vxl[ijidx*sz_rot:(ijidx+1)*sz_rot-1 ,kidx] = (EVALS[i_slice[ijidx],j_slice[ijidx],k_iter,*])[egvl_iter]
							ENDELSE
						ENDFOR
					ENDFOR

					FOR kidx = 0, (sz_k/2)-1 DO BEGIN
						IF kidx EQ 0 THEN BEGIN
							d_dp_slf_all = egvc_slf[*,kidx*2] * SQRT(egvl_slf[*,kidx*2]) * $
								egvc_slf[*,kidx*2+1] * SQRT(egvl_slf[*,kidx*2+1])

							d_dp_vxl_all = egvc_vxl[*,kidx*2] * SQRT(egvl_vxl[*,kidx*2]) * $
								egvc_vxl[*,kidx*2+1] * SQRT(egvl_vxl[*,kidx*2+1])
						ENDIF ELSE BEGIN
							d_dp_slf_all = [d_dp_slf_all, egvc_slf[*,kidx*2] * SQRT(egvl_slf[*,kidx*2]) * $
								egvc_slf[*,kidx*2+1] * SQRT(egvl_slf[*,kidx*2+1])]

							d_dp_vxl_all = [d_dp_vxl_all, egvc_vxl[*,kidx*2] * SQRT(egvl_vxl[*,kidx*2]) * $
								egvc_vxl[*,kidx*2+1] * SQRT(egvl_vxl[*,kidx*2+1])]
						ENDELSE
					ENDFOR

					; calculate the d_dp_r for each of the neighboring voxels with itself

					FOR idx = 0, num_3D_neighbors-1 DO BEGIN
						d_dp_slf_iter = 0
						FOR elmt_idx=num_rot_elements*idx,num_rot_elements*(idx+1)-1 DO BEGIN
							d_dp_slf_iter = d_dp_slf_iter + TOTAL(d_dp_slf_all[(idx_thr+(elmt_idx*3))])^2
						ENDFOR

						d_dp_vxl_iter = 0
						FOR elmt_idx=num_rot_elements*idx,num_rot_elements*(idx+1)-1 DO BEGIN
							d_dp_vxl_iter = d_dp_vxl_iter + TOTAL(d_dp_vxl_all[(idx_thr+(elmt_idx*3))])^2
						ENDFOR

						d_dp_vxl_ital = d_dp_vxl_iter - (1.0/3 * trace_ref * trace_mat[idx])

						kidx = FLOOR(idx / num_3D_neighbors) * 2

						; d_dp_ref_slf calculated above for 2D neighbors
						lattice_3D[i,j,idx]  = SQRT((3*d_dp_vxl_ital)/(8*d_dp_vxl_iter)) + (3.0/4)*(d_dp_vxl_ital/(SQRT(d_dp_ref_slf*d_dp_slf_iter)))
						lattice_wgt[i,j,idx] = SQRT(((i_slice[idx MOD (sz_k/2)]-i)*xres)^2 + $
							((j_slice[idx MOD (sz_k/2)]-j)*yres)^2 + ((k_slice[kidx]-k)*zres)^2)
					ENDFOR
				ENDIF
			ENDFOR
			progressBar -> Update, (float(j+1)/float(ydim))*100.0
		ENDFOR

		nan_idx = ~FINITE(lattice_n)
		IF TOTAL(nan_idx) NE 0 THEN	lattice_n[WHERE(nan_idx EQ 1)] = 0

		IF LA_3D_FLAG EQ 1 THEN BEGIN
			nan_idx = ~FINITE(lattice_3D)
			IF TOTAL(nan_idx) NE 0 THEN	lattice_3D[WHERE(nan_idx EQ 1)] = 0
		ENDIF

		FOR j=1, ydim-2 DO BEGIN
			FOR i=1, xdim-2 DO BEGIN
				i_p = 2*i+2
				j_p = 2*j+2

				norm_fac = MIN([xres,yres])
				diag_res = SQRT(xres^2 + yres^2)

				adj_lr = (lattice_n[i_p-1,j_p] + lattice_n[i_p+1,j_p]) * (norm_fac/xres)
				adj_tb = (lattice_n[i_p,j_p-1] + lattice_n[i_p,j_p+1]) * (norm_fac/yres)
				diags  = (lattice_n[i_p-1,j_p-1] + lattice_n[i_p-1,j_p+1] + $
					lattice_n[i_p+1,j_p-1] + lattice_n[i_p+1,j_p+1]) * (norm_fac/diag_res)

				IF LA_3D_FLAG EQ 1 THEN BEGIN
					iter_wghts = REPLICATE(norm_fac,num_3D_neighbors) / REFORM(lattice_wgt[i,j,*])
					lattice_3D_iter = TOTAL(REFORM(lattice_3D[i,j,*]) * iter_wghts)

					norm_total = 2 * (norm_fac/xres) + 2 * (norm_fac/yres) + 4 * (norm_fac/diag_res) + TOTAL(iter_wghts)

					LA[i,j,k] = (adj_lr + adj_tb + diags + lattice_3D_iter) / norm_total
				ENDIF ELSE BEGIN
					norm_total = 2 * (norm_fac/xres) + 2 * (norm_fac/yres) + 4 * (norm_fac/diag_res)
					LA[i,j,k] = (adj_lr + adj_tb + diags) / norm_total
				ENDELSE
			ENDFOR
			progressBar -> Update, (float(j+1)/float(ydim))*100.0
		ENDFOR
		progressbar -> Destroy
	ENDFOR

	project.dataArray[CI].latt_Ani = ptr_new(LA)
END


; Subroutine name: update_ADT_Slider_And_label
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.


pro update_ADT_Slider_And_label

    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets

    CI = project.CI
    display_type = project.procpramArray[CI].adt_display_type
    ADT_Start = project.procPramArray[CI].ADT_Start
    EigenVec_Start = project.procPramArray[CI].EigenVec_Start

    ;change the label above the ADT slider so that when
    ;adt type 0|1 is active then slider and label are
    ;sensitive
    ;else they are not sensitive
    if display_type eq 0 or display_type eq 1 then begin
       widget_control, Slider_ADT_Start_label , SENSITIVE=1
       widget_control, Slider_ADT_Start , SENSITIVE=1
    end
    if display_type eq 2 or display_type eq 3 or display_type eq 4 then begin
       widget_control, Slider_ADT_Start_label , SENSITIVE=0
       widget_control, Slider_ADT_Start , SENSITIVE=0
       widget_control, wid_ADT_MOVIE, sensitive=1
    end


    if display_type eq 0 then begin

       widget_control , Slider_ADT_Start , SET_VALUE = ADT_Start, SET_SLIDER_MAX=9

       if ADT_Start eq 0 then  widget_control, Slider_ADT_Start_label, SET_VALUE ='S0'
       if ADT_Start eq 1 then  widget_control, Slider_ADT_Start_label, SET_VALUE ='XX'
       if ADT_Start eq 2 then  widget_control, Slider_ADT_Start_label, SET_VALUE ='YY'
       if ADT_Start eq 3 then  widget_control, Slider_ADT_Start_label, SET_VALUE ='ZZ'
       if ADT_Start eq 4 then  widget_control, Slider_ADT_Start_label, SET_VALUE ='XY'
       if ADT_Start eq 5 then  widget_control, Slider_ADT_Start_label, SET_VALUE ='XZ'
       if ADT_Start eq 6 then  widget_control, Slider_ADT_Start_label, SET_VALUE ='YZ'
       if ADT_Start eq 7 then  widget_control, Slider_ADT_Start_label, SET_VALUE ='Color Trace'
       if ADT_Start eq 8 then  widget_control, Slider_ADT_Start_label, SET_VALUE ='Orientation'
       if ADT_Start eq 9 then  widget_control, Slider_ADT_Start_label, SET_VALUE ='Color Sphere'

       widget_control, wid_ADT_MOVIE, sensitive=(ADT_Start eq 9) ? 0 : 1

    end
    if display_type eq 1 then begin

       widget_control , Slider_ADT_Start , SET_VALUE = EigenVec_Start, SET_SLIDER_MAX=2
       if EigenVec_Start eq 0 then  widget_control, Slider_ADT_Start_label, SET_VALUE ='E1'
       if EigenVec_Start eq 1 then  widget_control, Slider_ADT_Start_label, SET_VALUE ='E2'
       if EigenVec_Start eq 2 then  widget_control, Slider_ADT_Start_label, SET_VALUE ='E3'

    end
end


; Subroutine name: redraw_ADT
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

PRO redraw_ADT
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    CI = project.CI


    IF N_ELEMENTS(ADT_WINDOW_BASE) EQ 1 THEN $
      IF WIDGET_INFO(ADT_WINDOW_BASE, /VALID_ID ) EQ 1 THEN BEGIN

        ;;if the current scan does not have a good b-matrix then make every thing insensitive
        SENSITIVE =  ptr_valid(project.imndArray[CI].b_matrix)

        title = STRING('Apparent Diffusion Tensor Parameters Scan #', STRTRIM(CI+1,2))
        WIDGET_CONTROL, ADT_WINDOW_BASE, BASE_SET_TITLE=title

        WIDGET_CONTROL, SLIDER_ADTThreshold, SENSITIVE=SENSITIVE

        update_ADT_Slider_And_label

        WIDGET_CONTROL, Slider_ADT_Start, SENSITIVE=SENSITIVE
        WIDGET_CONTROL, ADT_ROI, SENSITIVE=SENSITIVE
        WIDGET_CONTROL, ADT_FLT, SENSITIVE=SENSITIVE
        WIDGET_CONTROL, ADT_TIF, SENSITIVE=SENSITIVE
;			WIDGET_CONTROl, ADT_FIT, SENSITIVE=0;;SENSITIVE
        WIDGET_CONTROL, ADT_DISPLAY, SENSITIVE=SENSITIVE

        WIDGET_CONTROL, wid_ADT_MOVIE, SENSITIVE=SENSITIVE*((project.procPramArray[CI].ADT_Start eq 9 ? 0 : 1)+(project.procpramArray[CI].adt_display_type ne 0 ? 1 : 0))

        WIDGET_CONTROL, wid_ADT_SLICER, SENSITIVE=SENSITIVE

        ;; this will update the two sets of exclusive buttons so that they act like one set
        IF project.procPramArray[CI].adt_display_multi EQ 1 THEN BEGIN
            WIDGET_CONTROL, ADT_MULTI_display, SENSITIVE=SENSITIVE, SET_BUTTON=1
            WIDGET_CONTROL, ADT_Type0, SENSITIVE=SENSITIVE, SET_BUTTON=0
            WIDGET_CONTROL, ADT_Type4, SENSITIVE=SENSITIVE, SET_BUTTON=0
            WIDGET_CONTROL, ADT_Type2, SENSITIVE=SENSITIVE, SET_BUTTON=0
            WIDGET_CONTROL, ADT_Type3, SENSITIVE=SENSITIVE, SET_BUTTON=0
            WIDGET_CONTROL, ADT_Type5, SENSITIVE=SENSITIVE, SET_BUTTON=0
            WIDGET_CONTROL, ADT_Type6, SENSITIVE=SENSITIVE, SET_BUTTON=0
        ENDIF ELSE BEGIN
            IF project.procPramArray[CI].adt_display_type EQ 0 THEN BEGIN
                WIDGET_CONTROL, ADT_MULTI_display, SENSITIVE=SENSITIVE, SET_BUTTON=0
                WIDGET_CONTROL, ADT_Type0, SENSITIVE=SENSITIVE, SET_BUTTON=1
                WIDGET_CONTROL, ADT_Type4, SENSITIVE=SENSITIVE, SET_BUTTON=0
                WIDGET_CONTROL, ADT_Type2, SENSITIVE=SENSITIVE, SET_BUTTON=0
                WIDGET_CONTROL, ADT_Type3, SENSITIVE=SENSITIVE, SET_BUTTON=0
                WIDGET_CONTROL, ADT_Type5, SENSITIVE=SENSITIVE, SET_BUTTON=0
                WIDGET_CONTROL, ADT_Type6, SENSITIVE=SENSITIVE, SET_BUTTON=0
            ENDIF ELSE BEGIN
                WIDGET_CONTROL, ADT_MULTI_display, SENSITIVE=SENSITIVE, SET_BUTTON=0
                WIDGET_CONTROL, ADT_Type0, SENSITIVE=SENSITIVE, SET_BUTTON=0
                WIDGET_CONTROL, ADT_Type4, SENSITIVE=SENSITIVE, SET_BUTTON=0
                WIDGET_CONTROL, ADT_Type2, SENSITIVE=SENSITIVE, SET_BUTTON=0
                WIDGET_CONTROL, ADT_Type3, SENSITIVE=SENSITIVE, SET_BUTTON=0
                WIDGET_CONTROL, ADT_Type5, SENSITIVE=SENSITIVE, SET_BUTTON=0
                WIDGET_CONTROL, ADT_Type6, SENSITIVE=SENSITIVE, SET_BUTTON=0

                CASE project.procPramArray[CI].adt_display_type OF
                    4: WIDGET_CONTROL, ADT_Type4, SET_BUTTON=1
                    2: WIDGET_CONTROL, ADT_Type2, SET_BUTTON=1
                    3: WIDGET_CONTROL, ADT_Type3, SET_BUTTON=1
                    5: WIDGET_CONTROL, ADT_Type5, SET_BUTTON=1
                    6: WIDGET_CONTROL, ADT_Type6, SET_BUTTON=1
                    ELSE:
                ENDCASE
            ENDELSE
        ENDELSE

                                ;check to see if they are running in the virtual machine if they are then dont
                                ;allow mpeg creation
        WIDGET_CONTROL, wid_ADT_MPG ,SENSITIVE= (lmgr(/vm) gt 0) ? 0 : SENSITIVE

        WIDGET_CONTROL, ADT_scale_avg_d_max, SENSITIVE=SENSITIVE
        WIDGET_CONTROL, ADT_scale_avg_d_min, SENSITIVE=SENSITIVE
        WIDGET_CONTROL, ADT_scale_off_dia, SENSITIVE=SENSITIVE
        WIDGET_CONTROL, ADT_scale_dia_max, SENSITIVE=SENSITIVE
        WIDGET_CONTROL, ADT_scale_dia_min, SENSITIVE=SENSITIVE
        WIDGET_CONTROL, ADT_scale_fa_max, SENSITIVE=SENSITIVE

        WIDGET_CONTROL, ADT_scale_avg_d_max, SET_VALUE=project.procPramArray[CI].avg_d_max
        WIDGET_CONTROL, ADT_scale_avg_d_min, SET_VALUE=project.procPramArray[CI].avg_d_min
        WIDGET_CONTROL, ADT_scale_off_dia, SET_VALUE=project.procPramArray[CI].off_dia
        WIDGET_CONTROL, ADT_scale_dia_max, SET_VALUE=project.procPramArray[CI].dia_max
        WIDGET_CONTROL, ADT_scale_dia_min, SET_VALUE=project.procPramArray[CI].dia_min
        WIDGET_CONTROL, ADT_scale_fa_max, SET_VALUE=project.procPramArray[CI].fa_max
        WIDGET_CONTROL, SLIDER_ADTThreshold, SET_VALUE=project.procPramArray[CI].adtThreshold*100

        ;; Reupdate if necessary
;        if (sensitive ne 0) then begin
;            widget_control, wid_adt_movie, sensitive=(1-project.procpramarray[ci].adt_glyph_overlay)
;            widget_control, wid_adt_mpg, sensitive=(1-project.procpramarray[ci].adt_glyph_overlay)
;            widget_control, wid_adt_slicer, sensitive=(1-project.procpramarray[ci].adt_glyph_overlay)
;        endif



    endif

end


; Subroutine name: mas_adt_regress_event
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.
    ;Edited by CD 2007/06/28
    ;Added DOT GUI events
    ;Edited by MK 2014/11/28
    ;1) Added selective curve fit tool
    ;2) Zero pad DWI option included

PRO mas_adt_regress_event, Event
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets

    CI = project.CI
    case Event.id of
        SLIDER_ADTThreshold: begin
            project.procPramArray[CI].adtThreshold = float(float(event.value)/100.0)
            project.procpramArray[project.CI].adt_proccess_flag = 0
        end

        ADT_FLT: begin
          adt_save_flt
        end
        ADT_TIF: begin
          adt_save_tif
        end

 ;       ADT_FIT : begin
  ;         mas_smooth
   ;     end
        ADT_DISPLAY : begin
           display_adt
        end
        ADT_SELECT_DATA : select_data_gui, index = CI
        ADT_Type0:begin
           project.procpramArray[CI].adt_display_type = 0
           project.procPramArray[CI].adt_display_multi = 0
          update_ADT_Slider_And_label
         end
        ADT_Type4:begin
         project.procpramArray[CI].adt_display_type = 4
         project.procPramArray[CI].adt_display_multi = 0
           update_ADT_Slider_And_label
         end
        ADT_Type2:begin
           project.procpramArray[CI].adt_display_type = 2
           project.procPramArray[CI].adt_display_multi = 0
         update_ADT_Slider_And_label
         end
        ADT_Type3:begin
         project.procpramArray[CI].adt_display_type = 3
         project.procPramArray[CI].adt_display_multi = 0
         update_ADT_Slider_And_label
         end
         ADT_Type5:begin
             project.procpramArray[CI].adt_display_type = 5
             project.procPramArray[CI].adt_display_multi = 0
             update_ADT_Slider_And_label
         end

         ADT_Type6:begin
             project.procpramArray[CI].adt_display_type = 6
             project.procPramArray[CI].adt_display_multi = 0
             update_ADT_Slider_And_label
         end

         Slider_ADT_Start:begin
          if project.procpramArray[CI].adt_display_type eq 0 then $
              project.procPramArray[CI].ADT_Start = event.value

          if project.procpramArray[CI].adt_display_type eq 1 then $
              project.procPramArray[CI].EigenVec_Start = event.value
          update_ADT_Slider_And_label
         end

         ADT_ROI: adt_ROI

         ADT_MULTI_display: begin
            project.procPramArray[CI].adt_display_multi = 1
            project.procpramArray[CI].adt_display_type = -1
            ;if the users selectes mulit and they have alos slected diff-weighted then change
            ;diff weighted setting to display diff-Tensor.

            ;if project.procpramArray[CI].adt_display_type eq 4 then begin

            ;
             ;  update_ADT_Slider_And_label
            ; end
        end
         wid_ADT_MOVIE: ADT_MOVIE

         wid_ADT_MPG:ADT_MPEG_WIDGET

         wid_ADT_SLICER:write_slicer_diffusion_image

         ADT_scale_avg_d_max:begin
           if event.value le project.procPramArray[CI].avg_d_min then $
             project.procPramArray[CI].avg_d_max = project.procPramArray[CI].avg_d_min $
           else $
               project.procPramArray[CI].avg_d_max = event.value

            widget_control, ADT_scale_avg_d_max , SET_VALUE= project.procPramArray[CI].avg_d_max

         end
         ADT_scale_avg_d_min:begin
          if event.value ge project.procPramArray[CI].avg_d_max then $
              project.procPramArray[CI].avg_d_min = project.procPramArray[CI].avg_d_max $
          else $
              project.procPramArray[CI].avg_d_min = event.value

          widget_control, ADT_scale_avg_d_min , SET_VALUE= project.procPramArray[CI].avg_d_min

         end

         ADT_scale_off_dia:project.procPramArray[CI].off_dia = event.value

         ADT_scale_dia_max:begin
          if event.value le project.procPramArray[CI].dia_min then $
              project.procPramArray[CI].dia_max = project.procPramArray[CI].dia_min $
          else $
              project.procPramArray[CI].dia_max = event.value

          widget_control, ADT_scale_dia_max , SET_VALUE= project.procPramArray[CI].dia_max

         end

         ADT_scale_dia_min:begin
          if event.value ge project.procPramArray[CI].dia_max then $
              project.procPramArray[CI].dia_min = project.procPramArray[CI].dia_max $
          else $
              project.procPramArray[CI].dia_min = event.value

          widget_control, ADT_scale_dia_min , SET_VALUE= project.procPramArray[CI].dia_min
         end

         ADT_scale_fa_max:project.procPramArray[CI].fa_max = event.value

;         ADT_GLYPH_OVERLAY: begin
;             project.procpramarray[ci].adt_glyph_overlay = event.select
;         end

;         ADT_GLYPH_CONFIG: begin
;             mas_adt_glyph_configure
;         end

         ADT_REG_STANDARD: begin
             project.procpramarray[ci].adt_proccess_flag=0
             project.procpramarray[ci].adt_regress_method = 1
         end

         ADT_REG_MOW: begin
             project.procpramarray[ci].adt_proccess_flag=0
             project.procpramarray[ci].adt_regress_method = 0
         end

         else:

    endcase

	redraw_ADT
end


; Subroutine name: mas_adt_regress
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.
    ;Edited by CD 2007/06/28
    ;Added DOT GUI
    ;Edited by MK 2014/11/28
    ;1) Added selective curve fit tool
    ;2) Zero pad DWI option included

PRO mas_adt_regress

    HEAP_GC
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    CI = project.CI

    IF N_ELEMENTS(ADT_WINDOW_BASE) EQ 1 THEN $
      IF WIDGET_INFO(ADT_WINDOW_BASE, /VALID_ID ) EQ 1 THEN RETURN

    ;;if the current scan does not have a good b-matrix then make every thing insensitive
    SENSITIVE =  PTR_VALID(project.imndArray[CI].b_matrix)

    title = STRING('Apparent Diffusion Tensor Parameters Scan #', STRTRIM(CI+1,2))
    ADT_WINDOW_BASE = WIDGET_BASE(TITLE=title, UVALUE = ADT_WINDOW_BASE, ROW=2, XOFFSET=420)

    ;; outline basic structure

    ADT_BASE_ROW1 = WIDGET_BASE(ADT_WINDOW_BASE, /ROW)

    ADT_ROW1_COL1 = WIDGET_BASE(ADT_BASE_ROW1, /COLUMN, FRAME=1)
    ADT_ROW1_COL2 = WIDGET_BASE(ADT_BASE_ROW1, /COLUMN, FRAME=1)
    ADT_ROW1_COL3 = WIDGET_BASE(ADT_BASE_ROW1, /COLUMN, FRAME=1)

    ADT_BASE_ROW2 = WIDGET_BASE(ADT_WINDOW_BASE, /COLUMN, FRAME=1)

;    ADT_BASE_ROW2_t = WIDGET_TAB(ADT_WINDOW_BASE)
;    ADT_BASE_ROW2 = widget_base(adt_base_row2_t, /column, title='Display Scaling')

    trash = WIDGET_LABEL(ADT_BASE_ROW2, VALUE='ADT Display Scaling')
    ADT_ROW2_ROW1 = WIDGET_BASE(ADT_BASE_ROW2, /ROW)
    ADT_ROW2_ROW2 = WIDGET_BASE(ADT_BASE_ROW2, /ROW)

    ;; row 1, column 1

    trash = WIDGET_LABEL(ADT_ROW1_COL1 ,VALUE='Save Output')

    ADT_FLT = WIDGET_BUTTON(ADT_ROW1_COL1, UNAME='ADT_FLT', $
                            /ALIGN_LEFT ,VALUE='Save flt files', $
                            SENSITIVE=SENSITIVE)

    ADT_TIF = WIDGET_BUTTON(ADT_ROW1_COL1, UNAME='ADT_TIF' ,/ALIGN_LEFT ,VALUE='Save tif files', SENSITIVE=SENSITIVE)

    wid_ADT_MPG = WIDGET_BUTTON(ADT_ROW1_COL1,value='Save MPEG Movie', SENSITIVE=SENSITIVE)

    wid_ADT_SLICER = WIDGET_BUTTON(ADT_ROW1_COL1,value='Save To Slicer', SENSITIVE=SENSITIVE)

    trash = widget_label(adt_row1_col1, value="ROI Processing")
    ADT_ROI = WIDGET_BUTTON(ADT_ROW1_COL1,value='ROI Tools',SENSITIVE=SENSITIVE)

    ;; row 1, column 2

    trash = WIDGET_LABEL(ADT_ROW1_COL2 ,VALUE='Display Selection')

    ex_base = WIDGET_BASE(ADT_ROW1_COL2 ,COLUMN=1 ,/EXCLUSIVE, SENSITIVE=SENSITIVE)

    ADT_MULTI_display = WIDGET_BUTTON(ex_base,value='All DTI Parameters', SENSITIVE=SENSITIVE)
    IF project.procPramArray[CI].adt_display_multi EQ 1 THEN WIDGET_CONTROL, /SET_BUTTON, ADT_TIF

    ADT_Type2 = WIDGET_BUTTON(ex_base, /ALIGN_LEFT ,VALUE='Average Diffusion')
    if project.procPramArray[CI].adt_display_type eq 2 then WIDGET_CONTROL, /SET_BUTTON, ADT_Type2, SENSITIVE=SENSITIVE

    ADT_Type3 = WIDGET_BUTTON(ex_base, /ALIGN_LEFT ,VALUE='Fract Anisotropy')
    if project.procPramArray[CI].adt_display_type eq 3 then WIDGET_CONTROL, /SET_BUTTON, ADT_Type3, SENSITIVE=SENSITIVE

    ADT_Type5 = WIDGET_BUTTON(ex_base, /ALIGN_LEFT ,VALUE='2D Lattice Anisotropy')
    if project.procPramArray[CI].adt_display_type eq 5 then WIDGET_CONTROL, /SET_BUTTON, ADT_Type5, SENSITIVE=SENSITIVE

    ADT_Type6 = WIDGET_BUTTON(ex_base, /ALIGN_LEFT ,VALUE='3D Lattice Anisotropy')
    if project.procPramArray[CI].adt_display_type eq 6 then WIDGET_CONTROL, /SET_BUTTON, ADT_Type6, SENSITIVE=SENSITIVE

    ADT_Type4 = WIDGET_BUTTON(ex_base, /ALIGN_LEFT ,VALUE='Diffusion-Weighted')
    if project.procPramArray[CI].adt_display_type eq 4 then WIDGET_CONTROL, /SET_BUTTON, ADT_Type4, SENSITIVE=SENSITIVE

    ADT_Type0 = WIDGET_BUTTON(ex_base, /ALIGN_LEFT ,VALUE='Diffusion Tensor:', SENSITIVE=SENSITIVE)
    IF project.procPramArray[CI].adt_display_type EQ 0 THEN WIDGET_CONTROL, /SET_BUTTON, ADT_Type0

    Slider_ADT_Start_label = WIDGET_LABEL(ADT_ROW1_COL2 ,VALUE='S0', XSIZE=90, /ALIGN_CENTER, SENSITIVE=SENSITIVE)

    Slider_ADT_Start = WIDGET_SLIDER(ADT_ROW1_COL2, MINIMUM=0, FRAME=0, MAXIMUM=9,/SUPPRESS_VALUE,  $
                                     VALUE= project.procPramArray[CI].ADT_Start, SENSITIVE=SENSITIVE)

;    ADT_DISPLAY = WIDGET_BUTTON(ADT_ROW1_COL2,value='Display' ,UNAME = 'ADT_DISPLAY', SENSITIVE=SENSITIVE)

;    wid_ADT_MOVIE = WIDGET_BUTTON(ADT_ROW1_COL2,value='Movie', SENSITIVE=SENSITIVE)

;    ADT_ROI = WIDGET_BUTTON(ADT_ROW1_COL2,value='ROI',SENSITIVE=SENSITIVE)

    ;; row 1, column 3

    ex_base = widget_base(ADT_ROW1_COL3, /column, xpad=0, ypad=0)
    ex_base1 = widget_base(ex_base, /column, /nonexclusive, xpad=0, ypad=0)
    trash = WIDGET_LABEL(ex_base, VALUE='Threshold')

    SLIDER_ADTThreshold = CW_FSLIDER(ex_base, $ ;ADT_ROW1_COL3,         $
                                     MINIMUM = 0,                  $
                                     /edit,                        $
                                     FRAME = 0,                    $
                                     MAXIMUM = 20,                 $
                                     SCROLL=0.1, /DRAG,            $
                                     TITLE = '% Threshold',        $
                                     UNAME ='SLIDER_ADTThreshold', $
                                     VALUE = project.procPramArray[CI].adtThreshold * 100)

    WIDGET_CONTROL, SLIDER_ADTThreshold, SENSITIVE=SENSITIVE
    
    ADT_SELECT_DATA   = WIDGET_BUTTON(ex_base,value = 'Select Data to Fit', UNAME = 'ADT_SELECT_DATA', SENSITIVE=SENSITIVE)    
  ;;  ADT_FIT = WIDGET_BUTTON(ADT_ROW1_COL3,value='Choose Regression Type' ,UNAME = 'mas_smooth', SENSITIVE=0);SENSITIVE)


    ADT_DISPLAY = WIDGET_BUTTON(ex_base,value='Display' ,UNAME = 'ADT_DISPLAY', SENSITIVE=SENSITIVE)

    wid_ADT_MOVIE = WIDGET_BUTTON(ex_base,value='Movie', SENSITIVE=SENSITIVE)
;    ex_base = WIDGET_BASE(ADT_ROW1_COL3, UNAME='ex_base' ,COLUMN=1 ,/EXCLUSIVE, SENSITIVE=SENSITIVE)


;    ADT_MULTI_display = WIDGET_BUTTON(ex_base,value='All DTI Parameters', SENSITIVE=SENSITIVE)
;    IF project.procPramArray[CI].adt_display_multi EQ 1 THEN WIDGET_CONTROL, /SET_BUTTON, ADT_TIF

;     ADT_Type0 = WIDGET_BUTTON(ex_base, /ALIGN_LEFT ,VALUE='Diffusion Tensor', SENSITIVE=SENSITIVE)
;     IF project.procPramArray[CI].adt_display_type EQ 0 THEN WIDGET_CONTROL, /SET_BUTTON, ADT_Type0

;     Slider_ADT_Start_label = WIDGET_LABEL(ADT_ROW1_COL3 ,VALUE='S0', XSIZE=90, /ALIGN_CENTER, SENSITIVE=SENSITIVE)

;     Slider_ADT_Start = WIDGET_SLIDER(ADT_ROW1_COL3, MINIMUM=0, FRAME=0, MAXIMUM=8,/SUPPRESS_VALUE,  $
;                                      VALUE= project.procPramArray[CI].ADT_Start, SENSITIVE=SENSITIVE)

;     update_ADT_Slider_And_label

;     ADT_DISPLAY = WIDGET_BUTTON(ADT_ROW1_COL3,value='Display' ,UNAME = 'ADT_DISPLAY', SENSITIVE=SENSITIVE)

;     wid_ADT_MOVIE = WIDGET_BUTTON(ADT_ROW1_COL3,value='Movie', SENSITIVE=SENSITIVE)

    ;; row 2, row 1

    start_value =  project.procPramArray[CI].avg_d_min
    ADT_scale_avg_d_min = CW_FSLIDER(ADT_ROW2_ROW1, MINIMUM=0.0, FRAME=0, MAXIMUM=0.001, $
                                     TITLE = 'Avg Diff Min' ,  /EDIT , VALUE= start_value)

    start_value =  project.procPramArray[CI].avg_d_max
    ADT_scale_avg_d_max = CW_FSLIDER(ADT_ROW2_ROW1, MINIMUM=0.0, FRAME=0, MAXIMUM=0.01, $
                                     TITLE = 'Avg Diff Max' ,  /EDIT , VALUE=start_value)

    start_value =  project.procPramArray[CI].fa_max
    ADT_scale_fa_max = CW_FSLIDER(ADT_ROW2_ROW1, MINIMUM=0.0, FRAME=0, MAXIMUM=1.0, $
                                  TITLE = 'Fractional Ani Max' ,  /EDIT , VALUE= start_value)

    ;; row 2, row 2

    start_value =  project.procPramArray[CI].dia_min
    ADT_scale_dia_min = CW_FSLIDER(ADT_ROW2_ROW2, MINIMUM=0.0, FRAME=0, MAXIMUM=0.004, $
                                   TITLE = 'Diagonal Min' ,  /EDIT , VALUE=start_value)

    start_value =  project.procPramArray[CI].dia_max
    ADT_scale_dia_max = CW_FSLIDER(ADT_ROW2_ROW2, MINIMUM=0.0, FRAME=0, MAXIMUM=0.004, $
                                   TITLE = 'Diagonal Max' ,  /EDIT , VALUE= start_value)

    start_value =  project.procPramArray[CI].off_dia
    ADT_scale_off_dia = CW_FSLIDER(ADT_ROW2_ROW2, MINIMUM=0.0, FRAME=0, MAXIMUM=0.001, $
                                   TITLE = 'Off Diagonal' ,  /EDIT , VALUE= start_value)
                                           
    update_ADT_Slider_And_label


    WIDGET_CONTROL, ADT_WINDOW_BASE, /realize

    XMANAGER, 'mas_adt_regress', ADT_WINDOW_BASE, /no_block, GROUP_LEADER=WID_BASE_MAIN

END
