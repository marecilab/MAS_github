;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

; Subroutine name: mas_Trace
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/02
    ;Fix spelling mistakes and commenting.


pro mas_Trace , phase=phase
    compile_Opt idl2
    COMMON scan_data

    mas_load_state_1

    CI = project.CI
    sdim_start = project.procpramArray[CI].sdim_start
    adim_start = project.procpramArray[CI].adim_start
    pfov = project.imndArray[ci].p_fov
    ffov = project.imndArray[ci].f_fov

    if keyword_set(phase) then begin


       sz_sig = size( (*project.dataArray[CI].state1))


        fdim = sz_sig[1]
        pdim = sz_sig[2]
        sdim = sz_sig[3]
        adim = sz_sig[4]

       trace = fltarr(pdim);(*project.dataArray[CI].state1)[*,*,sdim_start,adim_start]

       for ii=0 ,pdim-1 do begin
         trace[ii] = mean((*project.dataArray[CI].state1)[*,ii,sdim_start,adim_start])
       end

       trace = reverse(trace)

       print, pfov

       xaxis = findgen(fdim)/(fdim-1)*pfov
       ;print, xaxis


       IPLOT, xaxis, trace , $
       title = 'Mean Freq intensity vs Phase', $
       ;IDENTIFIER=iPlot_indentifer, $
       NAME= roi_name, $
       ;DIMENSIONS = [pfov,16.0] ,$
       YRange=[0.0,16.0], $
       XTitle = 'cm', $
       YTitle = 'Intensity'



       display_stats, strtrim(trace,2) , 'Trace in the Phase direction'


    end
END


; Subroutine name: mas_algebra_event
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; This procedure will allow users to +-*/and other operations on 2 2d or 3d images that have the same dimensions

; Editing Information:
    ;Edited by HS 2006/10/02
    ;Fix spelling mistakes and commenting.

pro mas_algebra_event, event
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets

    catch,error_status
    if (error_status ne 0) then begin
        help, calls= trace_back
        
        dummy = DIALOG_MESSAGE(['Please send a copy of this error message to the developer of this program.', $
                                !ERROR_STATE.msg,trace_back], /ERROR, $
                               TITLE='Error in WID_BASE_MAIN_event')
        RETURN
    endif

    CI = project.CI
    sdim_start = project.procpramArray[CI].sdim_start
    adim_start = project.procpramArray[CI].adim_start

    case Event.id of
        alg_load1: begin
            mas_load_state_1
            project.algebra1 = ptr_new((*project.dataArray[CI].state1)[*,*,sdim_start,adim_start])
            print, 'current image loaded into X'
        end
        
        alg_load2: begin
            mas_load_state_1
            project.algebra2 = ptr_new((*project.dataArray[CI].state1)[*,*,sdim_start,adim_start])
            print, 'current image loaded into Y'
        end
        
        alg_add: begin
            if ptr_valid(project.algebra1)eq 0 or ptr_valid(project.algebra2) eq 0 then begin
                update_status_bar,'please load X and Y with data'
                return
            end
            widget_control, alg_text, GET_VALUE = K
            print, k
            display = *project.algebra1 + (K*(*project.algebra2))
            print, size(display)
            iimage ,display
            
        end
        alg_sub: begin
            if ptr_valid(project.algebra1)eq 0 or ptr_valid(project.algebra2) eq 0 then begin
                update_status_bar,'please load X and Y with data'
                return
            end
            widget_control, alg_text, GET_VALUE = K
            print,k
            display = *project.algebra1 - (K*(*project.algebra2))
            print, size(display)
            iimage ,display
            
        end
        alg_mult:begin
            if ptr_valid(project.algebra1)eq 0 or ptr_valid(project.algebra2) eq 0 then begin
                update_status_bar,'please load X and Y with data'
                return
            end
            
            display = *project.algebra1 * *project.algebra2
            print, size(display)
            iimage ,display
            
        end
        alg_div:begin
            if ptr_valid(project.algebra1)eq 0 or ptr_valid(project.algebra2) eq 0 then begin
                update_status_bar,'please load X and Y with data'
                return
            end
            
            display = *project.algebra1 / *project.algebra2
            print, size(display)
            iimage ,display
            
        end
        alg_abs:begin
            if ptr_valid(project.algebra1)eq 0 or ptr_valid(project.algebra2) eq 0 then begin
                update_status_bar,'please load X and Y with data'
                return
            end
            
            display = abs(*project.algebra1)
            print, size(display)
            iimage ,display
            
        end
        alg_inv:begin
            if ptr_valid(project.algebra1)eq 0 or ptr_valid(project.algebra2) eq 0 then begin
                update_status_bar,'please load X and Y with data'
                return
            end
            
            display = 1/ *project.algebra1
            print, size(display)
            iimage ,display
            
        end
        alg_text: print,'this is just a holder to make sure all events are handeled'
    endcase
end


; Subroutine name: mas_algebra
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/02
    ;Fix spelling mistakes and commenting.

pro mas_algebra
    HEAP_GC
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    CI = project.CI

    IF N_ELEMENTS(alg_window_base) EQ 1 THEN $
    IF     widget_info( alg_window_base,  /VALID_ID ) eq 1 THEN return


    title = 'ImageAlgebra'
    alg_window_base = widget_base(TITLE=title,XOFFSET=420 ,/col )


    w=widget_label(alg_window_base,value='Image Algebra')
    alg_text = CW_FIELD(alg_window_base, /FLOATING ,title='K=', value='1.0')
    alg_load1 = widget_button(alg_window_base,value='              Load X             ' )
    alg_load2 = widget_button(alg_window_base,value='Load Y' )
    alg_add = widget_button(alg_window_base,value='X + K*Y' )
    alg_sub = widget_button(alg_window_base,value='X - K*Y' )
    alg_mult = widget_button(alg_window_base,value='X * Y' )
    alg_div = widget_button(alg_window_base,value='X / Y' )
    alg_abs = widget_button(alg_window_base,value='Abs( X )' )
    alg_inv = widget_button(alg_window_base,value='1/X' )


    widget_control, alg_window_base, /realize

    xmanager, 'mas_algebra',alg_window_base,/no_block

    mas_load_state_1
end
