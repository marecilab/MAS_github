;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

;comments for this procedure can be found bellow


; Subroutine name: mas_subset_3d
; Created by
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_subset_3d , scan
    COMPILE_OPT IDL2
    common scan_data
    single_Multi_flag = project.procPramArray[project.ci].single_Multi_flag
    slice_axis     =  project.procPramArray[project.ci].slice_axis
    fdim_start     =  project.procPramArray[project.ci].fdim_start
    pdim_start     =  project.procPramArray[project.ci].pdim_start
    sdim_start     =  project.procPramArray[project.ci].sdim_start
    adim_start     =  project.procPramArray[project.ci].adim_start
    sort_dir       =  project.procPramArray[project.ci].sort_dir
    Freq_zoom   =   project.procPramArray[project.ci].Freq_zoom
    Phase_zoom   =   project.procPramArray[project.ci].Phase_zoom
    Slice_zoom   =   project.procPramArray[project.ci].Slice_zoom
    sz_scan         =    size(scan)



    if 1 eq 1 then begin

       if single_Multi_flag eq 0 then begin

         if slice_axis eq 0 then begin
          if Slice_zoom gt 1 then begin
              slice = float(sdim_start)/float(Slice_zoom)
              top_slice = ceil(slice)
              bot_slice = floor(slice)
                if top_slice eq bot_slice then begin
                 scan = reform(scan[*,*,top_slice,adim_start])
              end else begin
                 if top_slice ge sz_scan[3] then $
                   scan = reform(scan[*,*,bot_slice,adim_start]) $
                 else $
                   scan = reform(scan[*,*,bot_slice:top_slice,adim_start])
              end
          endif else scan = reform(scan[*,*,sdim_start,adim_start])
         endif;if slice_axis eq 0

         if slice_axis eq 1 then begin

          if Phase_zoom gt 1 then begin
              slice = float(pdim_start)/float(Phase_zoom)
              top_slice = ceil(slice)
              bot_slice = floor(slice)

                if top_slice eq bot_slice then begin
                 scan = reform(scan[*,top_slice,*,adim_start])
              end else begin

                 if top_slice ge sz_scan[3] then $
                   scan = reform(scan[*,bot_slice,*,adim_start]) $
                 else $
                   scan = reform(scan[*,bot_slice:top_slice,*,adim_start])
              end
          endif else scan = reform(scan[*,pdim_start,*,adim_start])
         endif

         if slice_axis eq 2 then begin

          if Freq_zoom gt 1 then begin
              slice = float(fdim_start)/float(Freq_zoom)
              top_slice = ceil(slice)
              bot_slice = floor(slice)

                if top_slice eq bot_slice then begin
                 scan = reform(scan[top_slice,*,*,adim_start])
              end else begin
                 if top_slice ge sz_scan[3] then $
                   scan = reform(scan[bot_slice,*,*,adim_start]) $
                 else $
                   scan = reform(scan[bot_slice:top_slice,*,*,adim_start])
              end
          endif else scan = reform(scan[fdim_start,*,*,adim_start])
         end
       endif;if single_Multi_flag eq 0

       if single_Multi_flag eq 1 then begin

          if sort_dir eq 0 then begin

              if slice_axis eq 0 then begin
                 scan = reform(scan[*,*,*,adim_start])
              endif
              if slice_axis eq 1 then begin
                 data_out = fltarr(sz_scan[1],sz_scan[3],sz_scan[2])
                 for ii=0, sz_scan[2]-1 do $
                   data_out[*,*,ii] = scan[*,ii,*,adim_Start]
                 scan = reform(data_out)
                 data_out = 0

              endif
              if slice_axis eq 2 then begin
                 data_out = fltarr(sz_scan[2],sz_scan[3],sz_scan[1])
                 for ii=0, sz_scan[1]-1 do $
                   data_out[*,*,ii] = scan[ii,*,*,adim_Start]
                 scan = reform(data_out)
                 data_out = 0
              end
          end;if sort_dir eq 0
          if sort_dir eq 1  then begin

              ;sort the data by the adim
                if slice_axis eq 0 then begin
                 scan = reform(scan[*,*,sdim_start,*])

              endif
              if slice_axis eq 1 then begin
                 data_out = fltarr(sz_scan[1],sz_scan[3],sz_scan[2])
                 for ii=0, sz_scan[2]-1 do $
                   data_out[*,*,ii] = scan[*,ii,sdim_start,*]
                 scan = reform(data_out)
                 data_out = 0

              endif
              if slice_axis eq 2 then begin
                 data_out = fltarr(sz_scan[2],sz_scan[3],sz_scan[1])
                 for ii=0, sz_scan[1]-1 do $
                   data_out[*,*,ii] = scan[ii,*,sdim_Start,*]
                 scan = reform(data_out)
                 data_out = 0
              end

          end;if sort_dir eq 1

       endif;if single_Multi_flag eq 1

    endif ;if zoom_3rd_dim_flag eq 1


end


; Subroutine name: mas_subset_2d
; Created by
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_subset_2d , scan
    COMPILE_OPT IDL2
    common scan_data
    single_Multi_flag = project.procPramArray[project.ci].single_Multi_flag
    slice_axis     =  project.procPramArray[project.ci].slice_axis
    fdim_start     =  project.procPramArray[project.ci].fdim_start
    pdim_start     =  project.procPramArray[project.ci].pdim_start
    sdim_start     =  project.procPramArray[project.ci].sdim_start
    adim_start     =  project.procPramArray[project.ci].adim_start
    sort_dir       =  project.procPramArray[project.ci].sort_dir
    sz_scan        =    size(scan)

    ;did the user select a single or multi image set.

    if single_Multi_flag eq 0 then begin
       ;the user wants a single image, which direction do they want to cut

       if slice_axis eq 0 then $ ;the user wants to slice the data on the freq phase axis
         scan = reform(scan[*,*,sdim_start,adim_start])

       if slice_axis eq 1 then $ ;the user want to slice the data on the freq slice axis
        scan = reform(scan[*,pdim_start,*,adim_start])

       if slice_axis eq 2 then $ ;the user want to slice the data on the phase slice axis
         scan = reform(scan[fdim_start,*,*,adim_start])

    endif; single_Multi_flag eq 0


    if single_Multi_flag eq 1 then begin
    ;the user wants a multi data set
    ;which direction should the data be sorted by

       if sort_dir eq 0 then begin
         ;sort the data by the sdim

         if slice_axis eq 0 then $ ;the user wants to slice the data on the freq phase axis
          scan = reform(scan[*,*,*,adim_start])

         if slice_axis eq 1 then begin
          ;the user wants to slice the data on the freq slice axis
          data_out = fltarr(sz_scan[1],sz_scan[3],sz_scan[2])
          for ii=0, sz_scan[2]-1 do $
              data_out[*,*,ii] = scan[*,ii,*,adim_Start]
          scan = reform(data_out)
          data_out = 0
         endif

         if slice_axis eq 2 then begin
          ;the user want to slice the data on the phase slice axis
          data_out = fltarr(sz_scan[2],sz_scan[3],sz_scan[1])
          for ii=0, sz_scan[1]-1 do $
              data_out[*,*,ii] = scan[ii,*,*,adim_Start]
          scan = reform(data_out)
          data_out = 0
         end

       end;if sort_dir eq 0

       if sort_dir eq 1  then begin

         ;sort the data by the adim
       if slice_axis eq 0 then $ ;the user wants to slice the data on the freq phase axis
          scan = reform(scan[*,*,sdim_start,*])

         if slice_axis eq 1 then begin
          ;the user want to slice the data on the freq slice axis
          data_out = fltarr(sz_scan[1],sz_scan[3],sz_scan[2])
          for ii=0, sz_scan[2]-1 do $
              data_out[*,*,ii] = scan[*,ii,sdim_start,*]
          scan = reform(data_out)
          data_out = 0
         endif

         if slice_axis eq 2 then begin
          ;the user want to slice the data on the phase slice axis
          data_out = fltarr(sz_scan[2],sz_scan[3],sz_scan[1])
          for ii=0, sz_scan[1]-1 do $
              data_out[*,*,ii] = scan[ii,*,sdim_Start,*]
          scan = reform(data_out)
          data_out = 0
         endif



       endif ;if sort_dir eq 1

    endif ;if single_Multi_flag eq 1
end


;*****************************************************************************************************
;
; NAME:
;   mas_subset
;
; PURPOSE:
;   take in the scan and give a smaller subset that is what the user specified.
;
; ARGUMENTS:
;
;
; MODIFICATION HISTORY:
; Edited by HS 2006/10/04.
; Fix spelling mistakes and commenting
;
;*****************************************************************************************************

function mas_subset , p_data
    COMPILE_OPT IDL2
    common scan_data

    single_Multi_flag = project.procPramArray[project.ci].single_Multi_flag
    slice_axis      =  project.procPramArray[project.ci].slice_axis
    fdim_start      =  project.procPramArray[project.ci].fdim_start
    pdim_start      =  project.procPramArray[project.ci].pdim_start
    sdim_start      =  project.procPramArray[project.ci].sdim_start
    adim_start      =  project.procPramArray[project.ci].adim_start
    sort_dir        =  project.procPramArray[project.ci].sort_dir
    ;;Freq_interp     =  project.procPramArray[project.ci].Freq_interp
    ;;Phase_interp    =  project.procPramArray[project.ci].Phase_interp
    ;;Slice_interp    =  project.procPramArray[project.ci].Slice_interp
    sz_p_data       =  size((*p_data))



    ;does the users want single or multiple slices?

    if single_Multi_flag eq 0 then begin
        ;the user wants a single image
        if slice_axis eq 0 then begin
            ;slice along the freq_phase
             slice = ptr_new((*p_data)[*,*,sdim_start,adim_start])

        end else if slice_axis eq 1 then begin
            ;slice along the freq_slice
            slice = ptr_new((*p_data)[*,pdim_start,*,adim_start])

        end else if slice_axis eq 2 then begin
            ;slice along the phase_slice
            slice = ptr_new((*p_data)[fdim_start,*,*,adim_start])

        end
        *slice = reform(*slice)
        return, slice

    end
    if single_Multi_flag eq 1 then begin
        ;the user wants a multiple image

        ;which dimension do they want to traverse
        if sort_dir eq 0 then begin
            ;sort the data by the sdim
            if slice_axis eq 0 then begin
                ;slice along the freq_phase
                ;we have to put the data together in the right direction

                p_multi_slices = ptr_new(fltarr(sz_p_data[1], sz_p_data[2], sz_p_data[3]))
                for ii=0, sz_p_data[3]-1 do $
                  (*p_multi_slices)[*,*,ii] = (*p_data)[*,*,ii,adim_start]

            end else if slice_axis eq 1 then begin
                ;slice along the freq_slice
                ;we have to put the data together in the right direction

                p_multi_slices = ptr_new(fltarr(sz_p_data[1], sz_p_data[3], sz_p_data[2]))
                for ii=0, sz_p_data[2]-1 do $
                  (*p_multi_slices)[*,*,ii] = reform((*p_data)[*,ii,*,adim_start])

            end else if slice_axis eq 2 then begin
                ;slice along the phase_slice
                ;we have to put the data together in the right direction

                p_multi_slices = ptr_new(fltarr(sz_p_data[2], sz_p_data[3], sz_p_data[1]))
                for ii=0, sz_p_data[1]-1 do $
                  (*p_multi_slices)[*,*,ii] = reform((*p_data)[ii,*,*,adim_start])
            end


        end else begin
            ;sort the data by the adim

            ;does the adim exist if it does not then there will only be one image
            ;else let them traverse the adim
            if project.imndArray[project.ci].adim eq 1 then sz_p_data[4] = 1



            if slice_axis eq 0 then begin
                ;slice along the freq_phase
                ;we have to put the data together in the right direction

                p_multi_slices = ptr_new(fltarr(sz_p_data[1], sz_p_data[2], sz_p_data[4]))
                for ii=0, sz_p_data[4]-1 do $
                  (*p_multi_slices)[*,*,ii] = (*p_data)[*,*,sdim_start,ii]

            end else if slice_axis eq 1 then begin
                ;slice along the freq_slice
                ;we have to put the data together in the right direction

                p_multi_slices = ptr_new(fltarr(sz_p_data[1], sz_p_data[3], sz_p_data[4]))
                for ii=0, sz_p_data[4]-1 do $
                  (*p_multi_slices)[*,*,ii] = reform((*p_data)[*,pdim_start,*,ii])

            end else if slice_axis eq 2 then begin
                ;slice along the phase_slice
                ;we have to put the data together in the right direction

                p_multi_slices = ptr_new(fltarr(sz_p_data[2], sz_p_data[3], sz_p_data[4]))
                for ii=0, sz_p_data[4]-1 do $
                  (*p_multi_slices)[*,*,ii] = reform((*p_data)[fdim_start,*,*,ii])
            end


        end


        return, p_multi_slices

    end


end
