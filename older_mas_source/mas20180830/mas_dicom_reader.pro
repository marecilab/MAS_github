;; $Id$

;;
;; Subroutine name: mas_dicom_get_dicom_obj
;; Created by: BT, 2012-02
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Returns a dicom object based on the type of machine IDL is running on.
;;
;; Editing Information:

pro mas_dicom_get_dicom_obj, dicom_obj=dicom_obj


    if (!VERSION.MEMORY_BITS eq 64 and !VERSION.OS eq 'darwin') then begin
        dicom_obj = obj_new('gdlffdicom')
    endif else if (!VERSION.MEMORY_BITS eq 64 and !VERSION.OS eq 'Win32') then begin
        dicom_obj = obj_new('gdlffdicom')
    endif else begin
        dicom_obj = obj_new('idlffdicom')
    endelse

end

;;
;; Subroutine name: mas_dicom_reader_checkval
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Checks that a value returned from idlffdicom->getValue() is a valid pointer
;;
;; Editing Information:

function mas_dicom_reader_checkval, val

    sz = size(val, /structure)
    
    if (sz.type_name eq 'POINTER') then begin
        
        for i = 0, n_elements(val)-1 do begin
            if (n_elements(*val[i]) ne 0) then begin
                return, 1
            endif else begin
                return, 0
            endelse
        endfor
    
    endif else begin
    
        return, 0
            
    endelse

end

;; $Id$
;;
;; Subroutine name: mas_dicom_group_acq
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Called recursively to group dicom array into acquisitions as specified by the user.
;;
;; Editing Information:

function mas_dicom_group_acq, dcm_array, indices, acq_indices, offset

    index = indices[0]
    
    reordering    = sort(dcm_array.(index))
    dcm_array_new = dcm_array[reordering]
    
    un = uniq(dcm_array_new.(index))
    
    if (n_elements(acq_indices) eq 0) then begin
        acq_indices = un + offset
    endif else begin
        acq_indices = [ acq_indices, un+offset ]
        acq_indices = acq_indices[sort(acq_indices)]
        acq_indices = acq_indices[uniq(acq_indices)]
    endelse

    ;; one field selected, we're done
    if (n_elements(indices) eq 1) then return, dcm_array_new
    
    start_ind = 0
    
    ;; now we recursively loop through the additional fields sorting on each one.
    for u = 0, n_elements(un)-1 do begin
        
        end_ind = un[u]
        
        dcm_array_new[start_ind:end_ind] = mas_dicom_group_acq(dcm_array_new[start_ind:end_ind], indices[1:*], acq_indices, start_ind)
        
        start_ind = un[u]+1

    endfor
    
    return, dcm_array_new
    
end

;; $Id$
;;
;; Subroutine name: mas_dicom_read_single_file
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Reads a single dicom file and parses the fields into the appropriate structure.
;;
;; Editing Information:

function mas_dicom_read_single_file, file, $
                                          reuse_object=reuse_object,$
                                          success=success, $
                                          images_only=images_only, $
                                          dcm_struct=dcm_struct
    
    destroy_on_return = 1
    
    if (keyword_set(reuse_object)) then begin

        if (not obj_valid(reuse_object)) then begin
            success = 0
            return, 0
            ;; ERROR
        endif else begin
            odicom = reuse_object
            odicom->Reset
            destroy_on_return = 0
        endelse
        
    endif else begin
    
        mas_dicom_get_dicom_obj, dicom_obj=odicom
        
    endelse

    read_ok = odicom->read(file)
    
    if (read_ok eq 0) then begin
        success = 0
        if (destroy_on_return eq 1) then begin
            obj_destroy, odicom
        endif
        
        return, read_ok
    endif else if (keyword_set(images_only)) then begin
        ;; only consider this if it has pixel data.
        tmp = odicom->getValue('7FE0'x, '0010'x, /no_copy)
        if (not mas_dicom_reader_checkval(tmp)) then begin
            print, "Skipping file: "+file+" (no pixel data)"
            return, 0
        endif
    endif
    
    ;; Manufacturer
    tmp = odicom->getValue('0008'x, '0070'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
        manufacturer = strtrim(*tmp[0],2)
    endif else begin
        manufacturer = 'UNKNOWN (empty)'
    endelse
    
    ; Center frequency
    temp = oDICOM->GetValue('0014'x,'401B'x,/NO_COPY)
    if temp eq -1 then begin
      temp = oDICOM->GetValue('0018'x,'0087'x,/NO_COPY)
      sf = *temp[0]*42.576
    endif else sf = *temp[0]
    
    str_name = 'MAS_DICOM_SINGLE_FILE'
    case manufacturer of 
    
        'Philips Medical Systems': str_name += '_PHILIPS'
        'SIEMENS':                 str_name += '_SIEMENS'
        'Bruker BioSpin MRI GmbH': str_name += '_BRUKER'
        'GE MEDICAL SYSTEMS':      str_name += '_GE'
        'MAS DICOM EXPORTER':      str_name += '_MAS'
        else: str_name += '_MAS'
        
    endcase
    
    dcm_single = create_struct(name=str_name)

    dcm_single.spect_freq = sf
    dcm_single.file_location = file
    dcm_single.parent_dir = file_dirname(file)
    dcm_single.manufacturer = manufacturer

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Required "Generic" DICOM tags                                                  ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    ;; Image Type
    tmp = odicom->getValue('0008'x, '0008'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       dcm_single.image_type = strtrim((*tmp[0]),2)
    endif
    
    ;; Study Date
    tmp = odicom->getValue('0008'x, '0020'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       dcm_single.study_date = long(*tmp[0])
    endif

    ;; Study Time
    tmp = odicom->getValue('0008'x, '0030'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       dcm_single.study_time = long(*tmp[0])
    endif

    ;; Image content time
    ;;tmp = odicom->getValue('0008'x, '0033'x, /no_copy)
    ;;if (mas_dicom_reader_checkval(tmp)) then begin
    ;;    print, *tmp[0]
    ;;endif
    
    ;; Protocol Name
    tmp = odicom->getValue('0018'x, '1030'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       dcm_single.protocol_name = *tmp[0]
    endif

    ;; Sequence Name
    tmp = odicom->getValue('0018'x, '0024'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       dcm_single.sequence_name = *tmp[0]
    endif
    
    ;; Patient Position (HFS/HFP/FFS/FFP/etc.)
    tmp = odicom->getValue('0018'x, '5100'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       pat_pos = strcompress(*tmp[0], /remove_all)
    endif else begin
       pat_pos = 'UNK'
    endelse
    dcm_single.patient_pos = pat_pos

    ;; pixel spacing (row,col)
    tmp = odicom->getValue('0028'x, '0030'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       vox = strsplit(*tmp[0], '\', /extract)
       dcm_single.pixel_spacing_x = float(vox[1])
       dcm_single.pixel_spacing_y = float(vox[0])
    endif

    ;; Image Number
    tmp = odicom->getValue('0020'x, '0013'x, /no_copy) 
    if (mas_dicom_reader_checkval(tmp)) then begin
        dcm_single.image_number = *tmp[n_elements(tmp)-1]
    endif

    ;; Acquisition Number
    tmp = odicom->getValue('0020'x, '0012'x, /no_copy) 
    if (mas_dicom_reader_checkval(tmp)) then begin
        dcm_single.acquisition_number = *tmp[0]
    endif

    ;; Series Number
    tmp = odicom->getValue('0020'x, '0011'x, /no_copy) 
    if (mas_dicom_reader_checkval(tmp)) then begin
        dcm_single.series_number = *tmp[0]
    endif

    ;; image position
    tmp = odicom->getValue('0020'x, '0032'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
        ;dcm_single.image_position = float(*tmp[0])
        dcm_single.image_position = strtrim(*tmp[0], 2)
    endif
    
    ;; image orientation
    tmp = odicom->getValue('0020'x, '0037'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
        cosines = strsplit(*tmp[0], '\', /extract)
        col1 = double(cosines[0:2]) & col2 = double(cosines[3:5])
        dcm_single.pat_orientation = strtrim(*tmp[0],2) 
    endif

    ;; slice thickness
    tmp = odicom->getValue('0018'x, '0050'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
        dcm_single.slice_thickness = float(*tmp[0])
    endif

    ;; slice spacing
    tmp = odicom->getValue('0018'x, '0088'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
        dcm_single.slice_spacing = float(*tmp[0])
    endif

    ;; slice location
    tmp = odicom->getValue('0020'x, '1041'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
        dcm_single.slice_location = float(*tmp[0])
    endif

    ;; samples per pixel => indicates color. should be 1 for b/w
    tmp = odicom->getValue('0028'x, '0002'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
        if (*tmp[0] ne 1) then begin
            print, "Detected multichannel imagery. Expect strange results."
        endif
    endif

    ;; number of pixel rows
    tmp = odicom->getValue('0028'x, '0010'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       dcm_single.y_dim = fix(*tmp[0])
    endif

    ;; number of pixel columns
    tmp = odicom->getValue('0028'x, '0011'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       dcm_single.x_dim = fix(*tmp[0])
    endif

    ;; smallest pixel value
    ;tmp = odicom->getValue('0028'x, '0106'x, /no_copy)
    ;if (mas_dicom_reader_checkval(tmp)) then begin
    ;endif

    ;; largest pixel value
    ;tmp = odicom->getValue('0028'x, '0107'x, /no_copy)
    ;if (mas_dicom_reader_checkval(tmp)) then begin
    ;endif

    ;; TR
    tmp = odicom->getValue('0018'x, '0080'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       dcm_single.TR = float(*tmp[0])
    endif

    ;; TE
    tmp = odicom->getValue('0018'x, '0081'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       dcm_single.TE = float(*tmp[0])
    endif

    ;; TI
    tmp = odicom->getValue('0018'x, '0082'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       dcm_single.TI = float(*tmp[0])
    endif

    ;; Number of averages
    tmp = odicom->getValue('0018'x, '0083'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       dcm_single.num_avgs = fix(*tmp[0])
    endif

    ;; Number of echos
    tmp = odicom->getValue('0018'x, '0086'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       dcm_single.num_echos = fix(*tmp[0])
    endif

    ;; Rescale slope
    tmp = odicom->getValue('0028'x, '1053'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       dcm_single.rescale_slope = float(*tmp[0])
    endif else begin
       dcm_single.rescale_slope = 1.0
    endelse

    ;; Rescale intercept
    tmp = odicom->getValue('0028'x, '1052'x, /no_copy)
    if (mas_dicom_reader_checkval(tmp)) then begin
       dcm_single.rescale_intercept = float(*tmp[0])
    endif else begin
       dcm_single.rescale_intercept = 0.0
    endelse

    ;; SOP instance id
    ;;tmp = odicom->getValue('0008'x, '0018'x, /no_copy)
    ;;if (mas_dicom_reader_checkval(tmp)) then begin
    ;;   dcm_single.sop_instance_id = strtrim(*tmp[0], 2)
    ;;endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Proprietary PHILIPS DICOM tags                                                 ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    if (str_name eq 'MAS_DICOM_SINGLE_FILE_PHILIPS') then begin
       ;; Slice Number    2001    100A
       tmp = odicom->getValue('2001'x, '100A'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
          dcm_single.slice_number = long(*tmp[0])
       endif
       
       ;; Diffusion Volume Index   2005    1413
       ;; This is computed from these two undocumented tags. They
       ;; embody the diffusion b-value number (1 ... n) and the
       ;; gradient direction for that b-value.
       tmp0 = odicom->getValue('2005'x, '1413'x, /no_copy)
       tmp1 = odicom->getValue('2005'x, '1412'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
          dcm_single.volume_index = long(string(*tmp0[0],*tmp1[0], format='(%"%d%d")'))
       endif

       ;; scan type (?)         0008    9209 
       tmp = odicom->getValue('0008'x, '9209'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
          dcm_single.scan_type = strtrim(*tmp[0],2)
       endif
       
       ;; bval                  2001    1003
       tmp = odicom->getValue('2001'x, '1003'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
          dcm_single.b_value = float((*tmp[0])[0])
       endif
       
       ;; # Slices              2001    1018
       
       ;; RLcos                 2005    10B0
       tmp = odicom->getValue('2005'x, '10B0'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
          dcm_single.bvec_X = float((*tmp[0])[0])
       endif
       
       ;; APcos                 2005    10B1
       tmp = odicom->getValue('2005'x, '10B1'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
          dcm_single.bvec_Y = float((*tmp[0])[0])
       endif
       
       ;; FHcos                 2005    10B2
       tmp = odicom->getValue('2005'x, '10B2'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
          dcm_single.bvec_Z = float((*tmp[0])[0])
       endif
       
       ;; scale slope           2005    100E
       
       ;; orientation           2001    100B (tra/cor/sag)
       tmp = odicom->getValue('2001'x, '100B'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
          dcm_single.orient_string = strtrim(*tmp[0],2)
       endif

    endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Proprietary SIEMENS DICOM tags                                                 ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    if (str_name eq 'MAS_DICOM_SINGLE_FILE_SIEMENS') then begin

       ;; Check for explicit entries for Diffusion info, otherwise resort to
       ;; trying to parse the CSA shadow header
        
       ;; bval   
       tmp = odicom->getValue('0019'x, '100C'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
          dcm_single.b_value = float(*tmp[0])
          
          tmp = odicom->getValue('0019'x,'100E'x, /no_copy)
          if (mas_dicom_reader_checkval(tmp)) then begin
              dcm_single.bvec_X = float((*tmp[0])[0])
              dcm_single.bvec_Y = float((*tmp[0])[1])
              dcm_single.bvec_Z = float((*tmp[0])[2])
          endif
          
       endif else begin
       
           ;; CSA Shadow Header will be converted to a single string, then searched. 
           ;;'(DiffusionGradientDirection).FD..M..M.(-?[0-9\.]+)..M.(-?[0-9\.]+)..M.(-?[0-9\.]+)'
           ;;'(B_value)(.{57})IS..M..M.([0-9]+)'
       
           tmp = odicom->getValue('0029'x, '1010'x, /no_copy)
           if (mas_dicom_reader_checkval(tmp)) then begin
           
               tmp_str = strarr(n_elements(*tmp[0]))
               for i = 0, n_elements(*tmp[0])-1 do begin
                   tmp_str[i] = string((*tmp[0])[i])
               endfor
               tmp_str = strjoin(tmp_str, '')

               reg = stregex(tmp_str, '(DiffusionGradientDirection)(.{1,64})FD..M..M.(-?[0-9\.]+)..M.(-?[0-9\.]+)..M.(-?[0-9\.]+)', $
                              /extract, /subexpr)
               dcm_single.bvec_x = float(reg[3])
               dcm_single.bvec_y = float(reg[4])
               dcm_single.bvec_z = float(reg[5])
               
               reg = stregex(tmp_str, '(B_value)(.{1,60})IS..M..M.([0-9]+)', $
                              /extract, /subexpr)
               dcm_single.b_value = float(reg[3])

; Not Used
;               reg = stregex(tmp_str, '(ProtocolSliceNumber)(.{44}).IS..M..M.([0-9]+)', $
;                              /extract, /subexpr)               
;               dcm_single.slice_number = long(reg[3])
;               print, long(reg[3])
           endif
           
       endelse
       
        ;; Orientation String 
       tmp = odicom->getValue('0051'x, '100E'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
;          dcm_single.orient_string = strtrim(*tmp[0], 2)
       endif
 
       ;;dcm_single.volume_index = dcm_single.acquisition_number
       
   endif
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Proprietary GE DICOM tags                                                      ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    if (str_name eq 'MAS_DICOM_SINGLE_FILE_GE') then begin

       ;; bval   
       tmp = odicom->getValue('0043'x, '1039'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
          temp = strsplit(*tmp[0], '\', /extract)
          dcm_single.b_value = float(temp[0])
       endif
;       
;       ;; # Slices              2001    1018
;       
       ;; RLcos  
       ;tmp = odicom->getValue('0021'x,'105A'x, /no_copy)
       ;if (mas_dicom_reader_checkval(tmp)) then begin
       ;   dir = float(*tmp[0])
       ;   dcm_single.bvec_X = float((*tmp[0])[0])
       ;   dcm_single.bvec_Y = float((*tmp[0])[1])
       ;   dcm_single.bvec_Z = float((*tmp[0])[2])
       ;endif


       tmp = odicom->getValue('0027'x, '1040'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
          dcm_single.orient_string = strtrim(string(*tmp[0]), 2)
       endif
 
       dcm_single.volume_index = dcm_single.acquisition_number
       
   endif
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Proprietary BRUKER DICOM tags                                                  ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Proprietary MAS DICOM tags                                                     ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    if (str_name eq 'MAS_DICOM_SINGLE_FILE_MAS') then begin
       ;; bval                  2001    1003
       tmp = odicom->getValue('2001'x, '1003'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
          dcm_single.b_value = float(fix(*tmp[0], type=7))
       endif
       
       ;; # Slices              2001    1018
       
       ;; RLcos                 2005    10B0
       tmp = odicom->getValue('2005'x, '10B0'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
          dcm_single.bvec_X = float(fix(*tmp[0], type=7))
       endif
       
       ;; APcos                 2005    10B1
       tmp = odicom->getValue('2005'x, '10B1'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
          dcm_single.bvec_Y = float(fix(*tmp[0], type=7))
       endif
       
       ;; FHcos                 2005    10B2
       tmp = odicom->getValue('2005'x, '10B2'x, /no_copy)
       if (mas_dicom_reader_checkval(tmp)) then begin
          dcm_single.bvec_Z = float(fix(*tmp[0], type=7))
       endif

    endif



    if (destroy_on_return eq 1) then begin
        obj_destroy, odicom
    endif
    
;    grad_vec = [dcm_single.bvec_Y, dcm_single.bvec_X, dcm_single.bvec_Z]
;    if (total(abs(grad_vec)) eq 0 and dcm_single.b_value gt 0) then begin
;        print, "Skipping: "+file+" (bvec = [0,0,0] & bval > 0)"
;        success = 0
;        return, success
;    endif
    
    success = 1
    
    dcm_struct = temporary(dcm_single)
    
    return, success

end

;; $Id$
;;
;; Subroutine name: mas_dicom_choose_fields_event
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Event handler for the add/remove acquisitions buttons.
;;
;; Editing Information:

pro mas_dicom_choose_fields_event, ev

    widget_control, ev.top, get_uvalue=state
    
    name = widget_info(ev.id, /uname)
    
    if (ptr_valid((*state).sel_active)) then begin
        sel_active = *((*state).sel_active)
        n_active = n_elements(sel_active)
    endif else begin
        n_active = 0
    endelse
    
    if (ptr_valid((*state).sel_inactive)) then begin
        sel_inactive = *((*state).sel_inactive)
        n_inactive = n_elements(sel_inactive)
    endif else begin
        n_inactive = 0
    endelse
    
    case name of
    
        'btn_add': begin
        
            sel_index = widget_info((*state).inact_list, /list_select)
            
            if (n_inactive eq 0) then return
            
            if (n_active eq 0) then begin
                sel_active = [ sel_inactive[sel_index] ]
            endif else begin
                sel_active = [ sel_active, sel_inactive[sel_index] ]
            endelse
            n_active++
            
            if (n_inactive eq 1) then begin
                ptr_free, (*state).sel_inactive
            endif else begin
                sel_inactive_new = bytarr(n_inactive-1)
                ind = 0
                for i = 0, n_inactive-1 do begin
                    if (i eq sel_index) then continue
                    sel_inactive_new[ind++] = sel_inactive[i]
                endfor
                sel_inactive = temporary(sel_inactive_new)
            endelse
            n_inactive--
            
        end
        
        'btn_remove': begin

            sel_index = widget_info((*state).act_list, /list_select)
            
            if (n_active eq 0) then return
            
            if (n_inactive eq 0) then begin
                sel_inactive = [ sel_active[sel_index] ]
            endif else begin
                sel_inactive = [ sel_inactive, sel_active[sel_index] ]
            endelse
            n_inactive--
            
            if (n_active eq 1) then begin
                ptr_free, (*state).sel_active
            endif else begin
                sel_active_new = bytarr(n_active-1)
                ind = 0
                for i = 0, n_active-1 do begin
                    if (i eq sel_index) then continue
                    sel_active_new[ind++] = sel_active[i]
                endfor
                sel_active = temporary(sel_active_new)
            endelse
            n_active--

        end
        
        'act_list': begin
            widget_control, (*state).inact_list, set_list_select=-1
            widget_control, (*state).act_list, set_list_select=ev.index
            return
        end
        
        'inact_list': begin
            widget_control, (*state).act_list, set_list_select=-1
            widget_control, (*state).inact_list, set_list_select=ev.index
            return
        end
        
        else: begin
            help, ev, /structure
            return
        end
        
    endcase
    
    values_inactive = ''
    values_active = ''
    
    pdcm_array = (*state).dcm_array
    
    if (n_inactive gt 0) then begin
        sel_inactive = sel_inactive[sort(sel_inactive)]
        values_inactive = (*state).sel_values[sel_inactive]
        ptr_free, (*state).sel_inactive
        (*state).sel_inactive = ptr_new(sel_inactive)
    endif
    
    if (n_active gt 0) then begin
        values_active = (*state).sel_values[sel_active]
        ptr_free, (*state).sel_active
        (*state).sel_active = ptr_new(sel_active)
        *pdcm_array = mas_dicom_group_acq(*pdcm_array, sel_active, acq_indices, 0)        
    endif else begin
        sel_active = [0]
        *pdcm_array = (*pdcm_array)[sort( (*pdcm_array).(0) )]
        acq_indices = [ n_elements(*pdcm_array) - 1]
    endelse
    
    widget_control, (*state).inact_list, set_value=values_inactive
    widget_control, (*state).act_list,   set_value=values_active
    
    (*state).acq_indices = ptr_new(acq_indices)
    ptr_free, (*state).acq_skip
    (*state).acq_skip    = ptr_new(bytarr(n_elements(acq_indices))+0)
    
    widget_control, (*state).table_id, set_value=*((*state).dcm_array)
    
end

;; $Id$
;;
;; Subroutine name: mas_dicom_import_tree_event
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; event handler for the dicom import tree
;;
;; Editing Information:

pro mas_dicom_import_tree_event, event

    widget_control, event.top, get_uvalue=tree_state
    
    name = widget_info(event.id, /uname)
    
    case name of
    
        'leaf_node': begin
            widget_control, event.id, get_uvalue=dcm
            gui_state = (*tree_state).gui_state
            col_names = (*gui_state).sel_values
            max_len = max(strlen(col_names))+1
            
            if (n_elements(dcm) ne 0) then begin
                col_names = tag_names(dcm)
                for i = 0, n_elements(col_names)-1 do begin
                
                    label = (*tree_state).label_ids[i]
                    
                    val_str = strtrim(strcompress(string(dcm.(i))), 2)
                    if (strlen(val_str) gt 60) then begin
                        val_str = '...'+strmid(val_str, strlen(val_str)-60, strlen(val_str))
                    endif else if (strlen(val_str) eq 0) then begin
                        val_str = '<empty>'
                    endif
                    
                    new_value = strmid(val_str, 0, 80)
                    
                    widget_control, label, set_value=new_value
                    
                endfor
            endif
        end
        
        'btn_import_ok': begin
            (*tree_state).import_ok = 1
            widget_control, event.top, /destroy
        end

        'btn_remove_selected': begin
        
            gui_state = (*tree_state).gui_state
            top_root = widget_info(event.top, find_by_uname='top_root')
            remove_node = widget_info(top_root, /tree_select)
            
            if (widget_info(remove_node, /valid_id) eq 0) then return
            
            ;; get the acq index of the selected node
            widget_control, remove_node, get_uvalue=selected_index
            if (n_elements(selected_index) eq 0) then return
            
            ;; Set it to be skipped
            (*((*gui_state).acq_skip))[selected_index] = 1
            
            ;; Remove the node from the tree
            widget_control, widget_info(top_root, /tree_select), /destroy
            
        end
        
        'btn_import_cancel': begin
            (*tree_state).import_ok = 0
            widget_control, event.top, /destroy
        end
        
        else: print, "mas_dicom_import_tree_event: recv'd event from unknown widget: "+name
        
    endcase
    
end

;; $Id$
;;
;; Subroutine name: mas_dicom_select_event
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; event handler for the dicom image sorting table
;;
;; Editing Information:

pro mas_dicom_select_event, ev

    name = widget_info(ev.id, /uname)
    
    widget_control, ev.top, get_uvalue=state
    
    case name of 
        
        ;; called when user click the checkbox to select volumes to sort.
        'btn_ident_vol': begin
        
            dl_id = widget_info(ev.top, find_by_uname='dl_volumes')
            dl_selected = widget_info(dl_id, /droplist_select)
            
            widget_control, dl_id, sensitive=ev.select
            
            ;; constructs a widget event which will be thrown by us. this
            ;; event triggers the re-sorting of the table based on what
            ;; ever is currently selected in the volumes drop list. 
            dl_event = create_struct(name='WIDGET_DROPLIST')
            dl_event.id = dl_id
            dl_event.top = ev.top
            dl_event.handler = ev.top
                            
            if (ev.select eq 0) then begin
                dl_event.index = 0
            endif else begin
                dl_event.index = dl_selected
            endelse

            ;; throw droplist event to trigger re-sorting
            mas_dicom_select_event, dl_event
            
        end
        
        'dl_acquisitions': begin
        
            pdcm_array = (*state).dcm_array

            *pdcm_array = mas_dicom_group_acq(*pdcm_array, *((*state).sel_active), acq_indices, 0)

            (*state).acq_indices = ptr_new(acq_indices)
                        
            dl_id = widget_info(ev.top, find_by_uname='dl_volumes')
            dl_event = create_struct(name='WIDGET_DROPLIST')
            dl_event.id = dl_id
            dl_event.top = ev.top
            dl_event.handler = ev.top
            dl_event.index = (*state).vol_key
             
            mas_dicom_select_event, dl_event
        end
        
        ;; called when user selects an dicom header field to 
        ;; group acquistions.
        'dl_volumes': begin

            pdcm_array = (*state).dcm_array
            
            index = ev.index
            
            acq_indices = *((*state).acq_indices)
            
            start = 0
            
            ;; This sorts by each field selected in the acquisition
            ;; selection box.
            for a = 0, n_elements(acq_indices)-1 do begin
            
                end_ind = acq_indices[a]
                
                splice = sort( (*pdcm_array)[start:end_ind].(index) ) + start
                (*pdcm_array)[start:end_ind] = (*pdcm_array)[splice]
                
                start = end_ind+1
            
            endfor
            
            (*state).vol_key = index

            dl_id = widget_info(ev.top, find_by_uname='dl_images')
            dl_event = create_struct(name='WIDGET_DROPLIST')
            dl_event.id = dl_id
            dl_event.top = ev.top
            dl_event.handler = ev.top
            dl_event.index = (*state).img_key
            
            mas_dicom_select_event, dl_event
            ;widget_control, (*state).table_id, set_value=*((*state).dcm_array)
        
        end
        
        'dl_images': begin
        
            pdcm_array = (*state).dcm_array
            
            index = ev.index
            
            acq_indices = *((*state).acq_indices)
            
            start_acq = 0
            
            ;; Loops through and breaks file array into chunks containing
            ;; the acquisitions and volumes
            for a = 0, n_elements(acq_indices)-1 do begin
            
                end_acq = acq_indices[a]
                
                ;; for the current acquisition, pull out unique instances of the
                ;; volume search key. This will give the number of volumes as well
                ;; as the indices at which they are broken up.
                if ((*state).vol_key ne 0) then begin
                    vol_indices = uniq( (*pdcm_array)[start_acq:end_acq].((*state).vol_key) )
                    n_volumes = n_elements(vol_indices)
                endif else begin
                    vol_indices = [end_acq-start_acq] ;;acq_indices
                    n_volumes = 1
                endelse
                
                start_vol = start_acq
                
                ;; For each volume in the acquisition, loop through and sort the 
                ;; images by the user-selected image key.
                for v = 0, n_volumes-1 do begin
                    
                    end_vol = vol_indices[v] + start_acq
                    
                    splice = sort( (*pdcm_array)[start_vol:end_vol].(index) ) + start_vol
                    (*pdcm_array)[start_vol:end_vol] = (*pdcm_array)[splice]
                    
                    start_vol = end_vol+1
                    
                endfor
                
                start_acq = end_acq + 1
                    
            endfor

            ;; update the table widget with new data
            widget_control, (*state).table_id, set_value=*((*state).dcm_array)
            
            (*state).img_key = index
            
        end
        
        'btn_import': begin
        
            ;; display an overview of the what will be imported. Note that this
            ;; blocks and allows user to manipulate the importing. 
            mas_dicom_show_import_tree, state, group_leader=ev.top, import_ok=import_ok
        
            ;; User decides to import?
            if (import_ok eq 1) then begin
                mas_dicom_do_import, state
                
                ptr_free, (*state).dcm_array
                ptr_free, (*state).sel_active
                ptr_free, (*state).sel_inactive
                ptr_free, (*state).acq_indices
                ptr_free, state
    
                widget_control, ev.top, /destroy
            endif
            
        end
        
        'btn_cancel': begin
        
           ptr_free, (*state).dcm_array
           ptr_free, (*state).sel_active
           ptr_free, (*state).sel_inactive
           ptr_free, (*state).acq_indices
           ptr_free, state
            
           widget_control, ev.top, /destroy
           
        end
        
        else: begin
            print, "mas_dicom_select_event: Unknown event: "+name
            help, ev, /structure
        end
        
    endcase
    
end

;; $Id$
;;
;; Subroutine name: mas_dicom_show_import_table
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Displays the import table, after reading a dicom directory
;;
;; Editing Information:

pro mas_dicom_show_import_table, dcm_array, group_leader=group_leader

    col_names = tag_names(dcm_array[0])

    if (not keyword_set(group_leader)) then group_leader = 0

    if (widget_info(group_leader, /valid_id)) then begin
        base = widget_base(title='DICOM File List', /column, group_leader=group_leader)
    endif else begin
        base = widget_base(title='DICOM File List', /column)
    endelse

    base_2col = widget_base(base, /row)
    
    ;; section for user to select acquisistions
    base_acq = widget_base(base_2col, /column, /frame)
    lbl = widget_label(base_acq, value="Choose the sort keys to identify acquisitions:")
    base_dbl = widget_base(base_acq, /row)
    b1 = widget_base(base_dbl)
    l1 = widget_list(b1, value=col_names, ysize=10, uname='inact_list', event_pro='mas_dicom_choose_fields_event')
    b2 = widget_base(base_dbl, /column)
    btn = widget_button(b2, value="Add >>", uname="btn_add", event_pro='mas_dicom_choose_fields_event')
    btn = widget_button(b2, value="<< Remove", uname='btn_remove', event_pro='mas_dicom_choose_fields_event')
    b3 = widget_base(base_dbl)
    l2 = widget_list(b3, value=col_names, ysize=10, uname='act_list', event_pro='mas_dicom_choose_fields_event')

    base_controls = widget_base(base_2col, /column, /frame)

    ;; section for selecting volumes
    b = widget_base(base_controls, /row)
    bex = widget_base(b, /nonexclusive)
    btn = widget_button(bex,  value="Choose a key to identify volumes:", uname='btn_ident_vol')
    dl = widget_droplist(b, value=col_names, uname='dl_volumes', sensitive=0)

    ; section for selecting images
    b = widget_base(base_controls, /row)
    bex = widget_base(b);, /nonexclusive)
    btn = widget_label(bex, value="Choose a key to order images:", uname='btn_ident_img')
    dl = widget_droplist(b, value=col_names, uname='dl_images')
    
    ;; cancel buttons
    base_btns = widget_base(base_controls, /row)
    btn_imp = widget_button(base_btns, value="Next >>", uname='btn_import')
    btn_can = widget_button(base_btns, value="Cancel", uname='btn_cancel')
    
    ;; create the table to display the files and their properties
    base_table = widget_base(base)
    table = widget_table(base_table, /scroll, /resizeable_columns, $
                         uname='tbl_images', $
                         value=dcm_array, $
                         $;;/all_events, /disjoint_selection, $
                         column_labels=col_names, $
                         scr_ysize=400, $
                         scr_xsize=900)

    ;; create some widget state information 
    state = { inact_list: l1, $
              act_list: l2, $
              use_vol_key: 0, $
              use_img_key: 1, $
              sel_values: col_names, $
              sel_active: ptr_new(), $
              sel_inactive: ptr_new(indgen(n_elements(col_names))), $
              sel_map: bytarr(n_elements(col_names)), $
              table_id: table, $
              acq_key: 0, $
              acq_indices: ptr_new([ n_elements(dcm_array)-1 ]), $
              acq_skip: ptr_new([0]), $
              vol_key: 0, $
              img_key:0, $
              dcm_array: ptr_new(dcm_array, /no_copy) }
    
    ;; enter event loop
    widget_control, base, set_uvalue=ptr_new(state, /no_copy), /realize
    widget_control, l2, set_value=''
    xmanager, 'mas_dicom_select', base, /no_block
    ;;;;;;;;


end

;; $Id$
;;
;; Subroutine name: mas_dicom_show_import_tree
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Displays an import review window, with a tree representing the
;; sorting that the user requested.
;;
;; Editing Information:

pro mas_dicom_show_import_tree, state, group_leader=group_leader, import_ok=import_ok

    pdcm_array  = (*state).dcm_array
    col_names   = tag_names(*pdcm_array)
    acq_indices = (*state).acq_indices
    sel_active  = (*state).sel_active
    vol_key     = (*state).vol_key
    img_key     = (*state).img_key
    import_ok   = 0
    
    ;;if not ptr_valid(sel_active) then return
    if (not keyword_set(group_leader)) then group_leader = 0
    
    if (widget_info(group_leader, /valid_id)) then begin
        base = widget_base(title='DICOM Import Tree', /column, group_leader=group_leader, /modal)
    endif else begin
        base = widget_base(title='DICOM Import Tree', /column)
    endelse
    
    top_base = widget_base(base, /row, /frame)
    bot_base = widget_base(base, /row)
    
    lbl = widget_label(top_base, value="Please confirm the import results")
    btn_remove = widget_button(top_base, value='Remove selected acquisition', /align_right, uname='btn_remove_selected')
    btn_ok = widget_button(top_base, value="Start Importing", /align_right, uname="btn_import_ok")
    btn_cancel = widget_button(top_base, value="Cancel", /align_right, uname="btn_import_cancel")
    
    tree = widget_tree(bot_base, xsize=300, uname='top_root') ;; , /context_events)
    info_base = widget_base(bot_base, /row, /frame)
    
    info_lbls = widget_base(info_base, /column, xpad=0, ypad=0, /base_align_top)
    info_txts = widget_base(info_base, /column, xpad=0, ypad=0, /base_align_top)
    
    max_len = max(strlen(col_names))+1
    labels = lonarr(n_elements(col_names))
    
    for i = 0, n_elements(col_names)-1 do begin
    
        col_name = col_names[i] + ': '

        lbl       = widget_label(info_lbls, value=col_name, /align_right, /dynamic_resize)
        labels[i] = widget_label(info_txts, value=strjoin(replicate(' ',20)), /align_left, /dynamic_resize)
        
    endfor

    progressBar = obj_new('progressbar', title="Operation in Progress", text="Building Import Tree", /fast_loop)
    progressBar->Start
    
    start_acq = 0
    item_ct = 0L
    
    for a = 0, n_elements(*acq_indices)-1 do begin
    
        node_title = ''
        acq = (*acq_indices)[a]
        
        if (ptr_valid(sel_active)) then begin
            for s = 0, n_elements(*sel_active)-1 do begin
                tag = (*sel_active)[s]
                node_title += col_names[tag] + ': ' 
                node_title += strcompress(string((*pdcm_array)[acq].(tag)))
                if (s ne n_elements(*sel_active)-1) then node_title += ', '
            endfor
        endif else begin
            node_title = 'ACQUISITION: Single'
        endelse
        
        root = widget_tree(tree, value=node_title, /folder, uname="root_node", uvalue=a)
        
        if (vol_key eq 0) then begin
            volume_indices = [acq-start_acq]
            n_volumes = 1
        endif else begin
            volume_indices = uniq( (*pdcm_array)[start_acq:acq].(vol_key) )
            n_volumes = n_elements(volume_indices)
;            sl_per_vol = volume_indices[*] - [-1, volume_indices[0:n_volumes-1]]
;            slice_layout = float(sl_per_vol) / max(sl_per_vol)
;            incomplete_volumes = where(abs(1.-slice_layout) gt 1e-4, n_incomplete)
;            print, slice_layout
        endelse
        
        start_vol = start_acq
        
        for v = 0, n_volumes-1 do begin
        
            vol = volume_indices[v] + start_acq
            
            if (vol_key eq 0) then begin
                node_title = 'VOLUME: Single'
            endif else begin
                node_title  = 'VOLUME: ' + string(v, format='(I0)') 
                node_title += ' (' + col_names[vol_key] + ': ' 
                node_title += strtrim(strcompress(string((*pdcm_array)[vol].(vol_key))),2) + ')'
            endelse
            
            branch = widget_tree(root, value=node_title, /folder, uname="branch_node")
            
            for im = start_vol, vol do begin

                node_title = col_names[img_key] + ': ' + strcompress(string((*pdcm_array)[im].(img_key)))
                leaf = widget_tree(branch, value=node_title, uvalue=(*pdcm_array)[im], uname="leaf_node")
                progressBar->Update, float(++item_ct)/float(n_elements(*pdcm_array)) * 100.0
                if (item_ct mod 10 eq 0 && progressBar->CheckCancel()) then begin
                    progressBar->Destroy
                    return
                endif
            endfor
        
            start_vol = vol + 1
        
        endfor
        
        start_acq = acq + 1
        
    endfor
    
    progressBar->Destroy
    
    tree_state = ptr_new({label_ids: labels, gui_state: state, import_ok: 0}, /no_copy)
    
    import_ok = 0
    
    widget_control, base, set_uvalue=tree_state, /realize
    xmanager, 'mas_dicom_import_tree', base
    
    import_ok = (*tree_state).import_ok
    
    ptr_free, tree_state

end

;; $Id$
;;
;; Subroutine name: mas_dicom_do_import
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Performs the import into MAS of the selected acquisitions.
;;
;; Editing Information:

pro mas_dicom_do_import, state

    common scan_data
    forward_function mnf_mat2qt
    
    pdcm_array = (*state).dcm_array
    
    acq_key = (*state).acq_key
    vol_key = (*state).vol_key
    img_key = (*state).img_key
    
    acq_skip    = *((*state).acq_skip)
    acq_indices = *((*state).acq_indices)
    
    start_acq = 0
    
    for a = 0, n_elements(acq_indices)-1 do begin

        end_acq = acq_indices[a]

        ;; this will be true if one of the acquisitions was marked
        ;; for removal in the import tree window by the user.
        if ( acq_skip[a] ne 0) then begin
            start_acq = end_acq + 1
            continue
        endif
        
        ni   = project.ni
        imnd = project.imndarray[ni]
        
        if (vol_key eq 0) then begin
            volume_indices = acq_indices
            n_volumes = 1
            sdim = (end_acq - start_acq)/n_volumes + 1
            imnd.DICOMslicelayout = ptr_new([sdim])
        endif else begin
            volume_indices = uniq( (*pdcm_array)[start_acq:end_acq].(vol_key) )
            n_volumes = n_elements(volume_indices)
            sl_per_vol = volume_indices[*] - [-1, volume_indices[0:n_volumes-1]]
            sdim = max(sl_per_vol)
            slice_layout = sl_per_vol
            imnd.DICOMslicelayout = ptr_new(slice_layout, /no_copy)
            volume_indices = volume_indices +  start_acq
        endelse
        
        imnd.file_Path = file_dirname( (*pdcm_array)[start_acq].file_location )
        imnd.display_name = imnd.file_path + ': ' + (*pdcm_array)[start_acq].protocol_name
        imnd.sdim = sdim
        imnd.slices = imnd.sdim
        imnd.spect_bf1 = ((*pdcm_array[0]).spect_freq)[0]
        imnd.adim = n_volumes
        imnd.multi_scan_file_array = ptr_new( (*pdcm_array)[start_acq:end_acq].file_location )
        imnd.DICOMfiles = imnd.multi_scan_file_Array
        imnd.DICOMtype = 0
        imnd.image_type = 18
        imnd.dimensions = 2
        
        ;; Study Date
        datestr = stregex(strtrim((*pdcm_array)[start_acq].study_date, 2), $
            '^([0-9]{4})([0-9]{2})([0-9]{2})$', /extract, /subexpr)
        months = ['???','Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', $
            'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        mon = fix(datestr[2]) > 0
        if (mon gt 12) then mon = 0
        imnd.scan_date = months[mon]+' '+datestr[3]+' '+datestr[1]
        
        ;; Study Time
        timestr = strtrim((*pdcm_array)[start_acq].study_time,2)
        timelen = strlen(timestr)
        if (timelen lt 6) then begin
            timestr = strjoin(replicate('0', 6-timelen))+timestr
        endif
        timestr = stregex(timestr, '^([0-9]{2})([0-9]{2})([0-9]{2}).*$',$
                           /extract, /subexpr)
            
        imnd.scan_date += ' '+timestr[1]+':'+timestr[2]+':'+timestr[3]
        
        ;; Protocol Name
        imnd.scan_name = (*pdcm_array)[start_acq].protocol_name
        
        ;; pixel spacing (x,y)
        imnd.f_voxsz = (*pdcm_array)[start_acq].pixel_spacing_x / 10.
        imnd.p_voxsz = (*pdcm_array)[start_acq].pixel_spacing_y / 10. 
        
        ;; image position
        qoffset_xyz = float(strsplit((*pdcm_array)[start_acq].image_position, '\', /extract))
        qoffset_xyz[0] = -qoffset_xyz[0]
        
        ;; image orientation
        cosines = strsplit((*pdcm_array)[start_acq].pat_orientation, '\', /extract)
        col1 = double(cosines[0:2]) & col2 = double(cosines[3:5])
        col3 = crossp(col1, col2)
        col_mult = float([1,1,1])
        
        orient = [ [ 1, 0, 0 ], $ ;; TRANSVERSE
                   [ 0, 1, 0 ], $
                   [ 1, 0, 0 ], $ ;; CORONAL
                   [ 0, 0,-1 ], $
                   [ 0, 1, 0 ], $ ;; SAGITTAL
                   [ 0, 0,-1 ] ]
            
        for j = 0,5,2 do begin
            ref_a = orient[*,j] & ref_b = orient[*,j+1]
            tst_a = total(abs(ref_a - round(col1)))
            tst_b = total(abs(ref_b - round(col2)))
            if (tst_a lt 1e-3 and tst_b lt 1e-3) then begin
                case j of
                    0: begin
                        imnd.orientation[0] = 'TRANSVERSE'
                        imnd.orientation[1] = 'F'
                        imnd.orientation[2] = 'H'
                    end
                    2: begin
                        imnd.orientation[0] = 'CORONAL'
                        imnd.orientation[1] = 'A'
                        imnd.orientation[2] = 'P'
                    end
                    4: begin
                        imnd.orientation[0] = 'SAGITTAL'
                        imnd.orientation[1] = 'L'
                        imnd.orientation[2] = 'R'
                        col_mult[0:1] = -1. ;; this is necessary
                    end
                endcase
                break
            endif
        endfor
        
        ;; Patient Position (HFS/HFP/FFS/FFP/etc.)
        pat_pos = (*pdcm_array)[start_acq].patient_pos
        
        imnd.orientation[0] +=  ' ('+pat_pos+')'
        mat = transpose([ [col_mult[0] * col1], $
                          [col_mult[1] * col2], $
                          [col_mult[2] * col3] ])
        
        axes_tx = float(round(mat))
        axes_tx[*,0] = -axes_tx[*,0] ;; flips right & left to accomodate our orientation scheme.
        print, axes_tx
        ijk = mas_orient_from_matrix(axes_tx)
        
        project.procPramArray[ni].orientation_Imin = ijk[0]
        project.procPramArray[ni].orientation_Jmin = ijk[1]
        project.procPramArray[ni].orientation_Kmin = ijk[2]
        project.procpramarray[ni].have_orientation = 1
        
        ;; slice thickness
        imnd.thick = (*pdcm_array)[start_acq].slice_thickness; / 10.
        imnd.s_voxsz = imnd.thick / 10.
        imnd.s_fov = (imnd.thick * imnd.sdim); * 10.
        
        ;; slice spacing
        if ((*pdcm_array)[start_acq].slice_spacing ne 0) then begin
            imnd.thick = (*pdcm_array)[start_acq].slice_spacing; / 10.
            imnd.s_voxsz = imnd.thick / 10.
            imnd.s_fov = (imnd.thick * imnd.sdim); * 10.
        endif
        
        ;; number of pixel rows
        imnd.fdim = (*pdcm_array)[start_acq].x_dim
        imnd.f_fov = (imnd.fdim * imnd.f_voxsz); / 10.
        
        ;; number of pixel columns
        imnd.pdim = (*pdcm_array)[start_acq].y_dim
        imnd.p_fov = (imnd.pdim * imnd.p_voxsz); / 10.

        ;; sform
        sform_mat= fltarr(4,4)
        sform_mat[0:2,0:2] = axes_tx # diag_matrix([imnd.f_voxsz, imnd.p_voxsz, imnd.s_voxsz]*10.0)
        sform_mat[3,0:2] = qoffset_xyz
        sform_mat[3,3] = 1
        imnd.qoffset_xyz = ptr_new(qoffset_xyz)
        imnd.sform_matrix = ptr_new(sform_mat)
        imnd.qfac = 1.0
        tmp = axes_tx
        if (determ(tmp) ne 1.0) then begin
            tmp[2,*] *= -1.0
            imnd.qfac = -1.0
        endif
        imnd.qform_bcd = ptr_new((mnf_mat2qt(tmp))[1:3])
            
        ;; TR
        imnd.recov_time = (*pdcm_array)[start_acq].TR
        if (n_volumes gt 1 and ((*state).sel_values)[vol_key] eq 'TR') then begin
            TR_array_ind = uniq((*pdcm_array)[start_acq:end_acq].TR)
            if (n_elements(TR_array_ind) eq n_volumes) then begin
                tr_array = ( (*pdcm_array)[start_acq:end_acq].TR )[TR_array_ind]
                imnd.rep_time_ptr = ptr_new(tr_array, /no_copy)
            endif
        endif
        
        ;; TE
        imnd.echo_time = (*pdcm_array)[start_acq].TE
        if (n_volumes gt 1 and ((*state).sel_values)[vol_key] eq 'TE') then begin
            TE_array_ind = uniq((*pdcm_array)[start_acq:end_acq].TE)
            if (n_elements(TE_array_ind) eq n_volumes) then begin
                te_array = ( (*pdcm_array)[start_acq:end_acq].TE )[TE_array_ind]
                imnd.echo_time_ptr = ptr_new(te_array, /no_copy)
            endif
        endif

        ;; TI
        if (n_volumes gt 1 and ((*state).sel_values)[vol_key] eq 'TI') then begin
            TI_array_ind = uniq((*pdcm_array)[start_acq:end_acq].TI)
            if (n_elements(TI_array_ind) eq n_volumes) then begin
                ti_array = ( (*pdcm_array)[start_acq:end_acq].TI )[TI_array_ind]
                imnd.inversion_time_ptr = ptr_new(ti_array, /no_copy)
            endif
        endif
        
        ;; Number of averages
        imnd.n_avg = (*pdcm_array)[start_acq].num_avgs
        imnd.pn_avg = ptr_new(imnd.n_avg)
        
        ;; Number of echos
        imnd.n_echo = (*pdcm_array)[start_acq].num_echos
        
        if (n_volumes gt 1) then begin
        
            bval_array  = fltarr(n_volumes)
            b_matrix    = fltarr(6,n_volumes)
            angle_theta = fltarr(n_volumes)
            angle_phi   = fltarr(n_volumes)
            
            ;; now we iterate through the volumes and extract
            ;; per-volume information, such as diffusion weighting
            ;; and variable parameters
                
            for vol = 0, n_volumes-1 do begin
            
                vol_index = volume_indices[vol]
                vol_rep   = (*pdcm_array)[vol_index]
            
                bval_array[vol] = vol_rep.b_value
                gr_x = vol_rep.bvec_X
                gr_y = vol_rep.bvec_Y
                gr_z = vol_rep.bvec_Z
                
                gr   = [gr_x, gr_y, gr_z]
                gr_x =  total(gr * col1)
                gr_y = -total(gr * col2)
                gr_z =  total(gr * col3)
                
                sph_grad = (cv_coord(from_rect=[gr_x, gr_y, gr_z], /to_sphere, /degrees))
                angle_theta[vol] = 90. - sph_grad[1]
                angle_phi[vol] = sph_grad[0]
                
                n_neg = where(angle_phi lt 0, n_n_neg)
                if (n_n_neg gt 0) then begin
                    angle_phi[n_neg] += 360.
                endif
                
                b_matrix[*,vol] = vol_rep.b_value * [ gr_x^2, gr_y^2, gr_z^2, 2*gr_x*gr_y, 2*gr_x*gr_z, 2*gr_y*gr_z ]
            
            endfor
            
            imnd.n_bvals     = n_elements(bval_array)
            imnd.bval_array  = ptr_new(bval_array, /no_copy)
            imnd.angle_theta = ptr_new(angle_theta, /no_copy)
            imnd.angle_phi   = ptr_new(angle_phi, /no_copy)
            imnd.b_matrix    = ptr_new(b_matrix, /no_copy)
            imnd.big_delta   = ptr_new(replicate(29.73, n_volumes), /no_copy)
            imnd.small_delta = ptr_new(replicate(16.91, n_volumes), /no_copy)
            
        endif
        
        imnd.state1_load_procedure = 'mas_dicom_load_state_1'
        
        project.current_path = imnd.file_Path
        project.procpramarray[ni].signal_type = 8
        project.scan_open_flag = 1
        project.imndarray[ni] = temporary(imnd)
        project.ci = project.ni
        project.ni++
        mas_add_current_scan_to_scan_selection
        mas_redraw_GUI
        
        start_acq = end_acq + 1
        
    endfor
    
end

;+
; :Description:
;    Given an array of dicom files and a "recipe" for how to sort them, this
;    procedure will read and sort them accordingly, and load them into the 
;    MAS window..
;
; :Params:
;    files  - a list of dicom files to sort
;    recipe - a recipe for sorting (see mas_dicom_import_recipe__define)
; 
; :Keywords:
;    progressbar - if set, a progressbar gui window will show file
;                  reading progress,
; :Author: wtriplett
;-
pro mas_dicom_import_by_recipe, files, recipe, progressbar=progressbar

    if (n_elements(files) eq 0) then return
    
    n_files = n_elements(files)
    
     mas_dicom_get_dicom_obj, dicom_obj=dicom_obj
    
    data = mas_dicom_read_single_file(files[0], reuse_object=dicom_obj, dcm_struct=dcm_tmp)
    
    dcm = replicate(dcm_tmp, n_files)
  
    n=0
    message, /info, 'Reading DICOM files...'
    if (keyword_set(progressbar)) then begin
        pbar = obj_new('PROGRESSBAR', title='Loading DICOMs.', text='Loading Files...')
        pbar->Start
    endif
    
    for i = 0, n_files-1 do begin 
       
       retval = mas_dicom_read_single_file(files[i], reuse_object=dicom_obj, dcm_struct=dcm_tmp)
       if (retval eq 0) then begin
           message, /info, "Skipping: "+files[i]+" (retval = 0)"
           continue
       endif
       
       dcm[n++] = dcm_tmp
       if (i mod 50 eq 0) then begin
           message, /info, 'Percent Complete: '+string(float(i)/n_files*100.0, format='(%"%0.1f")')
           if (obj_valid(pbar)) then pbar->update, float(i)/n_files*100.0
       endif
       
    endfor
  
    if (obj_valid(pbar)) then begin
        pbar->Destroy
    endif
    
    dcm = dcm[0:n-1]
  
    acq_tag  = strupcase(recipe.acq_tag)
    vol_tag  = strupcase(recipe.vol_tag)
    img_tag  = strupcase(recipe.img_tag)
    want_acq = strupcase(recipe.want_acq)
    
    if (acq_tag ne '') then begin
        acq_key = where(tag_names(dcm) eq acq_tag)
        dcm_sorted = mas_dicom_group_acq(temporary(dcm), acq_key, acq_indices, 0)
    endif else begin
        acq_key = 0L
        dcm_sorted = temporary(dcm)
        acq_indices = n_elements(dcm_sorted)-1
    endelse
    
    vol_key = where(tag_names(dcm_sorted) eq vol_tag)
    img_key = where(tag_names(dcm_sorted) eq img_tag)
  
    acq_key = acq_key[0]
    vol_key = vol_key[0]
    img_key = img_key[0]
  
    acq_start = 0
    
    skip = replicate(1,  n_elements(acq_indices))
    
    for a = 0, n_elements(acq_indices)-1 do begin
  
       acq_end   = acq_indices[a]
       
       dcm_sorted[acq_start:acq_end] = mas_dicom_group_acq(dcm_sorted[acq_start:acq_end], $
                                                           [vol_key, img_key], vol_indices, 0)
       if (strpos(dcm_sorted[acq_start].image_type, want_acq) ne -1) then begin
           skip[a] = 0
       endif
       
       acq_start = acq_end+1
       
    endfor
  
    state = ptr_new({ dcm_array: ptr_new(dcm_sorted), $
                      acq_key: acq_key, $
                      vol_key: vol_key, $
                      img_key: img_key, $
                      acq_skip: ptr_new(skip), $
                      acq_indices: ptr_new(acq_indices), $
                      sel_values: tag_names(dcm_sorted) })
  
    obj_destroy, dicom_obj
    
    mas_dicom_do_import, state
    
    ptr_free, (*state).dcm_array
    ptr_free, (*state).acq_skip
    ptr_free, (*state).acq_indices
    ptr_free, state
  
end

;+
; :Description:
;    Structure definition for a DICOM directory import recipe.
;
;    acq_tag: the tag name used to sort out the separate acquisitions
;    vol_tag: the tag name used to sort out volumes within an acquisition
;    img_tag: the tag name used to sort out images within a volume
;    want_acq: set this to a string that if matched to the DICOM "ImageType" element,
;              will be loaded into mas. Acquisitions that do not match will not
;              be loaded. If this is set to '' then all acquisitons will be loaded.
;    
;    See: mas_dicom_single_file__define for the tag names recognized by MAS.
;    
; :Author: wtriplett
;-
pro mas_dicom_import_recipe__define

    struct = { MAS_DICOM_IMPORT_RECIPE, $
               acq_tag: '', $
               vol_tag: '', $
               img_tag: '', $
               want_acq: '' }

end

;+
; :Description:
;    State 1 loader for generic dicom opening.
;
; :Author: wtriplett
;-
pro mas_dicom_load_state_1

    common scan_data, project
    
    ci = project.ci
    
    file_list = *project.imndarray[ci].DICOMfiles
    fdim = project.imndarray[ci].pdim
    pdim = project.imndarray[ci].fdim
    sdim = project.imndarray[ci].sdim
    adim = project.imndarray[ci].adim
    
    mas_dicom_get_dicom_obj, dicom_obj=odicom
    
    data_array = fltarr(pdim, fdim, sdim, adim)
    
    forward_function mas_dicom_reader_checkval
    
    prog = obj_new('progressbar', color='RED', /nocancel, $
        title="Loading DICOM", /fast_loop, $
        text='Reading DICOM files...')
        
    curr_volume = 0
    curr_slice = 0
    curr_invol_slice = 0
    slice_layout = project.imndarray[ci].DICOMslicelayout
    
    prog->Start
    
    rs = 1.0
    ri = 0.0
    
    ;; Loops through the files reading them in and applying any
    ;; scaling necessary.
    for f = 0, n_elements(file_list)-1 do begin
        result = odicom->read(file_list[f])
        if (result ne 0) then begin
            tmp = odicom->getValue('7FE0'x, '0010'x, /no_copy)
            
            ;; Rescale slope
            rs_tmp = odicom->getValue('0028'x, '1053'x, /no_copy)
            ;; Rescale intercept
            ri_tmp = odicom->getValue('0028'x, '1052'x, /no_copy)
            
            ;; Scale Slope -- Philips only, other vendors not known.
            ;; This value provides the scaling factor to bring the data
            ;; back to proportional -> "Abs MR Signal"
            ;; This is the "Scale Slope" parameter which exists in PAR/REC 
            ;; files.
            ss_tmp = odicom->getValue('2005'x, '100e'x, /no_copy)

            if (mas_dicom_reader_checkval(tmp)) then begin
            
                ntmp = n_elements(tmp)
                size_fits = 0
                                
                for t = 0, n_elements(tmp)-1 do begin
    
                    if (n_elements(rs_tmp) gt t && n_elements(ri_tmp) gt t) then begin
                        if (ptr_valid(rs_tmp[t])) then begin
                            rs = float(*rs_tmp[t])
                        endif else begin
                            rs = 1.0
                        endelse
        
                        if (ptr_valid(ri_tmp[t])) then begin
                            ri = float(*ri_tmp[t])
                        endif else begin
                            ri = 0.0
                        endelse
                        
                        if (ptr_valid(ss_tmp[t])) then begin
                            ss = (float(*ss_tmp[t]))[0]
                        endif else begin
                            ss = 1.0
                        endelse
                    endif else begin
                        rs = 1.0
                        ri = 0.0
                        ss = 1.0
                    endelse
                    
                    mat_size = size(*tmp[t], /dimensions)

                    if (mat_size[0] eq pdim and mat_size[1] eq fdim) then begin
                        data_array[*,*,curr_slice,curr_volume] = reverse((rs*float(*tmp[t]) + ri)/(ss*rs),2)
                        size_fits = 1
                        continue
                    endif
                endfor
                
                if (not size_fits) then print, "Problem reading "+file_list[f]+" (img data wrongly sized.)"
                
            endif else begin
                print, "Skipping (no image data): "+file_list[f]
            endelse
        endif
        
        curr_slice++
        if (curr_slice eq (*slice_layout)[curr_volume]) then begin
            ;if (curr_slice eq project.imndarray[ci].sdim) then begin
            curr_volume++
            curr_slice = 0
        endif
        
        if (f mod 4 eq 0) then begin
            prog->Update, float(f)/float(n_elements(file_list)) * 100.
        endif
        
    endfor
    
    prog->Destroy
    
    obj_destroy, odicom
    
    project.dataarray[ci].state1 = ptr_new(data_array, /no_copy)
   
    IF project.procPramArray[CI].mc_enable NE 0 THEN BEGIN
        mas_motion_correct
    ENDIF
    
    project.dataarray[ci].state1 = mas_interpolate(project.dataarray[ci].state1)
    
    ;;mas_signal_selection, project.dataarray[ci].state1
    mas_smoothing, project.dataarray[ci].state1
    
    project.procpramarray[ci].state_1 = 1
    
end

;; $Id$
;;
;; Subroutine name: mas_dicom_opendir
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Reads a directory containing DICOM files and prepares
;; them to be sorted by the user into acquisitions and volumes
;;
;; Editing Information:

pro mas_dicom_opendir

    common common_widgets
    common scan_data
    forward_function mas_get_text_dialog    
    
    src_dir = dialog_pickfile(/directory, /read, path=project.current_path)
    
    if (src_dir eq '') then return
    
    text = '' 

    if (text eq '') then begin
        file_prefix = '*'
    endif else begin
        file_prefix = text+'*'
    endelse
    
    prog = obj_new('progressbar', color='red', /fast_loop, $
                   title='Reading file list', $
                   text='Searching for DICOM files...')
    prog->start
    load_cancelled = 0
    
    flist = file_search(src_dir, file_prefix, count=n_files) 
    if (n_files eq 0) then begin
        void = dialog_message('No DICOM files found.', /center, /error)
        prog->Destroy
        return
    endif
    
    ;; obtain a dicom reader, either IDL's provided or GDLffDICOM depending 
    ;; on operating system.
    mas_dicom_get_dicom_obj, dicom_obj=odicom
    
    ;; tracks the files as they are loaded into array.
    load_index = 0

    ;; Find valid DICOM files and order them according to the DICOM image number.
    goodness_map = bytarr(n_files)
    
    for i = 0L, n_files-1 do begin
   
       ;; skip DICOMDIR files
       if (stregex(flist[i], 'DICOMDIR$') ne -1) then continue
       
       ;; Read DICOM file, if not possible then res == 0
       res = mas_dicom_read_single_file(flist[i], $
                                        reuse_object=odicom, $
                                        images_only=1, $
                                        dcm_struct=dcm)
       
       if (res eq 1) then begin
       
           ;; first time through the loop, create an array to hold the number of files
           ;; This is done inside the loop becasue we aren't sure what kind of dicoms
           ;; (philips, siemens, generic, etc) until we start reading them
           if (n_elements(dcm_array) eq 0) then begin
               dcm_array = replicate(create_struct(name=tag_names(dcm, /structure_name)) , n_files)
           endif
           
           ;; Mark file as good
           goodness_map[load_index] = 1B
           dcm.load_index = load_index
           dcm_array[load_index++] = temporary(dcm)

       endif
              
       if (i mod 20 eq 0) then begin
           prog->update, float(i)/float(n_files) * 100.0
           if (prog->CheckCancel()) then begin
               load_cancelled = 1
               break
           endif
       endif

    endfor

    prog->Destroy
    
    if (load_cancelled) then begin
        obj_destroy, odicom
        return
    endif
    
    ;; find the "good" dicom files, ie the ones that were readable.
    good_files = where(goodness_map ne 0, n_good)
   
    if (n_good ne 0) then begin
       dcm_array = dcm_array[good_files]
    endif else begin
       void = dialog_message('No usable DICOM files found in directory.', /error, /center)
       obj_destroy, odicom
       return
    endelse

    ;; displays import GUI table.
    mas_dicom_show_import_table, dcm_array, group_leader=WID_BASE_MAIN
    
end

;; $Id$
;;
;; Subroutine name: mas_dicom_single_file_philips__define
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Structure definition for dicom files from PHILIPS systems
;;
;; Editing Information:

pro mas_dicom_single_file_philips__define

    struct = { MAS_DICOM_SINGLE_FILE_PHILIPS, $
               inherits MAS_DICOM_SINGLE_FILE, $
               volume_index: long(0), $
               slice_number: long(0), $
               b_value: float(0), $
               bvec_X: float(0), $
               bvec_Y: float(0), $
               bvec_Z: float(0), $
               orient_string: '', $
               scan_type: '' }               

end

;; Subroutine name: mas_dicom_single_file_MAS__define
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Structure definition for dicom files from MAS systems
;;
;; Editing Information:

pro mas_dicom_single_file_MAS__define

    struct = { MAS_DICOM_SINGLE_FILE_MAS, $
               inherits MAS_DICOM_SINGLE_FILE, $
               volume_index: long(0), $
               slice_number: long(0), $
               b_value: float(0), $
               bvec_X: float(0), $
               bvec_Y: float(0), $
               bvec_Z: float(0) }               
               
end

;; $Id$
;;
;; Subroutine name: mas_dicom_single_file_ge__define
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Structure definition for dicom files from GE systems
;;
;; Editing Information:

pro mas_dicom_single_file_ge__define

    struct = { MAS_DICOM_SINGLE_FILE_GE, $
               inherits MAS_DICOM_SINGLE_FILE, $
               volume_index: long(0), $
               slice_number: long(0), $
               b_value: float(0), $
               bvec_X: float(0), $
               bvec_Y: float(0), $
               bvec_Z: float(0), $
               orient_string: '', $
               scan_type: '' }    

end

;; $Id$
;;
;; Subroutine name: mas_dicom_single_file_bruker__define
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Structure definition for dicom files from BRUKER systems
;;
;; Editing Information:

pro mas_dicom_single_file_bruker__define

    struct = { MAS_DICOM_SINGLE_FILE_BRUKER, $
               inherits MAS_DICOM_SINGLE_FILE  }              

end

;; $Id$
;;
;; Subroutine name: mas_dicom_single_file_siemens__define
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Structure definition for dicom files from SIEMENS systems
;;
;; Editing Information:

pro mas_dicom_single_file_siemens__define

    struct = { MAS_DICOM_SINGLE_FILE_SIEMENS, $
               inherits MAS_DICOM_SINGLE_FILE, $
               orient_string:'', $
               b_value: float(0), $
               bvec_X: float(0), $
               bvec_Y: float(0), $
               bvec_Z: float(0) }
               
end

;; $Id$
;;
;; Subroutine name: mas_dicom_single_file__define
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Structure definition for dicom files from any system. This
;; contains the basic tags from which other systems inherit
;;
;; Editing Information:

pro mas_dicom_single_file__define

    struct = { MAS_DICOM_SINGLE_FILE, $
               load_index: long(0), $
               image_type: '', $
               spect_freq:float(0),$
               manufacturer: '', $
               study_date: long(0), $
               study_time: long(0), $
               protocol_name: '', $
               sequence_name: '', $
               series_number: long(0), $
               acquisition_number: long(0), $
               image_number: long(0), $
               slice_location: float(0), $
               x_dim: long(0), $
               y_dim: long(0), $
               pixel_spacing_x: float(0), $
               pixel_spacing_y: float(0), $
               slice_thickness: float(0), $
               slice_spacing: float(0), $
               image_position: '', $
               pat_orientation: '', $
               file_location: '', $
               parent_dir: '', $
               patient_name: '', $
               patient_pos: '', $
               num_avgs: 0, $
               num_echos: 0, $
               rescale_slope: float(0), $
               rescale_intercept: float(0), $
               TR: float(0), $
               TE: float(0), $
               TI: float(0)  }
end

;; $Id$
;;
;; Subroutine name: mas_dicom_reader
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; See: mas_dicom_opendir
;;
;; Editing Information:

pro mas_dicom_reader

    
    mas_dicom_opendir

end

