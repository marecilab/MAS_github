; Subroutine name: select_data_gui
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Data selection GUI which allows which arrays to fit for analysis (T1,T2,ADC & ADT)
;
; Editing Information:

pro write_data_txt_gui, index = index
  common scan_data

  adim = project.imndarray[index].adim
  if (adim le 1) then begin
    mbox = dialog_message('Dataset is not arrayed',/error)
    return
  endif else project.procpramArray[index].array_select = ptr_new(indgen(adim))

  if xregistered('select_data_gui') ge 1 then return
  top_base = widget_base(title = 'Select Data : Scan #' + strtrim(index+1,2),tlb_frame_attr=1,/grid_layout)
  main_base = widget_base(top_base,row = 2)
  dat_buttons = strtrim(indgen(adim),2)
  dat_button_grp = cw_bgroup(main_base,row = sqrt(adim),dat_buttons, label_top = 'Array Indices', uvalue = 'dat_buttons', /nonexclusive, button_uvalue = dat_buttons, frame = 2, set_value=make_array(adim,value = 1),/no_release)
  misc_buttons_txt = ['Apply', 'Reset', 'Close']
  misc_button_grp = cw_bgroup(main_base,column = 3, misc_buttons_txt, uvalue = 'misc_buttons', button_uvalue = misc_buttons_txt, /no_release, space = 10*sqrt(adim))

  infoptr = ptr_new({dat_button_grp:dat_button_grp, misc_button_grp:misc_button_grp, index:index})
  widget_control, top_base, set_uvalue=infoptr
  widget_control, top_base, /realize
  xmanager, 'select_data_gui', top_base, cleanup='select_data_gui_cleanup',/no_block

end

; Subroutine name: select_data_gui_event
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Event handler for the data selection GUI
;
; Editing Information:

pro select_data_gui_event, event

  common scan_data

  widget_control, event.top, get_uvalue = infoptr
  info = *infoptr
  widget_control, event.id, get_uvalue = widget

  case widget of

    'misc_buttons' : begin
      ; Set ADT process flag to false
      project.procpramArray[info.index].adt_proccess_flag = 0
      case event.value of
        'Apply' :  begin
          widget_control, info.dat_button_grp, get_value = indices
          project.procpramArray[info.index].array_select = ptr_new(where(indices eq 1))
        end
        'Reset'   : widget_control, info.dat_button_grp, set_value= make_array(project.imndArray[info.index].adim,value = 1)
        'Close'   : widget_control, event.top,/destroy
      endcase
    end
    else          : return
  endcase

end

; Subroutine name: select_data_gui_cleanup
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Cleanup procedure before exiting the data selection GUI
;
; Editing Information:

pro select_data_gui_cleanup, top_base

  common scan_data

  widget_control, top_base, get_uvalue=state
  if ptr_valid(state) then begin
    if ptr_valid(project.procpramArray[(*state).index].array_select) then ptr_free,project.procpramArray[(*state).index].array_select
    ptr_free, state
  endif
  return

end