;; $Id$
;;

; IDL Widget Interface Procedures. This Code is automatically
;     generated and should not be modified.

;
; Generated on:	08/05/2007 19:05.33
;
pro WID_BASE_Frequency_Domain_lowpass_event, Event
	common scan_data
	CI = project.ci

  wTarget = (widget_info(Event.id,/NAME) eq 'TREE' ?  $
      widget_info(Event.id, /tree_root) : event.id)


  wWidget =  Event.top

  case wTarget of

    Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_apply'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin


		project.procPramArray[ci].filterflag = 'Frequency'
		project.procPramArray[ci].subtype = 'Lowpass'

		; Order
			id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_order')
				widget_control, id, GET_VALUE=(temp)
				temp = FIX(temp)
			if (temp eq 0) then temp = 1
			project.procPramArray[ci].order = temp

		; Pixelwidth
			id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_fstop')
				widget_control, id, GET_VALUE=(temp)
				temp = FIX(temp)
			;if (temp eq 0) then temp = long( project.imndArray[CI].fdim/2)
			project.procPramArray[ci].fstop = temp

		; Design
		id=widget_info(Event.top, FIND_BY_UNAME='WID_DROPLIST_type')
		selected=widget_info(id, /DROPLIST_SELECT)
			CASE selected of
			0: project.procPramArray[ci].design = 'Butterworth'
			1: project.procPramArray[ci].design = 'Hanning'
			2: project.procPramArray[ci].design = 'Hamming'
			3: project.procPramArray[ci].design = 'Kaiser-Bessel'
			4: project.procPramArray[ci].design = 'Blackman'
			ENDCASE


		; Clear state1
		project.procPramArray[ci].state_1 = 0
		project.procPramArray[ci].state_2 = 0

      endif

      end ;Widget_Info END

	Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_reset'): begin
  		if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin

			project.procPramArray[ci].filterflag = ''
			project.procPramArray[ci].subtype = ''
			project.procPramArray[ci].state_1 = 0
			project.procPramArray[ci].state_2 = 0

		endif
	end ;Widget_Info END


	Widget_Info(wWidget, FIND_BY_UNAME='WID_DROPLIST_type'): begin
		if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_DROPLIST' )then begin

			id=widget_info(Event.top, FIND_BY_UNAME='WID_DROPLIST_type')
			selected=widget_info(id, /DROPLIST_SELECT)
				If (selected eq 0) then begin

					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_fstop')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Cutoff Frequency (Pixel Width)'
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_fstop')
					widget_control, id, EDITABLE=1, SET_VALUE=''

					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_order')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Filter Order ='
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_order')
					widget_control, id, EDITABLE=1, SET_VALUE = ''

				endif else if (selected eq 1) then begin

					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_fstop')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Alpha = '
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_fstop')
					widget_control, id, EDITABLE=1, SET_VALUE=''

					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_order')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Filter Order ='
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_order')
					widget_control, id, EDITABLE=0

				endif else if (selected eq 2) then begin

					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_fstop')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Alpha = '
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_fstop')
					widget_control, id, EDITABLE=0, SET_VALUE='0.54'

					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_order')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Filter Order ='
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_order')
					widget_control, id, EDITABLE=0, SET_VALUE=''


				endif else if (selected eq 3) then begin

					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_fstop')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Alpha ='
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_fstop')
					widget_control, id, EDITABLE=1, SET_VALUE='4'

					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_order')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Filter Order ='
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_order')
					widget_control, id, EDITABLE=0,SET_VALUE=''


				endif else if (selected eq 4) then begin
					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_fstop')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Alpha zero ='
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_fstop')
					widget_control, id, EDITABLE=1, SET_VALUE='0.42'

					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_order')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Alpha one ='
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_order')
					widget_control, id, EDITABLE=1, SET_VALUE='0.5'

				endif

		endif
	end ;Widget_Info END

    ;else:
  endcase

end



pro WID_BASE_Frequency_Domain_highpass_event, Event
	common scan_data
	CI = project.ci

  wTarget = (widget_info(Event.id,/NAME) eq 'TREE' ?  $
      widget_info(Event.id, /tree_root) : event.id)


  wWidget =  Event.top

  case wTarget of

    Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_apply'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin


		project.procPramArray[ci].filterflag = 'Frequency'
		project.procPramArray[ci].subtype = 'Highpass'

		; Order
			id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_order')
				widget_control, id, GET_VALUE=(temp)
				temp = FIX(temp)
			if (temp eq 0) then temp = 1
			project.procPramArray[ci].order = temp

		; Pixelwidth
			id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_fstart')
				widget_control, id, GET_VALUE=(temp)
				temp = FIX(temp)
			if (temp eq 0) then temp = long( project.imndArray[CI].fdim/2)
			project.procPramArray[ci].fstart = temp

		; Design
		id=widget_info(Event.top, FIND_BY_UNAME='WID_DROPLIST_type')
		selected=widget_info(id, /DROPLIST_SELECT)
			CASE selected of
			0: project.procPramArray[ci].design = 'Butterworth'
			1: project.procPramArray[ci].design = 'Distance'
			ENDCASE



		; Clear state1
		project.procPramArray[ci].state_1 = 0
		project.procPramArray[ci].state_2 = 0

      endif

      end ;Widget_Info END

	Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_reset'): begin
  		if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin

			project.procPramArray[ci].filterflag = ''
			project.procPramArray[ci].subtype = ''
			project.procPramArray[ci].state_1 = 0
			project.procPramArray[ci].state_2 = 0

		endif
	end ;Widget_Info END


	Widget_Info(wWidget, FIND_BY_UNAME='WID_DROPLIST_type'): begin
		if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_DROPLIST' )then begin

			id=widget_info(Event.top, FIND_BY_UNAME='WID_DROPLIST_type')
			selected=widget_info(id, /DROPLIST_SELECT)
			If (selected eq 0) then begin

					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_fstart')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Start Frequency (Pixel Width)'
					widget_control, id, /DYNAMIC_RESIZE

					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_fstart')
					widget_control, id, EDITABLE=1, SET_VALUE=''

					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_order')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Filter Order ='
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_order')
					widget_control, id, EDITABLE=1, SET_VALUE = ''

				endif else if (selected eq 1) then begin

					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_fstart')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Start Frequency (Pixel Width)'
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_fstart')
					widget_control, id, EDITABLE=0, SET_VALUE=''

					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_order')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Filter Order ='
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_order')
					widget_control, id, EDITABLE=0, SET_VALUE=''

				endif


			;print, 'Inside Droplist type'

		endif
	end ;Widget_Info END


    ;else:
  endcase

end



pro WID_BASE_Frequency_Domain_lowpass, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_

  Resolve_Routine, 'mas_image_filters',/COMPILE_FULL_FILE  ; Load event callback routines

  WID_BASE_Frequency_Domain_lowpass = Widget_Base( GROUP_LEADER=wGroup,  $
      UNAME='WID_BASE_Frequency_Domain_lowpass' ,XOFFSET=5 ,YOFFSET=5  $
      ,SCR_XSIZE=300 ,SCR_YSIZE=200 ,TITLE='Frequency Domain'+ $
      ' Filtering' ,SPACE=3 ,XPAD=3 ,YPAD=3)

  WID_LABEL_type = Widget_Label(WID_BASE_Frequency_Domain_lowpass,  $
      UNAME='WID_LABEL_type' ,XOFFSET=15 ,YOFFSET=10 ,/ALIGN_LEFT  $
      ,VALUE='Type:')


  WID_DROPLIST_type = Widget_Droplist(WID_BASE_Frequency_Domain_lowpass,  $
      UNAME='WID_DROPLIST_type' ,XOFFSET=15 ,YOFFSET=30 ,VALUE=[  $
      'Butterworth', 'Hanning Window', 'Hamming Window', 'Kaiser-Bessel Window', 'Blackman Window'])


  WID_LABEL_fstop = Widget_Label(WID_BASE_Frequency_Domain_lowpass,  $
      UNAME='WID_LABEL_fstop' ,XOFFSET=15 ,YOFFSET=65  $
      ,/ALIGN_LEFT ,VALUE='Cutoff Frequency (Pixel Width)')


  WID_TEXT_fstop = Widget_Text(WID_BASE_Frequency_Domain_lowpass,  $
      UNAME='WID_TEXT_fstop' ,XOFFSET=15 ,YOFFSET=85  $
      ,SCR_XSIZE=60 ,SCR_YSIZE=26 ,/EDITABLE ,XSIZE=20 ,YSIZE=1)


  WID_BUTTON_reset = Widget_Button(WID_BASE_Frequency_Domain_lowpass,  $
      UNAME='WID_BUTTON_reset' ,XOFFSET=170 ,YOFFSET=110  $
      ,SCR_XSIZE=116 ,SCR_YSIZE=22 ,/ALIGN_CENTER ,VALUE='Reset')


  WID_BUTTON_apply = Widget_Button(WID_BASE_Frequency_Domain_lowpass,  $
      UNAME='WID_BUTTON_apply' ,XOFFSET=170 ,YOFFSET=140  $
      ,SCR_XSIZE=116 ,SCR_YSIZE=22 ,/ALIGN_CENTER ,VALUE='Apply')


  WID_LABEL_order = Widget_Label(WID_BASE_Frequency_Domain_lowpass,  $
      UNAME='WID_LABEL_order' ,XOFFSET=15 ,YOFFSET=115 ,/ALIGN_LEFT  $
      ,VALUE='Filter Order = ')


  WID_TEXT_order = Widget_Text(WID_BASE_Frequency_Domain_lowpass,  $
      UNAME='WID_TEXT_order' ,XOFFSET=15 ,YOFFSET=130 ,SCR_XSIZE=60  $
      ,SCR_YSIZE=26 ,/EDITABLE ,XSIZE=20 ,YSIZE=1)

  Widget_Control, /REALIZE, WID_BASE_Frequency_Domain_lowpass

  XManager, 'WID_BASE_Frequency_Domain_lowpass', WID_BASE_Frequency_Domain_lowpass, /NO_BLOCK

end



pro WID_BASE_Frequency_Domain_highpass, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_

  Resolve_Routine, 'mas_image_filters',/COMPILE_FULL_FILE  ; Load event callback routines

  WID_BASE_Frequency_Domain_highpass = Widget_Base( GROUP_LEADER=wGroup,  $
      UNAME='WID_BASE_Frequency_Domain_highpass' ,XOFFSET=5 ,YOFFSET=5  $
      ,SCR_XSIZE=300 ,SCR_YSIZE=200 ,TITLE='Frequency Domain'+ $
      ' Filtering' ,SPACE=3 ,XPAD=3 ,YPAD=3)


  WID_LABEL_type = Widget_Label(WID_BASE_Frequency_Domain_highpass,  $
      UNAME='WID_LABEL_type' ,XOFFSET=15 ,YOFFSET=10 ,/ALIGN_LEFT  $
      ,VALUE='Type:')


  WID_DROPLIST_type = Widget_Droplist(WID_BASE_Frequency_Domain_highpass,  $
      UNAME='WID_DROPLIST_type' ,XOFFSET=15 ,YOFFSET=30 ,VALUE=[  $
      'Butterworth', 'Distance Window'])


  WID_LABEL_fstart = Widget_Label(WID_BASE_Frequency_Domain_highpass,  $
      UNAME='WID_LABEL_fstart' ,XOFFSET=15 ,YOFFSET=65  $
      ,/ALIGN_LEFT ,VALUE='Start Frequency')

  WID_TEXT_fstart = Widget_Text(WID_BASE_Frequency_Domain_highpass,  $
      UNAME='WID_TEXT_fstart' ,XOFFSET=15 ,YOFFSET=85  $
      ,SCR_XSIZE=60, /EDITABLE ,SCR_YSIZE=26  ,XSIZE=20 ,YSIZE=1)


  WID_BUTTON_reset = Widget_Button(WID_BASE_Frequency_Domain_highpass,  $
      UNAME='WID_BUTTON_reset' ,XOFFSET=170 ,YOFFSET=110  $
      ,SCR_XSIZE=116 ,SCR_YSIZE=22 ,/ALIGN_CENTER ,VALUE='Reset')


  WID_BUTTON_apply = Widget_Button(WID_BASE_Frequency_Domain_highpass,  $
      UNAME='WID_BUTTON_apply' ,XOFFSET=170 ,YOFFSET=140  $
      ,SCR_XSIZE=116 ,SCR_YSIZE=22 ,/ALIGN_CENTER ,VALUE='Apply')


  WID_LABEL_order = Widget_Label(WID_BASE_Frequency_Domain_highpass,  $
      UNAME='WID_LABEL_order' ,XOFFSET=15 ,YOFFSET=115 ,/ALIGN_LEFT  $
      ,VALUE='Order')


  WID_TEXT_order = Widget_Text(WID_BASE_Frequency_Domain_highpass,  $
      UNAME='WID_TEXT_order' ,XOFFSET=16 ,YOFFSET=130 ,SCR_XSIZE=60  $
      ,SCR_YSIZE=26 ,/EDITABLE ,XSIZE=20 ,YSIZE=1)

  Widget_Control, /REALIZE, WID_BASE_Frequency_Domain_highpass

  XManager, 'WID_BASE_Frequency_Domain_highpass', WID_BASE_Frequency_Domain_highpass, /NO_BLOCK

end



; HS 20070810
; Handles the events from the Image Domain Filtering GUI
; Types: 'Lee Adaptive Filter','Wavelet Transform', 'Sharpening', 'Low Pass Filter', 'High Pass Filter'
; Parameters edited: image_domain_filterflag: '', image_domain_subtype:'',image_domain_parameter:1

pro WID_BASE_IMAGE_DOMAIN_FILTERING_event, Event
	common scan_data
	CI = project.ci

  wTarget = (widget_info(Event.id,/NAME) eq 'TREE' ?  $
      widget_info(Event.id, /tree_root) : event.id)


  wWidget =  Event.top

  case wTarget of

    Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_apply'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin

		project.procPramArray[ci].image_domain_filterflag = 'Image Domain'

		; Filter Parameter
		id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_parameter')
		widget_control, id, GET_VALUE=(temp)
			temp = FIX(temp)
			if (temp eq 0) then temp = 1
			project.procPramArray[ci].image_domain_parameter = temp

		; Design
		id=widget_info(Event.top, FIND_BY_UNAME='WID_DROPLIST_subtype')
		selected=widget_info(id, /DROPLIST_SELECT)
			CASE selected of
			0: project.procPramArray[ci].image_domain_subtype = 'Lee Adaptive Filter'
			1: project.procPramArray[ci].image_domain_subtype = 'Wavelet Transform'
			2: project.procPramArray[ci].image_domain_subtype = 'Sharpening'
			3: project.procPramArray[ci].image_domain_subtype = 'Low Pass Filter'
			4: project.procPramArray[ci].image_domain_subtype = 'High Pass Filter'
			ENDCASE


		; Clear state1
		project.procPramArray[ci].state_1 = 0
		project.procPramArray[ci].state_2 = 0

      endif

      end ;Widget_Info END

	Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_reset'): begin
  		if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then begin

			project.procPramArray[ci].image_domain_filterflag = ''
			project.procPramArray[ci].image_domain_subtype = ''
			project.procPramArray[ci].state_1 = 0
			project.procPramArray[ci].state_2 = 0

		endif
	end ;Widget_Info END

; Change Widget names to Image Domain appropriate one

	Widget_Info(wWidget, FIND_BY_UNAME='WID_DROPLIST_subtype'): begin
		if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_DROPLIST' )then begin

			id=widget_info(Event.top, FIND_BY_UNAME='WID_DROPLIST_subtype')
			selected=widget_info(id, /DROPLIST_SELECT)
			If (selected eq 0) then begin

					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_parameter')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Filter box size = '
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_parameter')
					widget_control, id, EDITABLE=1, SET_VALUE='1'


				endif else if (selected eq 1) then begin

					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_parameter')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Percent of Image Power to Keep = '
					widget_control, id, /DYNAMIC_RESIZE

					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_parameter')
					widget_control, id, EDITABLE=1, SET_VALUE='99'

				endif else if (selected eq 2) then begin
					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_parameter')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Sharpening Kernel = '
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_parameter')
					widget_control, id, EDITABLE=0, SET_VALUE='1'

				endif else if (selected eq 3) then begin
					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_parameter')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='Low Pass Kernel = '
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_parameter')
					widget_control, id, EDITABLE=0, SET_VALUE='1'

				endif else if (selected eq 4) then begin
					id=widget_info(Event.top, FIND_BY_UNAME='WID_LABEL_parameter')
					widget_control, id, /DYNAMIC_RESIZE
					widget_control, id, SET_VALUE ='High Pass Kernel = '
					widget_control, id, /DYNAMIC_RESIZE
					id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_parameter')
					widget_control, id, EDITABLE=0, SET_VALUE='1'

				endif


			;print, 'Inside Droplist type'

		endif
	end ;Widget_Info END


    ;else:
  endcase

end





; HS 20070810
; Launches the Image domain filtering GUI
PRO WID_BASE_IMAGE_DOMAIN_FILTERING, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_


  WID_BASE_IMAGE_DOMAIN_FILTERING = Widget_Base( GROUP_LEADER=wGroup,  $
      UNAME='WID_BASE_IMAGE_DOMAIN_FILTERING' ,XOFFSET=5 ,YOFFSET=5  $
      ,SCR_XSIZE=300 ,SCR_YSIZE=200 ,TITLE='Image Domain'+ $
      ' Filtering' ,SPACE=3 ,XPAD=3 ,YPAD=3)


  WID_LABEL_subtype = Widget_Label(WID_BASE_IMAGE_DOMAIN_FILTERING,  $
      UNAME='WID_LABEL_subtype' ,XOFFSET=15 ,YOFFSET=10 ,/ALIGN_LEFT  $
      ,VALUE='Type:')


  WID_DROPLIST_subtype = Widget_Droplist(WID_BASE_IMAGE_DOMAIN_FILTERING,  $
      UNAME='WID_DROPLIST_subtype' ,XOFFSET=15 ,YOFFSET=30 ,VALUE=[  $
      'Lee Adaptive Filter ','Wavelet Transform', 'Sharpening', 'Low Pass Filter', 'High Pass Filter'])


  WID_LABEL_parameter = Widget_Label(WID_BASE_IMAGE_DOMAIN_FILTERING,  $
      UNAME='WID_LABEL_parameter' ,XOFFSET=15 ,YOFFSET=65  $
      ,/ALIGN_LEFT ,VALUE='Filtering Parameter')

  WID_TEXT_parameter = Widget_Text(WID_BASE_IMAGE_DOMAIN_FILTERING,  $
      UNAME='WID_TEXT_parameter' ,XOFFSET=15 ,YOFFSET=85  $
      ,SCR_XSIZE=60, /EDITABLE  ,XSIZE=20 , VALUE='1')


  WID_BUTTON_reset = Widget_Button(WID_BASE_IMAGE_DOMAIN_FILTERING,  $
      UNAME='WID_BUTTON_reset' ,XOFFSET=170 ,YOFFSET=110  $
      ,SCR_XSIZE=116 ,SCR_YSIZE=22 ,/ALIGN_CENTER ,VALUE='Reset')


  WID_BUTTON_apply = Widget_Button(WID_BASE_IMAGE_DOMAIN_FILTERING,  $
      UNAME='WID_BUTTON_apply' ,XOFFSET=170 ,YOFFSET=140  $
      ,SCR_XSIZE=116 ,SCR_YSIZE=22 ,/ALIGN_CENTER ,VALUE='Apply')

;
;  WID_LABEL_order = Widget_Label(WID_BASE_IMAGE_DOMAIN_FILTERING,  $
;      UNAME='WID_LABEL_order' ,XOFFSET=15 ,YOFFSET=115 ,/ALIGN_LEFT  $
;      ,VALUE='Order')
;
;
;  WID_TEXT_order = Widget_Text(WID_BASE_IMAGE_DOMAIN_FILTERING,  $
;      UNAME='WID_TEXT_order' ,XOFFSET=16 ,YOFFSET=130 ,SCR_XSIZE=60  $
;      ,SCR_YSIZE=26 ,/EDITABLE ,XSIZE=20 ,YSIZE=1)

  Widget_Control, /REALIZE, WID_BASE_IMAGE_DOMAIN_FILTERING

  XManager, 'WID_BASE_IMAGE_DOMAIN_FILTERING', WID_BASE_IMAGE_DOMAIN_FILTERING, /NO_BLOCK

END

;
; Empty stub procedure used for autoloading.
;
pro mas_filtering, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_

end

