;; $id$
;; Copyright 2003 University of Florida. All Rights Reserved

;*****************************************************************************************************
;
; NAME:
;   mas_display
;
; PURPOSE:
;   to post the selected image to the screen and proceed how the user wants it.
;
; ARGUMENTS:
;
;
; MODIFICATION HISTORY:
;
;
;*****************************************************************************************************


; Subroutine name: Pad_Image
; Created by: CD 2007/05/09
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:

FUNCTION Pad_Image, curr_ptr, vw_xsize, vw_ysize

	img_sz = SIZE(*curr_ptr)
	idx = (img_sz[1] LT vw_xsize) + (img_sz[2] LT vw_ysize) * 2

	idx = 0

	CASE idx OF
		0: BEGIN
			xstart = 0
			ystart = 0
			img_ptr = curr_ptr
		END

		1: BEGIN
			new_img = INTARR(vw_xsize,img_sz[2])
			new_img[*,*] = 255
			xstart  = FLOOR((vw_xsize - img_sz[1]) / 2)
			new_img[xstart:xstart+img_sz[1]-1,*] = *curr_ptr
			img_ptr = PTR_NEW(new_img)
			ystart = 0
		END

		2: BEGIN
			new_img = INTARR(img_sz[1],vw_ysize)
			new_img[*,*] = 255
			ystart  = FLOOR((vw_ysize - img_sz[2]) / 2)
			new_img[*,ystart:ystart+img_sz[2]-1] = *curr_ptr
			img_ptr = PTR_NEW(new_img)
			xstart = 0
        END

		3: BEGIN
			new_img = INTARR(vw_xsize,vw_ysize)
			new_img[*,*] = 255
			ystart  = FLOOR((vw_ysize - img_sz[2]) / 2)
			xstart  = FLOOR((vw_xsize - img_sz[1]) / 2)

			new_img[xstart:xstart+img_sz[1]-1,ystart:ystart+img_sz[2]-1] = *curr_ptr
			img_ptr = PTR_NEW(new_img)
		END
	ENDCASE

	RETURN, {img_ptr:img_ptr, x:xstart, y:ystart}

END



; Subroutine name: Array_Zoom_2D
; Created by: CD 2006/02/26
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:

FUNCTION Array_Zoom_2D, data, currzoom, newzoom, usefloat

	sz_data = SIZE(data)
	xsize = sz_data[1]

	sz_scan = sz_data / currzoom

	unq_idx_x = INDGEN(sz_scan[1]) * currzoom
	unq_idx_y = INDGEN(sz_scan[2]) * currzoom

	IF usefloat EQ 1 THEN BEGIN
		out_idx_x = FLTARR(sz_scan[1] * newzoom)
		out_idx_y = FLTARR(sz_scan[2] * newzoom)
	ENDIF ELSE BEGIN
		out_idx_x = INTARR(sz_scan[1] * newzoom)
		out_idx_y = INTARR(sz_scan[2] * newzoom)
	ENDELSE

	tmp_idx_x = INDGEN(sz_scan[1]) * newzoom
	tmp_idx_y = INDGEN(sz_scan[2]) * newzoom

	FOR i=0, newzoom-1 DO BEGIN
		out_idx_x(tmp_idx_x+i) = unq_idx_x
		out_idx_y(tmp_idx_y+i) = unq_idx_y
	ENDFOR

	outsz_x = (SIZE(out_idx_x))[1]
	outsz_y = (SIZE(out_idx_y))[1]

	IF usefloat EQ 1 THEN retarr = FLTARR(outsz_x,outsz_y) $
	ELSE retarr = INTARR(outsz_x,outsz_y)

	tmp_idx = INDGEN(outsz_x)

	FOR i = 0, outsz_y-1 DO BEGIN
		retarr[(tmp_idx+(i*outsz_x))] = data[(out_idx_x + (xsize * out_idx_y(i)))]
	ENDFOR

	RETURN, retarr
END


; Subroutine name: Display_DrawRes_Resize
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; Resize Event Callback Procedure

; Editing Information:
    ;Edited by HS 2006/10/04
    ;Fix spelling mistakes and commenting.
    ;Edited by CD 2007/01/26
    ;Add zoom capability that maintains aspect ratio
    ;Edited so that pixels are copied/deleted, no interpolation

PRO Display_DrawRes_Resize, ev

	Widget_Control, ev.top, Get_UValue=info, /No_Copy
	COMMON display_data

  	CASE ev.id OF

    	Widget_Info(ev.top, FIND_BY_UNAME='zoom_button'): BEGIN
			display_mouseStatus = 1
    	END

		Widget_Info(ev.top, FIND_BY_UNAME='scroll_button'): BEGIN
			display_mouseStatus = 2
    	END

		Widget_Info(ev.top, FIND_BY_UNAME='intensity_button'): BEGIN
			display_mouseStatus = 3
    	END

		Widget_Info(ev.top, FIND_BY_UNAME='resize_interp_button'): BEGIN
			display_mouseStatus = 4
    	END

		Widget_Info(ev.top, FIND_BY_UNAME='resize_zoom_button'): BEGIN
			display_mouseStatus = 5
    	END

		Widget_Info(ev.top, FIND_BY_UNAME='rotate_zoom_button'): BEGIN
			display_mouseStatus = 6
		END

		Widget_Info(ev.top, FIND_BY_UNAME='apply_button'): BEGIN

    	END

		Widget_Info(ev.top, FIND_BY_UNAME='do nothing'): BEGIN

    	END

		ELSE: BEGIN

    	END
	ENDCASE

	Widget_Control, ev.top, Set_UValue=info, /No_Copy
END


; Subroutine name: Display_DrawRes_Cleanup
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; Clean up pointers procedure

; Editing Information:
    ;Edited by HS 2006/10/04
    ;Fix spelling mistakes and commenting.
    ;Edited by CD 2007/01/31
    ;Delete objects from adding ROI/Zoom

PRO Display_DrawRes_Cleanup, tlb

	COMMON display_data
	display_baseID = 0
	display_mouseStatus = 0

    Widget_Control, tlb, Get_UValue=info, /No_Copy
    IF N_Elements(info) EQ 0 THEN RETURN

    COMPILE_OPT IDL2
    common scan_data
    ci = project.ci

    Ptr_Free, info.image
    Ptr_Free, info.currimage
    OBJ_DESTROY, info.oImage
    OBJ_DESTROY, info.oView
    OBJ_DESTROY, info.oModel
END


PRO Create_IMG_Tab, tabID, image, ffov, pfov, FPVALS=fpvals, TABTITLE=tabtitle

	scf_dims = GET_SCREEN_SIZE()
	vw_xsize = scf_dims[0] - 100
	vw_ysize = scf_dims[1] - 200; 150

	IF KEYWORD_SET(fpvals) NE 1 THEN BEGIN
		fpvals = 0
		fpdisp = 0
	ENDIF ELSE fpdisp = 1

        sz_scan = size(image)

	vw_xloc  = 0
	vw_yloc  = 0

        offsets = Pad_Image(PTR_NEW(image), vw_xsize, vw_ysize)
        image   = offsets.img_ptr

	IF sz_scan[0] EQ 3 THEN BEGIN
		oImage = OBJ_NEW('IDLgrImage',*image)

		IF sz_scan[1] EQ 3 THEN oImage->SetProperty, INTERLEAVE=0 $
		ELSE oImage->SetProperty, INTERLEAVE=2

		interleave = 1
	ENDIF ELSE BEGIN
		oImage = OBJ_NEW('IDLgrImage', *image)
		interleave = 0
	ENDELSE

    drawID = Widget_Draw(tabID, XSIZE=vw_xsize, YSIZE=vw_ysize, $
    	BUTTON_EVENTS=1, EXPOSE_EVENTS=1, MOTION_EVENTS=1, GRAPHICS_LEVEL=2, $
    	/ALIGN_CENTER)
    Widget_Control, tabID, /Realize
    Widget_Control, drawID, Get_Value=oWindow

    oView = OBJ_NEW('IDLgrView', VIEWPLANE_RECT=[vw_xloc,vw_yloc,vw_xsize,vw_ysize])
    oModel = OBJ_NEW('IDLgrModel')
    oView->Add, oModel
	oModel->Add, oImage

	oWindow->Draw, oView

	xfov = pfov
	yfov = ffov

    info = {image:image, oWindow:oWindow, drawID:drawID, $
    	xsize:sz_scan[1], ysize:sz_scan[2], oView:oView, oModel:oModel, $
    	oROI:OBJ_NEW(), buttonXY: LONARR(2), buttonDown: 0B, $
    	oImage:oImage, currimage:image, initx:sz_scan[1], $
    	inity:sz_scan[2], roizfac:1, imgzfac:1, priorzfac:1, $
    	xfov:xfov, yfov:yfov, xfov_init:xfov, yfov_init:yfov, $
    	xfov_offset:0.0, yfov_offset:0.0, $
    	int_Cen:0, int_Max:0, int_Min:0, fpvals:PTR_NEW(fpvals), fpdisp:fpdisp, $
    	currfpvals:PTR_NEW(fpvals), interleave:interleave, $
    	vw_xloc:vw_xloc,vw_yloc:vw_yloc,vw_xloc_init:vw_xloc,vw_yloc_init:vw_yloc, $
    	vw_xsize:vw_xsize, vw_ysize:vw_ysize, padx:offsets.x, pady:offsets.y, $
    	padx_init:offsets.x, pady_init:offsets.y }

    Widget_Control, tabID, Set_UValue=info, /No_Copy

    XMANAGER, 'Display_DrawRes_Create', drawID, Event_Handler='Display_Zoom_Event', /NO_BLOCK

END

pro Display_DrawRes_Save, event

    help, event, /structure

end

; Subroutine name: Display_DrawRes_Create
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
; Also called in mas_adt_regress:adt_display (to allow ADT to use
; the tool enhanced image window

; Purpose of subroutine:
; Procedure to create resizeable draw widget for image display with
; constrast enhancement

; Editing Information:
    ;Edited by HS 2006/10/04
    ;Fix spelling mistakes and commenting.
    ;Edited by CD 2007/01/25
    ;Add zoom functionality

PRO Display_DrawRes_Create, image, display_title, ffov, pfov, FPVALS=fpvals, TAB_TITLE=tab_title

	COMMON display_data
        
        mas_display_multi, image
        return

        if ((size(image))[0] lt 2 or (size(image))[0] gt 3) then begin
            print, (size(image))[0]
            void = dialog_message('"image" is not a two dimensional image '+ $
                                  '@ Display_DrawRes_Create', /error, /center)
            return
        endif

	IF display_baseID EQ 0 THEN BEGIN

            if (not keyword_set(tab_title)) then begin
                tab_title = 'Image 1'
            endif else begin
                tab_title = strcompress(tab_title)
            endelse
            
            baseID = Widget_Base(Column=1, Title=display_title, TLB_Size_Events=1, mbar=mbar)

            btn_file = widget_button(mbar, value='File', /menu)
            btn_save = widget_button(btn_file, uname='btn_save', value='Save...', event_pro='Display_DrawRes_Save')

            display_baseID = baseID
            Widget_Control, baseID, /Realize
            
            tabID  = WIDGET_TAB(baseID, UNAME='do nothing')
            img_one_tab = WIDGET_BASE(tabID,title=tab_title, UNAME='do nothing')
            
            IF N_ELEMENTS(fpvals) EQ 0 THEN BEGIN
                Create_IMG_Tab, img_one_tab, image, ffov, pfov
            ENDIF ELSE BEGIN
                Create_IMG_Tab, img_one_tab, image, ffov, pfov, FPVALS=fpvals
            ENDELSE
            
            display_labelID = Widget_Label(baseID, SCR_XSIZE=200 ,SCR_YSIZE=30, $
                                           FRAME=0 ,/ALIGN_LEFT ,VALUE='Zoom Factor=1')
            
            butBaseID = WIDGET_BASE(baseID, /EXCLUSIVE, ROW=1, FRAME=2)
            butoneID  = WIDGET_BUTTON(butBaseID, VALUE = 'ROI', SCR_XSIZE=70, $
                                      /ALIGN_LEFT, UNAME='zoom_button') ; button_func 1
            buttwoID  = WIDGET_BUTTON(butBaseID, VALUE = 'Scroll', SCR_XSIZE=70, $
                                      /ALIGN_LEFT, UNAME='scroll_button') ; button_func 2
            butthreeID  = WIDGET_BUTTON(butBaseID, VALUE = 'Intensity', SCR_XSIZE=70, $
                                        /ALIGN_LEFT, UNAME='intensity_button') ; button_func 3
            butfourID  = WIDGET_BUTTON(butBaseID, VALUE = 'Interp RS', SCR_XSIZE=70, $
                                       /ALIGN_LEFT, UNAME='resize_interp_button') ; button_func 4
            butfiveID  = WIDGET_BUTTON(butBaseID, VALUE = 'Zoom', SCR_XSIZE=70, $
                                       /ALIGN_LEFT, UNAME='resize_zoom_button') ; button_func 5
            butsixID  = WIDGET_BUTTON(butBaseID, VALUE = 'Rotate', SCR_XSIZE=70, $
                                      /ALIGN_LEFT, UNAME='rotate_zoom_button') ; button_func 6
            
            butapply  = WIDGET_BUTTON(baseID, VALUE = 'APPLY SIZING TO ALL IMAGES', $
                                      /ALIGN_CENTER, UNAME='apply_button')
            
            XMANAGER, 'Display_DrawRes_Create', baseID, Event_Handler='Display_DrawRes_Resize', Cleanup='Display_DrawRes_Cleanup',/NO_BLOCK
        ENDIF ELSE BEGIN

            tabIDs = Widget_Info(display_baseID, /ALL_CHILDREN)
            tabID = tabIDs[1]

            numtabs = WIDGET_INFO(tabID, /TAB_NUMBER)
            
            if (not keyword_set(tab_title)) then begin
                tab_title = 'Image '+ STRTRIM(STRING(numtabs+1),2)
            endif else begin
                tab_title = strcompress(tab_title)
            endelse

            newtab = WIDGET_BASE(tabID,title=tab_title,  UNAME='do nothing')
            
            IF N_ELEMENTS(fpvals) EQ 0 THEN BEGIN
                Create_IMG_Tab, newtab, image, ffov, pfov
            ENDIF ELSE BEGIN
                Create_IMG_Tab, newtab, image, ffov, pfov, FPVALS=fpvals
            ENDELSE
            
            widget_control, tabID, set_tab_current=numtabs

        ENDELSE
    END
    
; Subroutine name: Display_Zoom_Event
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:


PRO Display_Zoom_Event, event

	COMMON scan_data
	COMMON display_data

	ZOOM_SENSITIVITY = 50
	ZOOM_MAX = 10
	SCROLL_SENSITIVITY = 1

        tabIDs   = Widget_Info(event.top, /ALL_CHILDREN)
        tabID = tabIDs[1]
	currTab = Widget_Info(tabID, /TAB_CURRENT)
	tabID   = Widget_Info(tabID, /CHILD)
	FOR ii=0, currTab-1 DO tabID = Widget_Info(tabID, /SIBLING)

	; reset image on double click
	IF event.clicks EQ 2 THEN BEGIN

		WIDGET_CONTROL, tabID, GET_UVALUE=info, /No_Copy

		sz_scan = size((*info.image))

		xsize = sz_scan[1]
		ysize = sz_scan[2]

		info.oView->SetProperty, VIEWPLANE_RECT=[info.vw_xloc_init, $
			info.vw_yloc_init, info.vw_xsize, info.vw_ysize]
		info.oImage->SetProperty, DATA=(*info.image)

		info.xsize = xsize
		info.ysize = ysize
		info.initx = xsize
		info.inity = ysize
		info.currimage = info.image
		info.currfpvals = info.fpvals
		info.roizfac = 1
		info.imgzfac = 1
		info.priorzfac = 1
		info.xfov_offset = 0
		info.yfov_offset = 0
		info.xfov = info.xfov_init
		info.yfov = info.yfov_init
		info.vw_xloc = info.vw_xloc_init
		info.vw_yloc = info.vw_yloc_init
		info.padx    = info.padx_init
		info.pady    = info.pady_init

		info.oWindow->Draw, info.oView

		WIDGET_CONTROL, display_labelID, SET_VALUE='Zoom Factor=1'
		WIDGET_CONTROL, tabID, SET_UVALUE=info, /No_Copy

		RETURN
	END

	CASE event.type OF
		; Button Press
		0: BEGIN
			WIDGET_CONTROL, tabID, GET_UVALUE=info, /No_Copy

			CASE display_mouseStatus OF
				1: BEGIN ;roi

					xcoord = event.x+info.vw_xloc-info.padx
					ycoord = event.y+info.vw_yloc-info.pady

    				oROI = OBJ_NEW('IDLgrROI',COLOR=BYTE([255, 0, 0]),STYLE=0)
					oROI->AppendData, [xcoord,ycoord,0]

					info.oModel->Add, oROI

    				info.oWindow->Draw, info.oView

					info.oROI = oROI
					info.buttonDown = 1b
            		info.buttonXY = [xcoord,ycoord]
            	END

				2: BEGIN ;scroll
				    info.buttonDown = 1b
            		info.buttonXY = [event.x, event.y]
				END

				3: BEGIN ;intensity
					; set the initial sliders
    				info.int_Cen = MAX(*info.currimage)
    				info.int_Max = info.int_Cen
    				info.int_Min = 0

            		info.buttonDown = 1b
            		info.buttonXY = [event.x, event.y]
				END

				5: BEGIN ;zoom
				    info.buttonDown = 1b
            		info.buttonXY = [event.x, event.y]
				END

				6: BEGIN ;rotate
				    info.buttonDown = 1b
            		info.buttonXY = [event.x, event.y]
				END

				ELSE: BEGIN
				END
			ENDCASE

            WIDGET_CONTROL, tabID, SET_UVALUE=info, /No_Copy
		END

        ; Button Release
        1: BEGIN
	        WIDGET_CONTROL, tabID, GET_UVALUE=info, /No_Copy

			CASE display_mouseStatus OF

				1: BEGIN ;roi
					IF (info.buttonDown NE 0) THEN BEGIN

						info.buttonDown=0b

            			oROI = info.oROI
						oROI->GetProperty, DATA=roiData

						; check to make sure rectangale has four sides
						IF ((N_ELEMENTS(roiData)/3) EQ 4) THEN BEGIN

							x1 = MIN(roiData[0,*]) - info.padx
							y1 = MIN(roiData[1,*]) - info.pady

							x2 = MAX(roiData[0,*]) - info.padx
							y2 = MAX(roiData[1,*]) - info.pady

							currx = info.xsize
							curry = info.ysize

							; check to make sure all coordinates of rectangle are on the image
							IF (x1 LT currx) AND (x1 GE 0) AND (x2 LT currx) AND (x2 GE 0) AND $
								(y1 LT curry) AND (y1 GE 0) AND (y2 LT curry) AND (y2 GE 0) $
								THEN BEGIN

								; check to see that that the selection includes all of
								; zoom voxel, if not adjust

								x1 = x1 - (x1 MOD info.imgzfac) + info.padx
								y1 = y1 - (y1 MOD info.imgzfac) + info.pady

								x2 = x2 + (info.imgzfac - (x2 MOD info.imgzfac)) - 1 + info.padx
								y2 = y2 + (info.imgzfac - (y2 MOD info.imgzfac)) - 1 + info.pady

								xsize = x2-x1+1
								ysize = y2-y1+1

								IF info.interleave EQ 1 THEN BEGIN
									new_image = BYTARR(xsize,ysize,3)
									new_image[*,*,0] = (*info.currimage)[x1:x2,y1:y2,0]
									new_image[*,*,1] = (*info.currimage)[x1:x2,y1:y2,1]
									new_image[*,*,2] = (*info.currimage)[x1:x2,y1:y2,2]
								ENDIF ELSE BEGIN
									new_image = (*info.currimage)[x1:x2,y1:y2]
								ENDELSE

								IF info.fpdisp EQ 1 THEN BEGIN
									info.currfpvals = PTR_NEW(((*info.currfpvals)[x1:x2,y1:y2]))
								ENDIF

								offsets = Pad_Image(PTR_NEW(new_image), info.vw_xsize, info.vw_ysize)
								info.currimage = offsets.img_ptr
								info.oImage->SetProperty, DATA=*info.currimage
								info.oView->SetProperty, VIEWPLANE_RECT=[0,0,info.vw_xsize,info.vw_ysize]

								WIDGET_CONTROL, display_labelID, SET_VALUE='Zoom Factor='+STRTRIM(STRING(info.imgzfac),2)

								info.vw_xloc = 0
								info.vw_yloc = 0
								info.padx    = offsets.x
								info.pady    = offsets.y

								info.xsize = xsize
								info.ysize = ysize
								info.initx = xsize
								info.inity = ysize
								info.priorzfac = info.imgzfac
								info.roizfac = 1

								; adjust field of view
								info.xfov_offset = x1 / currx * info.xfov
								info.xfov = xsize / currx * info.xfov

								info.yfov_offset = y1 / curry * info.yfov
								info.yfov = ysize / curry * info.yfov

							ENDIF
						ENDIF

						OBJ_DESTROY, oROI
						info.oROI = OBJ_NEW()

						info.oWindow->Draw, info.oView
					ENDIF
				END

				2: BEGIN ; scroll
					info.buttonDown=0b
				END

				3: BEGIN ; intensity
					info.buttonDown=0b
            	END

				5: BEGIN ; zoom

					IF (info.buttonDown NE 0) THEN BEGIN

						info.buttonDown=0b

						curr_zfac = info.roizfac

						x0 = info.buttonXY[0]
                		x1 = event.x

						dX = (x1 - x0) / ZOOM_SENSITIVITY + 1

						IF dX LT 0 THEN info.roizfac = -1 / CEIL(dX)  $
						ELSE info.roizfac = FLOOR(dX)

						info.roizfac = FLOOR(curr_zfac * info.roizfac)

						IF info.roizfac LT 1 THEN info.roizfac = 1
						IF info.roizfac GT ZOOM_MAX THEN info.roizfac = ZOOM_MAX

						; maintain aspect ratio, size of new image
						info.xsize = info.initx * info.roizfac
						info.ysize = info.inity * info.roizfac

						IF info.interleave EQ 0 THEN BEGIN
							new_image = Array_Zoom_2D((*info.currimage),curr_zfac,info.roizfac,0)
						ENDIF ELSE BEGIN
							sz_scan = size((*info.currimage))
							IF sz_scan[1] EQ 3 THEN BEGIN
								new_image = BYTARR(3,info.xsize, info.ysize)
								new_image[0,*,*] = Array_Zoom_2D(((*info.currimage)[0,*,*]),curr_zfac,info.roizfac,0)
								new_image[1,*,*] = Array_Zoom_2D(((*info.currimage)[1,*,*]),curr_zfac,info.roizfac,0)
								new_image[2,*,*] = Array_Zoom_2D(((*info.currimage)[2,*,*]),curr_zfac,info.roizfac,0)
							ENDIF ELSE BEGIN
								new_image = BYTARR(info.xsize, info.ysize, 3)
								new_image[*,*,0] = Array_Zoom_2D(((*info.currimage)[*,*,0]),curr_zfac,info.roizfac,0)
								new_image[*,*,1] = Array_Zoom_2D(((*info.currimage)[*,*,1]),curr_zfac,info.roizfac,0)
								new_image[*,*,2] = Array_Zoom_2D(((*info.currimage)[*,*,2]),curr_zfac,info.roizfac,0)

							ENDELSE
						ENDELSE

						info.imgzfac = info.priorzfac * info.roizfac

						IF info.fpdisp EQ 1 THEN BEGIN
							info.currfpvals = PTR_NEW(Array_Zoom_2D((*info.currfpvals),curr_zfac,info.roizfac,1))
						ENDIF

						offsets = Pad_Image(PTR_NEW(new_image), info.vw_xsize, info.vw_ysize)
						info.currimage = offsets.img_ptr
						info.oImage->SetProperty, DATA=*info.currimage
						info.oView->SetProperty, VIEWPLANE_RECT=[0,0,info.vw_xsize,info.vw_ysize]

						WIDGET_CONTROL, display_labelID, SET_VALUE='Zoom Factor='+STRTRIM(STRING(info.imgzfac),2)

						info.vw_xloc = 0
						info.vw_yloc = 0
						info.padx = offsets.x
						info.pady = offsets.y

						info.oWindow->Draw, info.oView

					ENDIF
				END

				ELSE: BEGIN
				END

			ENDCASE

			WIDGET_CONTROL, tabID, SET_UVALUE=info, /No_Copy
        END

        ; Motion
        2: BEGIN
			WIDGET_CONTROL, tabID, GET_UVALUE=info, /No_Copy

			CASE display_mouseStatus OF

				1: BEGIN ;roi
            		IF (info.buttonDown NE 0) THEN BEGIN

        				oROI = info.oROI

                		style_point = 0
                		style_line = 1
                		style_closed = 2

                		x0 = info.buttonXY[0]
                		y0 = info.buttonXY[1]
                		x1 = event.x + info.vw_xloc - info.padx
                		y1 = event.y + info.vw_xloc - info.pady

                		IF (x0 EQ x1) THEN BEGIN
                    		IF (y0 EQ y1) THEN BEGIN
                        		newBox = [[x0,y0,0.0]]
                        		style = style_point
                    		ENDIF ELSE BEGIN
                        		newBox = [[x0,y0,0.0], [x0,y1,0.0]]
                        		style = style_line
                    		ENDELSE
                		ENDIF ELSE IF (y0 EQ y1) THEN BEGIN
                    		newBox = [[x0,y0,0.0], [x1,y0,0.0]]
                    		style = style_line
                		ENDIF ELSE BEGIN
                    		newBox = [[x0,y0,0.0],[x1,y0,0.0],[x1,y1,0.0],[x0,y1,0.0]]
                    		style = style_closed
                		ENDELSE

                		oROI->GetProperty, N_VERTS=nVerts
               	 		oROI->ReplaceData, newBox, START=0, FINISH=nVerts-1
                		oROI->SetProperty, STYLE=style

						info.oWindow->Draw, info.oView
					ENDIF
				END

				2: BEGIN ;scroll
					IF (info.buttonDown NE 0) THEN BEGIN

                		x0 = info.buttonXY[0]
                		y0 = info.buttonXY[1]
                		x1 = event.x
                		y1 = event.y

						dX = (x1 - x0) / SCROLL_SENSITIVITY
						dY = (y1 - y0) / SCROLL_SENSITIVITY

						info.vw_xloc = info.vw_xloc - dX
						info.vw_yloc = info.vw_yloc - dY

						info.buttonXY = [x1, y1]

						info.oView->SetProperty, VIEWPLANE_RECT=[info.vw_xloc, info.vw_yloc, $
							info.vw_xsize, info.vw_ysize]
						info.oWindow->Draw, info.oView
					END
				END

				3: BEGIN ;intensity

					IF info.buttonDown NE 0 THEN BEGIN

						SENSITIVITY = 0.5

                		x0 = info.buttonXY[0]
                		y0 = info.buttonXY[1]
                		x1 = event.x
                		y1 = event.y

						dX = (x1 - x0) / SENSITIVITY
						dY = (y1 - y0) / SENSITIVITY

						val_Cen = info.int_Cen + dX
						IF val_Cen GT 255 THEN val_Cen = 255
						IF val_Cen LT 0 THEN val_Cen = 0

						val_Max = info.int_Max + dX - dY
						IF val_Max GT 255 THEN val_Max = 255
						IF val_Max LT 0 THEN val_Max = 0

						IF info.interleave EQ 0 THEN BEGIN
							scale_img = PTR_NEW(BYTSCL((*info.currimage), $
								MAX=val_Max, MIN=0, TOP=val_Cen))
						ENDIF ELSE BEGIN
							sz_scan = size((*info.currimage))
							scale_img = FLTARR(sz_scan[1],sz_scan[2],sz_scan[3])
							IF sz_scan[1] EQ 3 THEN BEGIN
								scale_img[0,*,*] = BYTSCL(((*info.currimage)[0,*,*]), MAX=val_Max, MIN=0, TOP=val_Cen)
								scale_img[1,*,*] = BYTSCL(((*info.currimage)[1,*,*]), MAX=val_Max, MIN=0, TOP=val_Cen)
								scale_img[2,*,*] = BYTSCL(((*info.currimage)[2,*,*]), MAX=val_Max, MIN=0, TOP=val_Cen)
							ENDIF ELSE BEGIN
								scale_img[*,*,0] = BYTSCL(((*info.currimage)[*,*,0]), MAX=val_Max, MIN=0, TOP=val_Cen)
								scale_img[*,*,1] = BYTSCL(((*info.currimage)[*,*,1]), MAX=val_Max, MIN=0, TOP=val_Cen)
								scale_img[*,*,2] = BYTSCL(((*info.currimage)[*,*,2]), MAX=val_Max, MIN=0, TOP=val_Cen)
							ENDELSE
							scale_img = PTR_NEW(scale_img)
						ENDELSE

						info.oImage->SetProperty, DATA=(*scale_img)
						info.oWindow->Draw, info.oView

					ENDIF
				END

				5: BEGIN ;zoom

					IF info.buttonDown EQ 1b THEN BEGIN

						curr_zfac = info.roizfac

						x0 = info.buttonXY[0]
                		x1 = event.x

						dX = (x1 - x0) / ZOOM_SENSITIVITY + 1

						IF dX LT 0 THEN new_zfac = -1 / CEIL(dX)  $
						ELSE new_zfac = FLOOR(dX)

						new_zfac = FLOOR(curr_zfac * new_zfac)

						IF new_zfac LT 1 THEN new_zfac = 1
						IF new_zfac GT ZOOM_MAX THEN new_zfac = ZOOM_MAX
					ENDIF ELSE new_zfac = info.roizfac
				END

				6: BEGIN ;rotate
					IF info.buttonDown EQ 1b THEN BEGIN

						curr_zfac = info.roizfac

						x0 = info.buttonXY[0]
                		x1 = event.x

						dX = (x1 - x0) / ZOOM_SENSITIVITY + 1

						IF dX LT 0 THEN new_zfac = -1 / CEIL(dX)  $
						ELSE new_zfac = FLOOR(dX)

						new_zfac = FLOOR(curr_zfac * new_zfac)

						IF new_zfac LT 1 THEN new_zfac = 1
						IF new_zfac GT ZOOM_MAX THEN new_zfac = ZOOM_MAX
					ENDIF ELSE new_zfac = info.roizfac
				END

				ELSE: BEGIN
				END

			ENDCASE

			;update status bar
			xcoor = event.x - info.padx + info.vw_xloc
			ycoor = event.y - info.pady + info.vw_yloc

			; check to make sure over image
			IF (xcoor LT info.xsize) AND (xcoor GE 0) AND (ycoor GE 0) AND $
				(ycoor LT info.ysize) THEN BEGIN
				; convert from pixels to mm
				xpos = info.xfov_offset + (FLOAT(xcoor+1) / info.xsize) * info.xfov
				ypos = info.yfov_offset + (FLOAT(ycoor+1) / info.ysize) * info.yfov

				; format output
				xstr = STRTRIM(STRING(FORMAT='(F0.3)', xpos),2)
				ystr = STRTRIM(STRING(FORMAT='(F0.3)', ypos),2)

				IF info.fpdisp EQ 1 THEN BEGIN
				    pstr = STRTRIM(STRING(FORMAT='(E0.3)', FLOAT((*info.currfpvals)[xcoor,ycoor])),2)
				    plab = ' FLOAT = '
				ENDIF ELSE BEGIN
					IF info.interleave EQ 1 THEN BEGIN
						pstr = ''
						plab = ''
					ENDIF ELSE BEGIN
						plab = ' PIXEL = '
						pstr = STRTRIM(STRING(FORMAT='(E0.3)', FLOAT((*info.currimage)[xcoor,ycoor])),2)
					ENDELSE
				END

				; if zooming update to new potential zoom
				IF display_mouseStatus EQ 5 THEN BEGIN
					zstr = STRTRIM(STRING(new_zfac),2)
				ENDIF ELSE zstr = STRTRIM(STRING(info.imgzfac),2)

				statusStr = 'XPOS =' + xstr + ' mm' + ' YPOS =' + ystr + ' mm' + $
					STRING(13b) + STRING(10b) + 'ZOOM =' + zstr + plab + pstr

				WIDGET_CONTROL, display_labelID, SET_VALUE=statusStr
			ENDIF

			WIDGET_CONTROL, tabID, SET_UVALUE=info, /No_Copy
        END

        ; Expose
        4: BEGIN
	        WIDGET_CONTROL, tabID, GET_UVALUE=info, /No_Copy
		    info.oWindow->Draw, info.oView
		    WIDGET_CONTROL, tabID, SET_UVALUE=info, /No_Copy
		END

		ELSE: BEGIN
		END

	ENDCASE
END


; Subroutine name: mas_display
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04
    ;Fix spelling mistakes and commenting.


pro mas_display
    COMPILE_OPT IDL2
    common scan_data
    CI = project.ci
    
    if (project.imndarray[ci].image_type eq 99) then begin
    
        mas_load_state_1
        
        adim_start = project.procPramArray[ci].adim_start
        ts_data = (*project.dataarray[ci].state1)[*,adim_start]
        
        title = '1D Spect: Spectrum #'+string(adim_start, format='(I0)')
        
        mas_display_1d_spect, ts_data, title=title
        
        return
     
    endif
    
    
    adim_start = project.procPramArray[ci].adim_start
    sdim_start = project.procPramArray[ci].sdim_start
    ffov = project.imndArray[ci].f_fov
    pfov = project.imndArray[ci].p_fov
    sfov = project.imndArray[ci].s_fov
    rotate_Direction = project.procpramArray[ci].rotate_Direction
    slice_axis = project.procpramArray[ci].slice_axis
    display_title = project.imndArray[CI].display_Name +' Slice #'+strtrim(sdim_start,2)

    ;before we display we have to make sure that the image is ready for display
    mas_load_state_2
    if (project.procpramarray[ci].state_2 eq 0) then return
    
    ;measure the size of the image are going to display
    displaySize = size(*project.dataArray[ci].state2)

    ;now that the image is loaded we have to set the size in the x and
    ;y direction correctly with respect to the rotation and
    ;the direction of the slice.
    if slice_axis eq 0 then begin
       x_fov = ffov
       y_fov = pfov
    end
    if slice_axis eq 1 then begin
       x_fov = ffov
       y_fov = sfov
    end
    if slice_axis eq 2 then begin
       x_fov = pfov
       y_fov = sfov
    end

    ;if the image is rotate 90 or 180 then swap the fov's
    if rotate_Direction eq 1 or rotate_Direction eq 3 then begin
       temp = x_fov
       x_fov = y_fov
       y_fov = temp
    end


    ;if the image is a 2d data set
    if displaySize[0] eq 2 then begin
        ;2dim data

       ;check to see if they want a simple display using window
       if project.iImage_flag eq 0 then begin
         image = (*project.dataArray[ci].state2)[*,*]
         fpvals = (*project.dataArray[ci].state2)[*,*]

         dataValueMin = project.procPramArray[ci].min_display_thresh
         dataValueMax = project.procPramArray[ci].max_display_thresh

         mas_windowing, image
         mas_display_multi, $
           image, fp_range=[dataValueMin, dataValueMax], $
           fov_x=x_fov, fov_y=y_fov, fov_units=1, $
           fp_vals=fpvals, tab_title=display_title, /standalone
;         Display_DrawResS_Create, image, display_title, ffov, pfov, FPVALS=fpvals

       end else begin
         ; if the user wants a complex display then do so.
         iImage, (*project.dataArray[ci].state2)$[*,*] $
                ,IMAGE_DIMENSIONS=[x_fov,y_fov] $
                ,xtitle='cm' $
                ,ytitle='cm' $
                ,title = display_title
        endelse

        update_status_bar, ' '
    endif

    if displaySize[0] eq 3 then begin

        totalFrames = displaySize[3]

        screen_dimensions = GET_SCREEN_SIZE(RESOLUTION=resolution)

        xSize = ceil(sqrt(totalFrames))
        ySize = ceil(float(totalFrames)/float(xSize))

        sheet = mas_make_sheet_display( project.dataArray[ci].state2,x_fov,y_fov )

        if project.iImage_flag eq 0 then begin

            mas_windowing , *sheet.image
            update_status_bar, ' '

            mas_display_multi, *sheet.image, tab_title='Multiple', $
                               fov_x=sheet.x_fov, fov_y=sheet.y_fov, $
                               fov_units=1, /standalone

;;            Display_DrawResS_Create, *sheet.image, 'Multiple', ffov, pfov

            ptr_free, sheet.image

        endif else begin

            iImage, *sheet.image $
                    ,IMAGE_DIMENSIONS=[sheet.x_fov,sheet.y_fov] $
                    ,xtitle='cm' $
                    ,ytitle='cm' $
                    ,title = display_title

            ptr_free, sheet.image

        end
    endif
end
