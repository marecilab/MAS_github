;; $Id$
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


; Subroutine name: Display_DrawResS_Resize
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

PRO Display_DrawResS_Resize, ev
    Widget_Control, ev.top, Get_UValue=info, /No_Copy

	curr_zfac = info.roizfac

	; use the window size to determine zoomfactor
;;	IF (ev.x GT info.xsize) AND (ev.y GT info.ysize) THEN BEGIN
	IF (ev.x - info.xsize gt 20) AND (ev.y - info.ysize gt 20) THEN BEGIN

		xval = CEIL(FLOAT(ev.x)/info.initx)
		yval = CEIL(FLOAT(ev.y)/info.inity)

		IF xval LT yval THEN info.roizfac = xval ELSE info.roizfac = yval

;;	ENDIF ELSE IF (ev.x LT info.xsize) AND (ev.y LT info.ysize) AND $
	ENDIF ELSE IF (info.xsize - ev.x gt 20) AND (info.ysize - ev.y gt 20) AND $
		(info.roizfac GT 1) THEN BEGIN

		xval = FLOOR(FLOAT(ev.x)/info.initx)
		yval = FLOOR(FLOAT(ev.y)/info.inity)

		IF xval GT yval THEN info.roizfac = xval ELSE info.roizfac = yval
		IF info.roizfac LT 1 THEN info.roizfac = 1
	END

	; maintain aspect ratio, size of new image
	info.xsize = info.initx * info.roizfac
	info.ysize = info.inity * info.roizfac

	IF info.interleave EQ 1 THEN BEGIN
		new_image = BYTARR(info.xsize, info.ysize, 3)
		new_image[*,*,0] = Array_Zoom_2D(((*info.currimage)[*,*,0]),curr_zfac,info.roizfac,0)
		new_image[*,*,1] = Array_Zoom_2D(((*info.currimage)[*,*,1]),curr_zfac,info.roizfac,0)
		new_image[*,*,2] = Array_Zoom_2D(((*info.currimage)[*,*,2]),curr_zfac,info.roizfac,0)
	ENDIF ELSE BEGIN
		new_image = Array_Zoom_2D((*info.currimage),curr_zfac,info.roizfac,0)
	ENDELSE

	Widget_Control, info.drawID, Draw_XSize=info.xsize, Draw_YSize=info.ysize

	info.oView->SetProperty, VIEWPLANE_RECT=[0,0,info.xsize,info.ysize]
	info.oImage->SetProperty, DATA=new_image

	info.currimage = PTR_NEW(new_image)
	info.imgzfac = info.priorzfac * info.roizfac

	IF info.fpdisp EQ 1 THEN BEGIN
		info.currfpvals = PTR_NEW(Array_Zoom_2D((*info.currfpvals),curr_zfac,info.roizfac,1))
	ENDIF

    Widget_Control, ev.top, Set_UValue=info, /No_Copy
END


; Subroutine name: Display_DrawResS_Cleanup
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

PRO Display_DrawResS_Cleanup, tlb
    Widget_Control, tlb, Get_UValue=info, /No_Copy
    IF N_Elements(info) EQ 0 THEN RETURN

    ;COMPILE_OPT IDL2
    ;common scan_data
    ;ci = project.ci

	; apply mask
	; project.dataArray [ni].state1 = ptr_new(imgdata)
	; add flag to mas_load_state_1 to act like DICOM load

;	project.procPramArray[ci].state_1 = 0
;	project.procPramArray[project.ci].state_2 = 0

;	mas_load_state_2

;    	project.imndArray[ni].fdim = dicom.fdim
;    	project.imndArray[ni].pdim = dicom.pdim
;    	project.imndArray[ni].sdim = dicom.sdim
;    	project.imndArray[ni].adim = dicom.adim

;    	project.imndArray[ni].f_fov = dicom.f_fov
;    	project.imndArray[ni].p_fov = dicom.p_fov
;    	project.imndArray[ni].s_fov = dicom.s_fov

;    	project.imndArray[ni].slices = dicom.slices

    Ptr_Free, info.image
    Ptr_Free, info.currimage
    OBJ_DESTROY, info.oImage
    OBJ_DESTROY, info.oView
    OBJ_DESTROY, info.oModel
END


; Subroutine name: Display_DrawResS_Create
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

PRO Display_DrawResS_Create, image, display_title, ffov, pfov, FPVALS=fpvals

;    mas_display_multi, image, /standalone
;    return

    IF KEYWORD_SET(fpvals) NE 1 THEN BEGIN
        fpvals = 0
        fpdisp = 0
        print, 'WARNING : FLOATING POINT DISPLAY NOT SET'
    ENDIF ELSE fpdisp = 1

    sz_scan = size(image)

    IF sz_scan[0] EQ 3 THEN BEGIN
        oImage = OBJ_NEW('IDLgrImage',image)
        oImage->SetProperty, INTERLEAVE=2
        interleave = 1
    ENDIF ELSE BEGIN
        oImage = OBJ_NEW('IDLgrImage', image)
        interleave = 0
    ENDELSE

    screen_dimensions = GET_SCREEN_SIZE(RESOLUTION=resolution)
    baseID = Widget_Base(Column=1, Title=display_title, TLB_Size_Events=1)
    
    drawID = Widget_Draw(baseID, XSIZE=sz_scan[1], YSIZE=sz_scan[2], $
                   ;;      x_scroll_size = min([sz_scan[1], fix(screen_dimensions[0]*.85)]), $
                   ;;      y_scroll_size = min([sz_scan[2], fix(screen_dimensions[1]*.85)]), $
                         BUTTON_EVENTS=1, EXPOSE_EVENTS=1, MOTION_EVENTS=1, GRAPHICS_LEVEL=2, $
                         /ALIGN_CENTER);, /scroll)

    Widget_Control, baseID, /Realize
    Widget_Control, drawID, Get_Value=oWindow
    
    oView  = OBJ_NEW('IDLgrView', VIEWPLANE_RECT=[0,0,sz_scan[1],sz_scan[2]])
    oModel = OBJ_NEW('IDLgrModel')

    oView->Add, oModel
    oModel->Add, oImage
    
    oStatus = Widget_Label(baseID, SCR_XSIZE=200 ,SCR_YSIZE=30, $
                           FRAME=0 ,/ALIGN_LEFT ,VALUE='Zoom Factor=1')
    
    oWindow->Draw, oView
    
    xfov = pfov
    yfov = ffov
    
    info = {image:Ptr_New(image), oWindow:oWindow, drawID:drawID, $
            xsize:sz_scan[1], ysize:sz_scan[2], oView:oView, oModel:oModel, $
            oROI:OBJ_NEW(), lbuttonXY: LONARR(2), lbButtonDown: 0B, $
            oImage:oImage, currimage:PTR_NEW(image), initx:sz_scan[1], $
            inity:sz_scan[2], roizfac:1, imgzfac:1, priorzfac:1, $
            oStatus:oStatus, xfov:xfov, yfov:yfov, xfov_init:xfov, yfov_init:yfov, $
            xfov_offset:0.0, yfov_offset:0.0, rbuttonXY:LONARR(2), rbButtonDown:0B, $
            int_Cen:MAX(image), int_Max:MAX(image), int_Min:0, fpvals:PTR_NEW(fpvals), fpdisp:fpdisp, $
            currfpvals:PTR_NEW(fpvals), interleave:interleave }
    
    Widget_Control, baseID, Set_UValue=info, /No_Copy
    
    XMANAGER, 'Display_DrawResS_Create', baseID, Event_Handler='Display_DrawResS_Resize', Cleanup='Display_DrawResS_Cleanup',/NO_BLOCK
    XMANAGER, 'Display_DrawResS_Create', drawID, Event_Handler='DisplayS_Zoom_Event', /NO_BLOCK
END

; Subroutine name: Display_Zoom_Event
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:


pro DisplayS_Zoom_Event, event

    common scan_data

    ;; reset image on double click
    IF event.clicks EQ 2 THEN BEGIN
        
        WIDGET_CONTROL, event.top, GET_UVALUE=info, /No_Copy
        
        info.int_Cen = MAX(*info.currimage)
        info.int_Max = info.int_Cen
        info.int_Min = 0
        
        sz_scan = size((*info.image))
        
        xsize = sz_scan[1]
        ysize = sz_scan[2]
        
        Widget_Control, info.drawID, Draw_XSize=xsize, Draw_YSize=ysize
        info.oView->SetProperty, VIEWPLANE_RECT=[0,0,xsize,ysize]
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
        
        info.oWindow->Draw, info.oView
        
        WIDGET_CONTROL, info.oStatus, SET_VALUE='Zoom Factor=1'
        WIDGET_CONTROL, event.top, SET_UVALUE=info, /No_Copy
        
        RETURN
    END

    CASE event.type OF
        ;; Button Press
        0: BEGIN
            WIDGET_CONTROL, event.top, GET_UVALUE=info, /No_Copy
            
            IF event.press EQ 1 THEN BEGIN
                ;; Process Left Mouse Button
                oROI = OBJ_NEW('IDLgrROI',COLOR=BYTE([255, 0, 0]),STYLE=0)
                oROI->AppendData, [event.x, event.y, 0]
                
                info.oModel->Add, oROI
                
                info.oWindow->Draw, info.oView
                
                info.oROI = oROI
                info.lbButtonDown = 1b
                info.lbuttonXY = [event.x, event.y]
                
            ENDIF ELSE BEGIN
                ;; Process Right Mouse Button
                IF info.interleave EQ 0 THEN BEGIN
                    ;; set the initial sliders
                    
                    info.rbButtonDown = 1b
                    info.rbuttonXY = [event.x, event.y]
                ENDIF
            ENDELSE
            
            WIDGET_CONTROL, event.top, SET_UVALUE=info, /No_Copy
        END
        
        ;; Button Release
        1: BEGIN
            WIDGET_CONTROL, event.top, GET_UVALUE=info, /No_Copy
            
            IF (info.lbButtonDown NE 0) THEN BEGIN
                
                info.lbButtonDown=0b
                
                oROI = info.oROI
                
                oROI->GetProperty, DATA=roiData
                
                ;; check to make sure rectangale has four sides
                IF ((N_ELEMENTS(roiData)/3) EQ 4) THEN BEGIN
                    
                    x1 = MIN(roiData[0,*])
                    y1 = MIN(roiData[1,*])
                    
                    x2 = MAX(roiData[0,*])
                    y2 = MAX(roiData[1,*])
                    
                    currx = info.xsize
                    curry = info.ysize
                    
                    ;; check to make sure all coordinates of rectangle are on the image
                    IF (x1 LT currx) AND (x1 GE 0) AND (x2 LT currx) AND (x2 GE 0) AND $
                      (y1 LT curry) AND (y1 GE 0) AND (y2 LT curry) AND (y2 GE 0) $
                      THEN BEGIN
                        
                        ;; check to see that that the selection includes all of
                        ;; zoom voxel, if not adjust
                        
                        x1 = x1 - (x1 MOD info.imgzfac)
                        y1 = y1 - (y1 MOD info.imgzfac)
                        
                        x2 = x2 + (info.imgzfac - (x2 MOD info.imgzfac)) - 1
                        y2 = y2 + (info.imgzfac - (y2 MOD info.imgzfac)) - 1
                        
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
                        
                        WIDGET_CONTROL, info.drawID, Draw_XSize=xsize, Draw_YSize=ysize
                        info.oView->SetProperty, VIEWPLANE_RECT=[0,0,xsize,ysize]
                        info.oImage->SetProperty, DATA=new_image
                        
                        WIDGET_CONTROL, info.oStatus, SET_VALUE='Zoom Factor='+STRTRIM(STRING(info.imgzfac),2)
                        
                        info.currimage = PTR_NEW(new_image)
                        info.xsize = xsize
                        info.ysize = ysize
                        info.initx = xsize
                        info.inity = ysize
                        info.priorzfac = info.imgzfac
                        info.roizfac = 1
                        
                        IF info.fpdisp EQ 1 THEN BEGIN
                            info.currfpvals = PTR_NEW(((*info.currfpvals)[x1:x2,y1:y2]))
                        ENDIF
                        
                        ;; adjust field of view
                        info.xfov_offset = x1 / currx * info.xfov
                        info.xfov = xsize / currx * info.xfov
                        
                        info.yfov_offset = y1 / curry * info.yfov
                        info.yfov = ysize / curry * info.yfov
                        
                    END
                END
                
                OBJ_DESTROY, oROI
                info.oROI = OBJ_NEW()
                
            ENDIF ELSE IF (info.rbButtonDown NE 0) AND (info.interleave EQ 0) THEN BEGIN
                info.rbButtonDown=0b
            ENDIF
            
            WIDGET_CONTROL, event.top, SET_UVALUE=info, /No_Copy
        END
        
        ;; Motion
        2: BEGIN
            WIDGET_CONTROL, event.top, GET_UVALUE=info, /No_Copy
            
            ;; check to make sure over image
            IF (event.x LT info.xsize) AND (event.x GE 0) AND (event.y GE 0) AND $
              (event.y LT info.ysize) THEN BEGIN
                ;; convert from pixels to mm
                ;; xpos = info.xfov_offset + (FLOAT(event.x+1) / info.xsize) * info.xfov
                ;; ypos = info.yfov_offset + (FLOAT(event.y+1) / info.ysize) * info.yfov
                
                xpos = event.x
                ypos = event.y
                
                ;; format output
                xstr = STRTRIM(STRING(FORMAT='(F0.3)', xpos),2)
                ystr = STRTRIM(STRING(FORMAT='(F0.3)', ypos),2)
                
                IF info.fpdisp EQ 1 THEN BEGIN
                    pstr = STRTRIM(STRING(FORMAT='(E0.3)', FLOAT((*info.currfpvals)[event.x,event.y])),2)
                    plab = ' FLOAT = '
                ENDIF ELSE BEGIN
                    IF info.interleave EQ 1 THEN BEGIN
                        pstr = ''
                        plab = ''
                    ENDIF ELSE BEGIN
                        plab = ' PIXEL = '
                        pstr = STRTRIM(STRING(FORMAT='(E0.3)', FLOAT((*info.currimage)[event.x,event.y])),2)
                    ENDELSE
                END
                
                zstr = STRTRIM(STRING(info.imgzfac),2)
                
                statusStr = 'XPOS =' + xstr + ' mm' + ' YPOS =' + ystr + ' mm' + $
                  STRING(13b) + STRING(10b) + 'ZOOM =' + zstr + plab + pstr
                
                WIDGET_CONTROL, info.oStatus, SET_VALUE=statusStr
            ENDIF
            
            IF (info.lbButtonDown NE 0) THEN BEGIN
                
                oROI = info.oROI
                
                style_point = 0
                style_line = 1
                style_closed = 2
                
                x0 = info.lbuttonXY[0]
                y0 = info.lbuttonXY[1]
                x1 = event.x
                y1 = event.y
                
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
                
            ENDIF ELSE IF (info.rbButtonDown NE 0) AND (info.interleave EQ 0) THEN BEGIN
                
                SENSITIVITY = 50
                
                x0 = info.rbuttonXY[0]
                y0 = info.rbuttonXY[1]
                x1 = event.x
                y1 = event.y
                
                dX = (x1 - x0) / SENSITIVITY
                dY = (y1 - y0) / SENSITIVITY
                
                val_Cen = info.int_Cen + dX
                IF val_Cen GT 255 THEN val_Cen = 255
                IF val_Cen LT 1   THEN val_Cen = 1
                
                val_Max = info.int_Max - dY
                IF val_Max GT 255 THEN val_Max = 255
                IF val_Max LT 1   THEN val_Max = 1
                
                IF val_Max NE info.int_Max OR val_Cen NE info.int_Cen THEN BEGIN
                    ;; update image
                    scale_img = PTR_NEW(BYTSCL((*info.currimage), $
                                               MAX=val_Max, MIN=0, TOP=val_Cen))
                    
                    info.int_Cen = val_Cen
                    info.int_Max = val_Max
                    
                    info.oImage->SetProperty, DATA=(*scale_img)
                    info.oWindow->Draw, info.oView
                ENDIF
                
            ENDIF
            
            WIDGET_CONTROL, event.top, SET_UVALUE=info, /No_Copy
        END
        
        ;; Expose
        4: BEGIN
            WIDGET_CONTROL, event.top, GET_UVALUE=info, /No_Copy
            info.oWindow->Draw, info.oView
            WIDGET_CONTROL, event.top, SET_UVALUE=info, /No_Copy
        END
        ELSE: BEGIN
        END
    ENDCASE
END
