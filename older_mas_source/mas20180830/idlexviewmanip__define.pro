; $Id: //depot/idl/IDL_64/idldir/lib/utilities/idlexviewmanip__define.pro#1 $
;
; Copyright (c) 1997-2007, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;+
; NAME:
;   IDLexViewManip
;
; PURPOSE:
;   This object serves to simplify dynamic manipulation of the
;   position and size of an IDLgrView object within a destination
;   object.
;
; CATEGORY:
;   Object graphics examples.
;
; CALLING SEQUENCE:
;
;   oVm = OBJ_NEW('IDLexViewManip')
;
;   IDLexViewManip::SetProperty
;   IDLexViewManip::GetProperty
;   IDLexViewManip::MouseUp, mouseXY, oDest
;   IDLexViewManip::MouseTrack, mouseXY, oDest
;   IDLexViewManip::MouseDown, mouseXY, oDest
;   IDLexViewManip::Reshape
;   IDLexViewManip::SetTarget, oView, oDest
;
; KEYWORD PARAMETERS:
;   IDLexViewManip::Init:
;
;
; MODIFICATION HISTORY:
;    Written by:   RJF, Nov 1996
;-
;----------------------------------------------------------------------------
FUNCTION IDLexViewManip::Init, $
    SELECTOR_COLOR=selector_color, $
    MANIPULATOR_COLOR=manipulator_color, $
    MANIP_SHOW_MANIP=manip_show_manip, $
    MANIP_SHOW_SEL=manip_show_sel, $
    MODE=mode, $
    _EXTRA=e

    IF self->IDLgrView::Init(transparent=1, _EXTRA=e) NE 1 THEN $
        RETURN, 0

    ; Create the local tree.
    self.oModel = OBJ_NEW('IDLgrModel')
    self->Add, self.oModel

    ; Create the manipulation and selection polyline objects
    self.oSelPoly = obj_new('IDLgrPolyline')
    self.oModel->Add, self.oSelPoly
    self.oManipPoly = obj_new('IDLgrPolyline')
    self.oModel->Add,self.oManipPoly

    ; Initialize default state.
    self.manip_show_manip = 1b
    self.manip_show_sel = 1b
    if N_ELEMENTS(selector_color) EQ 0 THEN $
        selector_color=[128,128,128]
    if N_ELEMENTS(manipulator_color) EQ 0 THEN $
        manipulator_color=[128,128,128]
    self->SetProperty, $
        SELECTOR_COLOR=selector_color, $
        MANIPULATOR_COLOR=manipulator_color, $
        MANIP_SHOW_MANIP=manip_show_manip, $
        MANIP_SHOW_SEL=manip_show_sel, $
        MODE=mode

    RETURN, 1
END

;----------------------------------------------------------------------------
;
;Create the selection visual polylines
;
PRO IDLexViewManip::Build_Visuals

    self->CalcRect

    self->GetProperty, ZCLIP=zclip
    tiny_amount = (zclip[0] - zclip[1]) * .001 ; Somewhat arbitrary.

    ; Selection visual first (very simple).
    v = DBLARR(3,4)
    v[0,0] = [self.target_rect[0], self.target_rect[1]]
    v[0,1] = [ $
        self.target_rect[0] + self.target_rect[2]-1, $
        self.target_rect[1] $
        ]
    v[0,2] = [ $
        self.target_rect[0] + self.target_rect[2]-1, $
        self.target_rect[1] + self.target_rect[3]-1 $
        ]
    v[0,3] = [ $
        self.target_rect[0], $
        self.target_rect[1] + self.target_rect[3]-1 $
        ]
    v[2,*] = zclip[0] - tiny_amount
    self.oSelPoly->SetProperty, DATA=v, POLYLINES=[5,0,1,2,3,0]

    ; Manipulation visual.
    inset = 2
    inset2 = 6
    CASE self.mode OF
        0: BEGIN   ;corners
            IF OBJ_ISA(self.oDest, 'IDLgrWindow') THEN $
                self.oDest->IDLgrWindow::SetCurrentCursor, "ORIGINAL"

            verts = FLTARR(3,12) ; Corners.
            verts[0,0] = v[0,0] + 2*inset2
            verts[1,0] = v[1,0]
            verts[0,1] = v[0,0] + 2*inset2
            verts[1,1] = v[1,0] + 2*inset2
            verts[0,2] = v[0,0]
            verts[1,2] = v[1,0] + 2*inset2

            verts[0,3] = v[0,1] - 2*inset2
            verts[1,3] = v[1,1]
            verts[0,4] = v[0,1] - 2*inset2
            verts[1,4] = v[1,1] + 2*inset2
            verts[0,5] = v[0,1]
            verts[1,5] = v[1,1] + 2*inset2

            verts[0,6] = v[0,2] - 2*inset2
            verts[1,6] = v[1,2]
            verts[0,7] = v[0,2] - 2*inset2
            verts[1,7] = v[1,2] - 2*inset2
            verts[0,8] = v[0,2]
            verts[1,8] = v[1,2] - 2*inset2

            verts[0,9] = v[0,3] + 2*inset2
            verts[1,9] = v[1,3]
            verts[0,10] = v[0,3] + 2*inset2
            verts[1,10] = v[1,3] - 2*inset2
            verts[0,11] = v[0,3]
            verts[1,11] = v[1,3] - 2*inset2

            conn = [3,0,1,2, 3,3,4,5, 3,6,7,8, 3,9,10,11]
           END
        1: BEGIN   ;translate
            IF OBJ_ISA(self.oDest, 'IDLgrWindow') THEN $
                self.oDest->IDLgrWindow::SetCurrentCursor, "MOVE"

            verts = fltarr(3,12) ; Corners.
            verts[0,0] = v[0,0] + 2*inset2
            verts[1,0] = v[1,0]
            verts[0,1] = v[0,0] + 2*inset2
            verts[1,1] = v[1,0] + 2*inset2
            verts[0,2] = v[0,0]
            verts[1,2] = v[1,0] + 2*inset2

            verts[0,3] = v[0,1] - 2*inset2
            verts[1,3] = v[1,1]
            verts[0,4] = v[0,1] - 2*inset2
            verts[1,4] = v[1,1] + 2*inset2
            verts[0,5] = v[0,1]
            verts[1,5] = v[1,1] + 2*inset2

            verts[0,6] = v[0,2] - 2*inset2
            verts[1,6] = v[1,2]
            verts[0,7] = v[0,2] - 2*inset2
            verts[1,7] = v[1,2] - 2*inset2
            verts[0,8] = v[0,2]
            verts[1,8] = v[1,2] - 2*inset2

            verts[0,9] = v[0,3] + 2*inset2
            verts[1,9] = v[1,3]
            verts[0,10] = v[0,3] + 2*inset2
            verts[1,10] = v[1,3] - 2*inset2
            verts[0,11] = v[0,3]
            verts[1,11] = v[1,3] - 2*inset2

            conn = [3,0,1,2, 3,3,4,5, 3,6,7,8, 3,9,10,11]
           END
        2: BEGIN   ;UpperLeft
            IF OBJ_ISA(self.oDest, 'IDLgrWindow') THEN $
                self.oDest->IDLgrWindow::SetCurrentCursor, "SIZE_NW"

            verts = fltarr(3,15) ; Corners.
            verts[0,0] = v[0,0] + 2*inset2
            verts[1,0] = v[1,0]
            verts[0,1] = v[0,0] + 2*inset2
            verts[1,1] = v[1,0] + 2*inset2
            verts[0,2] = v[0,0]
            verts[1,2] = v[1,0] + 2*inset2

            verts[0,3] = v[0,1] - 2*inset2
            verts[1,3] = v[1,1]
            verts[0,4] = v[0,1] - 2*inset2
            verts[1,4] = v[1,1] + 2*inset2
            verts[0,5] = v[0,1]
            verts[1,5] = v[1,1] + 2*inset2

            verts[0,6] = v[0,2] - 2*inset2
            verts[1,6] = v[1,2]
            verts[0,7] = v[0,2] - 2*inset2
            verts[1,7] = v[1,2] - 2*inset2
            verts[0,8] = v[0,2]
            verts[1,8] = v[1,2] - 2*inset2

            verts[0,9] = v[0,3] + inset
            verts[1,9] = v[1,3] - inset
            verts[0,10] = verts[0,9] + 2*inset2
            verts[1,10] = verts[1,9] - 2*inset2
            verts[0,11] = verts[0,9]
            verts[1,11] = verts[1,9] - inset2
            verts[0,12] = verts[0,9] + inset2
            verts[1,12] = verts[1,9]
            verts[0,13] = verts[0,10]
            verts[1,13] = verts[1,10] + inset2
            verts[0,14] = verts[0,10] - inset2
            verts[1,14] = verts[1,10]
            conn = [3,0,1,2, 3,3,4,5, 3,6,7,8, $
                    2,9,10, 3,11,9,12, 3,13,10,14]
           END
        3: BEGIN   ;UpperRight
            IF OBJ_ISA(self.oDest, 'IDLgrWindow') THEN $
                self.oDest->IDLgrWindow::SetCurrentCursor, "SIZE_NE"

            verts = fltarr(3,15) ; Corners.
            verts[0,0] = v[0,0] + 2*inset2
            verts[1,0] = v[1,0]
            verts[0,1] = v[0,0] + 2*inset2
            verts[1,1] = v[1,0] + 2*inset2
            verts[0,2] = v[0,0]
            verts[1,2] = v[1,0] + 2*inset2

            verts[0,3] = v[0,1] - 2*inset2
            verts[1,3] = v[1,1]
            verts[0,4] = v[0,1] - 2*inset2
            verts[1,4] = v[1,1] + 2*inset2
            verts[0,5] = v[0,1]
            verts[1,5] = v[1,1] + 2*inset2

            verts[0,6] = v[0,3] + 2*inset2
            verts[1,6] = v[1,3]
            verts[0,7] = v[0,3] + 2*inset2
            verts[1,7] = v[1,3] - 2*inset2
            verts[0,8] = v[0,3]
            verts[1,8] = v[1,3] - 2*inset2

            verts[0,9] = v[0,2] - inset
            verts[1,9] = v[1,2] - inset
            verts[0,10] = verts[0,9] - 2*inset2
            verts[1,10] = verts[1,9] - 2*inset2
            verts[0,11] = verts[0,9]
            verts[1,11] = verts[1,9] - inset2
            verts[0,12] = verts[0,9] - inset2
            verts[1,12] = verts[1,9]
            verts[0,13] = verts[0,10]
            verts[1,13] = verts[1,10] + inset2
            verts[0,14] = verts[0,10] + inset2
            verts[1,14] = verts[1,10]

            conn = [3,0,1,2, 3,3,4,5, 3,6,7,8, $
                    2,9,10, 3,11,9,12, 3,13,10,14]
           END
        4: BEGIN   ;LowerRight
            IF OBJ_ISA(self.oDest, 'IDLgrWindow') THEN $
                self.oDest->IDLgrWindow::SetCurrentCursor, "SIZE_SE"

            verts = fltarr(3,15) ; Corners.
            verts[0,0] = v[0,0] + 2*inset2
            verts[1,0] = v[1,0]
            verts[0,1] = v[0,0] + 2*inset2
            verts[1,1] = v[1,0] + 2*inset2
            verts[0,2] = v[0,0]
            verts[1,2] = v[1,0] + 2*inset2

            verts[0,3] = v[0,2] - 2*inset2
            verts[1,3] = v[1,2]
            verts[0,4] = v[0,2] - 2*inset2
            verts[1,4] = v[1,2] - 2*inset2
            verts[0,5] = v[0,2]
            verts[1,5] = v[1,2] - 2*inset2

            verts[0,6] = v[0,3] + 2*inset2
            verts[1,6] = v[1,3]
            verts[0,7] = v[0,3] + 2*inset2
            verts[1,7] = v[1,3] - 2*inset2
            verts[0,8] = v[0,3]
            verts[1,8] = v[1,3] - 2*inset2

            verts[0,9] = v[0,1] - inset
            verts[1,9] = v[1,1] + inset
            verts[0,10] = verts[0,9] - 2*inset2
            verts[1,10] = verts[1,9] + 2*inset2
            verts[0,11] = verts[0,9]
            verts[1,11] = verts[1,9] + inset2
            verts[0,12] = verts[0,9] - inset2
            verts[1,12] = verts[1,9]
            verts[0,13] = verts[0,10]
            verts[1,13] = verts[1,10] - inset2
            verts[0,14] = verts[0,10] + inset2
            verts[1,14] = verts[1,10]

            conn = [3,0,1,2, 3,3,4,5, 3,6,7,8, $
                    2,9,10, 3,11,9,12, 3,13,10,14]

           END
        5: BEGIN   ;LowerLeft
            IF OBJ_ISA(self.oDest, 'IDLgrWindow') THEN $
                self.oDest->IDLgrWindow::SetCurrentCursor, "SIZE_SW"

            verts = fltarr(3,15) ; Corners.

            verts[0,0] = v[0,1] - 2*inset2
            verts[1,0] = v[1,1]
            verts[0,1] = v[0,1] - 2*inset2
            verts[1,1] = v[1,1] + 2*inset2
            verts[0,2] = v[0,1]
            verts[1,2] = v[1,1] + 2*inset2

            verts[0,3] = v[0,2] - 2*inset2
            verts[1,3] = v[1,2]
            verts[0,4] = v[0,2] - 2*inset2
            verts[1,4] = v[1,2] - 2*inset2
            verts[0,5] = v[0,2]
            verts[1,5] = v[1,2] - 2*inset2

            verts[0,6] = v[0,3] + 2*inset2
            verts[1,6] = v[1,3]
            verts[0,7] = v[0,3] + 2*inset2
            verts[1,7] = v[1,3] - 2*inset2
            verts[0,8] = v[0,3]
            verts[1,8] = v[1,3] - 2*inset2

            verts[0,9] = v[0,0] + inset
            verts[1,9] = v[1,0] + inset
            verts[0,10] = verts[0,9] + 2*inset2
            verts[1,10] = verts[1,9] + 2*inset2
            verts[0,11] = verts[0,9]
            verts[1,11] = verts[1,9] + inset2
            verts[0,12] = verts[0,9] + inset2
            verts[1,12] = verts[1,9]
            verts[0,13] = verts[0,10]
            verts[1,13] = verts[1,10] - inset2
            verts[0,14] = verts[0,10] - inset2
            verts[1,14] = verts[1,10]

            conn = [3,0,1,2, 3,3,4,5, 3,6,7,8, $
                    2,9,10, 3,11,9,12, 3,13,10,14]
           END
        6: BEGIN    ;Rubber Line
            IF OBJ_ISA(self.oDest, 'IDLgrWindow') THEN $
                self.oDest->IDLgrWindow::SetCurrentCursor, "CROSSHAIR"
            verts = dblarr(3,2)
            verts[0:1,0] = self.target_rect[0:1] + self.target_rect[2:3] / 2
            verts[0:1,1] = self.mouseXY

            conn = [2,0,1]
           END
    ENDCASE
    verts[2,*] = zclip[0] - tiny_amount
    self.oManipPoly->SetProperty, DATA=verts, POLYLINES=conn
END

;----------------------------------------------------------------------------
;Destroy the object
;
PRO IDLexViewManip::Cleanup

    ; Remove from tree if needed
    IF (self.oTarget NE OBJ_NEW()) THEN BEGIN
        self->GetProperty, PARENT=myparent
    IF (OBJ_VALID(myparent)) THEN BEGIN
        list = myparent->Get(/ALL)
        IF (SIZE(list,/TYPE) EQ 11) THEN  $
                myparent->Remove, self
    END
        self.oTarget = OBJ_NEW()
        self.oDest = OBJ_NEW()
    ENDIF

    self->IDLgrView::Cleanup
END

;----------------------------------------------------------------------------
;Override some properties
;
PRO IDLexViewManip::SetProperty, $
    SELECTOR_COLOR=selector_color, $
    MANIPULATOR_COLOR=manipulator_color, $
    MANIP_SHOW_MANIP=manip_show_manip, $
    MANIP_SHOW_SEL=manip_show_sel, $
    mode=mode, $
    _EXTRA=e

    if n_elements(selector_color) gt 0 then $
        self.oSelPoly->SetProperty, color=selector_color

    if n_elements(manipulator_color) gt 0 then $
        self.oManipPoly->SetProperty, color=manipulator_color

    if n_elements(mode) gt 0 then begin
        self.mode = mode eq 1 ? 6 : 0
        end
    self.oManipPoly->SetProperty, hide=self.mode eq 6

    if N_ELEMENTS(manip_show_manip) gt 0 then begin
        self.oManipPoly->SetProperty, HIDE=keyword_set(manip_show_manip) eq 0
        self.manip_show_manip = manip_show_manip
        end

    if N_ELEMENTS(manip_show_sel) gt 0 then begin
        self.oSelPoly->SetProperty, HIDE=keyword_set(manip_show_sel) eq 0
        self.manip_show_sel = manip_show_sel
        end

    if keyword_set(manip_show_manip) or keyword_set(manip_show_sel) then $
        self->build_visuals

    self->IDLgrView::SetProperty, _EXTRA=e
END
;----------------------------------------------------------------------------
PRO IDLexViewManip::GetProperty, TARGET=target, ALL=all, _REF_EXTRA=e

    target = self.oTarget
    self->IDLgrView::GetProperty, _EXTRA=e

    IF ARG_PRESENT(all) NE 0 THEN BEGIN
        self->IDLgrView::GetProperty, ALL=all
        all = CREATE_STRUCT( $
            all, $
            'target', target $
            )
    END
END

;----------------------------------------------------------------------------
;Set the current target of manipulation
;
PRO IDLexViewManip::SetTarget, oView, oDest
    ON_ERROR, 2 ; Return to caller on error.

    IF N_ELEMENTS(oView) EQ 0 THEN $
        MESSAGE, 'First argument is undefined.'

    IF oView NE OBJ_NEW() THEN BEGIN
        IF NOT OBJ_ISA(oView, 'IDLgrView') THEN $
            MESSAGE, 'First argument must be a reference to ' + $
                'an IDLgrView or null.'
        IF N_ELEMENTS(oDest) LE 0 THEN begin
            MESSAGE, 'Second argument is undefined.
            end
        IF NOT OBJ_ISA(oDest, 'IDLgrSrcDest') THEN $
            MESSAGE, 'Second argument must be a reference to ' + $
                'an IDLgrSrcDest.'
    END

    ; remove self from old target
    IF (OBJ_VALID(self.oTarget)) THEN BEGIN
        self->GetProperty, PARENT=myparent
        myparent->Remove, self
        self.oTarget = OBJ_NEW()
        self.oDest = OBJ_NEW()
    ENDIF

    IF oView NE OBJ_NEW() THEN BEGIN
        self.oDest = oDest
        ; Paste us in...
        oView->GetProperty, PARENT=parent, UNITS=units
        parent->Add, self
        self.oTarget = oView
        self.targetUnits = units

        if self.manip_show_sel or self.manip_show_manip then begin
            self->build_visuals
            end
    END
END

;----------------------------------------------------------------------------
;Method to be called when the graphics destination changes size
;
PRO IDLexViewManip::Reshape

    IF (NOT OBJ_VALID(self.oTarget)) THEN RETURN

    ; Build selection visual
    self->build_visuals

END


;----------------------------------------------------------------------------
;Compute the size of the selection visual rectangle
;Also, adjust for any possible change in the size of the destination.
;
PRO IDLexViewManip::CalcRect

    IF (NOT OBJ_VALID(self.oTarget)) THEN $
        RETURN

    ; Resize myself to fill the destination at pixel res
    self.oDest->GetProperty, DIMENSIONS=dest_dim, RESOLUTION=res
    self->IDLgrView::SetProperty, LOCATION=[0,0], DIMENSIONS=dest_dim, $
                  VIEWPLANE_RECT=[0,0,dest_dim[0],dest_dim[1]], PROJECTION=1

    ; Get the target's Rect
    self.oTarget->GetProperty, LOCATION=view_loc,$
                                          DIMENSIONS=view_dim

    ; Convert to device units.
    CASE self.targetUnits OF
        0: BEGIN ; Device
               ; Do nothing.
           END
        1: BEGIN ; Inches
               view_loc = view_loc * 2.54 / resolution
               view_dim = view_dim * 2.54 / resolution
           END
        2: BEGIN ; Centimeters
               view_loc = view_loc / resolution
               view_dim = view_dim / resolution
           END
        3: BEGIN ; Normalized
               view_loc = view_loc * dest_dim
               view_dim = view_dim * dest_dim
           END
    ENDCASE

    IF ((view_dim[0] EQ 0) OR (view_dim[1] EQ 0)) THEN $
        view_dim = dest_dim

    ; save it off...
    self.target_rect[0:1] = view_loc
    self.target_rect[2:3] = view_dim
    self.target_cntr = [self.target_rect[0:1] + self.target_rect[2:3] / 2]
END


;----------------------------------------------------------------------------
;Start a dynamic manipulation by computing which corner the drag
;started in.
;
PRO IDLexViewManip::MouseDown, mouseXY, oDest

    self->CalcRect

    ; Check the hit location...
    inset2 = 12*2

    v = [[self.target_rect[0],self.target_rect[1]], $
         [self.target_rect[0]+self.target_rect[2],self.target_rect[1]], $
         [self.target_rect[0]+self.target_rect[2],$
          self.target_rect[1]+self.target_rect[3]], $
         [self.target_rect[0],self.target_rect[1]+self.target_rect[3]]]

    self.mouseXY = mouseXY
    self.downpt = mouseXY[0:1]
    self.orig_rect = self.target_rect

    self.oTarget->GetProperty, viewplane_rect=viewplane_rect
    self.down_viewplane_rect = viewplane_rect
    self.down_radius = sqrt(total((mouseXY - self.target_cntr)^2))

    case 1 of
        mouseXY[0] GT v[0,2] or $ ; Outside of self?
        mouseXY[0] LT v[0,0] or $
        mouseXY[1] GT v[1,2] or $
        mouseXY[1] LT v[1,0]: begin
            if self.mode ne 6 then begin
                self.mode = 0
                end
            end
        self.mode eq 6: $
            self.oManipPoly->SetProperty, HIDE=self.manip_show_manip eq 0
        ABS(mouseXY[0]-v[0,3]) LT inset2 AND $
        ABS(mouseXY[1]-v[1,3]) LT inset2: $
            self.mode = 2
        ABS(mouseXY[0]-v[0,2]) LT inset2 AND $
        ABS(mouseXY[1]-v[1,2]) LT inset2: $
            self.mode = 3
        ABS(mouseXY[0]-v[0,1]) LT inset2 AND $
        ABS(mouseXY[1]-v[1,1]) LT inset2: $
            self.mode = 4
        ABS(mouseXY[0]-v[0,0]) LT inset2 AND $
        ABS(mouseXY[1]-v[1,0]) LT inset2: $
            self.mode = 5
        else: begin
            if self.mode ne 6 then $
                self.mode = 1
            end
        endcase

    ; Show the world...
    self->build_visuals

END

;----------------------------------------------------------------------------
;Finish a drag and reset the current manipulation mode
;
PRO IDLexViewManip::MouseUp, mouseXY, oDest

    self->MouseTrack, mouseXY, oDest

    ; Revert to the default visuals
     if self.mode eq 6 then begin
        self.oManipPoly->SetProperty, /HIDE
        end $
    else begin
        self.mode = 0
        end

    self->build_visuals
END

;----------------------------------------------------------------------------
;Respond to a mouse move by resizing or dragging the target view
;
PRO IDLexViewManip::MouseTrack, mouseXY, oDest

    self.mouseXY = mouseXY
    inset2 = 12*2
    inset = 6

    rect = FLOAT(self.orig_rect)
    dx = FLOAT(mouseXY[0]-self.downpt[0])
    dy = FLOAT(mouseXY[1]-self.downpt[1])

    CASE self.mode OF
        0: BEGIN
            ; nothing hit...
           END
        1:BEGIN ; Translate
            rect[0] = rect[0] + dx
            rect[1] = rect[1] + dy
          END
        2:BEGIN ; UpperLeft (resize)
            rect[0] = rect[0] + dx
            rect[2] = rect[2] - dx
            rect[3] = rect[3] + dy
          END
        3:BEGIN ; UpperRight (resize)
            rect[2] = rect[2] + dx
            rect[3] = rect[3] + dy
          END
        4:BEGIN ; LowerRight (resize)
            rect[1] = rect[1] + dy
            rect[2] = rect[2] + dx
            rect[3] = rect[3] - dy
          END
        5:BEGIN ; LowerLeft (resize)
            rect[0] = rect[0] + dx
            rect[1] = rect[1] + dy
            rect[2] = rect[2] - dx
            rect[3] = rect[3] - dy
          END
        6:BEGIN ; Zoom
            radius = sqrt(total((mouseXY - self.target_cntr)^2))
            if radius eq 0 then $
                RETURN
            zoom = self.down_radius / radius
            if zoom eq 0 then $
                RETURN
            old_center = self.down_viewplane_rect[0:1] $
                + self.down_viewplane_rect[2:3] / 2
            new_x_size = self.down_viewplane_rect[2] * zoom
            new_y_size = self.down_viewplane_rect[3] * zoom

            viewplane_rect = [ $
                new_x_size / (-2) + old_center[0], $
                new_y_size / (-2) + old_center[1], $
                new_x_size, $
                new_y_size $
                ]
            self.oTarget->SetProperty, viewplane_rect=viewplane_rect
            self->build_visuals
            RETURN
          END
    ENDCASE

    ; Size cannot be negative/too small.
    IF (rect[2] LT inset2) THEN rect[2] = inset2
    IF (rect[3] LT inset2) THEN rect[3] = inset2

    ; Make sure it can be seen.
    oDest->GetProperty, DIMENSIONS=dims
    IF (rect[0] GT dims[0]-inset) THEN rect[0] = dims[0]-inset
    IF (rect[1] GT dims[1]-inset) THEN rect[1] = dims[1]-inset
    IF (rect[0]+rect[2] LT inset) THEN rect[2] = inset-rect[0]
    IF (rect[1]+rect[3] LT inset) THEN rect[3] = inset-rect[1]

    ; Convert to proper units.
    CASE self.targetUnits OF
        0: BEGIN ; Device
               ; Do nothing.
           END
        1: BEGIN ; Inches
               rect[0:1] = rect[0:1] * resolution / 2.54
               rect[2:3] = rect[2:3] * resolution / 2.54
           END
        2: BEGIN ; Centimeters
               rect[0:1] = rect[0:1] * resolution
               rect[2:3] = rect[2:3] * resolution
           END
        3: BEGIN ; Normalized
               self.oDest->GetProperty, DIMENSIONS=dest_dim
               rect[0:1] = rect[0:1] / dest_dim
               rect[2:3] = rect[2:3] / dest_dim
           END
    ENDCASE
    ; Make it so..
    self.oTarget->IDLgrView::SetProperty, LOCATION=rect[0:1],$
                                          DIMENSIONS=rect[2:3]
    self->build_visuals
END

;----------------------------------------------------------------------------
;Define the object
;
PRO IDLexViewManip__define
    struct = { IDLexViewManip, $
               INHERITS IDLgrView, $
               mode: 0, $
               manip_show_manip: 0b, $
               manip_show_sel: 0B, $
               target_rect: LONARR(4), $
               orig_rect: LONARR(4), $
               downpt: LONARR(2), $
               down_radius: 0.0, $
               center: FLTARR(2), $
               down_viewplane_rect: DBLARR(4), $
               oTarget: OBJ_NEW(), $
               oDest: OBJ_NEW(), $
               targetUnits: 0L, $
               target_cntr: FLTARR(2), $
               oModel: OBJ_NEW(), $
               oSelPoly: OBJ_NEW(), $
               mouseXY: LONARR(2), $
               oManipPoly: OBJ_NEW() $
             }
END

;----------------------------------------------------------------------------

