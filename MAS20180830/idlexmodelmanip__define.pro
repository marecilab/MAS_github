; $Id: //depot/idl/IDL_64/idldir/lib/utilities/idlexmodelmanip__define.pro#1 $
;
; Copyright (c) 1997-2007, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;+
; NAME:
;       IDLexModelManip
;
; PURPOSE:
;       This object serves to simplify dynamic manipulation of the
;       transformation matrix of an IDLgrModel object.
;
; CATEGORY:
;       Object graphics examples.
;
; CALLING SEQUENCE:
;
;       oMm = OBJ_NEW('IDLexModelManip')
;
;   IDLexModelManip::MouseTrack, pos, oDest
;   IDLexModelManip::MouseUp, pos, oDest
;   IDLexModelManip::MouseDown, [event.x, event.y], oDest
;   IDLexModelManip::SetTarget, oTarget, oDest, _extra=e
;   IDLexModelManip::Undo -- undo manipulation
;   IDLexModelManip::SetProperty
;   IDLexModelManip::GetProperty
;       TARGET - The IDLgrModel Object currently being manipulated (GET only)
;       MODE - The manipulation mode (0=translation,1=rotation,2=scale)
;       TWO_D - Set if the visuals should be 2D instead of 3D
;       MANIP_SHOW_SEL - True if the selection visual is to be displayed
;       MANIP_SHOW_MANIP - True if the manipulation visual is to be displayed
;       TRANSLATE - Set to a 3 element binary vector.  Element is true if
;           translation along that axis is to be allowed.  Axes are with
;           respect to self's orientation.  Note: the default is no
;           translation, i.e [0b,0b,0b]. Set or initialize this property
;           to a value other than the default to allow translation.
;       CONSTRAIN_REF_PT - Reference point for constrained translation
;       CONSTRAIN_TRANS - True if constrained translation is in effect
;       CONSTRAIN_XRANGE - Limits of translation in X [min,max]
;       CONSTRAIN_YRANGE - Limits of translation in Y [min,max]
;       CONSTRAIN_ZRANGE - Limits of translation in Z [min,max]
;
; KEYWORD PARAMETERS:
;       IDLexModelManip::Init:
;
;
; MODIFICATION HISTORY:
;       Written by:     RJF, Nov 1996
;-
;----------------------------------------------------------------------------
;Function to convert a rotation speecified in Quaternion form into a 4x4
;rotation matrix
;

FUNCTION IDLexModelManip__Quaternion_m3, q
compile_opt hidden

 x = q[0]
 y = q[1]
 z = q[2]
 w = q[3]

 m = [[ w^2+x^2-y^2-z^2, 2*(x*y-w*z), 2*(x*z+w*y), 0], $
      [ 2*(x*y+w*z), w^2-x^2+y^2-z^2, 2*(y*z-w*x), 0], $
      [ 2*(x*z-w*y), 2*(y*z+w*x), w^2-x^2-y^2+z^2, 0], $
      [ 0          , 0          , 0              , 1]]

 RETURN, m

END

;----------------------------------------------------------------------------
;Routine to return the bounding box of the graphics tree starting at OBJ
;with respect to TOP
;
;This routine traverses an object tree collecting the largest bounding
;box containing all the children of the object.
;
FUNCTION IDLexModelManip__get_bounds, obj, top
    compile_opt hidden

    FORWARD_FUNCTION IDLexModelManip__get_bounds

    outvalid = 0
    outbox = DBLARR(3,2)
    tmp = DBLARR(4,2)

    ; If a graphic, transform the xyzrange of the object and return it.
    IF OBJ_ISA(obj,'IDLgrGraphic') THEN BEGIN
        obj->IDLgrGraphic::GetProperty, XRANGE=xr, YRANGE=yr, ZRANGE=zr

        ; Because of rotation, we must consider all 8 points...
        p = DBLARR(4)
        tm = obj->GetCTM(TOP=top) ; CTM from graphic up to TOP

    ; cook up all 8 corner verticies
    FOR k=0,1 DO BEGIN
            FOR j=0,1 DO BEGIN
                FOR i=0,1 DO BEGIN

                    p[0] = xr[i]
                    p[1] = yr[j]
                    p[2] = zr[k]
                    p[3] = 1.0

                    ; Mult point by the transform to place in TOP space.
                    p = p # tm
                    ; Divide by W (just to be safe...)!!
                    IF (p[3] NE 0.0) THEN $
                        p[*] = p[*]/p[3]

                    ; first time through, just store the value
                    IF (outvalid EQ 0) THEN BEGIN
                        outbox[*,0] = p[0:2]
                        outbox[*,1] = p[0:2]
                        outvalid = 1
                    ; subsequent times collect min/max values
                    ENDIF ELSE BEGIN
                        outbox[0,0] = outbox[0,0] < p[0]
                        outbox[1,0] = outbox[1,0] < p[1]
                        outbox[2,0] = outbox[2,0] < p[2]
                        outbox[0,1] = outbox[0,1] > p[0]
                        outbox[1,1] = outbox[1,1] > p[1]
                        outbox[2,1] = outbox[2,1] > p[2]
                    ENDELSE
                ENDFOR
            ENDFOR
        ENDFOR

    ENDIF ELSE IF (OBJ_ISA(obj,'IDLgrModel')) THEN BEGIN

        ; Loop through all children of model.
        num = obj->Count()
        list = obj->Get(/ALL)
        IF (num GT 0) THEN BEGIN
            FOR i=0L,num-1 DO BEGIN

                ; Get the bounding box of the child
                box = IDLexModelManip__get_bounds(list[i],top)
                si = SIZE(box)
                IF (si[0] EQ 2) THEN BEGIN

                    ; first time through
                    IF (outvalid EQ 0) THEN BEGIN
                        outvalid = 1
                        outbox = box
                    ; collect the min/max values
                    ENDIF ELSE BEGIN
                        outbox[0,0]=MIN([outbox[0,0],box[0,0]])
                        outbox[1,0]=MIN([outbox[1,0],box[1,0]])
                        outbox[2,0]=MIN([outbox[2,0],box[2,0]])
                        outbox[0,1]=MAX([outbox[0,1],box[0,1]])
                        outbox[1,1]=MAX([outbox[1,1],box[1,1]])
                        outbox[2,1]=MAX([outbox[2,1],box[2,1]])
                    ENDELSE
                ENDIF

            ENDFOR

        ENDIF
    ENDIF

    IF (outvalid EQ 1) THEN BEGIN
        RETURN, outbox
    ENDIF ELSE $
        RETURN, -1
END

;----------------------------------------------------------------------------
;Initialize self.
;
FUNCTION IDLexModelManip::Init, $
    TRANSLATE=translate, $
    CONSTRAIN_REF_PT=constrain_ref_pt, $
    CONSTRAIN_TRANS=constrain_trans,$
    CONSTRAIN_XRANGE=constrain_xrange, $
    CONSTRAIN_YRANGE=constrain_yrange, $
    CONSTRAIN_ZRANGE=constrain_zrange, $
    TWO_D=two_d, $
    MANIP_SHOW_SEL=manip_show_sel,$
    MANIP_SHOW_MANIP=manip_show_manip,$
    MODE=mode, $
    SELECTOR_COLOR=selector_color, $
    MANIPULATOR_COLOR=manipulator_color, $
    ROT_RADIUS=rot_radius, $                ;IN: (opt) default .75
    _EXTRA=e

    compile_opt hidden

    CATCH, error_status
    IF error_status NE 0 THEN BEGIN
        CATCH, /CANCEL
        PRINT, !ERROR_STATE.MSG_PREFIX + !ERROR_STATE.MSG
        RETURN, 0
    END

    IF (self->IDLgrModel::init(_EXTRA=e) NE 1) THEN $
        MESSAGE, 'Failed to initialize IDLgrModel::'

    default_color = [128, 128, 128]
    IF (N_ELEMENTS(selector_color) EQ 0) THEN $
        selector_color = default_color
    IF N_ELEMENTS(manipulator_color) LE 0 THEN $
        manipulator_color = default_color
;
;   Create the local tree
;
    self.oManipScale = OBJ_NEW('IDLgrModel')
    self.oManipPolyline = OBJ_NEW('IDLgrPolyline',color=manipulator_color)
    self.oSelModel = OBJ_NEW('IDLgrModel')
    self.oSelPolyline = OBJ_NEW('IDLgrPolyline',color=selector_color)

    self->add,self.oSelModel
    self.oSelModel->add,self.oSelPolyline

    self.oManipScale->add,self.oManipPolyline
    self.oSelModel->add,self.oManipScale
;
;   Default state
;
    self.oTarget = OBJ_NEW()
    self.button_down = [-1,-1]
    self.mode = 0
    self.init_trans = IDENTITY(4)
    self.orig_xform = IDENTITY(4)
    self.two_d = 0
    self.translate = [0,0,0]
    self.constrain_ref_pt = [0.,0.,0.]
    self.constrain_xyzrange=[[0.,0.],[0.,0.],[0.,0.]]
    self.constrain_trans = 0b
    self.manip_show_sel = 1b
    self.manip_show_manip = 1b
    self.rot_radius = .75
;
;   Keywords
;
    IF (N_ELEMENTS(mode) EQ 1) THEN BEGIN
        IF ((mode GE 0) AND (mode LE 2)) THEN BEGIN
            self.mode = mode
        END
    END

    IF (KEYWORD_SET(two_d)) THEN $
        self.two_d = 1

    IF (N_ELEMENTS(translate) EQ 3) THEN $
        self.translate = translate

    n = N_ELEMENTS(constrain_ref_pt)
    IF (n eq 2) THEN BEGIN
        self.constrain_ref_pt[0:1] = constrain_ref_pt
        self.constrain_ref_pt[2] = 0
    ENDIF ELSE IF (n eq 3) THEN $
        self.constrain_ref_pt = constrain_ref_pt

    IF (N_ELEMENTS(constrain_xrange) EQ 2) THEN $
        self.constrain_xyzrange[*,0] = constrain_xrange
    IF (N_ELEMENTS(constrain_yrange) EQ 2) THEN $
        self.constrain_xyzrange[*,1] = constrain_yrange
    IF (N_ELEMENTS(constrain_zrange) EQ 2) THEN $
        self.constrain_xyzrange[*,2] = constrain_zrange

    IF (KEYWORD_SET(constrain_trans)) THEN $
        self.constrain_trans = 1

    IF (N_ELEMENTS(manip_show_sel) EQ 1) THEN $
       self.manip_show_sel = manip_show_sel
    IF (N_ELEMENTS(manip_show_manip) EQ 1) THEN $
       self.manip_show_manip = manip_show_manip

    IF (N_ELEMENTS(rot_radius) EQ 1) THEN $
       self.rot_radius  = rot_radius

    self->Build_Manip_Visual

    RETURN, 1 ; Success.
END

;----------------------------------------------------------------------------
;Create a polyline object representing the manipulation tool
;
PRO IDLexModelManip::Build_Manip_Visual
    compile_opt hidden

    ON_ERROR, 2 ; Return to caller on error.

    IF NOT OBJ_VALID(self.oTarget) THEN $
        RETURN

    self.oManipPolyline->SetProperty, HIDE=(self.manip_show_manip EQ 0)

    case self.mode of
        0 : BEGIN   ;translate
            verts = [[0,-1,-1],[0,-1,1],[0,1,1],[0,1,-1],$
                    [-1,0,-1],[-1,0,1],[1,0,1],[1,0,-1],$
                    [-1,-1,0],[-1,1,0],[1,1,0],[1,-1,0]]
            IF (self.two_d) THEN BEGIN
                conn = [5,8,9,10,11,8]
            ENDIF ELSE $
                conn = [5,0,1,2,3,0,5,4,5,6,7,4,5,8,9,10,11,8]
                self.oManipPolyline->SetProperty,DATA=verts,POLYLINES=conn
            END
        1 : BEGIN   ;rotate
            num = 21
            verts = DBLARR(3,3*num)
            conn = LONARR(3*(num+1))
            t = 0.0
            tinc = (2.*!PI)/FLOAT(num-1)
            FOR i=0,num-1 DO BEGIN
                verts[0,i] = COS(t)
                verts[1,i] = SIN(t)
                verts[2,i] = 0.0
                verts[1,i+num] = COS(t)
                verts[2,i+num] = SIN(t)
                verts[0,i+num] = 0.0
                verts[2,i+2*num] = COS(t)
                verts[0,i+2*num] = SIN(t)
                verts[1,i+2*num] = 0.0
                t = t + tinc
                conn[i+1] = i
                IF (NOT self.two_d) THEN BEGIN
                    conn[i+1+(num+1)] = i+num
                    conn[i+1+2*(num+1)] = i+(2*num)
                ENDIF
            END
            conn[0] = num
            IF (NOT self.two_d) THEN BEGIN
                conn[num+1] = num
                conn[2*(num+1)] = num
            ENDIF
            self.oManipPolyline->SetProperty,DATA=verts,POLYLINES=conn
            END
        2 : BEGIN   ;scale
            verts = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1],$
                     [.9,.1,0],[.9,-.1,0],[-.9,.1,0],[-.9,-.1,0],$
                     [0,.9,.1],[0,.9,-.1],[0,-.9,.1],[0,-.9,-.1],$
                     [.1,0,.9],[-.1,0,.9],[.1,0,-.9],[-.1,0,-.9]]
            conn = [2,0,1,2,2,3,2,4,5,3,6,0,7,3,8,1,9,$
                    3,10,2,11,3,12,3,13,3,14,4,15,3,16,5,17]
            self.oManipPolyline->SetProperty,data=verts,polylines=conn
            END
    END
;
;   Scale the manipulator.
;
    self.oManipScale->SetProperty, TRANSFORM=identity(4)
    d = MAX(ABS(self.sel_box))*sqrt(2)
    self.oManipScale->Scale,d,d,d
END

;----------------------------------------------------------------------------
;Destroy the object
;
PRO IDLexModelManip::Cleanup
    compile_opt hidden
;
;   Remove from tree if needed
;
    IF (OBJ_VALID(self.oTarget)) THEN BEGIN
    list = self.oTarget->Get(/ALL)
    IF (SIZE(list,/TYPE) EQ 11) THEN $
            self.oTarget->Remove, self
        self.oTarget = OBJ_NEW()
    ENDIF

    self->IDLgrModel::Cleanup
END

;----------------------------------------------------------------------------
;Change a property of the object
;
PRO IDLexModelManip::SetProperty, $
    TRANSLATE=translate, $
    CONSTRAIN_REF_PT=constrain_ref_pt, $
    CONSTRAIN_TRANS=constrain_trans,$
    CONSTRAIN_XRANGE=constrain_xrange, $
    CONSTRAIN_YRANGE=constrain_yrange, $
    CONSTRAIN_ZRANGE=constrain_zrange, $
    TWO_D=two_d, $
    MANIP_SHOW_SEL=manip_show_sel,$
    MANIP_SHOW_MANIP=manip_show_manip,$
    MODE=mode, $
    _EXTRA=e

    compile_opt hidden
    rebuild_visual = 0

    IF (N_ELEMENTS(mode) EQ 1) THEN BEGIN
        IF ((mode GE 0) AND (mode LE 2)) THEN BEGIN
            IF (self.mode NE mode) THEN BEGIN
                self.mode = mode
                rebuild_visual = 1
            END
        END
    END
    IF (N_ELEMENTS(two_d) EQ 1) THEN BEGIN
        self.two_d = two_d
        rebuild_visual = 1
    ENDIF
    IF (N_ELEMENTS(translate) EQ 3) THEN BEGIN
        self.translate = translate
        rebuild_visual = 1
    ENDIF

    n = N_ELEMENTS(constrain_ref_pt)
    IF (n eq 2) THEN BEGIN
        self.constrain_ref_pt[0:1] = constrain_ref_pt
        self.constrain_ref_pt[2] = 0
    ENDIF ELSE IF (n eq 3) THEN $
        self.constrain_ref_pt = constrain_ref_pt

    IF (N_ELEMENTS(constrain_xrange) EQ 2) THEN $
        self.constrain_xyzrange[*,0] = constrain_xrange
    IF (N_ELEMENTS(constrain_yrange) EQ 2) THEN $
        self.constrain_xyzrange[*,1] = constrain_yrange
    IF (N_ELEMENTS(constrain_zrange) EQ 2) THEN $
        self.constrain_xyzrange[*,2] = constrain_zrange

    IF (N_ELEMENTS(constrain_trans) EQ 1) THEN $
        self.constrain_trans = constrain_trans

    IF (N_ELEMENTS(manip_show_sel) EQ 1) THEN BEGIN
       self.manip_show_sel = manip_show_sel
       rebuild_visual = 1
    ENDIF
    IF (N_ELEMENTS(manip_show_manip) EQ 1) THEN BEGIN
       self.manip_show_manip = manip_show_manip
       rebuild_visual = 1
    ENDIF

    IF (rebuild_visual) THEN $
        self->Build_Manip_Visual

    self.oSelPolyline->SetProperty, HIDE=(self.manip_show_sel EQ 0)

    self->IDLgrModel::SetProperty, _EXTRA=e
END

;----------------------------------------------------------------------------
;Return properties of the object
;
PRO IDLexModelManip::GetProperty, $
    TARGET=target, $
    MODE=mode, $
    TWO_D=two_d, $
    MANIP_SHOW_SEL=manip_show_sel, $
    MANIP_SHOW_MANIP=manip_show_manip,$
    TRANSLATE=translate, $
    CONSTRAIN_REF_PT=constrain_ref_pt, $
    CONSTRAIN_TRANS=constrain_trans, $
    CONSTRAIN_XRANGE=constrain_xrange, $
    CONSTRAIN_YRANGE=constrain_yrange, $
    CONSTRAIN_ZRANGE=constrain_zrange, $
    ALL=ALL, $
    _REF_EXTRA=e

    compile_opt hidden

    target = self.oTarget
    mode = self.mode
    two_d = self.two_d
    manip_show_sel = self.manip_show_sel
    manip_show_manip = self.manip_show_manip
    translate = self.translate
    constrain_ref_pt = self.constrain_ref_pt
    constrain_trans  = self.constrain_trans
    constrain_xrange = self.constrain_xyzrange[*,0]
    constrain_yrange = self.constrain_xyzrange[*,1]
    constrain_zrange = self.constrain_xyzrange[*,2]
    self->IDLgrModel::GetProperty, _EXTRA=e

    IF ARG_PRESENT(all) NE 0 THEN BEGIN
        self->IDLgrModel::GetProperty, ALL=all
        all = CREATE_STRUCT( $
            all, $
            'target', target, $
            'mode', mode, $
            'two_d', two_d, $
            'manip_show_sel', manip_show_sel, $
            'manip_show_manip', manip_show_manip, $
            'translate', translate, $
            'constrain_ref_pt', constrain_ref_pt, $
            'constrain_trans',  constrain_trans, $
            'constrain_xrange', constrain_xrange, $
            'constrain_yrange', constrain_yrange, $
            'constrain_zrange', constrain_zrange $
            )
    END
END

;----------------------------------------------------------------------------
;Apply the original transform to the target
;
PRO IDLexModelManip::Undo
    compile_opt hidden

    ON_ERROR, 2 ; Return to caller on error.

    IF OBJ_VALID(self.oTarget) THEN $
        self.oTarget->SetProperty, TRANSFORM=self.orig_xform $
    ELSE $
        MESSAGE,'Can''t Undo. No valid target for manipulation.'
END

;----------------------------------------------------------------------------
;Point this manipulator at a new target object
;
PRO IDLexModelManip::SetTarget, oTarget, oDest, _EXTRA=e
    compile_opt hidden

    FORWARD_FUNCTION IDLexModelManip__get_bounds
    ON_ERROR, 2 ; Return to caller on error.

    IF N_ELEMENTS(oTarget) LE 0 THEN $
        MESSAGE, 'First argument is undefined.'

    IF NOT (oTarget EQ OBJ_NEW()) THEN $
        IF NOT OBJ_ISA(oTarget, 'IDLgrModel') THEN $
            MESSAGE, 'First argument must be a reference to an IDLgrModel ' + $
                     'or null object.'
;
;   Unhitch relationship with old target.
;
    IF OBJ_VALID(self.oTarget) THEN BEGIN
        self.oTarget->Remove, self
        self.oTarget = OBJ_NEW()
    ENDIF
;
;   Build selection visual.
;
    box = IDLexModelManip__get_bounds(oTarget, oTarget)

    IF (SIZE(box))[0] EQ 0 THEN BEGIN ; box is a scalar or undefined.
        self.oSelPolyline->SetProperty, /HIDE
        self.sel_box = [[-1,-1,-1],[1,1,1]]
    ENDIF ELSE BEGIN

        ; Create a box around the selected model heirarchy
        self.oSelPolyline->SetProperty, HIDE=0
        verts = [[box[0,0],box[1,0],box[2,0]], $
                 [box[0,1],box[1,0],box[2,0]], $
                 [box[0,1],box[1,1],box[2,0]], $
                 [box[0,0],box[1,1],box[2,0]], $
                 [box[0,0],box[1,0],box[2,1]], $
                 [box[0,1],box[1,0],box[2,1]], $
                 [box[0,1],box[1,1],box[2,1]], $
                 [box[0,0],box[1,1],box[2,1]]]

        conn = [5,0,1,2,3,0, 5,4,5,6,7,4,$
                2,0,4, 2,1,5, 2,2,6, 2,3,7]

        self.oSelPolyline->SetProperty, DATA=verts, POLYLINES=conn
        self.sel_box = box
    ENDELSE
;
;   Hitch up with new target.
;
    IF OBJ_VALID(oTarget) THEN BEGIN
        oTarget->GetProperty, TRANSFORM=tm
        self.orig_xform = tm
        oTarget->Add, self
    END ELSE $
        self.orig_xform = IDENTITY(4)
    self.oTarget = oTarget

    self->Build_Manip_Visual

    ; should we display the selection lines?
    self.oSelPolyline->SetProperty, HIDE=(self.manip_show_sel EQ 0)

    IF (self.button_down[0] NE -1) AND OBJ_VALID(oTarget) THEN BEGIN
        self->MouseDown, self.button_down, oDest
    ENDIF
END

;----------------------------------------------------------------------------
;Internal function to get the IDLgrView object above (self)
;
FUNCTION IDLexModelManip::Get_View
    compile_opt hidden
    tobj = self

    WHILE NOT OBJ_ISA(tobj,'IDLgrView') DO BEGIN
        tobj->IDLgrComponent::GetProperty, PARENT=parent
        IF NOT OBJ_VALID(parent) THEN $
            RETURN, parent
        tobj = parent
    ENDWHILE

    RETURN, tobj
END

;----------------------------------------------------------------------------
;Start a manipulation with a mouse down.
;
PRO IDLexModelManip::MouseDown, pos, dest
    compile_opt hidden

    on_error, 2 ; Return to caller on error.

    if N_ELEMENTS(pos) le 0 then $
        message, 'First argument is undefined.'

    if N_ELEMENTS(pos) lt 2 then $
        message, 'First argument should be a 2 element array.'

    if (SIZE(pos))[0] ne 1 then $
        message, 'First argument should be a 1-D array.'

    tname = SIZE(pos, /TNAME)
    if tname eq 'STRING' or $
       tname eq 'STRUCT' or $
       tname eq 'POINTER' or $
       tname eq 'OBJREF' then $
        message, 'First argument should be of numeric type.'

    if N_ELEMENTS(dest) le 0 then $
        message, 'Second argument is undefined.'

    if N_ELEMENTS(dest) gt 1 then $
        message, 'Second argument has more than one element.'

    if not OBJ_VALID(dest) then $
        message, 'Second argument is not a valid object.'

    if not OBJ_ISA(dest, 'IDLgrWindow') then $
        message, 'Second argument should be an IDLgrWindow.'
;
;   Each manipulation has a starting point, record that point now
;
    self.oTarget->GetProperty, TRANSFORM=tm
    self.init_trans=tm
    self.button_down = pos
;
;   Hide our selection lines?
;
    self.oSelPolyline->SetProperty, HIDE=(self.manip_show_sel EQ 0)
;
;   The manipulations need to be aware of their location in screen
;   space.  We do this here by obtaining the view we are in and
;   computing the necessary pixel/view transforms.
;
    view = self->Get_View()
    dest->IDLgrWindow::GetProperty, UNITS=units
    dest->IDLgrWindow::SetProperty, UNITS=0
    dest->IDLgrWindow::GetProperty, DIMENSIONS=dest_dim, RESOLUTION=resolution
    dest->IDLgrWindow::SetProperty, UNITS=TEMPORARY(units)
    IF OBJ_VALID(view) THEN BEGIN
        view->IDLgrView::GetProperty, LOCATION=view_loc, DIMENSIONS=view_dim, $
            UNITS=units
;
;       Convert to device units.
;
        CASE units OF
            0 : BEGIN ; Device
                ; Do nothing.
                END
            1 : BEGIN ; Inches
                view_loc = view_loc * 2.54 / resolution
                view_dim = view_dim * 2.54 / resolution
                END
            2 : BEGIN ; Centimeters
                view_loc = view_loc / resolution
                view_dim = view_dim / resolution
                END
            3 : BEGIN ; Normalized
                view_loc = view_loc * dest_dim
                view_dim = view_dim * dest_dim
                END
        ENDCASE
;
;       Handle the special case of a View of size (0,0) which is
;       sized to the destination.
;
        IF ((view_dim[0] EQ 0) OR (view_dim[1] EQ 0)) THEN BEGIN
            view_dim = dest_dim
            view_loc = [0,0]
        END
;
;       Get the current transformation matrix from our object to
;       screen space
;
        tm = self.oTarget->getCTM(DESTINATION=dest)
;
;       Transform the point (0,0,0) in object space to screen space
;
        pt = DBLARR(4)
        pt = [0,0,0,1]
        pt = pt # tm
        IF (pt[3] NE 0.0) THEN $
            pt = pt / pt[3]
;
;       Convert the point from normalized coords to pixels
;
        pt[0] = (pt[0] + 1.0)*0.5*view_dim[0]+view_loc[0]
        pt[1] = (pt[1] + 1.0)*0.5*view_dim[1]+view_loc[1]
        pt[2] = (pt[2] + 1.0)*0.5
        pt[3] = 1.0
;
;       Record the mouse down point
;
        self.point[0:1,0] = pos[0:1]
        self.point[2,0] = 0.0

;
;       Set up the transformation mode
;
        CASE self.mode OF
            0: BEGIN ; translation (object appears to move parallel to screen)

                ; Pick the axis to allow translation along
                doTrans = WHERE(self.translate NE 0, nTrans)
                IF (ntrans GT 0) THEN BEGIN
                    IF OBJ_ISA(dest, 'IDLgrWindow') THEN $
                        dest->IDLgrWindow::SetCurrentCursor, "MOVE"

                    ; Compute the motion in data space of pixel motion
                    ; in screen space.
                    ; Do this by computing the change in data space due
                    ; to a change in screen space (pixels) of one unit.

                    ; We need the inverse matrix to go from screen to
                    ; data units
                    inv_tm = INVERT(tm)
                    p = DBLARR(4)

                    ; Transform the point computed previously + 1 pixel in x
                    p[0] = ((pt[0]+1.0-view_loc[0])/FLOAT(view_dim[0]))*2.0-1.0
                    p[1] = ((pt[1]-view_loc[1])/FLOAT(view_dim[1]))*2.0-1.0
                    p[2] = ((pt[2])/1.0)*2.0-1.0
                    p[3] = 1.0

                    o = p # inv_tm
                    IF (o[3] NE 0.0) THEN $
                        o = o / o[3]

                    ; Store off the vector.  A unit change on the screen X
                    ; corresponds to adding this vector to our target model
                    self.point[0:2,1] = o[0:2]

                    ; Transform the point computed previously + 1 pixel in y
                    p[0] = ((pt[0]-view_loc[0])/FLOAT(view_dim[0]))*2.0-1.0
                    p[1] = ((pt[1]+1.0-view_loc[1])/FLOAT(view_dim[1]))*2.0-1.0
                    p[2] = ((pt[2])/1.0)*2.0-1.0
                    p[3] = 1.0

                    o = p # inv_tm
                    IF (o[3] NE 0.0) THEN $
                        o = o / o[3]

                    ; Store off the vector.  A unit change on the screen Y
                    ; corresponds to adding this vector to our target model
                    self.point[0:2,2] = o[0:2]

                ENDIF
               END
            1: BEGIN ; Rotation. Object rotates about its (local) center (0,0,0)
                    IF OBJ_ISA(dest, 'IDLgrWindow') THEN $
                        dest->IDLgrWindow::SetCurrentCursor, "CROSSHAIR"

                    ; Get the manipulation origin (in screen space).
                    ; This is the location where the target (0,0,0) location
                    ; will appear on the screen (computed previously).
                    self.point[0:1,1] = pt[0:1]

                    ; Mouse motion outside of a circle centered at the
                    ; object origin point will be converted into pure screen
                    ; Z axis rotation.
                    ; self.point[2,1] = (dest_dim[0] < dest_dim[1]) * .5 * .75
                    self.point[2,1] = $
                        (view_dim[0] < view_dim[1]) * .5 * self.rot_radius

                    ; Rotation is computed as follows:
                    ; Form a hemisphere at the object's center of rotation
                    ; Take the relative X,Y coords of the mouse point and
                    ; project these onto the hemisphere, forming a vector
                    ; from the origin to this point in 3D.  As the mouse moves,
                    ; project the new mouse point onto the hemisphere.  The
                    ; object is then rotated over a line perpendicular to
                    ; the hemisphere vectors formed by the initial and current
                    ; mouse points.  The angle of rotation is the angle between
                    ; the two vectors.

                    ; Project the initial point onto hemisphere.
                    xy = (self.point[0:1,1] - pos[0:1]) / self.point[2,1]
                    r = TOTAL(xy^2)
                    IF (r GT 1.0) THEN BEGIN
                        self.point[*,0] = [xy/SQRT(r), 0.0]
                    ENDIF ELSE BEGIN
                        self.point[*,0] = [xy, SQRT(1.0-r)]
                    ENDELSE

                    ; Transform the initial hemisphere vector to data space.
                    p0 = self->V_to_Data(dest, self.point[*,0])
                    self.point[*,0] = p0
               END
            2: BEGIN
                ; Isotropic scaling about the projected center of the
                ; object.  The object gets bigger as the mouse is dragged
                ; away from the projected center and smaller as it is dragged
                ; closer to the center (radii smaller than the initial point).

                ; Get the manipulation origin (in screen space).
                self.point[0:1,1] = pt[0:1]
                xy = pos[0:1] - self.point[0:1,1]

                ; Get the initial manipulation radius
                r = SQRT(TOTAL(xy^2))
                IF (r EQ 0) THEN r = 1
                self.point[2,1] = r

               END
        ENDCASE
    ENDIF
END

;----------------------------------------------------------------------------
;Finish the manipulation as the mouse is released
;
PRO IDLexModelManip::MouseUp, pos, dest

    compile_opt hidden
    tm = self.init_trans

    CASE self.mode OF
        0: BEGIN ;Translate. - compute a final translation
            doTrans = WHERE(self.translate NE 0, nTrans)
            IF (ntrans GT 0) THEN BEGIN
                IF OBJ_ISA(dest, 'IDLgrWindow') THEN $
                    dest->IDLgrWindow::SetCurrentCursor, "ORIGINAL"

                ; compute the screen space change in X and Y
                dx = FLOAT(pos[0]-self.point[0,0])
                dy = FLOAT(pos[1]-self.point[1,0])

                ; Motion is the screen space change times the data space vectors
                t = dx*self.point[*,1]
                t = t + dy*self.point[*,2]

                ; Only translate in requested directions.
                invalidT = WHERE(self.translate EQ 0, ninvalid)
                IF (ninvalid GT 0) THEN $
                    t[invalidT] = 0

                ; If appropriate, constrain to requested range.
                IF (self.constrain_trans) THEN BEGIN
                    xref_pt = [self.constrain_ref_pt,1.] # tm
                    xref_pt[0:2] = xref_pt[0:2]/xref_pt[3]
                    xform_pt = xref_pt[0:2] + t

                    ; clamp values which are too small
                    toosmall = WHERE(xform_pt LT self.constrain_xyzrange[0,*],$
                                     ntoosmall)
                    FOR i=0L,ntoosmall-1 DO BEGIN
                        xyz = toosmall[i]
                        IF (self.translate[xyz]) THEN BEGIN
                             t[xyz] = t[xyz] + $
                                      (self.constrain_xyzrange[0,xyz] - $
                                       xform_pt[xyz])
                        ENDIF
                    ENDFOR

                    ; clamp values which are too big
                    toobig = WHERE(xform_pt GT self.constrain_xyzrange[1,*], $
                                   ntoobig)
                    FOR i=0L,ntoobig-1 DO BEGIN
                        xyz = toobig[i]
                        IF (self.translate[xyz]) THEN BEGIN
                             t[xyz] = t[xyz] $
                                    - (xform_pt[xyz] - $
                                        self.constrain_xyzrange[1,xyz])
                        ENDIF
                    ENDFOR
                ENDIF

                ; Apply the final translation
                mult=[[1,0,0,t[0]],[0,1,0,t[1]],[0,0,1,t[2]],[0,0,0,1]]
                tm = mult # tm

            ENDIF
           END
        1: BEGIN ;Rotation. - compute the final rotation
                IF OBJ_ISA(dest, 'IDLgrWindow') THEN $
                    dest->IDLgrWindow::SetCurrentCursor, "ORIGINAL"

                ; Project the point onto the hemisphere
                xy = (self.point[0:1,1] - pos[0:1])/self.point[2,1]
                r = TOTAL(xy^2)
                IF (r GT 1.0) THEN BEGIN
                    self.point[*,2] = [xy/sqrt(r),0.0]
                ENDIF ELSE BEGIN
                    self.point[*,2] = [xy,sqrt(1.0-r)]
                ENDELSE

                ; Transform the hemisphere vector to data space.
                p0 = self->V_to_Data(dest,self.point[*,2])
                self.point[*,2] = p0

                ; Rotate about the axis in data space.
                ; The axis is the cross product of the two data space
                ; vectors.
                q = [crossp(self.point[*,0],self.point[*,2]),$
                     TOTAL(self.point[*,0]*self.point[*,2])]
                mult = IDLexModelManip__Quaternion_m3(q)
                tm = mult # tm
           END
        2: BEGIN ;Scale.
            ; find the vector from the centroid to the mouse
            xy = (pos[0:1] - self.point[0:1,1])
            sc = SQRT(TOTAL(xy^2))
            IF (sc EQ 0.0) THEN sc = 1.0

            ; compute the isotropic scale
            sc = sc/self.point[2,1]

            ; build the simple scaling matrix
            mult = [[sc,0,0,0],[0,sc,0,0],[0,0,sc,0],[0,0,0,1]]
            tm = mult # tm

           END
    ENDCASE

    self.oTarget->SetProperty, TRANSFORM=tm
    self.button_down = [-1,-1]

END

;----------------------------------------------------------------------------
;This function is called to update the object dynamically
;
PRO IDLexModelManip::MouseTrack, pos, dest

    compile_opt hidden
    tm = self.init_trans

    CASE self.mode OF
        0: BEGIN ;Translation.

            ; Compute the change in screen space
            dx = FLOAT(pos[0]-self.point[0,0])
            dy = FLOAT(pos[1]-self.point[1,0])

            ; Motion is the screen space change times the data space vectors
            t = dx*self.point[*,1]
            t = t + dy*self.point[*,2]

            ; Only translate in requested directions..
            invalidT = WHERE(self.translate EQ 0, ninvalid)
            IF (ninvalid GT 0) THEN $
                t[invalidT] = 0

            ; If appropriate, constrain to requested range.
            IF (self.constrain_trans) THEN BEGIN
                xref_pt = [self.constrain_ref_pt,1.] # tm
                xref_pt[0:2] = xref_pt[0:2]/xref_pt[3]
                xform_pt = xref_pt[0:2] + t

                ; clamp the case where the value is too small
                toosmall = WHERE(xform_pt LT self.constrain_xyzrange[0,*], $
                                 ntoosmall)
                FOR i=0L,ntoosmall-1 DO BEGIN
                    xyz = toosmall[i]
                    IF (self.translate[xyz]) THEN BEGIN
                        t[xyz] = t[xyz] + $
                                 (self.constrain_xyzrange[0,xyz] - $
                                  xform_pt[xyz])
                    ENDIF
                ENDFOR

                ; clamp the case where the value is too large
                toobig = WHERE(xform_pt GT self.constrain_xyzrange[1,*], $
                               ntoobig)
                FOR i=0L,ntoobig-1 DO BEGIN
                    xyz = toobig[i]
                    IF (self.translate[xyz]) THEN BEGIN
                        t[xyz] = t[xyz] - $
                                 (xform_pt[xyz] - $
                                  self.constrain_xyzrange[1,xyz])
                    ENDIF
                ENDFOR
            ENDIF

            ; Cook up the translation matrix
            mult=[[1,0,0,t[0]],[0,1,0,t[1]],[0,0,1,t[2]],[0,0,0,1]]
            tm = mult # tm

           END
        1: BEGIN ;Rotate.

                ; Project the point onto the hemisphere
                xy = (self.point[0:1,1] - pos[0:1])/self.point[2,1]
                r = TOTAL(xy^2)

                IF (r GT 1.0) THEN BEGIN
                    self.point[*,2] = [xy/sqrt(r),0.0]
                ENDIF ELSE BEGIN
                    self.point[*,2] = [xy,sqrt(1.0-r)]
                ENDELSE

                ; Transform new hemisphere vector to data space.
                p0 = self->V_to_Data(dest,self.point[*,2])
                self.point[*,2] = p0

                ; Rotate about the axis in data space.
                ; The axis is the cross product of the two data space vectors
                q = [crossp(self.point[*,0],self.point[*,2]),$
                    TOTAL(self.point[*,0]*self.point[*,2])]
                mult = IDLexModelManip__Quaternion_m3(q)
                tm = mult # tm
          END
       2: BEGIN ;Scale.
            ; compute the line from the centroid to the mouse pt
            xy = pos[0:1] - self.point[0:1,1]
            sc = SQRT(TOTAL(xy^2))
            IF (sc EQ 0.0) THEN sc = 1.0
            ; convert the radius into a scale
            sc = sc/self.point[2,1]
            ; build a matrix
            mult = [[sc,0,0,0],[0,sc,0,0],[0,0,sc,0],[0,0,0,1]]
            tm = mult # tm
          END
    ENDCASE

    self.oTarget->SetProperty, TRANSFORM=tm
END

;----------------------------------------------------------------------------
; routine to convert a vector from screen space into object data space
;
FUNCTION IDLexModelManip::V_to_Data,dest,v
    compile_opt hidden

    self.oTarget->SetProperty, TRANSFORM=self.init_trans
    tm = self.oTarget->getCTM(DESTINATION=dest)

    ; Project the point (0,0,0) in data space to the screen
    pt = [0.0,0.0,0.0,1.0]
    pt = pt # tm

    ; Scale PT to -1 to 1 domain.
    IF (pt[3] NE 0.0) THEN $
        pt = pt / pt[3]

    ; Compensate V for a non-equilateral view volume.
    oView = self->Get_View()
    oView->GetProperty, ZCLIP=zclip, VIEWPLANE_RECT=viewplane_rect
    eq_v = [ $
        v[0] / viewplane_rect[2], $
        v[1] / viewplane_rect[3], $
        v[2] / (zclip[0] - zclip[1]) $
        ]
    eq_v = eq_v / sqrt(total(eq_v^2.0))

    ; Add the vector to get the point (in screen space)
    pt[0:2] = pt[0:2] + eq_v

    ; transform PT back into data space
    inv_tm = invert(tm)
    pt = pt # inv_tm
    IF (pt[3] NE 0.0) THEN $
        pt = pt / pt[3]

    ; renormalize the resulting vector
    out = pt[0:2] / sqrt(total(pt[0:2]^2.0))
    RETURN,out
END

;----------------------------------------------------------------------------
;Routine to define the class structure
;
PRO IDLexModelManip__define
    compile_opt hidden
    struct_hide = { IDLexModelManip,  $
        INHERITS IDLgrModel, $
        oTarget:OBJ_NEW(), $
        button_down: LONARR(2), $ ; [event.x, event.y]
        sel_box:DBLARR(3,2), $
        mode:0, $
        point: DBLARR(3,3), $
        init_trans: DBLARR(4,4), $
        orig_xform: DBLARR(4,4), $
        two_d: 0b, $
        translate: bytarr(3), $
        constrain_ref_pt: DBLARR(3), $
        constrain_xyzrange: DBLARR(2,3), $
        constrain_trans: 0b, $
        manip_show_sel: 0b, $
        manip_show_manip: 0b, $
        oManipScale:OBJ_NEW(), $
        oManipPolyline:OBJ_NEW(), $
        oSelModel:OBJ_NEW(), $
        oSelPolyline:OBJ_NEW(), $
        rot_radius: 0.0 $
        }
END
;----------------------------------------------------------------------------
