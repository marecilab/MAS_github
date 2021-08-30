;$Id: //depot/idl/IDL_64/idldir/lib/utilities/xobjview.pro#1 $
;
;  Copyright (c) 1997-2007, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
;+
;  FILE:
;       xobjview.pro
;
;  CALLING SEQUENCE:
;       xobjview, oObj
;
;  INPUTS:
;       Argument oObj is a reference to an atomic graphics object, a
;       reference to an IDLgrModel, or an array of such references.
;
;  INPUT KEYWORDS:
;       XSIZE           Pixels.  Size of drawable.
;       YSIZE           Pixels.  Size of drawable.
;       STATIONARY      Non-moveable objects (e.g. lights) to be viewed.
;       GROUP           Group leader widget.
;       BLOCK           If set, block IDL command line.
;       JUST_REG        If set, register and return (XMANAGER)
;       MODAL           If set, block other IDL widgets from receiving
;                           events.
;       SCALE           Size factor for initial view.  Default: 1/sqrt(3).
;       TITLE           String to appear in Xobjview's title bar.
;       TEST            If set, don't require oObj arg.  Draw a blue
;                           sinusoidal surface instead.
;       DOUBLE_VIEW     Set this keyword for high precision (i.e.
;                           "double precision") graphics.
;       BACKGROUND      Three-element [r,g,b] color vector.
;       REFRESH         Set this keyword to a Top-level base of an XOBJVIEW
;                           to force a redraw of that XOBJVIEW's graphics.
;       XOFFSET         Pixels.  Horzontal postion of XOBJVIEW on the screen.
;       YOFFSET         Pixels.  Vertical postion of XOBJVIEW on the screen.
;       USE_INSTANCING  Set this keyword to enable "instancing."  XOBJVIEW
;                           will use IDL's CREATE_INSTANCE functionality when
;                           drawing graphics, and will use DRAW_INSTANCE if
;                           and when the current graphic is hidden and then
;                           exposed by the user.  See IDLLgrWindow::Draw
;                           documentation for definitions of CREATE_INSTANCE
;                           and DRAW_INSTANCE.
;       DEBUG           Set this keyword to disable error catching and
;                           recovery.  With this keyword set, execution
;                           of XOBJVIEW and/or its supporting routines
;                           will stop where and if an error occurs (unless
;                           other routines with error handling are overriding
;                           this behavior).  This stop-at-error behavior
;                           is often helpful for troubleshooting.
;
;  OUTPUT KEYWORD:
;       TLB             Top-level base widget.
;
;  PURPOSE:
;       Provide a quick way to see graphics objects on the screen.
;
;  CATEGORY:
;       Object Graphics.
;
;  REFERENCE: IDL Reference Guide, IDL User's Guide
;
;  A note about calling XOBJVIEW with the MODAL keyword set:
;       To be modal, XOBJVIEW does not require that its caller specify
;       a group leader.  This is unlike other IDL widget procedures
;       such as XLOADCT, which, to be modal, do require that their
;       caller specify a group leader.  In those other procedures,
;       the requirement exists to encourage the caller to create
;       a modal widget that will be well-behaved with respect to
;       Layering and Iconizing.  (See WIDGET_BASE in the IDL Reference
;       Guide for explanations of "layering" and "iconizing".)   For
;       the same reason the requirement exists in those procedures,
;       supplying an appropriate group leader when invoking
;       XOBJVIEW, /MODAL is good programming practice.  Sometimes,
;       however, it is desirable to invoke XOBJVIEW, /MODAL in a
;       program that uses no other IDL widgets. For this reason,
;       XOBJVIEW allows itself to be invoked MODAL with no group leader.
;
;
;  EXAMPLE 1:
;       IDL> oObj = obj_new('IDLgrSurface', dist(30))
;       IDL> xobjview, oObj, /bloc
;       IDL> obj_destroy, oObj
;
;  EXAMPLE 2:
;       oObj = obj_new('IDLgrSurface', dist(30))
;       xobjview, oObj, tlb=tlb
;       ...
;       oObj->SetProperty, color=[255, 0, 0]
;       xobjview, refresh=tlb
;       ...
;       obj_destroy, oObj
;
;  MODIFICATION HISTORY:
;       9/1999  PCS   - created.
;       11/1999 PCS   - added Update method and related changes.
;       01/2000 PCS   - added TITLE keyword and functionality.
;       07/2000 DMS   - added REFRESH, TLB and BACKGROUND keywords.
;-
;
;--------------------------------------------------------------------
pro xobjview__free, var
compile_opt idl2, hidden

on_error, 2

case size(var, /tname) of
    'OBJREF': obj_destroy, var
    'POINTER': ptr_free, var
    'LONG': widget_control, var, /destroy
    else:
    endcase

end
;--------------------------------------------------------------------
pro xobjview_cleanup, wID
compile_opt hidden

widget_control, wID, get_uvalue=pState

if (*pState).group_leader_is_fabricated then begin
    if widget_info(*(*pState).pGroupLeader, /valid_id) then begin
        widget_control, *(*pState).pGroupLeader, /destroy
        endif
    endif

obj_destroy, (*pState).oObjviewWid
obj_destroy, (*pState).oTestSurface
ptr_free, (*pState).pGroupLeader
ptr_free, pState
end
;--------------------------------------------------------------------
pro xobjview_event, event
compile_opt hidden
;
;Handle resize events.
;
if tag_names(event, /structure_name) eq 'WIDGET_BASE' then begin
    widget_control, event.top, get_uvalue=pState
    widget_control, /hourglass
    ;; >>>>>>>>>>>>>>>>>>>>>>>>
    extra_base = widget_info(event.top, find_by_uname='base_tlb')
    if (widget_info(extra_base, /valid_id)) then begin
       geom = widget_info(extra_base, /geometry)
    endif else geom = { ysize: 0 }
    pad_x = 4 ; Estimate.
    pad_y = 25
    xsize = event.x - pad_x
    ysize = event.y - (pad_y + geom.ysize)
    ;; <<<<<<<<<<<<<<<<<<<<<<<<<
    (*pState).oObjviewWid->SetProperty, xsize=xsize, ysize=ysize
    endif
end
;--------------------------------------------------------------------
pro xobjview, $
    oObj, $           ; IN
    xsize=xsize, $    ; IN: (opt) pixels.  Size of drawable.
    ysize=ysize, $    ; IN: (opt) pixels.  Size of drawable.
    stationary=oStationary, $
    group=group_leader, $   ; IN: (opt) group leader widget
    block=block, $
    just_reg=just_reg, $
    modal=modal, $
    scale=scale, $    ; IN: (opt) Scale size of initial view.
    test=test, $      ; IN: (opt) Bool.  Don't require oObj arg.
    title=title, $
    xoffset=xoffset, $; IN: (opt)
    yoffset=yoffset, $; IN: (opt)
    tlb=tlb, $        ;OUT: (opt), widget id of top level base
    refresh=refresh_tlb, $  ; IN: (opt) top-level base of existing XOBJVIEW.
    double_view=double_view, $  ; IN: (opt) set for high precision view.
    background=background, $    ; IN: (opt) view background color.
    renderer=renderer, $        ; IN: (opt) 1 = IDL's software renderer.
    use_instancing=use_instancing, $ ; IN: (opt) For expose events.
    debug=debug       ; IN: (opt)

on_error, 2 ; Return to caller on error.

catch, error_status
if error_status ne 0 then begin
    catch, /cancel
;
;   Clean up.
;
    xobjview__free, oTestSurface
    if keyword_set(group_leader_is_fabricated) then begin
        xobjview__free, group_leader
        ptr_free, ptr_new(group_leader, /no_copy)
        endif
    xobjview__free, oWindow
    xobjview__free, oScene
    if n_elements(oOriginalObj) gt 0 then $
        oObj = oOriginalObj
    if keyword_set(do_destroy_view) then begin
        xobjview__free, oCurrentView
        endif
;
;   Re-throw the error.
;
    message, !error_state.msg + ' ' + !error_state.sys_msg
    endif

if keyword_set(debug) then begin
    on_error, 0
    catch, /cancel
    endif

if n_elements(refresh_tlb) gt 0 then begin
    if not widget_info(refresh_tlb, /valid) then $
        message, 'Specified widget for REFRESH is invalid.', /noname

    widget_control, refresh_tlb, get_uvalue = pState
    if obj_valid((*pState).oObjviewWid) then $
        (*pState).oObjviewWid->Draw, /hourglass
    return
    end

if keyword_set(test) then begin
    if n_elements(oObj) gt 0 then $
        oOriginalObj = oObj
    oObj = obj_new('IDLgrSurface', $
        beselj(shift(dist(40), 20, 20) / 2,0) * 20, $
        color=[60, 60, 255], $
        style=2, $
        shading=1, $
        name='Test Surface'$
        )
    oTestSurface = oObj
    endif $
else begin
    oTestSurface = obj_new()
    endelse

if n_elements(oObj) eq 0 then begin
    if arg_present(oObj) then $
        message, 'Argument is undefined.', /noname $
    else $
        message, 'requires an argument.', /noname
    endif

if size(oObj[0], /tname) ne 'OBJREF' then $
    message, 'Argument must be of object reference type.', /noname

if n_elements(uniq(oObj, sort(oObj))) ne n_elements(oObj) then $
    message, 'Array argument must contain unique values.', /noname

if obj_isa(oObj[0], 'IDLgrView') then $
    message, 'Note: View arguments are an undocumented feature.', /inform

if n_elements(oStationary) gt 0 then begin
    if size(oStationary, /tname) ne 'OBJREF' then $
        message, $
            'Keyword STATIONARY must be of object reference type.', $
            /noname
    end

if n_elements(scale) eq 3 then begin
;
;   The restriction imposed by the following message is not technically
;   necessary, but it serves to simplify XOBJVIEW's command interface.
;   It leaves one way of doing independent x, y, z scaling, instead of
;   two ways that would interact with each other.
;
    message, 'SCALE must be a 1-element value.  (If you desire ' + $
        'different ammounts of scale in x, y & z, use an IDLgrModel.)'
    end

if n_elements(group_leader) ne 0 then begin
    if not widget_info(group_leader, /valid_id) then begin
        message, 'Specified Group Leader is not valid.', /noname
        endif
    endif $
else begin
    if keyword_set(modal) then begin
;
;       Modal widgets require a group leader.  A group leader was not
;       specified, so fabricate an invisible one.
;
        group_leader = widget_base(map=0)
        group_leader_is_fabricated = 1b
        endif
    endelse
;
;Create widgets.
;
if keyword_set(modal) then begin
    tlb = widget_base( $
        /column, $
        xpad=0, $
        ypad=0, $
        xoffset=xoffset, $
        yoffset=yoffset, $
        title=n_elements(title) eq 0 ? 'Xobjview' : title, $
        /tlb_size_events, $
        /modal, $
        group_leader=group_leader $
        )
    endif $
else begin
    tlb = widget_base( $
        /column, $
        xpad=0, $
        ypad=0, $
        xoffset=xoffset, $
        yoffset=yoffset, $
        title=n_elements(title) eq 0 ? 'Xobjview' : title, $
        /tlb_size_events, $
        group_leader=group_leader, $
        mbar=mbar $
        )
    endelse

oObjviewWid = obj_new('IDLexObjviewWid', $
    tlb, $
    oObj, $
    menu_parent=mbar, $
    scale=scale, $
    draw_xsize=xsize, $
    draw_ysize=ysize, $
    background=background, $
    stationary=oStationary, $
    double_view=double_view, $
    renderer=renderer, $
    use_instancing=use_instancing, $
    /include_refresh_button, $
    /include_full_reset_button, $
    debug=debug $
    )
if not obj_valid(oObjviewWid) then begin
    message, !error_state.msg + ' ' + !error_state.sys_msg, /noname
    end

widget_control, tlb, set_uvalue=ptr_new({ $
    oObjviewWid: oObjviewWid, $
    oTestSurface: obj_valid(oTestSurface) ? oTestSurface : obj_new(), $
    pGroupLeader: ptr_new(group_leader), $
    group_leader_is_fabricated: keyword_set(group_leader_is_fabricated) $
    })
widget_control, tlb, /realize
xmanager, $
    "xobjview", $
    tlb, $
    just_reg=keyword_set(just_reg),$
    no_block=keyword_set(block) eq 0, $
    cleanup='xobjview_cleanup'

if keyword_set(group_leader_is_fabricated) then begin
;
;   Leave GROUP_LEADER parameter like we found it: undefined.
;
    ptr_free, ptr_new(group_leader, /no_copy)
    endif
;
;Leave oObj argument unchanged.
;
if n_elements(oOriginalObj) gt 0 then $
    oObj = oOriginalObj

;;>>>>>>>>>>>>>>>>>>>
if (n_elements(mbar) ne 0 and widget_info(mbar, /valid_id)) then begin
   widget_control, mbar, set_uname='xobjview:mbar'
endif
;;<<<<<<<<<<<<<<<<<<<

end

