;$Id: //depot/idl/IDL_64/idldir/lib/utilities/idlexobjviewwid__define.pro#1 $
;
;  Copyright (c) 1997-2007, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
;
function IDLexObjviewWid::DefaultViewType
compile_opt hidden

return, 'IDLexObjView'
end
;--------------------------------------------------------------------
function IDLexObjviewWid::Rewrite, event
compile_opt hidden

if not self.own_menu then begin
    if event.handler eq self.wViewButton $
    or event.handler eq self.wFileButton $
    or event.handler eq self.wEditButton then begin
        return, event
        end
    end

return, self->IDLexWidget::Rewrite(event)
end
;--------------------------------------------------------------------
function IDLexObjviewWid::widget_button, $
    value=value, $
    menu=menu, $    ; IN: (opt) if set, this button will have children.
    trunk=trunk     ; IN: (opt) 'file', 'edit', 'view'
compile_opt hidden

if n_elements(menu) gt 0 then begin
    case strlowcase(trunk) of
        'view': wParent = self.wViewButton
        'file': wParent = self.wFileButton
        'edit': wParent = self.wEditButton
        else: wParent = self.wMenuBase
        endcase
    end

result = widget_button(wParent, value=value, menu=menu)
return, result
end
;--------------------------------------------------------------------
pro IDLexObjviewWid::Draw, $
    oDestination, $         ; IN: (opt)
    oPicture, $             ; IN: (opt)
    in_motion=in_motion_, $ ; IN: (opt) If set, don't create_instance.
    vector=vector, $
    hourglass=hourglass
compile_opt hidden

if n_elements(oDestination) eq 0 and obj_valid(self.oWindow) eq 0 then $
    return

oDest = n_elements(oDestination) eq 0 ? self.oWindow : oDestination
oPict = n_elements(oPicture) eq 0 ? self.oViewgroup : oPicture

in_motion = self.use_instancing eq 0
if n_elements(in_motion_) gt 0 then begin
    in_motion = in_motion_
    end
;
;Flush and print any accumulated math errors
;
void = check_math(/print)
;
;On some platforms, when IDLgrWindow::Draw is invoked, math
;errors (e.g. "% Program caused arithmetic error: Floating illegal
;operand") are printed.   Silently accumulate any such math errors.
;
orig_except = !except
!except = 0
;
;Draw.
;
if keyword_set(hourglass) then begin
    widget_control, /hourglass
    end

if n_elements(vector) gt 0 then begin
    oDest->Draw, oPict, vector=vector
    endif $
else begin
    if in_motion then begin
        oDest->Draw, oPict
        self.instance_exists = 0b
        end $
    else begin
        oDest->Draw, oPict, /create_instance
        oDest->Draw, self.oTransparentView, /draw_instance
        self.instance_exists = 1b
        end
    end
;
;Silently flush any accumulated math errors.
;
void = check_math()
;
;Restore original math error behavior.
;
!except = orig_except
end
;--------------------------------------------------------------------
pro IDLexObjviewWid::WriteImage, $
    filename, $             ; IN
    format, $               ; IN: (e.g. "bmp", "jpeg", etcetera.)
    dimensions=dimensions   ; IN: (opt) Use this pixel size. [x, y].

compile_opt hidden

self.oViewGroup->WriteImage, $
    filename, $
    format, $
    self.oWindow, $
    dimensions=dimensions

end
;--------------------------------------------------------------------
pro IDLexObjviewWid::Rotate, $
    axis, $     ; IN: three-element vector about which to rotate.
    angle, $    ; IN: The amount (measured in degrees) of rotation
    premultiply=premultiply ; IN: (opt) if set, do "data-centric" rotation.

on_error, 2

self.oViewgroup->CallAll, $
    'Rotate', axis, angle, premultiply=premultiply, isa='IDLexObjview'
self->Draw

end
;--------------------------------------------------------------------
pro IDLexObjviewWid::SetProperty, $
    xsize=xsize, $ ; Approximate size of self's widgets.
    ysize=ysize, $ ; Approximate size of self's widgets.
    draw_xsize=draw_xsize, $
    draw_ysize=draw_ysize, $
    background=background, $
    double_view=double_view, $
    drag_quality=drag_quality, $
    _extra=e

compile_opt hidden

pad = 3 ; Estimate.
if n_elements(ysize) gt 0 then begin
    ys = ysize ; Local copy so as not to change parameter.
    if self.own_toolbar then begin
        toolbar_geometry = widget_info(self.wToolbarBase, /geometry)
        ys = ys - toolbar_geometry.scr_ysize
        endif
    if self.own_menu then begin
        menu_geometry = widget_info(self.wMenuBase, /geometry)
        ys = ys - menu_geometry.scr_ysize
        endif
    self.instance_exists = 0b
    widget_control, self.wDraw, ysize=ys - pad
    endif

if n_elements(xsize) gt 0 then begin
    self.instance_exists = 0b
    widget_control, self.wDraw, xsize=xsize - pad
    endif

if n_elements(double_view) gt 0 then begin
    self.oViewgroup->CallAll, 'SetProperty', double=double_view
    end

case n_elements(background) of
    0: ; do nothing
    3: begin
        self.oViewgroup->CallAll, 'SetProperty', color=background
        self->Draw, /hourglass

        widget_control, self.wBlackButton, set_button=0
        widget_control, self.wCharcoalButton, set_button=0
        widget_control, self.wWhiteButton, set_button=0
        widget_control, self.wGrayButton, set_button=0
        end
    1: begin
        if size(background, /tname) ne 'STRING' $
        then $
            message, 'COLOR must be a string or a 3-element [r,g,b] array.'

        case strlowcase(background) of
            'black': clr = [0, 0, 0]
            'charcoal': clr = [80, 80, 80]
            'gray': clr = [127, 127, 127]
            'white': clr = [255, 255, 255]
            endcase

        self.oViewgroup->CallAll, 'SetProperty', color=clr
        self->Draw, /hourglass

        widget_control, $
            self.wBlackButton, $
            set_button=strlowcase(background) eq 'black'
        widget_control, $
            self.wCharcoalButton, $
            set_button=strlowcase(background) eq 'charcoal'
        widget_control, $
            self.wWhiteButton, $
            set_button=strlowcase(background) eq 'white'
        widget_control, $
            self.wGrayButton,  $
            set_button=strlowcase(background) eq 'gray'
        end
    else: $
        message, 'COLOR must be a 1-element (string) or 3-element ' + $
            '(numerical) value.'
    endcase

if n_elements(draw_xsize) ne 0 then begin
    self.instance_exists = 0b
    widget_control, self.wDraw, xsize=draw_xsize
    end
if n_elements(draw_ysize) ne 0 then begin
    self.instance_exists = 0b
    widget_control, self.wDraw, ysize=draw_ysize
    end

if n_elements(drag_quality) gt 0 then begin
    self.drag_quality = drag_quality
    self->UpdateDragQualMenu
    end

self->IDLexWidget::SetProperty, _extra=e
end
;--------------------------------------------------------------------
pro IDLexObjviewWid::UpdateDragQualMenu
compile_opt hidden

widget_control, self.wDragQualLow,    set_button=self.drag_quality eq 0
widget_control, self.wDragQualMedium, set_button=self.drag_quality eq 1
widget_control, self.wDragQualHigh,   set_button=self.drag_quality eq 2

end
;--------------------------------------------------------------------
pro IDLexObjViewWid::OnRealize
compile_opt hidden

widget_control, self.wDraw, get_value=oWindow
self.oWindow = oWindow
self.oViewgroup->CallAll, 'Reset', /full, oWindow
end
;--------------------------------------------------------------------
pro IDLexObjViewWid::OnExpose
compile_opt hidden

if self.use_instancing and self.instance_exists then begin
    self.oWindow->Draw, self.oTransparentView, /draw_instance
    end $
else begin
    ; Do not use /hourglass because it could mess up button events.
    self->Draw
    end

end
;--------------------------------------------------------------------
pro IDLexObjViewWid::OnMouseDown, event
compile_opt hidden

case event.press of
    4: begin ; Right mouse-button.
        end
    2: begin ; Middle mouse-button.
        end
    1: begin ; Left mouse button.
        self.oViewgroup->SetCurrent, event
        oCurrent = self.oViewgroup->Get(/current)
        if obj_valid(oCurrent) then begin
            if obj_isa(oCurrent, 'IDLexObjview') then $
                void = oCurrent->Update(event)
            end
        if self.mode eq 'select' then begin
            self.oViewgroup->CallCurrent, $
                'GetProperty', $
                isa='IDLexObjview', $
                selected=oSelected

            if (obj_valid(oSelected))[0] then begin
                oSelected[0]->GetProperty, name=name
                if name eq '' then $
                    name = obj_class(oSelected[0])
                endif $
            else begin
                name = ' '
                endelse

            widget_control, $
                self.wLabel, $
                set_value=name
            endif $
        else begin
            self.oWindow->SetProperty, $
                quality=self.drag_quality
            self->Draw, /in_motion
            endelse
        end
    else:
    endcase
end
;--------------------------------------------------------------------
pro IDLexObjViewWid::OnMouseMove, event
compile_opt hidden

oCurrent = self.oViewgroup->Get(/current)
if obj_valid(oCurrent) then begin
    if obj_isa(oCurrent, 'IDLexObjview') then begin
        if oCurrent->Update(event) then begin
            self->Draw, /in_motion
            end
        end
    end

end
;--------------------------------------------------------------------
pro IDLexObjViewWid::OnMouseUp, event
compile_opt hidden

case event.release of
    4: begin ; Right mouse-button.
        end
    2: begin ; Middle mouse-button.
        end
    1: begin ; Left mouse button.
        self.oWindow->setproperty, quality=2
        oCurrent = self.oViewgroup->Get(/current)
        if obj_valid(oCurrent) then begin
            if obj_isa(oCurrent, 'IDLexObjview') then begin
                if oCurrent->Update(event) then begin
                    if self.mode ne 'select' then begin
                        self->Draw, /hourglass
                        end
                    end
                end
            end
        end
    else:
    endcase

end
;--------------------------------------------------------------------
function IDLexObjViewWid::HandleEvent, event
compile_opt hidden

on_error, 2 ; Return to caller on error.

catch, error_status
if error_status ne 0 then begin
    catch, /cancel
    void = dialog_message( $
        dialog_parent=event.top, $
        title='Error', $
        /error, $
        !error_state.msg + ' ' + !error_state.sys_msg $
        )
    return, self->Rewrite(event)
    endif

if keyword_set(self.debug) then begin
    on_error, 0
    catch, /cancel
    endif

if not obj_valid(self.oWindow) then begin
    self->OnRealize
    end

case event.id of

    self.wReset: begin
        self.oViewgroup->CallAll, 'Reset', $
            self.oWindow, $
            isa='IDLexObjview'
        self->Draw, /hourglass
        end

    self.wRotate: begin
        ; Set the button state if called manually.
        if (WIDGET_INFO(self.wRotate, /BUTTON_SET) ne event.select) then $
            WIDGET_CONTROL, self.wRotate, SET_BUTTON=event.select
        widget_control, self.wLabel, set_value=' '
        self.mode = 'rotate'
        self.oViewgroup->CallAll, $
            'SetProperty', $
            mode=self.mode, $
            ISA='IDLexObjview'
        end

    self.wZoom: begin
        ; Set the button state if called manually.
        if (WIDGET_INFO(self.wZoom, /BUTTON_SET) ne event.select) then $
            WIDGET_CONTROL, self.wZoom, SET_BUTTON=event.select
        widget_control, self.wLabel, set_value=' '
        self.mode = 'zoom'
        self.oViewgroup->CallAll, $
            'SetProperty', $
            mode=self.mode, $
            isa='IDLexObjview'
        end

    self.wPan: begin
        ; Set the button state if called manually.
        if (WIDGET_INFO(self.wPan, /BUTTON_SET) ne event.select) then $
            WIDGET_CONTROL, self.wPan, SET_BUTTON=event.select
        widget_control, self.wLabel, set_value=' '
        self.mode = 'translate'
        self.oViewgroup->CallAll, $
            'SetProperty', $
            mode=self.mode, $
            isa='IDLexObjview'
        end

    self.wSelect: begin

        ; Set the button state if called manually.
        if (WIDGET_INFO(self.wSelect, /BUTTON_SET) ne event.select) then $
            WIDGET_CONTROL, self.wSelect, SET_BUTTON=event.select
        self.mode = 'select'
        self.oViewgroup->CallAll, $
            'SetProperty', $
            mode=self.mode, $
            isa='IDLexObjview'
        end

    self.wQuitButton: begin
        evnt = self->Rewrite(event)
        widget_control, event.top, /destroy
        return, evnt
        end

    self.wVRMLButton: begin
        void = self.oViewgroup->DialogWriteVRML( $
            self.oWindow, $
            dialog_parent=event.top $
            )
        end

    self.wExportImageButton: begin
        void = self.oViewgroup->DialogWriteImage( $
            self.oWindow, $
            dialog_parent=event.top $
            )
        end

    self.wClipboardButton: begin
        self.oViewgroup->CopyToClipboard, self.oWindow
        end

    self.wClipboardVectorButton: begin
        self.oViewgroup->CopyToClipboard, self.oWindow, /vector
        end

    self.wPrintButton: begin
        self.oViewgroup->DialogPrint, $
            self.oWindow, $
            /hourglass, $
            dialog_parent=event.top
        end

    self.wDragQualLow: begin
        self.drag_quality = 0
        self->UpdateDragQualMenu
        end

    self.wDragQualMedium: begin
        self.drag_quality = 1
        self->UpdateDragQualMenu
        end

    self.wDragQualHigh: begin
        self.drag_quality = 2
        self->UpdateDragQualMenu
        end

    self.wBlackButton: begin
        self->SetProperty, background='black'
        end

    self.wCharcoalButton: begin
        self->SetProperty, background='charcoal'
        end

    self.wWhiteButton: begin
        self->SetProperty, background='white'
        end

    self.wGrayButton: begin
        self->SetProperty, background='gray'
        end

    self.wRefreshButton: begin
        self->Draw, /hourglass
        end

    self.wFullResetButton: begin
        self.oViewgroup->CallAll, 'Reset', $
            /full, $
            self.oWindow, $
            isa='IDLexObjview'
        self->Draw, /hourglass
        end

    self.wDraw: begin
        case event.type of
            4: self->OnExpose
            0: self->OnMouseDown, event
            2: self->OnMouseMove, event
            1: self->OnMouseUp, event
            else:
            endcase
        end

    else:

    endcase

    return, self->Rewrite(event)
end
;--------------------------------------------------------------------
pro IDLexObjviewWid::CleanupNoninheritedMembers
compile_opt hidden

if obj_valid(self.oViewgroup) then begin
    if self.oViewgroup->Count() gt 0 then begin
        if not self.do_destroy_view then begin
            self.oViewgroup->Remove, /all
            end
        end
    end

obj_destroy, self.oTransparentview
obj_destroy, self.oViewgroup
obj_destroy, self.oWindow
ptr_free, self.pBitmapPath

end
;--------------------------------------------------------------------
pro IDLexObjviewWid::cleanup
compile_opt hidden

self->IDLexObjviewWid::CleanupNoninheritedMembers
self->IDLexWidget::cleanup

end
;--------------------------------------------------------------------
function IDLexObjviewWid::WidgetDrawBase
compile_opt hidden

return, widget_base(self.wBase)
end
;--------------------------------------------------------------------
function IDLexObjViewWid::init, $
    wParent, $                  ; IN: widget base
    oSubjects, $                ; IN: models, views or atomic graphics objects.
    menu_parent=wMenuBase, $    ; IN: (opt) button or base widget.
    toolbar_parent=wToolbarBase, $  ; IN: (opt) base widget.
    stationary=oStationary, $   ;
    double_view=double_view, $  ; IN: (opt) If set, use high precision view.
    scale=scale, $, $
    draw_xsize=draw_xsize, $    ; IN: (opt) pixel size of draw widget.
    draw_ysize=draw_ysize, $    ; IN: (opt) pixel size of draw widget.
    renderer=renderer, $        ; IN: (opt) 0==OpenGL (default). 1==software.
    background=background, $
    include_refresh_button=include_refresh_button, $
    include_full_reset_button=include_full_reset_button, $
    use_instancing=use_instancing, $ ; IN: (opt) For expose events.
    debug=debug, $              ; IN: (opt)
    _extra=e                    ; e.g. FRAME

compile_opt hidden

on_error, ([2, 0])[keyword_set(debug)]

catch, error_stat
if error_stat ne 0 then begin
    catch, /cancel
    return, 0
    end
if keyword_set(debug) then begin
    catch, /cancel
    end

self.debug = keyword_set(debug)

if not self->IDLexWidget::init(wParent, _extra=e) then $
    message, 'Failed to init IDLexWidget part of self.'

catch, error_stat
if error_stat ne 0 then begin
    catch, /cancel
    self->IDLexObjviewWid::CleanupNoninheritedMembers
    self->IDLexWidget::cleanup
    return, 0
    end
if keyword_set(debug) then begin
    catch, /cancel
    end

self.use_instancing = keyword_set(use_instancing)

if n_elements(scale) eq 0 then begin
    scale = 0
    endif
if n_elements(scale) ne 1 and n_elements(scale) ne 3 then $
    message, 'SCALE must be scalar, 1-element or 3-element array.'
if min(scale) lt 0 then $
    message, 'SCALE must be positive.'
;
;Create widgets.
;
self.wBase = widget_base( $
    self.wIDBase, $
    /column, $
    event_func='IDLexWidget__HandleEvent', $
    uvalue=self, $
    notify_realize='IDLexWidget__OnRealize' $
    )

if n_elements(wMenuBase) gt 0 then begin
    if size(wMenuBase, /n_dimensions) gt 0 then $
        message, 'MENU_PARENT must be a scalar.'
    if size(wMenuBase, /tname) ne 'LONG' then $
        message, 'MENU_PARENT must be of type LONG.'
    if not widget_info(wMenuBase, /valid) then $
        message, 'The given MENU_PARENT is not a valid widget.'
    if widget_info(wMenuBase, /type) ne 0 $
    and widget_info(wMenuBase, /type) ne 1 then $
        message, 'The given MENU_PARENT must be a button or base.'
    self.wMenuBase = wMenuBase
    end $
else begin
    self.wMenuBase = widget_base(self.wBase, /row, /frame)
    self.own_menu = 1b
    end

if n_elements(wToolbarBase) gt 0 then begin
    if size(wToolbarBase, /n_dimensions) gt 0 then $
        message, 'TOOLBAR_PARENT must be a scalar.'
    if size(wToolbarBase, /tname) ne 'LONG' then $
        message, 'TOOLBAR_PARENT must be of type LONG.'
    if not widget_info(wToolbarBase, /valid) then $
        message, 'The given TOOLBAR_PARENT is not a valid widget.'
    if widget_info(wToolbarBase, /type) ne 0 $
    and widget_info(wToolbarBase, /type) ne 1 then $
        message, 'The given TOOLBAR_PARENT must be a base widget.'
    self.wToolbarBase = wToolbarBase
    end $
else begin
    self.wToolbarBase = widget_base(self.wBase, /row, /frame)
    self.own_toolbar = 1b
    end
;
;Widget UNAMES assigned below are used internally to help
;automate testing of this code.  A prefix is given to the UNAMES
;to help ensure that they are unique from UNAMES in other programs.
;
prefix = 'xobjview:'
;
;Create the menu bar.
;
self.wFileButton = widget_button( $
    self.wMenuBase, $
    value='File', $
    /menu $
    )
if not self.own_menu then begin
    widget_control, $
        self.wFileButton, $
        event_func='IDLexWidget__HandleEvent', $
        set_uvalue=self
    end

    self.wExportImageButton = widget_button( $
        self.wFileButton, $
        uname=prefix + 'Export Image', $
        value='Export Image...' $
        )
    self.wVRMLButton = widget_button( $
        self.wFileButton, $
        uname=prefix + 'Export VRML', $
        value="Export VRML..." $
        )
    self.wPrintButton = widget_button( $
        self.wFileButton, $
        uname=prefix + 'Print', $
        value='Print...' $
        )
    self.wQuitButton = widget_button( $
        self.wFileButton, $
        /separator, $
        uname=prefix + 'Quit', $
        value='Quit' $
        )

self.wEditButton = widget_button( $
    self.wMenuBase, $
    /menu, $
    value='Edit' $
    )
if not self.own_menu then begin
    widget_control, $
        self.wEditButton, $
        event_func='IDLexWidget__HandleEvent', $
        set_uvalue=self
    end
    self.wClipboardVectorButton = widget_button( $
        self.wEditButton, $
        uname=prefix + 'Vectors to Clipboard', $
        value='Copy Scene to Clipboard (as Vectors)' $
        )
    self.wClipboardButton = widget_button( $
        self.wEditButton, $
        uname=prefix + 'Raster to Clipboard', $
        value='Copy Scene to Clipboard (as Raster)' $
        )

self.wViewButton = widget_button( $
    self.wMenuBase, $
    value='View', $
    /menu $
    )
if not self.own_menu then begin
    widget_control, $
        self.wViewButton, $
        event_func='IDLexWidget__HandleEvent', $
        set_uvalue=self
    end

    wDragQualButton = widget_button( $
        self.wViewButton, $
        value="Drag Quality", $
        /menu $
        )
        self.wDragQualLow = widget_button( $
            wDragQualButton, /check,$
            uname=prefix + 'Low', $
            value='Low' $
            )
        self.wDragQualMedium = widget_button( $
            wDragQualButton, /check,$
            uname=prefix + 'Medium', $
            value='Medium' $
            )
        self.wDragQualHigh = widget_button( $
            wDragQualButton, /check,$
            uname=prefix + 'High', $
            value='High' $
            )
    wBackgroundButton = widget_button( $
        self.wViewButton, $
        value='Background', $
        /menu $
        )
        self.wBlackButton = widget_button( $
            wBackgroundButton, /check,$
            uname=prefix + 'Black', $
            value='Black' $
            )
        self.wCharcoalButton = widget_button( $
            wBackgroundButton, /check,$
            uname=prefix + 'Charcoal', $
            value='Charcoal' $
            )
        self.wGrayButton = widget_button( $
            wBackgroundButton, /check,$
            uname=prefix + 'Gray', $
            value='Gray' $
            )
        self.wWhiteButton = widget_button( $
            wBackgroundButton, /check,$
            uname=prefix + 'White', $
            value='White' $
            )
    widget_control, self.wWhiteButton, set_button=1
    if keyword_set(include_full_reset_button) then $
        self.wFullResetButton = widget_button( $
            self.wViewButton, $
            uname=prefix + 'Full Reset', $
            value='Full Reset' $
            )
    if keyword_set(include_refresh_button) then $
        self.wRefreshButton = widget_button( $
            self.wViewButton, $
            uname=prefix + 'Refresh', $
            value='Refresh Display' $
            )

if lmgr(/demo) then begin
    widget_control, self.wPrintButton, sensitive=0
    widget_control, self.wExportImageButton, sensitive=0
    widget_control, self.wClipboardButton, sensitive=0
    widget_control, self.wVRMLButton, sensitive=0
    endif
;
;Create the Toolbar.
;
self.pBitmapPath = ptr_new(['resource', 'bitmaps'])

wResetBase = WIDGET_BASE(self.wToolbarBase, uname='xobjview:_resetbase', /ROW, /TOOLBAR)
;;wResetBase = WIDGET_BASE(self.wToolbarBase, /TOOLBAR)

self.wReset = widget_button( $
    wResetBase, $
    uname=prefix + 'reset', $
    value=filepath('reset.bmp', subdirectory=*self.pBitmapPath), $
    /bitmap, TOOLTIP='Reset' $
    )
if not self.own_toolbar then begin
    widget_control, $
        self.wReset, $
        event_func='IDLexWidget__HandleEvent', $
        set_uvalue=self
    end

void = widget_base(self.wToolbarBase, /row, scr_xsize=8)

wExcButtonBase = WIDGET_BASE(self.wToolbarBase, $
    /ROW, /TOOLBAR, /EXCLUSIVE, SPACE=0)

self.mode = 'rotate'
self.wRotate = widget_button( $
    wExcButtonBase, $
    uname=prefix + 'rotate', $
    value=filepath('rotate.bmp', subdirectory=*self.pBitmapPath), $
    /bitmap, TOOLTIP='Rotate' $
    )
WIDGET_CONTROL, self.wRotate, /SET_BUTTON
if not self.own_toolbar then begin
    widget_control, $
        self.wRotate, $
        event_func='IDLexWidget__HandleEvent', $
        set_uvalue=self
    end

self.wZoom = widget_button( $
    wExcButtonBase, $
    uname=prefix + 'zoom', $
    value=filepath('zoom.bmp', subdirectory=*self.pBitmapPath), $
    /bitmap, TOOLTIP='Zoom' $
    )
if not self.own_toolbar then begin
    widget_control, $
        self.wZoom, $
        event_func='IDLexWidget__HandleEvent', $
        set_uvalue=self
    end

self.wPan = widget_button( $
    wExcButtonBase, $
    uname=prefix + 'pan', $
    value=filepath('pan.bmp', subdirectory=*self.pBitmapPath), $
    /bitmap, TOOLTIP='Pan' $
    )
if not self.own_toolbar then begin
    widget_control, $
        self.wPan, $
        event_func='IDLexWidget__HandleEvent', $
        set_uvalue=self
    end

self.wSelect = widget_button( $
    wExcButtonBase, $
    uname=prefix + 'select', $
    value=filepath('select.bmp', subdirectory=*self.pBitmapPath), $
    /bitmap, TOOLTIP='Select' $
    )
if not self.own_toolbar then begin
    widget_control, $
        self.wSelect, $
        event_func='IDLexWidget__HandleEvent', $
        set_uvalue=self
    end

self.wLabel = widget_label(self.wToolbarBase, value=' ', /dynamic_resize)
;
self.wDraw = widget_draw( $
    self->WidgetDrawBase(), $
    xsize=n_elements(draw_xsize) gt 0 ? draw_xsize : 400, $
    ysize=n_elements(draw_ysize) gt 0 ? draw_ysize : 400, $
    /button_events, $
    /motion_events, $
    retain=0, $
    renderer=renderer, $
    /expose_events, $
    uname=prefix + 'draw', $
    graphics_level=2 $
    )

self.oViewgroup = obj_new('IDLexViewgroup', debug=debug)
;
;Obtain one or more views.
;
if obj_isa(oSubjects[0], 'IDLgrView') then begin
    found_exobjview = 0b
    for i=0L,n_elements(oSubjects)-1 do begin
        if not obj_isa(oSubjects[i], 'IDLgrView') then $
            message, 'Array argument containing a' + $ ; Too strict?
                'reference to a view must ' + $
                'contain only references to views.'
        found_exobjview = $
            found_exobjview or obj_isa(oSubjects[i], 'IDLexObjview')
        endfor
    if not found_exobjview then $  ; Too strict?
        message, 'At least one view must be an IDlexObjview.'

    self.oViewgroup->Add, oSubjects
    self.oViewgroup->SetProperty, current=0
    self.do_destroy_view = 0b
    endif $
else begin
    self.oViewgroup->Add, $
        obj_new(self->DefaultViewType(), debug=debug), $
        /make_current
    self.oViewgroup->CallCurrent, 'Add', oSubjects, /alias
    self.do_destroy_view = 1b

    if n_elements(oStationary) gt 0 then begin
        self.oViewgroup->CallCurrent, $
            'Add', $
            oStationary, $
            /stationary, $
            /alias
        end
    endelse

self.oViewgroup->CallAll, 'SetProperty', scale=scale

self.oTransparentView = obj_new('IDLgrView', /transparent)
self->SetProperty, background=background, double_view=double_view
self.drag_quality = 1
widget_control, self.wDragQualMedium, set_button=1
return, 1 ; Success.
end
;--------------------------------------------------------------------
pro IDLexObjViewWid__define

compile_opt hidden
on_error, 2

struct_hide, {IDLexObjViewWid, $
    inherits IDLexWidget, $
    wBase: 0L, $
    wExportImageButton: 0L, $
    wVRMLButton: 0L, $
    wPrintButton: 0L, $
    wQuitButton: 0L, $
    wClipboardVectorButton: 0L, $
    wClipboardButton: 0L, $
    wDragQualLow: 0L, $
    wDragQualMedium: 0L, $
    wDragQualHigh: 0L, $
    wBlackButton: 0L, $
    wCharcoalButton: 0L, $
    wWhiteButton: 0L, $
    wGrayButton: 0L, $
    wRefreshButton: 0L, $
    wFullResetButton: 0L, $
    wToolbarBase: 0L, $
    wMenuBase: 0L, $
    wFileButton: 0L, $
    wViewButton: 0L, $
    wEditButton: 0L, $
    wReset: 0L, $
    wRotate: 0L, $
    wZoom: 0L, $
    wPan: 0L, $
    wSelect: 0L, $
    wLabel: 0L, $
    wDraw: 0L, $
    pBitmapPath: ptr_new(), $
    debug: 0b, $
    mode: '', $
    oTransparentView: obj_new(), $
    use_instancing: 0b, $
    instance_exists: 0b, $
    drag_quality: 0, $
    do_destroy_view: 0b, $
    own_menu: 0b, $
    own_toolbar: 0b, $
    oWindow: obj_new(), $
    oViewgroup: obj_new() $
    }

end
