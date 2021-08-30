;$Id: //depot/idl/IDL_64/idldir/lib/utilities/idlexviewgroup__define.pro#1 $
;
;  Copyright (c) 1997-2007, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
;+
; NAME:
;   IDLexViewGroup
;
; PURPOSE:
;   Provide a subclass of IDLgrViewgroup that maintains the notion of a
;   "current" child, provides methods for operating on that child,
;   and has other convenience methods for printing, exporting to VRML, etc.
;
; CREATION:
;   oViewGroup = obj_new('IDLexViewGroup')
;
;-
pro IDLexViewgroup::Add, $
    oObj, $
    make_current=make_current, $
    position=position

compile_opt hidden

oCurrent = self->Get(/current)
;
;Perform a regular Add.
;
if n_elements(position) gt 0 then begin
    self->IDLgrViewgroup::Add, oObj, position=position
    end $
else begin
    self->IDLgrViewgroup::Add, oObj
    end
;
;Update self.current.
;
if keyword_set(make_current) then begin
    if n_elements(position) eq 0 then $
        position = self->Count()-1
    self->SetProperty, current=position
    end $
else begin
    if obj_valid(oCurrent) then begin
        void = self->IsContained(oCurrent, position=position)
        self->SetProperty, current=position
        end
    end

end
;--------------------------------------------------------------------
pro IDLexViewGroup::Remove, $
    oObj, $     ; IN: (opt)
    _extra=e

compile_opt hidden

oCurrent = self->Get(/current)
;
;Perform a regular Remove.
;
if n_params() eq 1 then begin
    self->IDLgrViewGroup::Remove, oObj, _extra=e
    end $
else begin
    self->IDLgrViewGroup::Remove, _extra=e
    end
;
;Update self.current.
;
if obj_valid(oCurrent) then begin
    void = self->IsContained(oCurrent, position=position)
    self.current = position
    end

end
;--------------------------------------------------------------------
function IDLexViewGroup::Get, $
    current=current, $      ; IN: (opt) if set, get the current object.
    _extra=e, $             ; e.g. ALL
    count=count


compile_opt hidden

if keyword_set(current) then begin
    if self.current eq -1 then $
        return, -1
    return, self->IDLgrViewGroup::Get(position=self.current, _extra=e, count=count)
    end $
else begin
    return, self->IDLgrViewGroup::Get(_extra=e, count=count)
    end

end
;--------------------------------------------------------------------
pro IDLexViewGroup::CallAll, $
    method_name, $
    isa=isa, $
    p1, p2, p3, p4, p5, p6, p7, $
    _ref_extra=e

compile_opt hidden

oObj = self->Get(/all, count=count, isa=isa)
for i=0L,count-1 do begin
    if n_elements(e) gt 0 then begin
        case n_params()-1 of
            0: call_method, method_name, oObj[i], _extra=e
            1: call_method, method_name, oObj[i], p1, _extra=e
            2: call_method, method_name, oObj[i], p1, p2, _extra=e
            3: call_method, method_name, oObj[i], p1, p2, p3, _extra=e
            4: call_method, method_name, oObj[i], p1, p2, p3, p4, _extra=e
            5: call_method, method_name, oObj[i], p1, p2, p3, p4, p5, $
                _extra=e
            6: call_method, method_name, oObj[i], p1, p2, p3, p4, p5, p6, $
                _extra=e
            7: call_method, method_name, oObj[i], p1, p2, p3, p4, p5, p6, $
                p7, _extra=e
            else: message, 'limit of 7 parameters.' ; Arbitrary.
            endcase
        endif $
    else begin
        case n_params()-1 of
            0: call_method, method_name, oObj[i]
            1: call_method, method_name, oObj[i], p1
            2: call_method, method_name, oObj[i], p1, p2
            3: call_method, method_name, oObj[i], p1, p2, p3
            4: call_method, method_name, oObj[i], p1, p2, p3, p4
            5: call_method, method_name, oObj[i], p1, p2, p3, p4, p5
            6: call_method, method_name, oObj[i], p1, p2, p3, p4, p5, p6
            7: call_method, method_name, oObj[i], p1, p2, p3, p4, p5, p6, p7
            else: message, 'limit of 7 parameters.' ; Arbitrary.
            endcase
        endelse
    end
end
;--------------------------------------------------------------------
pro IDLexViewGroup::CallCurrent, $
    method_name, $
    p1, p2, p3, p4, p5, p6, p7, $
    isa=isa, $
    _ref_extra=e
;
;Purpose: Call a given method on self's current child.
;
compile_opt hidden

oCurrent = self->Get(/current)
if obj_valid(oCurrent) then begin

    if n_elements(isa) gt 0 then begin
        if not obj_isa(oCurrent, isa) then $
            return
        end

    if n_elements(e) gt 0 then begin
        case n_params()-1 of
            0: call_method, method_name, oCurrent, _extra=e
            1: call_method, method_name, oCurrent, p1, _extra=e
            2: call_method, method_name, oCurrent, p1, p2, _extra=e
            3: call_method, method_name, oCurrent, p1, p2, p3, _extra=e
            4: call_method, method_name, oCurrent, p1, p2, p3, p4, _extra=e
            5: call_method, method_name, oCurrent, p1, p2, p3, p4, p5, $
                _extra=e
            6: call_method, method_name, oCurrent, p1, p2, p3, p4, p5, p6, $
                _extra=e
            7: call_method, method_name, oCurrent, p1, p2, p3, p4, p5, p6, $
                p7, _extra=e
            else: message, 'limit of 7 parameters.' ; Arbitrary.
            endcase
        endif $
    else begin
        case n_params()-1 of
            0: call_method, method_name, oCurrent
            1: call_method, method_name, oCurrent, p1
            2: call_method, method_name, oCurrent, p1, p2
            3: call_method, method_name, oCurrent, p1, p2, p3
            4: call_method, method_name, oCurrent, p1, p2, p3, p4
            5: call_method, method_name, oCurrent, p1, p2, p3, p4, p5
            6: call_method, method_name, oCurrent, p1, p2, p3, p4, p5, p6
            7: call_method, method_name, oCurrent, p1, p2, p3, p4, p5, p6, p7
            else: message, 'limit of 7 parameters.' ; Arbitrary.
            endcase
        endelse
    end
end
;--------------------------------------------------------------------
function IDLexViewGroup::CallCurrent, $
    method_name, $
    p1, p2, p3, p4, p5, p6, p7, $
    _ref_extra=e

compile_opt hidden

oCurrent = self->Get(/current)
if obj_valid(oCurrent) then begin
    if n_elements(e) gt 0 then begin
        case n_params()-1 of
            0: return, $
                call_method(method_name, oCurrent, _extra=e)
            1: return, $
                call_method(method_name, oCurrent, p1, _extra=e)
            2: return, $
                call_method(method_name, oCurrent, p1, p2, _extra=e)
            3: return, $
                call_method(method_name, oCurrent, p1, p2, p3, _extra=e)
            4: return, $
                call_method(method_name, oCurrent, p1, p2, p3, p4, _extra=e)
            5: return, $
                call_method(method_name, oCurrent, p1, p2, p3, p4, p5, _extra=e)
            6: return, $
                call_method(method_name, oCurrent, p1, p2, p3, p4, p5, p6, _extra=e)
            7: return, $
                call_method(method_name, oCurrent, p1, p2, p3, p4, p5, p6, p7, _extra=e)
            else: message, 'limit of 7 parameters.' ; Arbitrary.
            endcase
        endif $
    else begin
        case n_params()-1 of
            0: return, $
                call_method(method_name, oCurrent)
            1: return, $
                call_method(method_name, oCurrent, p1)
            2: return, $
                call_method(method_name, oCurrent, p1, p2)
            3: return, $
                call_method(method_name, oCurrent, p1, p2, p3)
            4: return, $
                call_method(method_name, oCurrent, p1, p2, p3, p4)
            5: return, $
                call_method(method_name, oCurrent, p1, p2, p3, p4, p5)
            6: return, $
                call_method(method_name, oCurrent, p1, p2, p3, p4, p5, p6)
            7: return, $
                call_method(method_name, oCurrent, p1, p2, p3, p4, p5, p6, p7)
            else: message, 'limit of 7 parameters.' ; Arbitrary.
            endcase
        endelse
    endif
end
;--------------------------------------------------------------------
function IDLexViewGroup::DialogWriteVRML, $
    oWindow, $  ; IN (opt) VRML will use same dimensions, etc. as oWindow.
    dialog_parent=dialog_parent
;
;Purpose: Output self's current child to VRML.  The current child must be
;   a view, else an error is thrown.
;
compile_opt hidden

oCurrent = self->Get(/current)

if not obj_valid(oCurrent) then $
    message, 'Current object is not valid.'
if not obj_isa(oCurrent, 'IDLgrView') then $
    message, 'Current object is not an IDLgrView'

filename = dialog_pickfile( $
    /write, $
    file='untitled.wrl', $
    dialog_parent=dialog_parent, $
    filter='*.wrl' $
    )
if filename ne '' then begin
    if obj_valid(oWindow) then begin
        oWindow->GetProperty, $
            dimensions=dimensions, $
            resolution=resolution,$
            units=units, $
            color_model=color_model, $
            n_colors=n_colors
        oVRML = obj_new('IDLgrVRML', $
            dimensions=dimensions, $
            resolution=resolution, $
            units=units, $
            color_model=color_model, $
            n_colors=n_colors $
            )
        end $
    else begin
        oVRML = obj_new('IDLgrVRML') ; OK?
        end
    oVRML->SetProperty, filename=filename
    oVRML->Draw, oCurrent

    obj_destroy, oVRML
    return, 1
    end $
else begin
    return, 0
    end

end
;--------------------------------------------------------------------
pro IDLexViewGroup::SetCurrent, event, _extra=e
;
;Purpose: select and store a current view indicated by a given
;   WIDGET_DRAW event.
;
compile_opt hidden

if tag_names(event, /structure_name) ne 'WIDGET_DRAW' then $
    return

widget_control, event.id, get_value=oWindow

if not obj_valid(oWindow) then $
    message, 'Event is not from a valid object.'

oSelectedView = oWindow->Select(self, [event.x, event.y], _extra=e)

void = self->IsContained(oSelectedView[0], position=position)
self.current = position

end
;--------------------------------------------------------------------
pro IDLexViewGroup::SetProperty, $
    current=current, $ ; IN: (opt) an IDL_Container POSITION in self.
    _extra=e

compile_opt hidden

if n_elements(current) ne 0 then begin
    if current ge self->Count() or current lt 0 then $
        message, 'specified CURRENT position is out of range.'
    self.current = current
    end

self->IDLgrViewGroup::SetProperty, _extra=e
end
;--------------------------------------------------------------------
pro IDLexViewGroup::GetProperty, $
    current=current ; OUT: (opt) an IDL_Container POSITION in self.
    _ref_extra=e

compile_opt hidden

current = self.current
self->IDLgrViewGroup::GetProperty, _extra=e
end
;--------------------------------------------------------------------
function IDLexViewGroup::DialogWriteImage, $
    oDevice, $               ; IN: (opt) Use this device's pixel size.
    dimensions=dimensions, $ ; IN: (opt) Use this pixel size. [x, y].
    _ref_extra=e             ; e.g. DIALOG_PARENT (IN), OPTIONS (OUT)
;
;Purpose: Export self's graphics to an image.  The resulting image
;   will be DIMENSIONS pixel size if that keyword is supplied, else
;   it will have oDevice's pixel size if that argument is supplied.
;   oDevice is typically an IDLgrWindow.
;
compile_opt hidden

if n_elements(dimensions) eq 0 and obj_valid(oDevice) then begin
;
;   Obtain DIMENSIONS from oDevice.
;
    oDevice->GetProperty, units=orig_units
    oDevice->SetProperty, units=0
    oDevice->GetProperty, dimensions=dimensions
    oDevice->SetProperty, units=orig_units
    end
;
;;;;;;;;;;;;;;;;;;;;;;;;; Added by BT
;print, "Fmt: ", format
;ext = (stregex('\.(tiff?|bmp|jpe?g|png|gif|srf|ppm)$', filename, /extract))[1]
;print, "ext"
;if (ext eq '') then begin
;   case format of
;      'bmp': ext = '.bmp'
;      'jpg': ext = '.jpg'
;      'tif': ext = '.tif'
;      'gif': ext = '.gif'
;      'ppm': ext = '.ppm'
;      'png': ext = '.png'
;      'srf': ext = '.srf'
;      else:
;   endcase
;   filename += ext
;endif
;;;;;;;;;;;;;;;;;;;;;;;;;; End of added by BT

oBuff = obj_new('IDLgrBuffer', dimensions=dimensions)
oBuff->Draw, self
oBuff->GetProperty, image_data=image_data
result = dialog_write_image(image_data, _extra=e)
obj_destroy, oBuff

return, result
end
;--------------------------------------------------------------------
pro IDLexViewGroup::WriteImage, $
    filename, $             ; IN
    format, $               ; IN
    oDevice, $              ; IN: (opt) Use this device's pixel size.
    dimensions=dimensions   ; IN: (opt) Use this pixel size. [x, y].

compile_opt hidden

if n_elements(dimensions) eq 0 and obj_valid(oDevice) then begin
;
;   Obtain DIMENSIONS from oDevice.
;
    oDevice->GetProperty, units=orig_units
    oDevice->SetProperty, units=0
    oDevice->GetProperty, dimensions=dimensions
    oDevice->SetProperty, units=orig_units
    end
;
oBuff = obj_new('IDLgrBuffer', dimensions=dimensions)
oBuff->Draw, self
oBuff->GetProperty, image_data=image_data

write_image, filename, format, image_data
obj_destroy, oBuff
end
;--------------------------------------------------------------------
pro IDLexViewGroup::CopyToClipboard, $
    oDevice, $  ; IN: (opt) Clippboard takes on properties of this device.
    dimensions=dimensions, $    ; IN: (opt)
    resolution=resolution, $    ; IN: (opt)
    color_model=color_model, $  ; IN: (opt)
    palette=palette, $          ; IN: (opt)
    units=units, $              ; IN: (opt)
    n_colors=n_colors, $        ; IN: (opt)
    vector=vector, $            ; IN: (opt)
    postscript=postcript        ; IN: (opt)
;
;Purpose: Export self's graphics to the operating system's clipboard.
;   The clipbaord will take on DIMENSIONS if that keyword is supplied,
;   else it will take on the dimensions of oDevice if that argument
;   is supplied.  The same behavior applies for RESOULTION, etcetera.
;   oDevice is typically an IDLgrWindow.
;
compile_opt hidden

if n_elements(dimensions) eq 0 and obj_valid(oDevice) then begin
    oDevice->GetProperty, dimensions=dimensions
    end
if n_elements(resolution) eq 0 and obj_valid(oDevice) then begin
    oDevice->GetProperty, resolution=resolution
    end
if n_elements(color_model) eq 0 and obj_valid(oDevice) then begin
    oDevice->GetProperty, color_model=color_model
    end
if n_elements(palette) eq 0 and obj_valid(oDevice) then begin
    oDevice->GetProperty, palette=palette
    end
if n_elements(units) eq 0 and obj_valid(oDevice) then begin
    oDevice->GetProperty, units=units
    end
if n_elements(n_colors) eq 0 and obj_valid(oDevice) then begin
    oDevice->GetProperty, n_colors=n_colors
    end

oClipboard = obj_new('IDLgrClipboard', $
    dimensions=dimensions, $
    resolution=resolution, $
    color_model=color_model, $
    palette=palette, $
    units=units, $
    n_colors=n_colors $
    )
oClipboard->Draw, self, vector=vector, postscript=postscript
obj_destroy, oClipboard
end
;--------------------------------------------------------------------
pro IDLexViewGroup::DialogPrint, $
    oDevice, $              ; IN: (opt) Typically an IDLgrWindow.
    hourglass=hourglass, $  ; IN: (opt) If set, set the cursor to hourglass.
    _extra=e                ; e.g. DIALOG_PARENT, etc.
;
;Purpose: Export self's graphics to hardcopy.  If argument oDevice
;   is supplied, make the printout appear approximately the same size
;   as the graphic would (or does) appear on oDevice.
;
;Bug: Sizing of the graphic output to match the appearance on oDevice
;   is applied only to views that are, or are subclasses of,
;   class IDLexInscribingView.
;
compile_opt hidden

oPrinter = obj_new('IDLgrPrinter')
if dialog_printersetup(oPrinter, _extra=e) then begin
    if dialog_printjob(oPrinter, _extra=e) then begin
        if keyword_set(hourglass) then $
            widget_control, /hourglass
;
;       Convert views to inches.
;
        oViews = self->Get(/all, count=count, isa='IDLexInscribingView')

        default_arr = bytarr(count > 1)
        unit_arr = bytarr(count > 1)
        dim_arr = fltarr(2, count > 1)
        loc_arr = fltarr(2, count > 1)

        for i=0,count-1 do begin
            oViews[i]->GetProperty, $
                units=units, $
                dimensions=dim, $
                location=loc
            unit_arr[i] = units
            dim_arr[0, i] = dim
            loc_arr[0, i] = loc

            dim = oViews[i]->GetViewportDimensions( $
                oDevice, $
                location=loc, $
                /inches, $
                defaulted=defaulted $
                )
            default_arr[i] = defaulted

            oViews[i]->SetProperty,$
                location=loc, $
                dimensions=dim, $
                units=1 ; Inches.
            endfor
;
;       Print.
;
        oPrinter->Draw, self
        oPrinter->NewDocument
;
;       Restore views to their original units.
;
        for i=0,count-1 do begin
            oViews[i]->SetProperty,$
                location=default_arr[i] ? [0, 0] : loc_arr[*, i], $
                dim=default_arr[i] ? [0, 0] : dim_arr[*, i], $
                units=unit_arr[i]

            endfor
        endif
    endif
obj_destroy,oPrinter
end
;--------------------------------------------------------------------
function IDLexViewGroup::init, debug=debug, _extra=e

compile_opt hidden

catch, error_stat
if error_stat ne 0 then begin
    catch, /cancel
    message, $
        !error_state.msg + ' ' + !error_state.sys_msg, $
        /informational
    return, 0
    end
if keyword_set(debug) then $
   catch, /cancel

if not self->IDLgrViewGroup::init(_extra=e) then $
    message, 'Failed to init IDLgrViewGroup part of self.', /noname

self.debug = keyword_set(debug)
self.current = -1
return, 1 ; Success.
end
;--------------------------------------------------------------------
pro IDLexViewGroup__define

compile_opt hidden
on_error, 2

struct_hide, {IDLexViewGroup, $
    inherits IDLgrViewGroup, $
    debug: 0b, $
    current: 0L $  ; IDL_Container POSITION of current child.
    }
end
