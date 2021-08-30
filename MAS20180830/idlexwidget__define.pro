;$Id: //depot/idl/IDL_64/idldir/lib/utilities/idlexwidget__define.pro#1 $
;
;  Copyright (c) 1997-2007, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
; NAME:
;   IDLexWidget
;
; PURPOSE:
;   Provide a base class for defining objects that function
;       as IDL compound widgets.
;
; CATEGORY:
;   Widgets.  Object Oriented Programming.
;
; CREATION:
;   It is intended that subclasses of class IDLexWidget, rather than
;   class "IDLexWidget" itself, get instantiated.  It is possible to
;   create a "IDLexWidget" object directly...
;
;       oWidget = OBJ_NEW('IDLexWidget', wParent)
;
;   ... but such objects would not have any buttons, sliders, etc.
;   in them, and would have a propensity to throw errors.
;
; METHODS:
;   IDLexWidget::Rewrite()
;       Return a structure that is appropriate for "rewriting" a given
;       event.  (Rewriting an event is described in IDL 5.x documentation
;       under "Event Processing and Callbacks.")  Typically,
;       self->Rewrite(event) is called by self's HandleEvent method.
;
;   IDLexWidget::GetProperty, $
;       WIDGET_ID=wID
;
;       Note: IDL widgets usually represent or control some value.  To
;       provide the ability to get this value from an IDLexWidget,
;       override IDLexWidget::GetProperty in your subclass so as to
;       return a value via keyword VALUE.  Then users of your
;       IDLexWidget subclass will be able to do this...
;
;           oWidget->GetProperty, WIDGET_ID=wID
;           WIDGET_CONTROL, wID, GET_VALUE=result
;
;       ...or this...
;
;           oWidget->GetProperty, VALUE=value
;
;       ...both of which are synonymous.
;
;   IDLexWidget::SetProperty
;
;       Note: IDL widgets usually represent or control some value.  To
;       provide the ability to set this value for an IDLexWidget,
;       override IDLexWidget::GetProperty in your subclass so as to
;       set a value via keyword VALUE.  Then users of your
;       IDLexWidget subclass will be able to do this...
;
;           oWidget->GetProperty, WIDGET_ID=wID
;           WIDGET_CONTROL, wID, SET_VALUE=value
;
;       ...or this...
;
;           oWidget->SetProperty, VALUE=value
;
;       ...both of which are synonymous.
;
;   IDLexWidget::HandleEvent()
;       This method is a no-op.  To provide the ability for your
;       IDLexWidget to handle GUI events, override function
;       IDLexWidget::HandleEvent(event) in a subclass.  Typically
;       your definition will return self->Rewrite(event) after
;       performing the desired custom operations.
;
; LIFECYCLE METHODS:
;   IDLexWidget::init()
;   IDLexWidget::cleanup
;
;--------------------------------------------------------------------
function IDLexWidget__HandleEvent, event
compile_opt hidden
;
;Purpose: provide a non-method version self's event handler.
;Use this function to attach self's event handler to a widget.
;
;Example use:
;   function Subclass_of_IDLexWidget::init, wParent...
;       .
;       .
;       .
;   widget_control, $
;       wid, $ ; Some widget ID, i.e. a child of self.wIDBase
;       event_func='IDLexWidget__HandleEvent, $
;       set_uvalue=self
;       .
;       .
;       .
;   end
;
;In the above example, we attach the event handler to a widget that
;is a child of self.wIDBase, rather than to self.wIDBase itself.
;This is so that a client of self (an "outside user") can attach an
;event handler to self's widget ID, if desired, without breaking
;the attachment that we have made.
;
on_error, 2 ; Return to caller on error.

widget_control, event.handler, get_uvalue=oObj
return, oObj->HandleEvent(event)
end
;----------------------------------------------------------------------
pro IDLexWidget__OnRealize, wid
compile_opt hidden

widget_control, wid, get_uvalue=oObj
oObj->OnRealize
end
;----------------------------------------------------------------------
pro IDLexWidget::cleanup
compile_opt hidden
end
;----------------------------------------------------------------------
function IDLexWidget__get_value, wid
compile_opt hidden

widget_control, widget_info(wid, /child), get_uvalue=oObj
oObj->GetProperty, value=value
return, value
end
;----------------------------------------------------------------------
pro IDLexWidget__set_value, wid, value
compile_opt hidden

widget_control, widget_info(wid, /child), get_uvalue=oObj
oObj->SetProperty, value=value
end
;----------------------------------------------------------------------
pro IDLexWidget::GetProperty, widget_id=widget_id
compile_opt hidden

widget_id=self.wIDBase
end
;----------------------------------------------------------------------
pro IDLexWidget::SetProperty, value=value ; No-op. See comments above.
compile_opt hidden
end
;----------------------------------------------------------------------
pro IDLexWidget::OnRealize
compile_opt hidden
end
;----------------------------------------------------------------------
function IDLexWidget::HandleEvent, event
compile_opt hidden
end
;----------------------------------------------------------------------
function IDLexWidget::Rewrite, event
compile_opt hidden
;
;Purpose: Fabricate an event suitable for return from self's event
;handler. The implementation provided here is rather generic, and thus
;will typically be overriden in subclasses.
;
return, {id: self.wIDBase, top: event.top, handler:0L}
end
;----------------------------------------------------------------------
function IDLexWidget::init, wParent, frame=frame
compile_opt hidden

catch, error_stat
if error_stat ne 0 then begin
    catch, /cancel
    MESSAGE, /INFO, !ERROR_STATE.MSG
    return, 0
    end

self.wIDBase = widget_base( $
    wParent, $
    frame=frame, $
    pro_set_value='IDLexWidget__set_value', $
    func_get_value='IDLexWidget__get_value' $
    )

return, 1 ; Success
end
;----------------------------------------------------------------------
pro IDLexWidget__define
compile_opt hidden
struct_hide, {IDLexWidget, $
;
;   IDL compound widgets are "black boxes" known to the outside world
;   by a widget ID.  If we think of self as a compound widget, "wIDBase"
;   is the widget ID by which self will be known to the outside world.
;
    wIDBase: 0L $
    }
end

