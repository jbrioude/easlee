' ===========================================================================
'  EASLEE - one double-click launcher (Windows), no typing required.
'  Double-click EASLEE.vbs:  the model starts and your browser opens at
'  http://localhost:8080 .  The server runs in a MINIMISED window on the
'  taskbar (not a terminal in your face) - click it and close it to stop EASLEE.
'
'  Tip: right-click EASLEE.vbs -> Send to -> Desktop (create shortcut), then
'  rename the shortcut "EASLEE" and give it an icon for a real app feel.
' ===========================================================================
Dim shell, fso, here
Set shell = CreateObject("WScript.Shell")
Set fso   = CreateObject("Scripting.FileSystemObject")

' Run from the folder this script lives in.
here = fso.GetParentFolderName(WScript.ScriptFullName)
shell.CurrentDirectory = here

' Window style 7 = minimised (taskbar), not activated. Use 0 for fully hidden.
' (Kept minimised, not hidden, so you still have a way to stop the server.)
shell.Run "cmd /c """ & here & "\run_gui.bat""", 7, False
