# Auto-generated by EclipseNSIS Script Wizard
# Sep 28, 2006 5:54:28 PM

# Modified by Ugo Varetto
#
# This script expects to be run in the build directory
# after Molekel has been built.
#
#   setup.nsi
#   dist
#      bin
#        Molekel.exe
#      resources
#        ...
#      doc
#        ...
#



;-------------------------------
; Test if Visual Studio Redistributables 2005+ SP1 installed
; Returns -1 if there is no VC redistributables intstalled
Function CheckVCRedist
   Push $R0
   ClearErrors
   ReadRegDword $R0 HKLM "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\{3C3D696B-0DB7-3C6D-A356-3DB8CE541918}" "Version"

   ; if VS 2005+ redist SP1 not installed, install it
   IfErrors 0 VSRedistInstalled
   StrCpy $R0 "-1"

VSRedistInstalled:
   Exch $R0
FunctionEnd



Name Molekel
# Defines
!define REGKEY "SOFTWARE\$(^Name)"
!define VERSION 5.4.0.8
!define COMPANY "Swiss National Supercomputing Centre"
!define URL http://www.cscs.ch

# Included files
!include Sections.nsh

# Reserved Files
ReserveFile "${NSISDIR}\Plugins\StartMenu.dll"

# Variables
Var StartMenuGroup

# Installer pages
Page license
Page directory
Page custom StartMenuGroupSelect "" ": Start Menu Folder"
Page instfiles

# Installer attributes
OutFile setup.exe
InstallDir $PROGRAMFILES\Molekel
CRCCheck on
XPStyle on
Icon Molekel.ico 
#"${NSISDIR}\Contrib\Graphics\Icons\classic-install.ico"
ShowInstDetails show
AutoCloseWindow false
LicenseData dist\license.txt
# Better without full-screen gradient
BGGradient 000080 000000 FFFFFF
VIProductVersion 5.4.0.8
VIAddVersionKey ProductName Molekel
VIAddVersionKey ProductVersion "${VERSION}"
VIAddVersionKey CompanyName "${COMPANY}"
VIAddVersionKey CompanyWebsite "${URL}"
VIAddVersionKey FileVersion ""
VIAddVersionKey FileDescription ""
VIAddVersionKey LegalCopyright "Copyright (c) 2006, 2007, 2008, 2009 Swiss National Supercomputing Centre"
InstallDirRegKey HKLM "${REGKEY}" Path
UninstallIcon "${NSISDIR}\Contrib\Graphics\Icons\classic-uninstall.ico"
ShowUninstDetails show

# Installer sections
Section -Main SEC0000
    SetOutPath $INSTDIR
    SetOverWrite on
    File dist\notes.txt
    File dist\license.txt
    SetOutPath $INSTDIR\bin
    SetOverwrite on
    File /r dist\bin\*
    SetOutPath $INSTDIR\resources
    SetOverwrite on
    File /r dist\resources\*
 #   SetOutPath $INSTDIR\doc
 #   SetOverwrite on
 #   File /r dist\doc\*
    SetOutPath $INSTDIR\data
    SetOverwrite on
    File /r dist\data\*
    SetOutPath $INSTDIR\shaders
    SetOverwrite on
    File /r dist\shaders\*
    WriteRegStr HKLM "${REGKEY}\Components" Main 1
SectionEnd

Section -post SEC0001
    WriteRegStr HKLM "${REGKEY}" Path $INSTDIR
    WriteRegStr HKLM "${REGKEY}" StartMenuGroup $StartMenuGroup
    WriteUninstaller $INSTDIR\uninstall.exe
    SetOutPath $SMPROGRAMS\$StartMenuGroup
    CreateShortcut "$SMPROGRAMS\$StartMenuGroup\Uninstall $(^Name).lnk" $INSTDIR\uninstall.exe
    WriteRegStr HKLM "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\$(^Name)" DisplayName "$(^Name)"
    WriteRegStr HKLM "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\$(^Name)" DisplayVersion "${VERSION}"
    WriteRegStr HKLM "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\$(^Name)" Publisher "${COMPANY}"
    WriteRegStr HKLM "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\$(^Name)" URLInfoAbout "${URL}"
    WriteRegStr HKLM "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\$(^Name)" DisplayIcon $INSTDIR\uninstall.exe
    WriteRegStr HKLM "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\$(^Name)" UninstallString $INSTDIR\uninstall.exe
    WriteRegDWORD HKLM "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\$(^Name)" NoModify 1
    WriteRegDWORD HKLM "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\$(^Name)" NoRepair 1
    #this is required to have the shortcut start the application in the bin folder
    SetOutPath $INSTDIR\bin
    CreateShortcut "$SMPROGRAMS\$StartMenuGroup\$(^Name).lnk" $INSTDIR\bin\Molekel.exe "" "$INSTDIR\resources\Molekel.ico"
    # VC++ redistributable
    Call CheckVCRedist
    Pop $R0
    StrCmp $R0 "-1" "" DontInstallVCRedist
    ExecWait '"$INSTDIR\bin\vc9redist_x86.exe" /q:a /c:"VC9RED~1.EXE /q:a /c:""msiexec /i vc_red.msi /qb!"" "'
DontInstallVCRedist:    
SectionEnd


# Macro for selecting uninstaller sections
!macro SELECT_UNSECTION SECTION_NAME UNSECTION_ID
    Push $R0
    ReadRegStr $R0 HKLM "${REGKEY}\Components" "${SECTION_NAME}"
    StrCmp $R0 1 0 next${UNSECTION_ID}
    !insertmacro SelectSection "${UNSECTION_ID}"
    GoTo done${UNSECTION_ID}
next${UNSECTION_ID}:
    !insertmacro UnselectSection "${UNSECTION_ID}"
done${UNSECTION_ID}:
    Pop $R0
!macroend

# Uninstaller sections
Section /o un.Main UNSEC0000
    RmDir /r /REBOOTOK $INSTDIR\resources
    RmDir /r /REBOOTOK $INSTDIR\bin
    RmDir /r /REBOOTOK $INSTDIR\doc
    RmDir /r /REBOOTOK $INSTDIR\data
    RmDir /r /REBOOTOK $INSTDIR\shaders
    DeleteRegValue HKLM "${REGKEY}\Components" Main
SectionEnd

Section un.post UNSEC0001
    DeleteRegKey HKLM "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\$(^Name)"
    Delete /REBOOTOK "$SMPROGRAMS\$StartMenuGroup\Uninstall $(^Name).lnk"
    Delete /REBOOTOK "$SMPROGRAMS\$StartMenuGroup\$(^Name).lnk"
    Delete /REBOOTOK $INSTDIR\uninstall.exe
    Delete /REBOOTOK $INSTDIR\bin\Molekel.exe
    Delete /REBOOTOK $INSTDIR\notes.txt
    Delete /REBOOTOK $INSTDIR\license.txt
    DeleteRegValue HKLM "${REGKEY}" StartMenuGroup
    DeleteRegValue HKLM "${REGKEY}" Path
    DeleteRegKey /IfEmpty HKLM "${REGKEY}\Components"
    DeleteRegKey /IfEmpty HKLM "${REGKEY}"
    RmDir /REBOOTOK $SMPROGRAMS\$StartMenuGroup
    RmDir /REBOOTOK $INSTDIR
SectionEnd

# Installer functions
Function StartMenuGroupSelect
    Push $R1
    StartMenu::Select /autoadd /text "Select the Start Menu folder in which you would like to create the program's shortcuts:" /lastused $StartMenuGroup Molekel
    Pop $R1
    StrCmp $R1 success success
    StrCmp $R1 cancel done
    MessageBox MB_OK $R1
    Goto done
success:
    Pop $StartMenuGroup
done:
    Pop $R1
FunctionEnd

Function .onInstSuccess
    ExecShell open $INSTDIR\Notes.txt
FunctionEnd
Function .onInit
    InitPluginsDir
FunctionEnd

# Uninstaller functions
Function un.onInit
    ReadRegStr $INSTDIR HKLM "${REGKEY}" Path
    ReadRegStr $StartMenuGroup HKLM "${REGKEY}" StartMenuGroup
    !insertmacro SELECT_UNSECTION Main ${UNSEC0000}
FunctionEnd


