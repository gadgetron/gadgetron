set INSTALL_DIR=%~dp0

: set the path of gadgetron
setx PATH "%PATH%;%INSTALL_DIR%\..\lib;%INSTALL_DIR%\..\bin"
: copy the gadgetron.xml file
copy /Y %INSTALL_DIR%\..\config\gadgetron.xml.example %INSTALL_DIR%\..\config\gadgetron.xml