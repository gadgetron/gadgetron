/** \file gpuxray_export.h
    \brief Required definitions for Windows, importing/exporting dll symbols 
*/

#ifndef GPUXRAY_EXPORT_H_
#define GPUXRAY_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GPUXRAY__) || defined (gpuxray_EXPORTS)
#define EXPORTGPUXRAY __declspec(dllexport)
#else
#define EXPORTGPUXRAY __declspec(dllimport)
#endif
#else
#define EXPORTGPUXRAY
#endif


#endif /* GPUXRAY_EXPORT_H_ */
