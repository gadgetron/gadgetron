/** \file gadgettools_export.h
    \brief Required definitions for Windows, importing/exporting dll symbols 
*/

#ifndef GADGETTOOLS_EXPORT_H_
#define GADGETTOOLS_EXPORT_H_


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GADGETTOOLS__) || defined (gadgetron_toolbox_gadgettools_EXPORTS)
#define EXPORTGADGETTOOLS __declspec(dllexport)
#else
#define EXPORTGADGETTOOLS __declspec(dllimport)
#endif
#else
#define EXPORTGADGETTOOLS
#endif


#endif /* GADGETTOOLS_EXPORT_H_ */
