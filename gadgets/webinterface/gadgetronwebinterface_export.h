/*
 * gadgetroncore_export.h
 *
 *  Created on: Nov 18, 2011
 *      Author: Michael S. Hansen
 */

#ifndef GADGETRONCORE_EXPORT_H_
#define GADGETRONCORE_EXPORT_H_


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_WEBINTERFACE__) || defined (gadgetronwebinterface_EXPORTS)
#define EXPORTGADGETSWEBINTERFACE __declspec(dllexport)
#else
#define EXPORTGADGETSWEBINTERFACE __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSWEBINTERFACE
#endif


#endif /* GADGETRONCORE_EXPORT_H_ */
