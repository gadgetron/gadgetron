/*
 * gadgetroncore_export.h
 *
 *  Created on: Nov 18, 2011
 *      Author: Michael S. Hansen
 */

#ifndef GADGETRONCORE_EXPORT_H_
#define GADGETRONCORE_EXPORT_H_


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_CORE__) || defined (gadgetroncore_EXPORTS)
#define EXPORTGADGETSCORE __declspec(dllexport)
#else
#define EXPORTGADGETSCORE __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSCORE
#endif


#endif /* GADGETRONCORE_EXPORT_H_ */
