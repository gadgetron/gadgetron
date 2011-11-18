/*
 * gadgetrongrappa_export.h
 *
 *  Created on: Nov 18, 2011
 *      Author: Michael S. Hansen
 */

#ifndef GADGETRONGRAPPA_EXPORT_H_
#define GADGETRONGRAPPA_EXPORT_H_


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GRAPPA__) || defined (gadgetrongrappa_EXPORTS)
#define EXPORTGADGETSGRAPPA __declspec(dllexport)
#else
#define EXPORTGADGETSGRAPPA __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSGRAPPA
#endif

#endif /* GADGETRONGRAPPA_EXPORT_H_ */
