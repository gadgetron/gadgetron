/*
 * gadgetroncore_export.h
 *
 *  Created on: Jan 28, 2013
 *      Author: Michael S. Hansen
 */

#ifndef GADGETRONOCTAVE_EXPORT_H_
#define GADGETRONOCTAVE_EXPORT_H_


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_OCTAVE__) || defined (gadgetronoctave_EXPORTS)
#define EXPORTGADGETSOCTAVE __declspec(dllexport)
#else
#define EXPORTGADGETSOCTAVE __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSOCTAVE
#endif


#endif /* GADGETRONOCTAVE_EXPORT_H_ */
