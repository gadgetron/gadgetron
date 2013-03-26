/*
 * gadgetroncore_export.h
 *
 *  Created on: Jan 28, 2013
 *      Author: Michael S. Hansen
 */

#ifndef GADGETRONMATLAB_EXPORT_H_
#define GADGETRONMATLAB_EXPORT_H_


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_MATLAB__) || defined (gadgetronmatlab_EXPORTS)
#define EXPORTGADGETSMATLAB __declspec(dllexport)
#else
#define EXPORTGADGETSMATLAB __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSMATLAB
#endif


#endif /* GADGETRONMATLAB_EXPORT_H_ */
