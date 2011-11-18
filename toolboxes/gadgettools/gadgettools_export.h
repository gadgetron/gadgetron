/*
 * gadgettools_export.h
 *
 *  Created on: Nov 18, 2011
 *      Author: Michael S. Hansen
 */

#ifndef GADGETTOOLS_EXPORT_H_
#define GADGETTOOLS_EXPORT_H_


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GADGETTOOLS__) || defined (gadgettools_EXPORTS)
#define EXPORTGADGETTOOLS __declspec(dllexport)
#else
#define EXPORTGADGETTOOLS __declspec(dllimport)
#endif
#else
#define EXPORTGADGETTOOLS
#endif


#endif /* GADGETTOOLS_EXPORT_H_ */
