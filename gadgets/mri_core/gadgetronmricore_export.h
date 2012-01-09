/*
 * gadgetronmricore_export.h
 *
 *  Created on: Dec 9, 2011
 *      Author: Michael S. Hansen
 */

#ifndef GADGETRONMRICORE_EXPORT_H_
#define GADGETRONMRICORE_EXPORT_H_


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_MRICORE__) || defined (gadgetronmricore_EXPORTS)
#define EXPORTGADGETSMRICORE __declspec(dllexport)
#else
#define EXPORTGADGETSMRICORE __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSMRICORE
#endif


#endif /* GADGETRONMRICORE_EXPORT_H_ */
