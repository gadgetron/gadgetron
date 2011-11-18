/*
 * gadgetroncgsense_export.h
 *
 *  Created on: Nov 18, 2011
 *      Author: Michael S. Hansen
 */

#ifndef GADGETRONCGSENSE_EXPORT_H_
#define GADGETRONCGSENSE_EXPORT_H_


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_CGSENSE__) || defined (gadgetroncgsense_EXPORTS)
#define EXPORTGADGETSCGSENSE __declspec(dllexport)
#else
#define EXPORTGADGETSCGSENSE __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSCGSENSE
#endif


#endif /* GADGETRONCGSENSE_EXPORT_H_ */
