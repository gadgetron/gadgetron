/*
 * hostutils_export.h
 *
 *  Created on: Nov 18, 2011
 *      Author: Michael S. Hansen
 */

#ifndef HOSTUTILS_EXPORT_H_
#define HOSTUTILS_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_HOSTUTILS__) || defined (gadgetron_toolbox_hostutils_EXPORTS)
#define EXPORTHOSTUTILS __declspec(dllexport)
#else
#define EXPORTHOSTUTILS __declspec(dllimport)
#endif
#else
#define EXPORTHOSTUTILS
#endif


#endif /* HOSTUTILS_EXPORT_H_ */
