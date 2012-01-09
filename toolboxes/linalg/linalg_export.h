/*
 * hostutils_export.h
 *
 *  Created on: Dec 9, 2011
 *      Author: Michael S. Hansen
 */

#ifndef LINALG_EXPORT_H_
#define LINALG_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_LINALG__) || defined (linalg_EXPORTS)
#define EXPORTLINALG __declspec(dllexport)
#else
#define EXPORTLINALG __declspec(dllimport)
#endif
#else
#define EXPORTLINALG
#endif


#endif /* LINALG_EXPORT_H_ */
