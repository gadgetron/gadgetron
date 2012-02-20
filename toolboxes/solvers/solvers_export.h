/*
 * solvers_export.h
 *
 *  Created on: Nov 18, 2011
 *      Author: Michael S. Hansen
 */

#ifndef SOLVERS_EXPORT_H_
#define SOLVERS_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_SOLVERS__) || defined (solvers_EXPORTS)
#define EXPORTSOLVERS __declspec(dllexport)
#else
#define EXPORTSOLVERS __declspec(dllimport)
#endif
#else
#define EXPORTSOLVERS
#endif


#endif /* SOLVERS_EXPORT_H_ */
