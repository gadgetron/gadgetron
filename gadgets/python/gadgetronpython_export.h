/*
 * gadgetronpython_export.h
 *
 *  Created on: Nov 18, 2011
 *      Author: Michael S. Hansen
 */

#ifndef GADGETRONPYTHON_EXPORT_H_
#define GADGETRONPYTHON_EXPORT_H_


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_PYTHON__) || defined (gadgetronpython_EXPORTS)
#define EXPORTGADGETSPYTHON __declspec(dllexport)
#else
#define EXPORTGADGETSPYTHON __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSPYTHON
#endif



#endif /* GADGETRONPYTHON_EXPORT_H_ */
