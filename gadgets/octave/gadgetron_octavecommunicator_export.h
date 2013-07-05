/*
 * gadgetroncore_export.h
 *
 *  Created on: Jan 28, 2013
 *      Author: Michael S. Hansen
 */

#ifndef GADGETRONOCTAVECOMMUNICATOR_EXPORT_H_
#define GADGETRONOCTAVECOMMUNICATOR_EXPORT_H_


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_OCTAVECOMMUNICATOR__) || defined (gadgetronOctaveCommunicator_EXPORTS)
#define EXPORTGADGETSOCTAVECOMMUNICATOR __declspec(dllexport)
#else
#define EXPORTGADGETSOCTAVECOMMUNICATOR __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSOCTAVECOMMUNICATOR
#endif


#endif /* GADGETRONOCTAVE_EXPORT_H_ */
