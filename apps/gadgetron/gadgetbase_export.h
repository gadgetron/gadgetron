#ifndef GADGETBASE_EXPORT_H_
#define GADGETBASE_EXPORT_H_


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GADGETBASE__) || defined (gadgetron_gadgetbase_EXPORTS)
#define EXPORTGADGETBASE __declspec(dllexport)
#else
#define EXPORTGADGETBASE __declspec(dllimport)
#endif
#else
#define EXPORTGADGETBASE
#endif


#endif /* GADGETBASE_EXPORT_H_ */
