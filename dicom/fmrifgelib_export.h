#ifndef FMRIFGE_EXPORT_H_
#define FMRIFGE_EXPORT_H_


#if defined (WIN32)
#if defined (gadgetronfmrifgelib_EXPORTS)
#define EXPORTGADGETSFMRIFGE __declspec(dllexport)
#else
#define EXPORTGADGETSFMRIFGE __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSFMRIFGE
#endif

#endif /* FMRIFGE_EXPORT_H_ */
