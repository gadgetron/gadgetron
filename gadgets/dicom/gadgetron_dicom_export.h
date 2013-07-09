#ifndef DICOM_EXPORT_H_
#define DICOM_EXPORT_H_


#if defined (WIN32)
#if defined (gadgetron_dicom_EXPORTS)
#define EXPORTGADGETSDICOM __declspec(dllexport)
#else
#define EXPORTGADGETSDICOM __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSDICOM
#endif

#endif /* DICOM_EXPORT_H_ */
