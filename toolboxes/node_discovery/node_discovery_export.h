#ifndef NODEDISCOVERY_EXPORT_H_
#define NODEDISCOVERY_EXPORT_H_

#if defined (WIN32)
   #if defined (__BUILD_GADGETRON_NODE_DISCOVERY__) || defined (gadgetron_toolbox_node_discovery__EXPORTS)
      #define EXPORTGADGETRONNODEDISCOVERY __declspec(dllexport)
   #else
      #define EXPORTGADGETRONNODEDISCOVERY __declspec(dllimport)
   #endif
#else
   #define EXPORTGADGETRONNODEDISCOVERY
#endif

#endif /* NODEDISCOVERY_EXPORT_H_ */
