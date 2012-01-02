#include "GadgetronConnector.h"
#include "ui_imageClient.h"

#include "ace/OS_NS_string.h"
#include "ace/Log_Msg.h"
#include "ace/Get_Opt.h"


class MainWindow : public QMainWindow, public Ui::MainWindowBase
{
	Q_OBJECT

public:
	MainWindow( QWidget *parent = 0x0 );

	int setup_gadgetron_connection();
	void set_config(ACE_TCHAR* filename);
	inline void set_hostname(ACE_TCHAR* val){ hostname_ = val; }
	inline void set_port_no(ACE_TCHAR* val){ port_no_ = val; }
	inline void set_num_repetitions(unsigned int val){ num_repetitions_ = val; }

public slots:
   void addFiles();
   void clearFiles();
   void loadConfig();
   void saveConfig();
   void go();

private:
   ACE_TCHAR *hostname_;
   ACE_TCHAR *port_no_;
   unsigned int num_repetitions_;
   GadgetronConnector con_;
};
